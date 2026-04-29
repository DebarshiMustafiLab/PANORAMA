import re
import math
import base64
import io
import os
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.subplots as sp
import requests
import dash
from itertools import zip_longest
import html as html_lib
from dash import html
from dash import Dash, Input, Output, State, dcc, html

GENE_NAME = "RB1"
CHROMOSOME = "chr13"


BASE_DIR = Path(__file__).resolve().parent

_data_dir_env = os.getenv("PANORAMA_DATA_DIR")
DATA_DIR = Path(_data_dir_env) if _data_dir_env else (BASE_DIR / "data")
if not DATA_DIR.is_absolute():
    DATA_DIR = (BASE_DIR / DATA_DIR).resolve()

_excel_env = os.getenv("PANORAMA_EXCEL_FILE")
_excel_basename = os.getenv("PANORAMA_EXCEL_BASENAME", "retinoblastoma_data_v1.xlsx")

# Candidate locations in priority order.
_candidates: list[Path] = []
if _excel_env:
    _candidates.append(Path(_excel_env))

# Expected repo location: ./data/<basename>
_candidates.append(DATA_DIR / _excel_basename)

# Common fallback: repo root
_candidates.append(BASE_DIR / _excel_basename)

# If the user didn't explicitly set a file, try to auto-detect a single .xlsx in DATA_DIR.
if not _excel_env:
    try:
        _xlsx = sorted(DATA_DIR.glob("*.xlsx"))
        if len(_xlsx) == 1:
            _candidates.append(_xlsx[0])
        elif len(_xlsx) > 1:
            preferred = [p for p in _xlsx if "retinoblastoma" in p.name.lower()]
            if preferred:
                _candidates.append(preferred[0])
    except Exception:
        pass

EXCEL_FILE: Path | None = None
for _p in _candidates:
    if _p is None:
        continue
    _p = Path(_p)
    if not _p.is_absolute():
        _p = (BASE_DIR / _p).resolve()
    if _p.exists():
        EXCEL_FILE = _p
        break

# Backward-compatibility: if you still have the legacy Windows path locally, use it automatically.
LEGACY_EXCEL_FILE = Path(r"~\VariantViewer\retinoblastoma_data_v1.xlsx")
if EXCEL_FILE is None and LEGACY_EXCEL_FILE.exists():
    EXCEL_FILE = LEGACY_EXCEL_FILE

INPUT_SHEET = os.getenv("PANORAMA_INPUT_SHEET", "all_fields")

if EXCEL_FILE is None:
    _dir_exists = DATA_DIR.exists()
    try:
        _dir_files = [p.name for p in sorted(DATA_DIR.glob("*"))][:50] if _dir_exists else []
    except Exception:
        _dir_files = []

    raise FileNotFoundError(
        "Curated data file not found.\n"
        f"Looked for: {', '.join(str(p) for p in _candidates if p)}\n"
        f"Resolved DATA_DIR: {DATA_DIR} (exists={_dir_exists}). Top files: {_dir_files}\n"
        "Fix options:\n"
        f"  1) Commit your Excel to the repo at data/{_excel_basename}\n"
        "  2) Or set PANORAMA_EXCEL_FILE to the correct path\n"
        "  3) Or set PANORAMA_DATA_DIR/PANORAMA_EXCEL_BASENAME accordingly"
    )

ENSEMBL_CANONICAL = "ENST00000267163"
RB1_GRCH38_START = 48_303_744
RB1_GRCH38_END = 48_599_436
HAVE_PREMAPPED_HG38 = True
UNIPROT_ACCESSION = "P06400"
AA_MAP_TSV = None
ENSEMBL_REST = "https://rest.ensembl.org"
PROMOTER_START = 48_303_707
PROMOTER_END = None

DMR_START = 48_318_500
DMR_END = 48_319_721
DMR_COLOR = "#9b59b6"  # purple
DMR_LABEL = "CpG85"

EXON_COLOR = "#3B516B"
INTRON_COLOR = "#D3D3D3"
PROMOTER_COLOR = "#BDC89E"

USER_VARIANT_COLOR = "#D14B25"

HYPER_METHYL_COLOR = "#d62728"
HYPO_METHYL_COLOR = "#1f77b4"
INHERITANCE_GRAY = '#9aa0a6'
USER_HAPLO_COLOR = '#111111'

GENOMIC_TICKFORMAT = ",d"

APP_BRAND = "PANORAMA"  # Top-right branding text

FIXED_SOURCE_COLORS = {
    "LOVD": "#3D6C88",
    "Publication": "#7B9041",
    "COSMIC": "#E4BC33",
    "Unknown": "#7f7f7f",
    "Multiple": "#141414",
}

FILTER_PANEL_STYLE = {
    "display": "flex",
    "flexWrap": "nowrap",
    "gap": "10px",
    "marginBottom": "16px",
    "alignItems": "flex-start",
    "justifyContent": "flex-start",
}

FILTER_BOX_STYLE = {
    "border": "1px solid #e1e5eb",
    "borderRadius": "8px",
    "padding": "6px 10px",
    "backgroundColor": "#f9fafb",
    "boxShadow": "0 1px 2px rgba(15, 23, 42, 0.06)",
    "minWidth": "150px",
}

FILTER_LABEL_STYLE = {
    "fontWeight": "600",
    "fontSize": "13px",
    "color": "#111827",
    "marginBottom": "4px",
    "display": "block",
}

CHECKLIST_INLINE_STYLE = {
    "display": "flex",
    "flexWrap": "wrap",
    "gap": "6px 12px",
    "marginTop": "4px",
}

MCONSEQ_LABEL_MAP = {
    "3_prime_UTR": "3' UTR",
    "5_prime_UTR": "5' UTR",
    "inframe_deletion": "inframe deletion",
    "inframe_indel": "inframe indel",
    "inframe_insertion": "inframe insertion",
    "splice_acceptor": "splice acceptor",
    "splice_donor": "splice donor",
    "stop_gained": "nonsense",
}


def mc_pretty_label(v: str) -> str:
    v = str(v)
    if v in MCONSEQ_LABEL_MAP:
        return MCONSEQ_LABEL_MAP[v]
    return v.replace("_", " ")


def mc_sort_key(v: str):
    v_low = str(v).lower()
    if v_low.startswith("3_prime"):
        group = 0
    elif v_low.startswith("5_prime"):
        group = 1
    elif v_low == "frameshift":
        group = 2
    elif v_low.startswith("inframe"):
        group = 3
    elif v_low in {"missense", "synonymous"}:
        group = 4
    elif v_low.startswith("splice"):
        group = 5
    elif v_low in {"start_lost", "stop_gained", "stop_lost"}:
        group = 6
    else:
        group = 7
    return (group, v_low)


def with_all_option(options, all_label="All"):
    return [{"label": all_label, "value": "__ALL__"}] + options


def fetch_cds_segments(transcript_id=ENSEMBL_CANONICAL):
    url = f"{ENSEMBL_REST}/overlap/id/{transcript_id}"
    params = {"feature": "cds"}
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, params=params, headers=headers, timeout=30)
    r.raise_for_status()
    cds_segments = r.json()

    records = []
    for seg in cds_segments:
        parent = seg.get("Parent")
        if parent is None or transcript_id not in str(parent):
            continue
        records.append(
            {
                "chrom": seg["seq_region_name"],
                "start": int(seg["start"]),
                "end": int(seg["end"]),
                "strand": int(seg["strand"]),
            }
        )

    df = pd.DataFrame(records)
    if df.empty:
        raise RuntimeError(f"No CDS segments returned for {transcript_id}")

    strand = df["strand"].value_counts().idxmax()
    df = df[df["strand"] == strand].copy()
    return df


def _split_multi_field(val):
    if val is None:
        return []
    # handle NaNs
    try:
        if math.isnan(val):
            return []
    except TypeError:
        pass

    parts = [p.strip() for p in str(val).split(";")]
    return [p for p in parts if p]


def _format_cdna_list(cdna_raw, max_len=3):
    items = _split_multi_field(cdna_raw)
    formatted = []
    for s in items:
        s = s.strip()
        if len(s) > max_len:
            formatted.append(html_lib.escape(s[:max_len] + "..."))
        else:
            formatted.append(html_lib.escape(s))
    return "<br>".join(formatted)


def _format_genomic_list(genomic_raw):
    items = _split_multi_field(genomic_raw)
    # Normalize chromosome like "chr13.0" -> "chr13" (common after Excel/Pandas dtype coercion)
    fixed_items = []
    for s in items:
        s = (s or "").strip()
        m = re.match(r'^(?P<chr>(?:chr)?[0-9A-Za-z\.]+)\s*:(?P<rest>.*)$', s)
        if m:
            ch = canonicalize_chrom(m.group("chr"))
            s = f"{ch}:{m.group('rest').lstrip()}"
        fixed_items.append(html_lib.escape(s))
    return "<br>".join(fixed_items)


def get_rb1_exons_from_ensembl(transcript_id=ENSEMBL_CANONICAL, use_cds=True):
    if use_cds:
        cds_df = fetch_cds_segments(transcript_id)

        # sort by transcript order (respect strand)
        strand = cds_df["strand"].value_counts().idxmax()
        cds_df = cds_df.sort_values("start", ascending=(strand == 1))

        exons = [
            (int(s), int(e))
            for s, e in zip(cds_df["start"].tolist(), cds_df["end"].tolist())
        ]
        return exons

    url = f"{ENSEMBL_REST}/lookup/id/{transcript_id}"
    params = {"expand": 1}
    headers = {"Content-Type": "application/json"}
    resp = requests.get(url, headers=headers, params=params, timeout=20)
    resp.raise_for_status()
    data = resp.json()
    exons = [(e["start"], e["end"]) for e in data.get("Exon", [])]

    strand = data.get("strand", 1)
    exons.sort(key=lambda e: e[0], reverse=(strand == -1))
    return exons


def build_protein_to_genome_map(cds_df, transcript_id=ENSEMBL_CANONICAL):
    strand = cds_df["strand"].iloc[0]

    if strand == 1:
        cds_df = cds_df.sort_values("start")
    else:
        cds_df = cds_df.sort_values("start", ascending=False)

    bases = []
    for _, row in cds_df.iterrows():
        start, end = int(row["start"]), int(row["end"])
        if strand == 1:
            coords = range(start, end + 1)
        else:
            coords = range(end, start - 1, -1)
        bases.extend(coords)

    chrom = cds_df["chrom"].iloc[0]
    num_codons = int(math.ceil(len(bases) / 3))
    rows = []

    for aa in range(1, num_codons + 1):
        i1 = (aa - 1) * 3
        i2 = i1 + 1
        i3 = i1 + 2
        base1 = bases[i1] if i1 < len(bases) else None
        base2 = bases[i2] if i2 < len(bases) else None
        base3 = bases[i3] if i3 < len(bases) else None

        rows.append(
            {
                "transcript_id": transcript_id,
                "aa_position": aa,
                "chrom": chrom,
                "strand": strand,
                "codon_base1_genomic_pos": base1,
                "codon_base2_genomic_pos": base2,
                "codon_base3_genomic_pos": base3,
            }
        )

    return pd.DataFrame(rows)


def load_or_build_aa_map():
    if AA_MAP_TSV:
        aa_map = pd.read_csv(AA_MAP_TSV, sep="\t")
    else:
        cds_df = fetch_cds_segments(ENSEMBL_CANONICAL)
        aa_map = build_protein_to_genome_map(cds_df, ENSEMBL_CANONICAL)

    return aa_map


def get_initiator_codon_bounds(aa_map):
    if aa_map is None or aa_map.empty:
        return None, None
    first = aa_map[aa_map["aa_position"] == 1]
    if first.empty:
        return None, None
    bases = first[
        ["codon_base1_genomic_pos", "codon_base2_genomic_pos", "codon_base3_genomic_pos"]
    ].iloc[0]
    coords = [int(b) for b in bases if not pd.isna(b)]
    if not coords:
        return None, None
    return min(coords), max(coords)


def build_protein_axis_ticks(aa_map, n_ticks=5):
    if aa_map is None or aa_map.empty:
        return [], []

    aa_min = int(aa_map["aa_position"].min())
    aa_max = int(aa_map["aa_position"].max())
    if aa_max <= aa_min:
        return [], []

    aa_positions = np.linspace(aa_min, aa_max, num=n_ticks)
    aa_positions = sorted({int(round(a)) for a in aa_positions})

    tickvals = []
    ticktext = []

    for aa in aa_positions:
        row = aa_map[aa_map["aa_position"] == aa]
        if row.empty:
            continue

        bases = row[
            [
                "codon_base2_genomic_pos",
                "codon_base1_genomic_pos",
                "codon_base3_genomic_pos",
            ]
        ].iloc[0]

        base = None
        for b in bases:
            if not pd.isna(b):
                base = int(b)
                break

        if base is None:
            continue

        tickvals.append(base)
        ticktext.append(str(aa))

    return tickvals, ticktext


def build_domain_edge_ticks(genomic_domains_df, aa_map):
    if (
            genomic_domains_df is None
            or genomic_domains_df.empty
            or aa_map is None
            or aa_map.empty
    ):
        return [], []

    aa_positions = set()
    for _, row in genomic_domains_df.iterrows():
        aa_positions.add(int(row["aa_start"]))
        aa_positions.add(int(row["aa_end"]))

    tickvals = []
    ticktext = []

    for aa in sorted(aa_positions):
        row = aa_map[aa_map["aa_position"] == aa]
        if row.empty:
            continue
        bases = row[
            [
                "codon_base2_genomic_pos",
                "codon_base1_genomic_pos",
                "codon_base3_genomic_pos",
            ]
        ].iloc[0]
        base = None
        for b in bases:
            if not pd.isna(b):
                base = int(b)
                break
        if base is None:
            continue
        tickvals.append(base)
        ticktext.append(str(aa))

    return tickvals, ticktext


# Protein Viewer
HGVS_PROT_POS_RE = re.compile(r"p\.\(?[A-Za-z]{1,3}(\d+)", re.IGNORECASE)


def extract_protein_aa_pairs(protein_field):
    """
    Given a Protein_Label field (string or list), extract [(aa_pos, protein_label)] pairs
    by parsing the AA position directly from the HGVS protein string, e.g.:
      p.Arg830Serfs*9 -> 830
      p.Ser829*       -> 829
      p.Pro822fs      -> 822
      p.Thr821Thr     -> 821

    Returns [] if nothing can be parsed.
    """
    out = []
    if protein_field is None:
        return out

    # Normalize input into a list of raw strings
    if isinstance(protein_field, list):
        raw_items = protein_field
    else:
        raw_items = [protein_field]

    for item in raw_items:
        if item is None:
            continue
        try:
            if isinstance(item, float) and math.isnan(item):
                continue
        except Exception:
            pass

        s = str(item).strip()
        if not s:
            continue

        # Split multi-valued strings
        parts = [p.strip() for p in re.split(r"[;\|,]+", s) if p and str(p).strip()]
        if not parts:
            parts = [s]

        for p_str in parts:
            p_clean = str(p_str).strip()
            if not p_clean:
                continue

            # Ensure "p." prefix for consistent parsing (but keep original label for hover)
            p_norm = p_clean.replace(" ", "")
            if not p_norm.lower().startswith("p."):
                if p_norm.lower().startswith("p"):
                    p_norm = "p." + p_norm[1:]
                else:
                    p_norm = "p." + p_norm

            m = HGVS_PROT_POS_RE.search(p_norm)
            if not m:
                continue
            try:
                aa_pos = int(m.group(1))
            except Exception:
                continue

            out.append((aa_pos, p_clean))

    return out


DOMAIN_PALETTE = [
    "#FDE68A", "#A7F3D0", "#BFDBFE", "#FBCFE8", "#DDD6FE", "#FED7AA",
    "#C7D2FE", "#BBF7D0", "#FECACA", "#BAE6FD", "#E9D5FF", "#FDE2E2",
]


def build_domain_color_map(uniprot_df: pd.DataFrame) -> dict:
    if uniprot_df is None or uniprot_df.empty or "description" not in uniprot_df.columns:
        return {}
    descs = (
        uniprot_df["description"]
        .astype(str)
        .fillna("")
        .str.strip()
        .replace({"nan": ""})
    )
    unique_descs = [d for d in sorted(set(descs)) if d]
    return {d: DOMAIN_PALETTE[i % len(DOMAIN_PALETTE)] for i, d in enumerate(unique_descs)}


def build_pos_to_aa_offset_map(aa_map: pd.DataFrame) -> dict:
    if aa_map is None or aa_map.empty:
        return {}
    m = {}
    for _, r in aa_map.iterrows():
        try:
            aa = int(r["aa_position"])
        except Exception:
            continue
        b1 = r.get("codon_base1_genomic_pos")
        b2 = r.get("codon_base2_genomic_pos")
        b3 = r.get("codon_base3_genomic_pos")
        if b1 is not None and not pd.isna(b1):
            m[int(b1)] = (aa, -0.22)
        if b2 is not None and not pd.isna(b2):
            m[int(b2)] = (aa, 0.0)
        if b3 is not None and not pd.isna(b3):
            m[int(b3)] = (aa, 0.22)
    return m


def compute_cds_exon_aa_ranges(exon_coords, aa_map: pd.DataFrame):
    if not exon_coords or aa_map is None or aa_map.empty:
        return []
    # Ensure numeric columns for between() checks
    b1 = pd.to_numeric(aa_map["codon_base1_genomic_pos"], errors="coerce")
    b2 = pd.to_numeric(aa_map["codon_base2_genomic_pos"], errors="coerce")
    b3 = pd.to_numeric(aa_map["codon_base3_genomic_pos"], errors="coerce")

    ranges = []
    for i, (s, e) in enumerate(exon_coords, start=1):
        s0, e0 = (int(s), int(e)) if int(s) <= int(e) else (int(e), int(s))
        mask = b1.between(s0, e0) | b2.between(s0, e0) | b3.between(s0, e0)
        sub = aa_map.loc[mask, "aa_position"]
        if sub.empty:
            continue
        aa0 = int(pd.to_numeric(sub, errors="coerce").min())
        aa1 = int(pd.to_numeric(sub, errors="coerce").max())
        ranges.append(
            {"exon_idx": i, "exon_display": i, "aa_start": aa0, "aa_end": aa1, "g_start": s0, "g_end": e0}
        )
    # Sort in protein order, but keep exon_display consistent with the gene-view exon numbering
    # (exon_display == exon_idx).
    ranges = sorted(ranges, key=lambda d: d["aa_start"])
    return ranges


def fetch_uniprot_domains_and_pockets(accession=UNIPROT_ACCESSION):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    headers = {"Accept": "application/json"}

    try:
        resp = requests.get(
            url,
            headers=headers,
            params={"fields": "ft_domain,ft_region"},
            timeout=20,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"[UniProt] Failed to fetch regions/domains for {accession}: {e}")
        return pd.DataFrame(columns=["type", "description", "start", "end"])

    feats = data.get("features") or []
    if not feats:
        return pd.DataFrame(columns=["type", "description", "start", "end"])

    rows = []

    for f in feats:
        ftype = (f.get("type") or "")
        if ftype not in {"Region", "Domain"}:
            continue

        loc = f.get("location") or {}
        start = end = None
        if isinstance(loc.get("start"), dict) and isinstance(loc.get("end"), dict):
            start = loc["start"].get("value")
            end = loc["end"].get("value")
        elif isinstance(loc.get("position"), dict):
            start = loc["position"].get("value")
            end = start

        try:
            start = int(start) if start is not None else None
            end = int(end) if end is not None else None
        except (TypeError, ValueError):
            start, end = None, None

        if start is None or end is None:
            continue

        desc = f.get("description") or ftype

        if desc.strip() == "Pocket; binds T and E1A":
            continue

        desc_low = desc.lower()
        is_domain_like = ("domain" in desc_low) or (ftype == "Domain")
        is_pocket_like = ("pocket" in desc_low) or ("binds" in desc_low)
        is_spacer_like = ("spacer" in desc_low)

        if not (is_domain_like or is_pocket_like or is_spacer_like):
            continue

        rows.append(
            {
                "type": ftype,
                "description": desc,
                "start": start,
                "end": end,
            }
        )

    df_feat = pd.DataFrame(rows)
    if df_feat.empty:
        print(f"[UniProt] No Domain/Pocket/Spacer features kept for {accession}")
    else:
        print(
            f"[UniProt] Loaded {len(df_feat)} Domain/Pocket/Spacer features for {accession} "
            f"(AA {df_feat['start'].min()}–{df_feat['end'].max()})"
        )
    return df_feat


def aa_interval_to_genomic_segments(aa_start, aa_end, aa_map):
    sub = aa_map[
        (aa_map["aa_position"] >= aa_start) & (aa_map["aa_position"] <= aa_end)
        ]
    if sub.empty:
        return []

    strand = int(sub["strand"].iloc[0])
    step = 1 if strand == 1 else -1

    sub = sub.sort_values("aa_position")

    segments = []
    current_start = None
    prev = None

    for _, row in sub.iterrows():
        for col in [
            "codon_base1_genomic_pos",
            "codon_base2_genomic_pos",
            "codon_base3_genomic_pos",
        ]:
            pos = row[col]
            if pd.isna(pos):
                continue
            pos = int(pos)

            if current_start is None:
                current_start = pos
                prev = pos
            else:
                if pos == prev + step:
                    prev = pos
                else:
                    seg_start = min(current_start, prev)
                    seg_end = max(current_start, prev)
                    segments.append((seg_start, seg_end))
                    current_start = pos
                    prev = pos

    if current_start is not None:
        seg_start = min(current_start, prev)
        seg_end = max(current_start, prev)
        segments.append((seg_start, seg_end))

    return segments


def convert_domains_to_genomic(domains_df, aa_map):
    if (
            domains_df is None
            or domains_df.empty
            or aa_map is None
            or aa_map.empty
    ):
        return pd.DataFrame(
            columns=["type", "description", "g_start", "g_end", "aa_start", "aa_end"]
        )

    genomic_domains = []
    for _, dom in domains_df.iterrows():
        aa_start = int(dom["start"])
        aa_end = int(dom["end"])

        segments = aa_interval_to_genomic_segments(aa_start, aa_end, aa_map)
        if not segments:
            continue

        for g_start, g_end in segments:
            if RB1_GRCH38_START and RB1_GRCH38_END:
                g_start = max(g_start, RB1_GRCH38_START)
                g_end = min(g_end, RB1_GRCH38_END)

            genomic_domains.append(
                {
                    "type": dom["type"],
                    "description": dom["description"],
                    "g_start": int(g_start),
                    "g_end": int(g_end),
                    "aa_start": aa_start,
                    "aa_end": aa_end,
                }
            )

    result = pd.DataFrame(genomic_domains)
    if result.empty:
        print("[Domains] Domain/Pocket/Spacer features did not map onto genomic coords.")
    else:
        result = (
            result.drop_duplicates(
                subset=["type", "description", "g_start", "g_end"]
            )
            .sort_values(["g_start", "g_end", "description"])
            .reset_index(drop=True)
        )
    return result


def add_domain_traces(fig, genomic_domains_df, row, col):
    if genomic_domains_df is None or genomic_domains_df.empty:
        fig.add_annotation(
            x=(RB1_GRCH38_START + RB1_GRCH38_END) / 2.0,
            y=0.5,
            text="No protein domain/pocket info",
            showarrow=False,
            row=row,
            col=col,
        )
        fig.update_yaxes(visible=False, showticklabels=False, row=row, col=col)
        return

    seg_df = genomic_domains_df.sort_values(["g_start", "g_end"]).copy()

    def overlaps(a_start, a_end, b_start, b_end):
        return not (a_end <= b_start or b_end <= a_start)

    seg_df["has_overlap"] = False
    idx_list = list(seg_df.index)
    for i, idx_i in enumerate(idx_list):
        s1 = int(seg_df.at[idx_i, "g_start"])
        e1 = int(seg_df.at[idx_i, "g_end"])
        if s1 == e1:
            e1 = s1 + 1
        for j, idx_j in enumerate(idx_list):
            if i == j:
                continue
            s2 = int(seg_df.at[idx_j, "g_start"])
            e2 = int(seg_df.at[idx_j, "g_end"])
            if s2 == e2:
                e2 = s2 + 1
            if overlaps(s1, e1, s2, e2):
                seg_df.at[idx_i, "has_overlap"] = True
                break

    lanes = [[] for _ in range(2)]
    seg_df["lane"] = None

    for idx, row_dom in seg_df.iterrows():
        if not row_dom["has_overlap"]:
            continue
        g_start = int(row_dom["g_start"])
        g_end = int(row_dom["g_end"])
        if g_start == g_end:
            g_end = g_start + 1

        assigned_lane = None
        for lane_idx in range(2):
            lane_segments = lanes[lane_idx]
            if not any(overlaps(g_start, g_end, s, e) for (s, e) in lane_segments):
                assigned_lane = lane_idx
                lane_segments.append((g_start, g_end))
                break

        if assigned_lane is None:
            assigned_lane = 0
            lanes[0].append((g_start, g_end))
        seg_df.at[idx] = row_dom
        seg_df.at[idx, "lane"] = assigned_lane

    palette = [
        "#636EFA",
        "#EF553B",
        "#00CC96",
        "#AB63FA",
        "#FFA15A",
        "#19D3F3",
        "#FF6692",
        "#B6E880",
    ]
    color_map = {}
    color_idx = 0

    for _, seg in seg_df.iterrows():
        g_start = int(seg["g_start"])
        g_end = int(seg["g_end"])
        if g_start == g_end:
            g_end = g_start + 1

        aa_start = int(seg["aa_start"])
        aa_end = int(seg["aa_end"])
        desc = str(seg["description"])

        if desc not in color_map:
            color_map[desc] = palette[color_idx % len(palette)]
            color_idx += 1
        color = color_map[desc]

        lane = seg["lane"]
        if lane is None:
            y0, y1 = 0.0, 1.0
        elif int(lane) == 0:
            y0, y1 = 0.0, 0.5
        else:
            y0, y1 = 0.5, 1.0

        fig.add_trace(
            go.Scatter(
                x=[g_start, g_end, g_end, g_start, g_start],
                y=[y0, y0, y1, y1, y0],
                mode="lines",
                line=dict(width=0, color=color),
                fill="toself",
                hoverinfo="text",
                hoveron="fills",
                text=(f"<b>{desc}</b><br>AA: {aa_start}–{aa_end}"),
                showlegend=False,
            ),
            row=row,
            col=col,
        )

    fig.update_yaxes(
        visible=False,
        showticklabels=False,
        range=[0, 1],
        row=row,
        col=col,
    )


def classify_position(pos, exons):
    for i, (start, end) in enumerate(exons):
        if start <= pos <= end:
            return f"Exon {i + 1}"
    for i in range(len(exons) - 1):
        if exons[i][1] <= pos < exons[i + 1][0]:
            return f"Intron {i + 1}"
    return "None"


def normalize_c_label(s):
    if s is None:
        return None
    s = str(s).strip()
    if not s or s.lower() in {"nan", "none", "n/a"}:
        return None
    s = re.sub(r"^\s*c\.\s*", "", s, flags=re.I).strip()
    return f"c. {s}"


def normalize_g_label(s):
    if s is None:
        return None
    s = str(s).strip()
    if not s or s.lower() in {"nan", "none", "n/a"}:
        return None
    s = re.sub(r"^\s*g\.\s*", "", s, flags=re.I).strip()
    return f"g. {s}"


def normalize_label_for_hover(lbl):
    s = str(lbl).strip()
    if not s or s.lower() in {"nan", "none", "n/a"}:
        return None
    s = re.sub(
        r"^\s*([CGPNMRcgpnmr])\.\s*",
        lambda m: m.group(1).lower() + ".",
        s,
    )
    s = re.sub(r"^\s*g\.\s*g\.\s*", "g.", s, flags=re.I)
    s = re.sub(r"^\s*c\.\s*c\.\s*", "c.", s, flags=re.I)
    s = re.sub(r"([acgt])", lambda m: m.group(1).lower(), s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def _clean_str(x):
    if x is None:
        return ""
    if isinstance(x, float) and math.isnan(x):
        return ""
    return str(x).strip()


def build_per_variant_genomic_map_for_row(row, per_variant_map):
    pos = row.get("GRCh38 coordinates")
    try:
        pos = int(pos)
    except Exception:
        pos = None

    mut_list = row["Mutation_Label"]
    if not isinstance(mut_list, list):
        mut_list = [mut_list]

    result = []
    if pos is None or pos not in per_variant_map:
        return [None] * len(mut_list)

    variants_at_pos = per_variant_map[pos]

    for raw_mut in mut_list:
        cdna_norm = normalize_label_for_hover(raw_mut)
        chosen = None

        if cdna_norm is not None:
            cdna_norm_low = cdna_norm.lower()
            for entry in variants_at_pos:
                ev = entry.get("cdna_norm")
                if ev and ev.lower() == cdna_norm_low:
                    chosen = entry
                    break

        if chosen is None and variants_at_pos:
            # Fallback: just take the first variant at that position
            chosen = variants_at_pos[0]

        if chosen is not None:
            result.append(chosen.get("g_str"))
        else:
            result.append(None)

    return result


def compact_ref_alt_in_label(label: str, max_ref_len: int = 3) -> str:
    if label is None:
        return ""
    s = str(label)

    pattern = re.compile(r"([ACGTNacgtn]+)>([ACGTNacgtn]+)")

    def _repl(m):
        ref = m.group(1)
        alt = m.group(2)

        if len(ref) <= max_ref_len:
            ref_comp = ref
        else:
            ref_comp = ref[:max_ref_len] + "..."

        if len(alt) <= max_ref_len:
            alt_comp = alt
        else:
            alt_comp = alt[:max_ref_len] + "..."

        return f"{ref_comp}>{alt_comp}"

    return pattern.sub(_repl, s, count=1)


def format_genomic_change_for_row(row):
    dna_combined = None
    if "dna_hg38_combined" in row.index and not pd.isna(row["dna_hg38_combined"]):
        dna_combined = _clean_str(row["dna_hg38_combined"])

    chrom = CHROMOSOME
    chrom_val = None

    if "chrom" in row.index and not pd.isna(row["chrom"]):
        chrom_val = row["chrom"]

    if chrom_val is not None:
        try:
            chrom = canonicalize_chrom(chrom_val)
        except Exception:
            chrom = CHROMOSOME

    pos = None
    for col in ("hg38_genomic_start", "GRCh38 coordinates"):
        if col in row.index and not pd.isna(row[col]):
            try:
                pos = int(row[col])
                break
            except Exception:
                continue

    ref = None
    alt = None

    if "ref_38" in row.index and not pd.isna(row["ref_38"]):
        ref = row["ref_38"]
    if "alt_38" in row.index and not pd.isna(row["alt_38"]):
        alt = row["alt_38"]

    if ref is None:
        for col in ["Ref", "REF", "Reference", "Ref_GRCh38", "Reference_GRCh38"]:
            if col in row.index and not pd.isna(row[col]):
                ref = row[col]
                break

    if alt is None:
        for col in ["Alt", "ALT", "Alternate", "Alt_GRCh38", "Alternate_GRCh38", "Variant"]:
            if col in row.index and not pd.isna(row[col]):
                alt = row[col]
                break

    ref = _clean_str(ref)
    alt = _clean_str(alt)

    if not ref and not alt and dna_combined:
        return dna_combined

    if not ref:
        ref = "?"
    if not alt:
        alt = "?"

    if len(ref) <= 3:
        ref_disp = ref
    else:
        ref_disp = ref[:3] + "..."

    if len(alt) <= 3:
        alt_disp = alt
    else:
        alt_disp = alt[:3] + "..."

    change_part = f"{ref_disp}>{alt_disp}"

    if pos is None:
        return f"{chrom} {change_part}"
    return f"{chrom}:{pos} {change_part}"


def format_genomic_change_for_row_full(row):
    dna_combined = None
    if "dna_hg38_combined" in row.index and not pd.isna(row["dna_hg38_combined"]):
        dna_combined = _clean_str(row["dna_hg38_combined"])

    chrom = CHROMOSOME
    chrom_val = None

    if "chrom" in row.index and not pd.isna(row["chrom"]):
        chrom_val = row["chrom"]

    if chrom_val is not None:
        try:
            chrom = canonicalize_chrom(chrom_val)
        except Exception:
            chrom = CHROMOSOME

    pos = None
    for col in ("hg38_genomic_start", "GRCh38 coordinates"):
        if col in row.index and not pd.isna(row[col]):
            try:
                pos = int(row[col])
                break
            except Exception:
                continue

    ref = None
    alt = None

    if "ref_38" in row.index and not pd.isna(row["ref_38"]):
        ref = row["ref_38"]
    if "alt_38" in row.index and not pd.isna(row["alt_38"]):
        alt = row["alt_38"]

    if ref is None:
        for col in ["Ref", "REF", "Reference", "Ref_GRCh38", "Reference_GRCh38"]:
            if col in row.index and not pd.isna(row[col]):
                ref = row[col]
                break

    if alt is None:
        for col in ["Alt", "ALT", "Alternate", "Alt_GRCh38", "Alternate_GRCh38", "Variant"]:
            if col in row.index and not pd.isna(row[col]):
                alt = row[col]
                break

    ref = _clean_str(ref)
    alt = _clean_str(alt)

    # If we only have the combined string, just use that
    if not ref and not alt and dna_combined:
        return dna_combined

    if not ref:
        ref = "?"
    if not alt:
        alt = "?"

    change_part = f"{ref}>{alt}"

    if pos is None:
        return f"{chrom} {change_part}"
    return f"{chrom}:{pos} {change_part}"


def _clean_category_value(v):
    if v is None:
        return None
    if isinstance(v, float):
        try:
            if math.isnan(v):
                return None
        except Exception:
            pass
    s = str(v).strip()
    if not s:
        return None
    if s.lower() in {"nan", "none", "n/a"}:
        return None
    return s


def _clean_gnomad_frequency_value(v, sig_digits=12):
    """Normalize gnomAD AF values so tiny float noise doesn't create duplicate entries.

    Example problem: 6.23713590719141e-07 vs 6.23713590719142e-07 (difference ~1e-20)
    should be treated as the same value for display/counting.
    """
    s = _clean_category_value(v)
    if s is None:
        return None
    # If the value already contains extra annotation, keep it as-is.
    if "(" in s and ")" in s:
        return s
    try:
        f = float(s)
        if not math.isfinite(f):
            return None
        # Use a stable significant-digit representation.
        return format(f, f".{int(sig_digits)}g")
    except Exception:
        return s


def get_combo_counts_for_row(
        row,
        per_variant_map=None,
        use_full_genomic=False,
        return_meta=False,
):
    if isinstance(row, dict):
        row = pd.Series(row)
    elif not isinstance(row, pd.Series):
        row = pd.Series(dict(row))

    # Fallback genomic string
    if use_full_genomic:
        row_fallback = format_genomic_change_for_row_full(row)
    else:
        row_fallback = format_genomic_change_for_row(row)

    if row_fallback is None:
        row_fallback = ""
    row_fallback = str(row_fallback)

    # Per-variant genomic strings, if we have them
    per_variant_g_list = None
    if per_variant_map is not None:
        per_variant_g_list = build_per_variant_genomic_map_for_row(
            row, per_variant_map
        )

    # Core lists
    mut_list = row.get("Mutation_Label", [])
    if not isinstance(mut_list, list):
        mut_list = [mut_list]

    freq_list = row.get("Frequency", [])
    if not isinstance(freq_list, list):
        freq_list = [freq_list]

    if "Protein_Label" in row.index:
        prot_list = row["Protein_Label"]
        if not isinstance(prot_list, list):
            prot_list = [prot_list]
    else:
        prot_list = []

    if "cDNA_Label_raw" in row.index:
        cdna_list = row["cDNA_Label_raw"]
        if not isinstance(cdna_list, list):
            cdna_list = [cdna_list]
    else:
        cdna_list = []

    if "gDNA_Label_raw" in row.index:
        g_list = row["gDNA_Label_raw"]
        if not isinstance(g_list, list):
            g_list = [g_list]
    else:
        g_list = []

    # Category lists (parallel to Mutation_Label where possible)
    if "genetic_origin" in row.index:
        origin_list = row["genetic_origin"]
        if not isinstance(origin_list, list):
            origin_list = [origin_list]
    else:
        origin_list = []

    if "laterality" in row.index:
        laterality_list = row["laterality"]
        if not isinstance(laterality_list, list):
            laterality_list = [laterality_list]
    else:
        laterality_list = []

    # gnomAD frequency (gnomad_frequency) and inheritance (inheritance)
    if "gnomad_frequency" in row.index:
        raw_af = row["gnomad_frequency"]
        if isinstance(raw_af, list):
            af_list = [_clean_gnomad_frequency_value(v) for v in raw_af]
        else:
            # broadcast single value if needed
            af_list = [_clean_gnomad_frequency_value(raw_af)] * max(1, len(mut_list))
    else:
        af_list = []

    if "inheritance" in row.index:
        raw_inh = row["inheritance"]
        if isinstance(raw_inh, list):
            inh_list = [_clean_category_value(v) for v in raw_inh]
        else:
            inh_list = [_clean_category_value(raw_inh)] * max(1, len(mut_list))
    else:
        inh_list = []

    # Variant-level scores (parallel to Mutation_Label where possible)
    # CADD (cadd_score), SpliceAI (spliceai_max), AlphaMissense (alpha_missense_pathogenicity_score)
    def _broadcast_score_list(colname):
        if colname not in row.index:
            return []
        raw = row[colname]
        if isinstance(raw, list):
            return [_clean_category_value(v) for v in raw]
        # broadcast a single value if needed
        return [_clean_category_value(raw)] * max(1, len(mut_list))

    cadd_list = _broadcast_score_list("cadd_score")
    spliceai_list = _broadcast_score_list("spliceai_max")
    alphamissense_list = _broadcast_score_list("alpha_missense_pathogenicity_score")

    # Pathogenicity
    path_list = []
    if "ClinVar_Category" in row.index:
        raw_pc = row["ClinVar_Category"]
        if isinstance(raw_pc, list):
            path_list = [_clean_category_value(v) for v in raw_pc]
        else:
            path_list = [_clean_category_value(raw_pc)] * max(1, len(mut_list))

    # Molecular consequence
    mc_list = []
    if "Molecular_Consequence" in row.index:
        raw_mc = row["Molecular_Consequence"]
        if isinstance(raw_mc, list):
            mc_list = [_clean_category_value(v) for v in raw_mc]
        else:
            mc_list = [_clean_category_value(raw_mc)] * max(1, len(mut_list))

    # Source fields
    source_val = None
    if "Source_Category" in row.index:
        source_val = _clean_category_value(row["Source_Category"])
    if (not source_val) and ("Display_Source" in row.index):
        source_val = _clean_category_value(row["Display_Source"])

    source_list = None
    if "Source_Atomic_List" in row.index:
        source_list = row["Source_Atomic_List"]
        if not isinstance(source_list, list):
            source_list = [source_list]

    label_entries = []

    for idx, raw_mut in enumerate(mut_list):
        # Coding label
        coding = normalize_label_for_hover(raw_mut)

        if idx < len(cdna_list):
            raw_cdna = cdna_list[idx]
            cdna_norm = normalize_label_for_hover(raw_cdna)
            if cdna_norm is not None:
                coding = cdna_norm

        if coding is None:
            continue

        # Genomic label
        g_lbl = None
        if idx < len(g_list):
            raw_g = g_list[idx]
            if raw_g is not None:
                g_lbl = normalize_label_for_hover(raw_g)

        if (not g_lbl) and (per_variant_g_list is not None) and (idx < len(per_variant_g_list)):
            g_lbl = per_variant_g_list[idx]

        if not g_lbl:
            g_lbl = row_fallback

        if g_lbl is None:
            g_lbl = ""
        g_lbl = str(g_lbl)

        # Count
        if idx < len(freq_list):
            try:
                ac = int(freq_list[idx])
            except Exception:
                ac = 1
        else:
            ac = 1

        # Protein
        prot = None
        if idx < len(prot_list):
            raw_prot = prot_list[idx]
            if raw_prot is not None and not (
                    isinstance(raw_prot, float) and math.isnan(raw_prot)
            ):
                prot_clean = str(raw_prot).strip()
                if prot_clean and prot_clean.lower() not in {"nan", "none", "n/a"}:
                    prot = prot_clean

        key = (g_lbl, coding, prot)

        if not return_meta:
            label_entries.append((key, ac))
        else:
            origin = _clean_category_value(origin_list[idx]) if idx < len(origin_list) else None
            later = _clean_category_value(laterality_list[idx]) if idx < len(laterality_list) else None

            path_class = None
            if path_list:
                path_class = path_list[idx] if idx < len(path_list) else path_list[0]

            mc = None
            if mc_list:
                mc = mc_list[idx] if idx < len(mc_list) else mc_list[-1]

            # Per-entry source
            if source_list is not None and idx < len(source_list):
                src = _clean_category_value(source_list[idx])
            else:
                src = source_val

            af_val = af_list[idx] if idx < len(af_list) else None
            inh_val = inh_list[idx] if idx < len(inh_list) else None

            cadd_val = cadd_list[idx] if idx < len(cadd_list) else None
            spliceai_val = spliceai_list[idx] if idx < len(spliceai_list) else None
            alphamissense_val = (
                alphamissense_list[idx] if idx < len(alphamissense_list) else None
            )
            label_entries.append(
                {
                    "key": key,
                    "count": ac,
                    "origin": origin,
                    "laterality": later,
                    "path_class": path_class,
                    "molconseq": mc,
                    "source": src,
                    "af": af_val,
                    "inheritance": inh_val,
                    "cadd": cadd_val,
                    "spliceai": spliceai_val,
                    "alphamissense": alphamissense_val,
                }
            )

    if not label_entries:
        return []

    if not return_meta:
        combo_counts = {}
        for key, cnt in label_entries:
            combo_counts[key] = combo_counts.get(key, 0) + cnt

        items = sorted(
            combo_counts.items(),
            key=lambda kv: (-kv[1], kv[0][1] or "", kv[0][2] or ""),
        )
        return items

    # Aggregate meta info per (genomic, coding, protein)
    combo_meta = {}
    for entry in label_entries:
        key = entry["key"]
        cnt = entry["count"]
        origin = entry["origin"]
        later = entry["laterality"]
        path_class = entry["path_class"]
        mc = entry["molconseq"]
        src = entry["source"]
        af = entry["af"]
        inh = entry["inheritance"]

        cadd = entry.get("cadd")
        spliceai = entry.get("spliceai")
        alphamissense = entry.get("alphamissense")
        cm = combo_meta.setdefault(
            key,
            {
                "genomic": key[0],
                "coding": key[1],
                "protein": key[2],
                "total_count": 0,
                "origin_counts": {},
                "laterality_counts": {},
                "path_class_counts": {},
                "molconseq_counts": {},
                "source_counts": {},
                # NEW:
                "cadd_counts": {},
                "spliceai_counts": {},
                "alphamissense_counts": {},
                "af_counts": {},
                "inheritance_counts": {},
            },
        )
        cm["total_count"] += cnt

        if origin:
            cm["origin_counts"][origin] = cm["origin_counts"].get(origin, 0) + cnt
        if later:
            cm["laterality_counts"][later] = cm["laterality_counts"].get(later, 0) + cnt
        if path_class:
            cm["path_class_counts"][path_class] = cm["path_class_counts"].get(path_class, 0) + cnt
        if mc:
            cm["molconseq_counts"][mc] = cm["molconseq_counts"].get(mc, 0) + cnt
        if src:
            cm["source_counts"][src] = cm["source_counts"].get(src, 0) + cnt
        if cadd:
            cm["cadd_counts"][cadd] = cm["cadd_counts"].get(cadd, 0) + 1
        if spliceai:
            cm["spliceai_counts"][spliceai] = cm["spliceai_counts"].get(spliceai, 0) + 1
        if alphamissense:
            cm["alphamissense_counts"][alphamissense] = cm["alphamissense_counts"].get(alphamissense, 0) + 1
        if af:
            cm["af_counts"][af] = cm["af_counts"].get(af, 0) + 1
        if inh:
            cm["inheritance_counts"][inh] = cm["inheritance_counts"].get(inh, 0) + cnt

    items = list(combo_meta.values())
    items.sort(
        key=lambda d: (-d["total_count"], d["coding"] or "", d["protein"] or "")
    )
    return items


def build_hover_text(row, per_variant_map=None):
    row_fallback = format_genomic_change_for_row(row)
    if row_fallback is None:
        row_fallback = ""
    row_fallback = str(row_fallback)
    max_genomic_len = 40
    if len(row_fallback) > max_genomic_len:
        row_fallback = row_fallback[: max_genomic_len - 3] + "..."

    per_variant_g_list = None
    if per_variant_map is not None:
        per_variant_g_list = build_per_variant_genomic_map_for_row(row, per_variant_map)

    lines = [
        f"<b>Region:</b> {row['Region_Label']}",
    ]

    source_info = None
    if "Source_List" in row.index:
        src_list = row["Source_List"]
        if isinstance(src_list, list) and len(src_list) > 0:
            unique_src = sorted({str(s) for s in src_list})
            source_info = ", ".join(unique_src)
    elif "Source_Category" in row.index:
        cat = row["Source_Category"]
        if isinstance(cat, str) and cat:
            source_info = cat

    if source_info:
        lines.append(f"<b>Source:</b> {source_info}")

    # Extra context for transcript/protein views (safe if maps are not initialized yet)
    gpos_i = None
    try:
        for _k in ("GRCh38 coordinates", "hg38_genomic_start", "Genomic_position", "Genomic Position"):
            if _k in row.index and pd.notnull(row.get(_k)):
                gpos_i = int(row.get(_k))
                break
    except Exception:
        gpos_i = None

    if gpos_i is not None:
        lines.append(f"<b>Genomic position (GRCh38):</b> {gpos_i:,}")

        _pos_to_aa = globals().get("POS_TO_AA_OFFSET", {}) or {}
        _mapped = _pos_to_aa.get(gpos_i)
        aa_pos_i = None
        if _mapped:
            try:
                aa_pos_i = int(_mapped[0])
            except Exception:
                aa_pos_i = None

        if aa_pos_i is not None:
            # lines.append(f"<b>Protein position (AA):</b> {aa_pos_i}")

            # CDS exon (protein-order exon numbering)
            _cds_exons = globals().get("CDS_EXON_AA_RANGES", []) or []
            exon_disp = None
            for _ex in _cds_exons:
                try:
                    if int(_ex.get("aa_start", -1)) <= aa_pos_i <= int(_ex.get("aa_end", -1)):
                        exon_disp = _ex.get("exon_display", None)
                        break
                except Exception:
                    continue
            # if exon_disp is not None:
            # lines.append(f"<b>CDS exon:</b> {exon_disp}")

            # Protein domains (UniProt features) at this AA position
            _udf = globals().get("uniprot_dom_df", None)
            if _udf is not None and hasattr(_udf, "empty") and not _udf.empty:
                try:
                    _sub = _udf[(_udf["start"] <= aa_pos_i) & (_udf["end"] >= aa_pos_i)]
                    _doms = sorted(set(_sub.get("description", pd.Series([], dtype=str)).astype(str).tolist()))
                except Exception:
                    _doms = []
                if _doms:
                    show = _doms[:4]
                    dom_txt = ", ".join(show)
                    if len(_doms) > 4:
                        dom_txt += f", … (+{len(_doms) - 4} more)"
                    # lines.append(f"<b>Protein domain:</b> {dom_txt}")

    # Core lists
    mut_list = row["Mutation_Label"]
    if not isinstance(mut_list, list):
        mut_list = [mut_list]

    freq_list = row["Frequency"] if isinstance(row["Frequency"], list) else [row["Frequency"]]

    prot_list = row["Protein_Label"] if ("Protein_Label" in row.index) else None
    if not isinstance(prot_list, list):
        prot_list = []

    cdna_list = row["cDNA_Label_raw"] if ("cDNA_Label_raw" in row.index) else None
    if not isinstance(cdna_list, list):
        cdna_list = []

    g_list = row["gDNA_Label_raw"] if ("gDNA_Label_raw" in row.index) else None
    if not isinstance(g_list, list):
        g_list = []

    label_quads = []

    for idx, raw_mut in enumerate(mut_list):

        coding = normalize_label_for_hover(raw_mut)

        if idx < len(cdna_list):
            raw_cdna = cdna_list[idx]
            cdna_norm = normalize_label_for_hover(raw_cdna)
            if cdna_norm is not None:
                coding = cdna_norm

        if coding is None:
            continue

        g_lbl = None

        if idx < len(g_list):
            raw_g = g_list[idx]
            g_lbl = normalize_label_for_hover(raw_g)

        if (not g_lbl) and (per_variant_g_list is not None) and (idx < len(per_variant_g_list)):
            g_lbl = per_variant_g_list[idx]

        if not g_lbl:
            g_lbl = row_fallback

        if g_lbl is None:
            g_lbl = ""
        g_lbl = str(g_lbl)

        g_lbl = compact_ref_alt_in_label(g_lbl)

        if len(g_lbl) > max_genomic_len:
            g_lbl = g_lbl[: max_genomic_len - 3] + "..."

        if idx < len(freq_list):
            try:
                ac = int(freq_list[idx])
            except Exception:
                ac = 1
        else:
            ac = 1

        prot = None
        if idx < len(prot_list):
            raw_prot = prot_list[idx]
            if raw_prot is not None and not (
                    isinstance(raw_prot, float) and math.isnan(raw_prot)
            ):
                prot_clean = str(raw_prot).strip()
                if prot_clean and prot_clean.lower() not in {"nan", "none", "n/a"}:
                    prot = prot_clean

        label_quads.append((g_lbl, coding, prot, ac))

    if not label_quads:
        lines.append("No labeled variants at this site")
        return "<br>".join(lines)

    combo_counts = {}
    for g_lbl, cdna, prot, cnt in label_quads:
        key = (g_lbl, cdna, prot)
        combo_counts[key] = combo_counts.get(key, 0) + cnt

    items = sorted(
        combo_counts.items(),
        key=lambda kv: (-kv[1], kv[0][1], kv[0][2] or "")
    )

    max_lines = 55
    display_items = items[:max_lines]
    hidden = max(0, len(items) - len(display_items))

    col0_w = 28  # Genomic
    col1_w = 26  # cDNA/coding
    col2_w = 18  # Protein

    header = (
        f"{'Genomic':<{col0_w}} "
        f"{'coding':<{col1_w}} "
        f"{'protein':<{col2_w}} "
        f"Count"
    )

    sep_char = "─"
    sep = sep_char * len(header)

    # Decide column widths from both headers and data
    col0_w = max(len("Genomic"), max(len(g) for (g, _, _), _ in display_items))
    col1_w = max(len("coding"), max(len(c) for (_, c, _), _ in display_items))
    col2_w = max(len("protein"), max(len((p or "-")) for (_, _, p), _ in display_items))
    col3_w = max(len("Count"), max(len(str(cnt)) for _, cnt in display_items))

    # Header row (names right above the line)
    header_text = (
        f"{'Genomic':<{col0_w}}  "
        f"{'Coding':<{col1_w}}  "
        f"{'Protein':<{col2_w}}  "
        f"{'Count':>{col3_w}}"
    )

    sep_char = "─"
    sep_text = sep_char * len(header_text)

    # Data rows
    table_lines = [header_text, sep_text]
    for (g_col, cdna, prot), cnt in display_items:
        prot_disp = prot or "-"
        row = (
            f"{g_col:<{col0_w}}  "
            f"{cdna:<{col1_w}}  "
            f"{prot_disp:<{col2_w}}  "
            f"{cnt:>{col3_w}}"
        )
        table_lines.append(row)

    table_block = (
            "<span style='font-family:monospace; font-size:11px; "
            "white-space:pre; line-height:1.1;'>"
            + "<br>".join(html_lib.escape(line) for line in table_lines)
            + "</span>"
    )

    lines.append(table_block)

    if hidden > 0:
        lines.append(
            "<span style='font-size:11px;'>… and "
            f"{hidden} more</span>"
        )

    return "<br>".join(lines)


def build_hover_text_subset(row, keep_indices, per_variant_map=None):
    """Build hover text for only the selected per-variant indices.

    The curated dataset is aggregated per genomic coordinate, with columns like
    Mutation_Label / Frequency / cDNA_Label_raw / Protein_Label stored as lists.
    For the methylation/parent-of-origin panels, we often need the hover to
    reflect only the variants that belong to that specific lollipop (e.g.
    methylated vs. unmethylated), rather than all variants at the coordinate.
    """

    if keep_indices is None:
        return build_hover_text(row, per_variant_map=per_variant_map)

    try:
        keep_set = set(int(i) for i in keep_indices if i is not None)
    except Exception:
        keep_set = set()

    # If we cannot interpret indices, fall back to full hover.
    if not keep_set:
        return build_hover_text(row, per_variant_map=per_variant_map)

    row2 = row.copy()

    def _subset_list(val):
        if isinstance(val, list):
            out = []
            for i in sorted(keep_set):
                if 0 <= i < len(val):
                    out.append(val[i])
            return out
        # Scalar value: keep only if index 0 is requested.
        return [val] if 0 in keep_set else []

    for col in [
        "Mutation_Label",
        "Frequency",
        "Protein_Label",
        "cDNA_Label_raw",
        "gDNA_Label_raw",
    ]:
        if col in row2.index:
            row2[col] = _subset_list(row2[col])

    return build_hover_text(row2, per_variant_map=per_variant_map)


def build_exon_hover_map(exons, genomic_domains_df):
    hover_map = {}
    if not exons:
        return hover_map

    for idx, (start, end) in enumerate(exons):
        exon_num = idx + 1

        exon_label_html = (
            f"<span style='color:{EXON_COLOR}; font-weight:bold;'>"
            f"Exon {exon_num}</span>"
        )

        lines = [
            exon_label_html,
            f"<b>Genomic:</b> {start}–{end}",
        ]

        if genomic_domains_df is not None and not genomic_domains_df.empty:
            overlap = genomic_domains_df[
                (genomic_domains_df["g_end"] >= start)
                & (genomic_domains_df["g_start"] <= end)
                ]
            if not overlap.empty:
                lines.append("<b>Protein domain:</b>")
                grouped = (
                    overlap.groupby("description")[["aa_start", "aa_end"]]
                    .agg({"aa_start": "min", "aa_end": "max"})
                    .reset_index()
                )
                for _, row in grouped.iterrows():
                    desc = str(row["description"])
                    aa_s = int(row["aa_start"])
                    aa_e = int(row["aa_end"])
                    lines.append(f"{desc} (AA {aa_s}–{aa_e})")
            else:
                lines.append("No overlapping protein domain")

        hover_map[exon_num] = "<br>".join(lines)

    return hover_map


def canonicalize_chrom(chrom_str):
    """
    Normalize chromosome identifiers to 'chrN' / 'chrX' / 'chrY' / 'chrM'.

    Handles common messy inputs like 13.0 (float), '13.0' (string), 'chr13', 'ChrX', etc.
    """
    if chrom_str is None:
        return None

    # Handle numeric chromosome columns (often read as floats like 13.0)
    try:
        if isinstance(chrom_str, (int, np.integer)):
            s = str(int(chrom_str))
        elif isinstance(chrom_str, float):
            if math.isnan(chrom_str):
                return None
            s = str(int(chrom_str)) if float(chrom_str).is_integer() else str(chrom_str)
        else:
            s = str(chrom_str).strip()
    except Exception:
        s = str(chrom_str).strip()

    s = str(s).strip()
    if not s or s.lower() in {"nan", "none"}:
        return None

    s = re.sub(r"^chr", "", s, flags=re.I).strip()

    # Remove trailing .0 / .00 etc (e.g., '13.0' -> '13')
    s = re.sub(r"\.0+$", "", s)

    # Map numeric encodings sometimes used for sex/mitochondrial chromosomes
    if s == "23":
        s = "X"
    elif s == "24":
        s = "Y"
    elif s == "25":
        s = "M"

    s = s.upper()

    # Harmonize mitochondrial naming
    if s == "MT":
        s = "M"

    return "chr" + s


def parse_variant_string(val):
    s = str(val).strip()
    if not s or s.lower() in {"nan", "none"}:
        return None, None

    m = re.search(r"(?i)\b(?P<chr>(?:chr)?[0-9xyzm]+):(?P<pos>\d+)", s)
    if m:
        chrom = canonicalize_chrom(m.group("chr"))
        try:
            pos = int(m.group("pos"))
        except ValueError:
            pos = None
        if pos is not None:
            return chrom, pos

    m = re.search(r"(?i)\b(?P<chr>(?:chr)?[0-9xyzm]+)-(?P<pos>\d+)", s)
    if m:
        chrom = canonicalize_chrom(m.group("chr"))
        try:
            pos = int(m.group("pos"))
        except ValueError:
            pos = None
        if pos is not None:
            return chrom, pos

    m = re.search(r"(chr[^:]*):(\d+)", s)
    if m:
        chrom = canonicalize_chrom(m.group(1))
        try:
            pos = int(m.group(2))
        except ValueError:
            pos = None
        if pos is not None:
            return chrom, pos

    return None, None


def parse_user_variant_upload(contents, filename=None, has_header=False, delimiter=","):
    empty_cols = [
        "GRCh38 coordinates",
        "Region_Label",
        "Total_AC",
        # Keep raw per-report labels so methylation-panel hovers can be filtered
        # to show only methylated vs. unmethylated reports at a coordinate.
        "Labels",
        "HaploList",
        "Hyper_Count",
        "Hypo_Count",
        "Hover_Text",
        "Jittered_Coord",
    ]
    if contents is None:
        return pd.DataFrame(columns=empty_cols)

    try:
        content_type, content_string = contents.split(",", 1)
    except ValueError:
        return pd.DataFrame(columns=empty_cols)

    try:
        data = base64.b64decode(content_string)
    except Exception:
        return pd.DataFrame(columns=empty_cols)

    is_excel = False
    if filename and filename.lower().endswith((".xls", ".xlsx")):
        is_excel = True
    elif data[:2] == b"PK":
        is_excel = True

    try:
        if is_excel:
            if has_header:
                user_df_raw = pd.read_excel(io.BytesIO(data))
            else:
                user_df_raw = pd.read_excel(io.BytesIO(data), header=None)
        else:
            try:
                text = data.decode("utf-8")
            except Exception:
                text = data.decode("latin-1")

            if delimiter == "tab":
                sep = "\t"
            elif delimiter == ",":
                sep = ","
            else:
                sep = None

            header_param = 0 if has_header else None
            user_df_raw = pd.read_csv(
                io.StringIO(text),
                header=header_param,
                sep=sep,
            )
    except Exception:
        return pd.DataFrame(columns=empty_cols)

    if user_df_raw.empty:
        return pd.DataFrame(columns=empty_cols)

    col_map = {str(c).strip().lower(): c for c in user_df_raw.columns}

    haplo_col = None

    for lower, orig in col_map.items():
        if "haplo" in lower or "methyl" in lower:
            haplo_col = orig
            break

    if haplo_col is None:

        for col in user_df_raw.columns:
            try:
                series = user_df_raw[col].astype(str).str.lower()
            except Exception:
                continue

            if (
                    series.str.contains("hyper", regex=False).any()
                    or series.str.contains("hypo", regex=False).any()
                    or series.str.contains("methyl", regex=False).any()
                    or series.str.contains("paternal", regex=False).any()
                    or series.str.contains("maternal", regex=False).any()
            ):
                haplo_col = col
                break

    chr_col = None
    pos_col = None
    ref_col = None
    alt_col = None
    for lower, orig in col_map.items():
        if lower in {"chr", "chrom", "chromosome"}:
            chr_col = orig
        elif lower in {"hg38_genomic_start", "pos", "position", "start"}:
            pos_col = orig
        elif lower in {"hg38_ref", "ref", "ref_allele"}:
            ref_col = orig
        elif lower in {"hg38_alt", "alt", "alt_allele"}:
            alt_col = orig

    def classify_haplotype(value):

        if value is None:
            return None
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass

        s = str(value).strip()
        if not s:
            return None

        s_low = s.lower()

        s_words = re.sub(r"[^a-z]+", " ", s_low)
        tokens = s_words.split()

        # Hypo / unmethylated / paternal : "hypo"
        if any(t.startswith("hypo") for t in tokens):
            return "hypo"
        if "unmethylated" in tokens or "unmeth" in tokens:
            return "hypo"
        if "paternal" in tokens:
            return "hypo"

        # Hyper / maternal : "hyper"
        if any(t.startswith("hyper") for t in tokens):
            return "hyper"
        if "maternal" in tokens:
            return "hyper"

        # Generic methylated/methylation : hyper
        if "methylated" in tokens or "methylation" in tokens or "methyl" in tokens:
            return "hyper"

        if "hypo" in s_low or "unmethyl" in s_low:
            return "hypo"
        if "hyper" in s_low:
            return "hyper"
        if "paternal" in s_low:
            return "hypo"
        if "maternal" in s_low:
            return "hyper"

        return None

    def canonicalize_chrom_local(chrom_str):
        if chrom_str is None:
            return None
        s = str(chrom_str).strip()
        if not s:
            return None
        s = re.sub(r"^chr", "", s, flags=re.I)
        s = s.upper()
        return "chr" + s

    def parse_variant_string_local(val):

        s = str(val).strip()
        if not s or s.lower() in {"nan", "none"}:
            return None, None

        m = re.search(r"(?i)\b(?P<chr>(?:chr)?[0-9xyzm]+):(?P<pos>\d+)", s)
        if m:
            chrom = canonicalize_chrom_local(m.group("chr"))
            try:
                pos = int(m.group("pos"))
            except ValueError:
                pos = None
            if pos is not None:
                return chrom, pos

        m = re.search(r"(?i)\b(?P<chr>(?:chr)?[0-9xyzm]+)-(?P<pos>\d+)", s)
        if m:
            chrom = canonicalize_chrom_local(m.group("chr"))
            try:
                pos = int(m.group("pos"))

            except ValueError:
                pos = None
            if pos is not None:
                return chrom, pos

        m = re.search(r"(chr[^:]*):(\d+)", s)
        if m:
            chrom = canonicalize_chrom_local(m.group(1))
            try:
                pos = int(m.group(2))
            except ValueError:
                pos = None
            if pos is not None:
                return chrom, pos

        return None, None

    variants = []

    if chr_col is not None and pos_col is not None:
        for _, row in user_df_raw.iterrows():
            chrom_val = row[chr_col]
            pos_val = row[pos_col]
            if pd.isna(chrom_val) or pd.isna(pos_val):
                continue
            try:
                pos = int(pos_val)
            except (TypeError, ValueError):
                continue

            canonical_chr = canonicalize_chrom_local(chrom_val)
            if canonical_chr != CHROMOSOME:
                continue

            if ref_col is not None and alt_col is not None:
                ref_val = row[ref_col]
                alt_val = row[alt_col]
                if not pd.isna(ref_val) and not pd.isna(alt_val):
                    label = (
                        f"{canonical_chr}:{pos} "
                        f"{str(ref_val).strip()}>{str(alt_val).strip()}"
                    )
                else:
                    label = f"{canonical_chr}:{pos}"
            else:
                label = f"{canonical_chr}:{pos}"

            haplo_type = None
            if haplo_col is not None:
                haplo_type = classify_haplotype(row[haplo_col])

            variants.append(
                {"label": label, "pos": pos, "haplo_type": haplo_type}
            )

    if not variants:

        for _, row in user_df_raw.iterrows():
            haplo_type = None
            if haplo_col is not None:
                haplo_type = classify_haplotype(row[haplo_col])

            row_vals = list(row.values)
            found = False
            for val in row_vals:
                if pd.isna(val):
                    continue
                chrom, pos = parse_variant_string_local(val)
                if chrom is None or pos is None:
                    continue
                if chrom != CHROMOSOME:
                    continue

                label = str(val).strip()
                variants.append(
                    {"label": label, "pos": pos, "haplo_type": haplo_type}
                )
                found = True
                break

            if not found:
                continue

    if not variants:
        return pd.DataFrame(columns=empty_cols)

    tmp = pd.DataFrame(variants)

    # Group by position
    grouped = (
        tmp.groupby("pos")
        .agg(Labels=("label", list), HaploList=("haplo_type", lambda x: list(x)))
        .reset_index()
        .rename(columns={"pos": "GRCh38 coordinates"})
    )

    grouped["Region_Label"] = grouped["GRCh38 coordinates"].apply(
        lambda x: classify_position(int(x), rb1_exons)
    )

    grouped["Total_AC"] = grouped["Labels"].apply(len)

    grouped["Hyper_Count"] = grouped["HaploList"].apply(
        lambda lst: sum(1 for h in lst if h == "hyper")
    )
    grouped["Hypo_Count"] = grouped["HaploList"].apply(
        lambda lst: sum(1 for h in lst if h == "hypo")
    )

    def build_user_hover(row):
        lines = [
            "<b>User variant(s)</b>",
            f"<b>Genomic Position:</b> {int(row['GRCh38 coordinates'])}",
            f"<b>Region:</b> {row['Region_Label']}",
        ]

        num_variants = int(row.get("Total_AC", 0) or 0)
        lines.append(f"<b>Report Count:</b> {num_variants}")

        labels = row.get("Labels")
        if isinstance(labels, list) and labels:
            label_counts = {}
            for lbl in labels:
                if lbl is None or (isinstance(lbl, float) and pd.isna(lbl)):
                    continue
                s = str(lbl).strip()
                if not s:
                    continue
                label_counts[s] = label_counts.get(s, 0) + 1

            if label_counts:
                lines.append("<b>Variant label(s):</b>")

                for lbl, cnt in sorted(
                        label_counts.items(),
                        key=lambda kv: (-kv[1], kv[0]),
                ):
                    lines.append(f"{lbl}, Report Count: {cnt}")

        hyper = int(row.get("Hyper_Count", 0) or 0)
        hypo = int(row.get("Hypo_Count", 0) or 0)
        if hyper or hypo:
            parts = []
            if hyper:
                parts.append(f"Hypermethylated: {hyper}")
            if hypo:
                parts.append(f"Hypomethylated: {hypo}")
            lines.append("<b>Methylation:</b> " + "; ".join(parts))

        return "<br>".join(lines)

    grouped["Hover_Text"] = grouped.apply(build_user_hover, axis=1)

    grouped["Jittered_Coord"] = grouped["GRCh38 coordinates"].astype(float)

    return grouped[empty_cols]


def search_in_curated_records(query, records):
    q = (query or "").strip()
    if not q or not records:
        return []

    results = []

    # 1) Genomic coordinate
    chrom, pos = parse_variant_string(q)
    if pos is None:
        if re.fullmatch(r"\d+", q):
            pos = int(q)
            chrom = CHROMOSOME

    if chrom == CHROMOSOME and pos is not None:
        for rec in records:
            try:
                coord_val = int(rec.get("GRCh38 coordinates"))
            except (TypeError, ValueError):
                continue
            if coord_val == pos and rec not in results:
                results.append(rec)

    # 2) cDNA
    norm_q = normalize_label_for_hover(q)
    if norm_q:
        nq = norm_q.lower()
        for rec in records:
            labels = rec.get("Mutation_Label")
            if not isinstance(labels, list):
                continue
            for lbl in labels:
                lbl_norm = normalize_label_for_hover(lbl)
                if not lbl_norm:
                    continue
                if lbl_norm.lower() == nq and rec not in results:
                    results.append(rec)
                    break

    # 3) Protein
    q_prot = q.strip()
    if q_prot:
        q_norm = q_prot.lower().replace(" ", "")
        if not q_norm.startswith("p."):
            if q_norm.startswith("p"):
                q_norm = "p." + q_norm[1:]
            else:
                q_norm = "p." + q_norm

        for rec in records:
            prot_list = rec.get("Protein_Label")
            if not isinstance(prot_list, list):
                continue
            for prot in prot_list:
                if (
                        prot is None
                        or (isinstance(prot, float) and math.isnan(prot))
                ):
                    continue
                p_norm = str(prot).strip().lower().replace(" ", "")
                if not p_norm.startswith("p."):
                    if p_norm.startswith("p"):
                        p_norm = "p." + p_norm[1:]
                    else:
                        p_norm = "p." + p_norm
                if p_norm == q_norm and rec not in results:
                    results.append(rec)
                    break

    return results


rb1_exons = get_rb1_exons_from_ensembl(ENSEMBL_CANONICAL)
aa_map_df = load_or_build_aa_map()

init_start, init_end = get_initiator_codon_bounds(aa_map_df)
if init_start is not None:
    PROMOTER_END = init_start
else:
    PROMOTER_END = PROMOTER_START

uniprot_dom_df = fetch_uniprot_domains_and_pockets()

n_term_df = pd.DataFrame(
    [
        {"type": "Domain", "description": "N-terminal domain", "start": 1, "end": 373}
    ]
)
c_term_df = pd.DataFrame(
    [
        {"type": "Domain", "description": "C-terminal domain", "start": 768, "end": 927}
    ]
)

if uniprot_dom_df is None or uniprot_dom_df.empty:
    uniprot_dom_df = pd.concat([n_term_df, c_term_df], ignore_index=True)
else:
    uniprot_dom_df = pd.concat(
        [uniprot_dom_df, n_term_df, c_term_df], ignore_index=True
    )

genomic_domains_df = convert_domains_to_genomic(uniprot_dom_df, aa_map_df)

# Precompute mappings for the protein / mRNA transcript track (used in the main figure callback)
DOMAIN_COLOR_MAP = build_domain_color_map(uniprot_dom_df)
POS_TO_AA_OFFSET = build_pos_to_aa_offset_map(aa_map_df)
CDS_EXON_AA_RANGES = compute_cds_exon_aa_ranges(rb1_exons, aa_map_df)

EXON_HOVER_MAP = build_exon_hover_map(rb1_exons, genomic_domains_df)

# --- Fast domain lookup (avoid Pandas boolean indexing inside callbacks) ---
DOM_INTERVALS = []
if isinstance(uniprot_dom_df, pd.DataFrame) and (not uniprot_dom_df.empty):
    try:
        DOM_INTERVALS = [
            (int(s), int(e), str(d))
            for s, e, d in zip(
                uniprot_dom_df.get('start', []),
                uniprot_dom_df.get('end', []),
                uniprot_dom_df.get('description', []),
            )
            if pd.notna(s) and pd.notna(e)
        ]
    except Exception:
        DOM_INTERVALS = []

def _domains_for_aa(aa_pos):
    """Return list of UniProt domain labels covering a given AA position."""
    try:
        aa = int(aa_pos)
    except Exception:
        return []
    # DOM_INTERVALS is tiny (usually a handful), so a simple loop is fastest
    return [d for (s, e, d) in DOM_INTERVALS if s <= aa <= e]

def _domains_for_range(aa0, aa1):
    """Return sorted unique UniProt domain labels overlapping [aa0, aa1]."""
    try:
        a0 = int(aa0)
        a1 = int(aa1)
    except Exception:
        return []
    lo, hi = (a0, a1) if a0 <= a1 else (a1, a0)
    return sorted({d for (s, e, d) in DOM_INTERVALS if s <= hi and e >= lo})

# --- Cache expensive per-row meta computation used in the main callback ---
ROW_META_CACHE = {}

def _needs_combo_parse(rec: dict) -> bool:
    """Heuristic: only run expensive combo parsing if row likely encodes multiple variants."""
    try:
        for k in ('Mutation_Label', 'cDNA_Label_raw', 'Protein_Label', 'gDNA_Label_raw'):
            v = rec.get(k, '')
            if isinstance(v, (list, tuple)):
                return True
            s = str(v)
            if ';' in s or '|' in s:
                return True
        return False
    except Exception:
        return True

def _get_variant_meta_cached(rec: dict):
    """Memoized wrapper around get_combo_counts_for_row (fast for single-variant rows)."""
    rid = rec.get('_row_id', None)
    if rid is not None and rid in ROW_META_CACHE:
        return ROW_META_CACHE[rid]
    meta = []
    if _needs_combo_parse(rec):
        try:
            meta = get_combo_counts_for_row(
                rec,
                per_variant_map=per_variant_map if isinstance(per_variant_map, dict) else None,
                use_full_genomic=True,
                return_meta=True,
            )
        except Exception:
            meta = []
    if rid is not None:
        ROW_META_CACHE[rid] = meta
    return meta


df = pd.read_excel(str(EXCEL_FILE), sheet_name=INPUT_SHEET)
df.columns = df.columns.str.strip()

SOURCE_COLUMN = None
for cand in ["source", "Source", "SOURCE"]:
    if cand in df.columns:
        SOURCE_COLUMN = cand
        break

if "hg38_genomic_start" not in df.columns and "hg38" in df.columns:
    tmp = df["hg38"].astype(str).str.extract(r":(\d+)")[0]
    df["hg38_genomic_start"] = pd.to_numeric(tmp, errors="coerce")

if "hg38_genomic_start" in df.columns:
    _lo, _hi = 48_300_000, 48_600_000
    oob = df[
        (df["hg38_genomic_start"] < _lo) | (df["hg38_genomic_start"] > _hi)
        ]
    if not oob.empty:
        print(
            f"[RB1] Dropping {oob.shape[0]} out-of-bounds variants outside "
            f"{_lo}-{_hi}"
        )
    df = df[
        (df["hg38_genomic_start"] >= _lo) & (df["hg38_genomic_start"] <= _hi)
        ].copy()

df = df.dropna(subset=["hg38_genomic_start"]).copy()
df["GRCh38 coordinates"] = df["hg38_genomic_start"].astype(int)

# Unique per-row identifier (used for memoization in callbacks)
df['_row_id'] = np.arange(df.shape[0], dtype=int)


if "# Reported" in df.columns:
    df["Frequency"] = (
        pd.to_numeric(df["# Reported"], errors="coerce").fillna(1).astype(int)
    )
elif "report_count" in df.columns:
    df["Frequency"] = (
        pd.to_numeric(df["report_count"], errors="coerce").fillna(1).astype(int)
    )
else:
    df["Frequency"] = 1

if "clinvar_variant_classification" in df.columns:
    base_cls = df["clinvar_variant_classification"]
elif "clinical_classification_combined" in df.columns:
    base_cls = df["clinical_classification_combined"]
else:
    base_cls = None

if base_cls is not None:
    df["ClinVar_Category"] = (
        base_cls.astype(str)
        .str.strip()
        .replace({"nan": "Unspecified", "": "Unspecified"})
    )
else:
    df["ClinVar_Category"] = "Unspecified"

if "lovd_classification" in df.columns:
    lovd_base = df["lovd_classification"]
    df["LOVD_Category"] = (
        lovd_base.astype(str)
        .str.strip()
        .replace({"nan": "Unspecified", "": "Unspecified"})
    )
else:
    df["LOVD_Category"] = "Unspecified"

GENETIC_ORIGIN_FILTER_OPTIONS = []
if "genetic_origin" in df.columns:
    go_vals = (
        df["genetic_origin"]
        .astype(str)
        .str.strip()
        .replace({"nan": "", "": np.nan})
        .dropna()
        .unique()
    )
    go_vals = sorted(go_vals)
    GENETIC_ORIGIN_FILTER_OPTIONS = [{"label": v, "value": v} for v in go_vals]
GENETIC_ORIGIN_FILTER_OPTIONS = with_all_option(GENETIC_ORIGIN_FILTER_OPTIONS)

MCONSEQ_FILTER_OPTIONS = []
if "molecular_consequence" in df.columns:
    mc_vals = (
        df["molecular_consequence"]
        .astype(str)
        .str.strip()
        .replace({"nan": "", "": np.nan})
        .dropna()
        .unique()
    )
    mc_vals = sorted(mc_vals, key=mc_sort_key)
    MCONSEQ_FILTER_OPTIONS = [
        {"label": mc_pretty_label(v), "value": v} for v in mc_vals
    ]
MCONSEQ_FILTER_OPTIONS = with_all_option(MCONSEQ_FILTER_OPTIONS)

LATERALITY_FILTER_OPTIONS = []
if "laterality" in df.columns:
    lat_vals = (
        df["laterality"]
        .astype(str)
        .str.strip()
        .replace({"nan": "", "": np.nan})
        .dropna()
        .unique()
    )
    lat_vals = sorted(lat_vals)
    LATERALITY_FILTER_OPTIONS = [{"label": v, "value": v} for v in lat_vals]
LATERALITY_FILTER_OPTIONS = with_all_option(LATERALITY_FILTER_OPTIONS)

LOVD_CLASS_FILTER_OPTIONS = []
if "LOVD_Category" in df.columns:
    lovd_vals = (
        df["LOVD_Category"]
        .astype(str)
        .str.strip()
        .replace({"nan": "", "": np.nan})
        .dropna()
        .unique()
    )
    lovd_vals = list(lovd_vals)
    if "Unspecified" in lovd_vals:
        core = sorted([v for v in lovd_vals if v != "Unspecified"])
        lovd_vals = core + ["Unspecified"]
    else:
        lovd_vals = sorted(lovd_vals)
    LOVD_CLASS_FILTER_OPTIONS = [{"label": v, "value": v} for v in lovd_vals]

CLINVAR_CLASS_FILTER_OPTIONS = []
category_order = [
    "Pathogenic",
    "Pathogenic/Likely pathogenic",
    "Likely pathogenic",
    "Likely benign",
    "Benign/Likely benign",
    "Benign",
    "Conflicting classifications of pathogenicity",
    "VUS",
    "Uncertain significance",
]
clinvar_vals = (
    df["ClinVar_Category"]
    .astype(str)
    .str.strip()
    .replace({"nan": "", "": np.nan})
    .dropna()
    .unique()
)
clinvar_set = set(clinvar_vals)
cats_present = [c for c in category_order if c in clinvar_set]
other_cats = sorted(list(clinvar_set - set(cats_present)))
clinvar_ordered = cats_present + other_cats
CLINVAR_CLASS_FILTER_OPTIONS = [
    {"label": c, "value": c} for c in clinvar_ordered
]
CLINVAR_CLASS_FILTER_OPTIONS = with_all_option(
    CLINVAR_CLASS_FILTER_OPTIONS
)

c_src = None
g_src = None

# Case-insensitive column picking (some sheets use "Protein"/"cDNA"/etc.)
_lower_to_col = {str(c).strip().lower(): c for c in df.columns}


def _pick_col(candidates):
    for key in candidates:
        key_l = str(key).strip().lower()
        if key_l in _lower_to_col:
            return _lower_to_col[key_l]
    return None


_cdna_col = _pick_col(["cdna", "hgvs_c", "cdna_hgvs", "c.hgvs"])
if _cdna_col is not None:
    c_src = df[_cdna_col].apply(normalize_c_label)

_g_col = _pick_col(["hg38_g_hgvs", "g-position", "hgvs_g", "g.hgvs"])
if _g_col is not None:
    g_src = df[_g_col].apply(normalize_g_label)
# keep raw per-row labels so we can aggregate lists later
df["cDNA_Label_raw"] = c_src
df["gDNA_Label_raw"] = g_src

if c_src is not None and g_src is not None:
    df["Mutation_Label"] = c_src.where(c_src.notna(), g_src).fillna("N/A")
elif c_src is not None:
    df["Mutation_Label"] = c_src.fillna("N/A")
elif g_src is not None:
    df["Mutation_Label"] = g_src.fillna("N/A")
else:
    df["Mutation_Label"] = "N/A"

# Protein label detection
PROT_VALUE_RE = re.compile(
    r"(?:^|[\s;,|])p\.\(?[A-Za-z]{1,3}\d+"
    r"|(?:^|[\s;,|])[A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|Ter|\*|=|fs|del|dup|ins)",
    flags=re.IGNORECASE,
)


def _score_prot_col(series: pd.Series) -> float:
    try:
        vals = series.dropna()
    except Exception:
        return 0.0
    if vals.empty:
        return 0.0
    vals = vals.astype(str).str.strip()
    vals = vals[vals != ""]
    if vals.empty:
        return 0.0
    vals = vals.head(2000)
    try:
        return float(vals.str.contains(PROT_VALUE_RE, na=False).mean())
    except Exception:
        return 0.0


def _guess_prot_col():
    # Prefer explicitly named columns, but verify by content.
    name_keys = [
        "protein", "hgvs_p", "protein_hgvs", "protein_label", "protein label", "p.hgvs",
        "protein_change", "protein change", "hgvs protein", "aa_change", "aa change",
        "amino acid change", "amino_acid_change", "protein_combined"
    ]
    candidates = []
    for k in name_keys:
        c = _pick_col([k])
        if c is not None and c not in candidates:
            candidates.append(c)

    # Add any column whose name suggests protein content
    for c in df.columns:
        cl = str(c).strip().lower()
        if ("protein" in cl) or ("protein_combined" in cl) or ("hgvs_p" in cl) or ("amino" in cl and "acid" in cl):
            if c not in candidates:
                candidates.append(c)

    best_col, best_score = None, -1.0
    for c in candidates:
        sc = _score_prot_col(df[c])
        if sc > best_score:
            best_col, best_score = c, sc

    # Fallback: scan object columns for p.HGVS-like content
    if best_score < 0.02:
        for c in df.columns:
            if c in candidates:
                continue
            if df[c].dtype != object:
                continue
            sc = _score_prot_col(df[c])
            if sc > best_score:
                best_col, best_score = c, sc

    return best_col, best_score


_prot_col, _prot_score = _guess_prot_col()
if _prot_col is not None:
    df["Protein_Label"] = (
        df[_prot_col]
        .astype(str)
        .str.strip()
        .replace(
            {"nan": "", "": np.nan, "None": np.nan, "none": np.nan, "N/A": np.nan, "n/a": np.nan}
        )
    )

    if df["Protein_Label"].notna().sum() == 0:

        def _extract_p_hgvs(x):
            if x is None or (isinstance(x, float) and math.isnan(x)):
                return np.nan
            s = str(x)
            m = re.search(r"p\.\(?[A-Za-z0-9*=_?]+", s)
            return m.group(0) if m else np.nan


        for _col in df.columns:
            if _col == _prot_col or df[_col].dtype != object:
                continue
            miss = df["Protein_Label"].isna()
            if not miss.any():
                break
            extracted = df.loc[miss, _col].apply(_extract_p_hgvs)
            if extracted.notna().any():
                df.loc[miss, "Protein_Label"] = extracted
else:
    df["Protein_Label"] = np.nan

try:
    _nn = int(df["Protein_Label"].notna().sum())
    # print(f"[Protein] Using column: {_prot_col!r} (score={_prot_score:.2f}); non-null Protein_Label={_nn}")
except Exception:
    pass

df["Allele_Freq"] = np.nan

outside_mask = (
        (df["GRCh38 coordinates"] < 48_300_000)
        | (df["GRCh38 coordinates"] > 48_600_000)
)
if outside_mask.any() and not HAVE_PREMAPPED_HG38:
    pass
else:
    if outside_mask.any():
        df = df[~outside_mask].copy()

before_drop = len(df)
df = df.dropna(subset=["GRCh38 coordinates"])

# Quick sanity check: do we have the expected per-variant score columns?
for _sc in ["cadd_score", "spliceai_max", "alpha_missense_pathogenicity_score"]:
    if _sc in df.columns:
        try:
            _nn = int(df[_sc].notna().sum())
            # print(f"[Scores] Found column {_sc!r}: non-null values = {_nn}")
        except Exception:
            pass

after_drop = len(df)

plot_cols = [
    "GRCh38 coordinates",
    "Frequency",
    "Allele_Freq",
    "Mutation_Label",
    "ClinVar_Category",
    "LOVD_Category",
]

if "Protein_Label" in df.columns:
    plot_cols.append("Protein_Label")

# also keep per-row cDNA/gDNA labels so we can build per-variant genomic strings
for extra in ["cDNA_Label_raw", "gDNA_Label_raw"]:
    if extra in df.columns and extra not in plot_cols:
        plot_cols.append(extra)

if "molecular_consequence" in df.columns:
    plot_cols.append("molecular_consequence")
if SOURCE_COLUMN is not None:
    plot_cols.append(SOURCE_COLUMN)
if "genetic_origin" in df.columns:
    plot_cols.append("genetic_origin")
if "laterality" in df.columns:
    plot_cols.append("laterality")

for extra_col in ["gnomad_frequency", "inheritance", "methylation"]:
    if extra_col in df.columns and extra_col not in plot_cols:
        plot_cols.append(extra_col)

# Optional per-variant pathogenicity / splicing scores (from the Excel sheet)
for score_col in ["cadd_score", "spliceai_max", "alpha_missense_pathogenicity_score"]:
    if score_col in df.columns and score_col not in plot_cols:
        plot_cols.append(score_col)

for extra_col in ["chrom", "hg38_genomic_start", "ref_38", "alt_38", "dna_hg38_combined"]:
    if extra_col in df.columns and extra_col not in plot_cols:
        plot_cols.append(extra_col)

plot_df = df[plot_cols].copy()

if "molecular_consequence" in plot_df.columns:
    plot_df["Molecular_Consequence"] = (
        plot_df["molecular_consequence"]
        .astype(str)
        .str.strip()
        .replace({"nan": "Unspecified", "": "Unspecified"})
    )
    plot_df.drop(columns=["molecular_consequence"], inplace=True)
else:
    plot_df["Molecular_Consequence"] = "Unspecified"

plot_df["Region_Label"] = plot_df["GRCh38 coordinates"].apply(
    lambda x: classify_position(int(x), rb1_exons)
)

agg_dict = {
    "Frequency": list,
    "Allele_Freq": list,
    "Mutation_Label": list,
    "Region_Label": "first",
    "ClinVar_Category": lambda x: x.mode()[0]
    if not x.mode().empty
    else "Uncertain significance",
    "LOVD_Category": lambda x: x.mode()[0]
    if not x.mode().empty
    else "Unspecified",
    "Molecular_Consequence": lambda x: list(x),
}
if "Protein_Label" in plot_df.columns:
    agg_dict["Protein_Label"] = list
if "cDNA_Label_raw" in plot_df.columns:
    agg_dict["cDNA_Label_raw"] = list
if "gDNA_Label_raw" in plot_df.columns:
    agg_dict["gDNA_Label_raw"] = list
if "gnomad_frequency" in plot_df.columns:
    agg_dict["gnomad_frequency"] = lambda x: list(x)

if "inheritance" in plot_df.columns:
    agg_dict["inheritance"] = lambda x: list(x)

if "methylation" in plot_df.columns:
    agg_dict["methylation"] = lambda x: list(x)

if any(col in plot_df.columns for col in ["cadd_score", "spliceai_max", "alpha_missense_pathogenicity_score"]):
    # Keep per-variant score lists aligned with Mutation_Label / cDNA_Label / Protein_Label lists
    for score_col in ["cadd_score", "spliceai_max", "alpha_missense_pathogenicity_score"]:
        if score_col in plot_df.columns:
            agg_dict[score_col] = list

if "chrom" in plot_df.columns:
    agg_dict["chrom"] = "first"
if "hg38_genomic_start" in plot_df.columns:
    agg_dict["hg38_genomic_start"] = "first"
if "ref_38" in plot_df.columns:
    agg_dict["ref_38"] = "first"
if "alt_38" in plot_df.columns:
    agg_dict["alt_38"] = "first"
if "dna_hg38_combined" in plot_df.columns:
    agg_dict["dna_hg38_combined"] = "first"

if "genetic_origin" in plot_df.columns:
    agg_dict["genetic_origin"] = lambda x: list(x)
if "laterality" in plot_df.columns:
    agg_dict["laterality"] = lambda x: list(x)

if SOURCE_COLUMN is not None and SOURCE_COLUMN in plot_df.columns:
    agg_dict[SOURCE_COLUMN] = list

grouped_df = (
    plot_df.groupby("GRCh38 coordinates")
    .agg(agg_dict)
    .reset_index()
)

per_variant_map = {}

if {"GRCh38 coordinates", "ref_38", "alt_38"}.issubset(df.columns):
    # Normalize cDNA from original df - same logic as for Mutation_Label
    if "cdna" in df.columns:
        c_src_full = df["cdna"].apply(normalize_c_label)
    elif "cdna " in df.columns:
        c_src_full = df["cdna "].apply(normalize_c_label)
    else:
        c_src_full = None

    for _, r in df.iterrows():
        try:
            pos0 = int(r["GRCh38 coordinates"])
        except Exception:
            continue

        chrom0 = r["chrom"] if "chrom" in r.index else CHROMOSOME
        try:
            chrom0 = canonicalize_chrom(chrom0)
        except Exception:
            chrom0 = CHROMOSOME

        ref0 = _clean_str(r.get("ref_38"))
        alt0 = _clean_str(r.get("alt_38"))

        if not ref0 or not alt0:
            continue

        g_str0 = f"{chrom0}:{pos0} {ref0}>{alt0}"

        cdna_raw0 = None
        if c_src_full is not None:
            cdna_raw0 = c_src_full.get(r.name, None)
        cdna_norm0 = normalize_label_for_hover(cdna_raw0)

        per_variant_map.setdefault(pos0, []).append(
            {
                "cdna_norm": cdna_norm0,
                "g_str": g_str0,
                "ref": ref0,
                "alt": alt0,
            }
        )

if SOURCE_COLUMN is not None and SOURCE_COLUMN in grouped_df.columns:

    def _clean_sources(src_list):
        if not isinstance(src_list, list):
            return []
        out = []
        for v in src_list:
            if pd.isna(v):
                continue
            s = str(v).strip()
            if not s or s.lower() == "nan":
                continue
            out.append(s)
        return sorted(set(out))


    grouped_df["Source_List"] = grouped_df[SOURCE_COLUMN].apply(_clean_sources)


    def _main_source(src_list):
        if not src_list:
            return "Unknown"
        if len(src_list) == 1:
            return src_list[0]
        return "Multiple"


    grouped_df["Source_Category"] = grouped_df["Source_List"].apply(
        _main_source
    )
    grouped_df.drop(columns=[SOURCE_COLUMN], inplace=True)
else:
    grouped_df["Source_List"] = [[] for _ in len(grouped_df) * [None]]
    grouped_df["Source_Category"] = "Unknown"

grouped_df["Hover_Text"] = grouped_df.apply(
    lambda row: build_hover_text(row, per_variant_map),
    axis=1,
)

grouped_df["Total_AC"] = grouped_df["Frequency"].apply(sum)

grouped_df["Jittered_Coord"] = grouped_df["GRCh38 coordinates"].astype(
    float
)

BASE_SOURCE_PALETTE = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]
SOURCE_COLOR_MAP = {}
LEGEND_SOURCE_ORDER = []
LEGEND_ORDER_MAP = {}

if "Source_Category" in grouped_df.columns:
    unique_sources = list(sorted(grouped_df["Source_Category"].unique()))
    if "Multiple" in unique_sources:
        unique_sources = [s for s in unique_sources if s != "Multiple"] + [
            "Multiple"
        ]

    LEGEND_SOURCE_ORDER = unique_sources
    LEGEND_ORDER_MAP = {src: i for i, src in enumerate(LEGEND_SOURCE_ORDER)}

    next_color_idx = 0
    for src in unique_sources:
        if src in FIXED_SOURCE_COLORS and FIXED_SOURCE_COLORS[src]:
            SOURCE_COLOR_MAP[src] = FIXED_SOURCE_COLORS[src]
        else:
            SOURCE_COLOR_MAP[src] = BASE_SOURCE_PALETTE[
                next_color_idx % len(BASE_SOURCE_PALETTE)
                ]
            next_color_idx += 1

if "Source_List" in grouped_df.columns:
    atomic_source_set = set()
    for lst in grouped_df["Source_List"]:
        if isinstance(lst, list):
            for s in lst:
                if s is None:
                    continue
                s_str = str(s).strip()
                if s_str and s_str.lower() != "nan":
                    atomic_source_set.add(s_str)
else:
    atomic_source_set = set()

SOURCE_FILTER_ATOMIC = sorted(atomic_source_set)

CURATED_SOURCE_OPTIONS = [
    {"label": src, "value": src} for src in SOURCE_FILTER_ATOMIC
]
CURATED_SOURCE_OPTIONS_WITH_ALL = with_all_option(CURATED_SOURCE_OPTIONS)
SOURCE_FILTER_OPTIONS = CURATED_SOURCE_OPTIONS_WITH_ALL + [
    {"label": "User data", "value": "User data"}
]

DEFAULT_SOURCE_CURATED_ONLY = ["__ALL__"] + SOURCE_FILTER_ATOMIC
DEFAULT_SOURCE_CURATED_PLUS_USER = DEFAULT_SOURCE_CURATED_ONLY + ["User data"]


def apply_all_logic(new_values, prev_values, option_values, special_exclude=None):
    all_val = "__ALL__"
    if special_exclude is None:
        special_exclude = []

    option_values = list(option_values)
    option_set = set(option_values)

    new_values = [v for v in (new_values or []) if v in option_set]
    prev_values = [v for v in (prev_values or []) if v in option_set]

    indiv_vals = [v for v in option_values if v not in ([all_val] + special_exclude)]
    special_vals = [v for v in option_values if v in special_exclude]

    option_index = {v: i for i, v in enumerate(option_values)}

    new_set = set(new_values)
    prev_set = set(prev_values)

    had_all = all_val in prev_set
    has_all = all_val in new_set

    # 1) All just turned on
    if has_all and not had_all:
        new_set = set(indiv_vals)
        for sv in special_vals:
            if sv in new_values:
                new_set.add(sv)
        new_set.add(all_val)
        return sorted(new_set, key=lambda v: option_index[v])

    # 2) All just turned off
    if not has_all and had_all:
        # Clear all regular options but keep any specials
        new_set = new_set.intersection(set(special_vals))
        return sorted(new_set, key=lambda v: option_index[v])

    # 3) All still selected drop it if selection is partial
    if has_all:
        new_indivs = new_set.intersection(set(indiv_vals))
        if len(new_indivs) < len(indiv_vals):
            new_set.discard(all_val)
        return sorted(new_set, key=lambda v: option_index[v])

    # 4) No 'All'
    new_indivs = new_set.intersection(set(indiv_vals))
    if indiv_vals and len(new_indivs) == len(indiv_vals):
        new_set.add(all_val)
    return sorted(new_set, key=lambda v: option_index[v])


# Genomic Range
_xs = []
if not grouped_df.empty:
    _xs.append(grouped_df["GRCh38 coordinates"].min())
    _xs.append(grouped_df["GRCh38 coordinates"].max())
if rb1_exons:
    _xs.append(min(s for s, _ in rb1_exons))
    _xs.append(max(e for _, e in rb1_exons))
if genomic_domains_df is not None and not genomic_domains_df.empty:
    _xs.append(genomic_domains_df["g_start"].min())
    _xs.append(genomic_domains_df["g_end"].max())

_xs.append(PROMOTER_START)
_xs.append(PROMOTER_END)

if DMR_START is not None and DMR_END is not None:
    _xs.append(DMR_START)
    _xs.append(DMR_END)

GENOMIC_MIN_RAW = min(_xs)
GENOMIC_MAX_RAW = max(_xs)

GENOMIC_MIN = (int(GENOMIC_MIN_RAW) // 100_000) * 100_000
GENOMIC_MAX = ((int(GENOMIC_MAX_RAW) + 100_000 - 1) // 100_000) * 100_000

app = Dash(__name__)
server = app.server

app.layout = html.Div(
    [
        dcc.Store(id="page-state", data="upload"),
        dcc.Store(id="pinned-variants", data=[]),
        dcc.Store(id="curated-grouped-data"),
        # dcc.Store(id="base-figure-store"),  # store for the base figure

        dcc.Store(
            id="clinvar-filter-store",
            data=[opt["value"] for opt in CLINVAR_CLASS_FILTER_OPTIONS],
        ),
        dcc.Store(
            id="mc-filter-store",
            data=[opt["value"] for opt in MCONSEQ_FILTER_OPTIONS],
        ),
        dcc.Store(
            id="genetic-origin-filter-store",
            data=[opt["value"] for opt in GENETIC_ORIGIN_FILTER_OPTIONS],
        ),
        dcc.Store(
            id="laterality-filter-store",
            data=[opt["value"] for opt in LATERALITY_FILTER_OPTIONS],
        ),
        dcc.Store(
            id="source-filter-store",
            data=DEFAULT_SOURCE_CURATED_ONLY,
        ),
        html.Div(
            [
                html.H2(
                    f"{GENE_NAME}",
                    style={
                        "textAlign": "left",
                        "fontFamily": "Helvetica, Arial, sans-serif",
                        "fontWeight": "bold",
                        "fontSize": "28px",
                        "color": "#2c3e50",
                        "margin": "0",
                    },
                ),
                html.Div(
                    APP_BRAND,
                    style={
                        "fontFamily": "Helvetica, Arial, sans-serif",
                        "fontWeight": "800",
                        "fontSize": "32px",
                        "color": "#2c3e50",
                        "letterSpacing": "1px",
                        "marginLeft": "auto",
                    },
                ),
            ],
            style={
                "display": "flex",
                "alignItems": "center",
                "justifyContent": "space-between",
                "marginTop": "10px",
                "marginBottom": "5px",
            },
        ),
        # Upload Page
        html.Div(
            id="upload-page",
            children=[
                html.Div(
                    [
                        html.P(
                            "This tool allows users to explore RB1 variants in retinoblastoma across the genomic region of the RB1 gene."
                            " RB1 variants were compiled from the Leiden Open Variation Database (LOVD), "
                            "11 cohort studies published between 2014 and 2025, and the Catalogue of Somatic "
                            "Mutations in Cancer (COSMIC) database. Variant case information was recovered, "
                            "reannotated in HGVS notation on transcript NM_000321.3, and mapped to the GRCh38 "
                            "reference genome. Users can also upload their own variant files to visualize their "
                            "data locally within the viewer, alongside the curated RB1 dataset.",
                            style={
                                "fontSize": "17px",
                                "color": "#111827",
                                "marginTop": "0px",
                                "marginBottom": "10px",
                                "lineHeight": "1.6",
                            },
                        ),
                    ],
                    style={"marginBottom": "10px"},
                ),
                html.Hr(style={"marginTop": "0px", "marginBottom": "15px"}),
                # First: view existing variants
                html.H3(
                    "Explore RB1 variants",
                    style={
                        "fontWeight": "bold",
                        "marginBottom": "10px",
                        "marginTop": "0px",
                    },
                ),
                html.Div(
                    [
                        html.Button(
                            "View",
                            id="view-lovd-only",
                            n_clicks=0,
                            style={
                                "padding": "6px 16px",
                                "backgroundColor": "#2c3e50",
                                "color": "white",
                                "border": "none",
                                "borderRadius": "4px",
                                "cursor": "pointer",
                            },
                        ),
                    ],
                    style={"marginBottom": "15px"},
                ),
                html.Hr(style={"marginTop": "10px", "marginBottom": "15px"}),
                # Second: optional upload section
                html.H3(
                    "View your own dataset (optional)",
                    style={
                        "fontWeight": "bold",
                        "marginBottom": "10px",
                        "marginTop": "0px",
                    },
                ),
                html.Div(
                    [
                        html.P(
                            "Upload a CSV, TSV or excel file with RB1 variants with the option to include case details. "
                            "The file must include a header row "
                            "and be formatted in one of the following formats.",
                            style={
                                "fontSize": "17px",
                                "color": "#555",
                                "marginBottom": "6px",
                            },
                        ),

                        # Two-column "template" block
                        html.Div(
                            [
                                # Required columns
                                html.Div(
                                    [
                                        html.Strong("Required columns", style={"fontSize": "16px"}),
                                        html.Ul(
                                            [
                                                html.Li(
                                                    [
                                                        html.Span("Option A (1 column): "),
                                                        html.Span("Genomic_positions "),
                                                        html.Span("e.g.,  "),
                                                        html.Span("chr13:48303995 C>G"),
                                                    ]
                                                ),
                                                html.Li(
                                                    [
                                                        html.Span("Option B (4 columns): "),
                                                        html.Span("CHR"),
                                                        html.Span(" / "),
                                                        html.Span("POS"),
                                                        html.Span(" / "),
                                                        html.Span("REF"),
                                                        html.Span(" / "),
                                                        html.Span("ALT"),
                                                        html.Span(" (e.g.,: "),
                                                        html.Span("chr13 | 48303995 | C | G"),
                                                        html.Span(")"),
                                                    ]
                                                ),
                                            ],
                                            style={
                                                "fontSize": "15px",
                                                "color": "#555",
                                                "marginTop": "4px",
                                                "marginBottom": "6px",
                                                "paddingLeft": "18px",
                                            },
                                        ),
                                    ],
                                    style={"flex": "1", "paddingRight": "8px"},
                                ),

                                # Optional columns
                                html.Div(
                                    [
                                        html.Strong("Optional columns", style={"fontSize": "16px"}),
                                        html.Ul(
                                            [
                                                html.Li(
                                                    [
                                                        html.Span("Methylation"),
                                                        html.Span(" (e.g. methylated, unmethylated)"),
                                                    ]
                                                ),
                                                html.Li(
                                                    [
                                                        html.Span("Parent_of_origin"),
                                                        html.Span(" (e.g. maternal, paternal, unknown)"),
                                                    ]
                                                ),
                                                html.Li(
                                                    [
                                                        html.Span("Laterality"),
                                                        html.Span(
                                                            " (e.g. bilateral, unilateral)"),
                                                    ]
                                                ),
                                                html.Li(
                                                    [
                                                        html.Span("Genetic_origin"),
                                                        html.Span(" (e.g. germline, somatic)"),
                                                    ]
                                                ),
                                            ],
                                            style={
                                                "fontSize": "15px",
                                                "color": "#555",
                                                "marginTop": "4px",
                                                "marginBottom": "6px",
                                                "paddingLeft": "18px",
                                            },
                                        ),
                                    ],
                                    style={"flex": "1", "paddingLeft": "8px"},
                                ),
                            ],
                            style={
                                "display": "flex",
                                "flexWrap": "wrap",
                                "border": "1px solid #e5e7eb",
                                "borderRadius": "6px",
                                "padding": "8px 10px",
                                "backgroundColor": "#f9fafb",
                                "marginBottom": "10px",
                            },
                        ),

                    ]
                ),

                dcc.Upload(
                    id="upload-data",
                    children=html.Div(["Drag and Drop or ", html.A("Browse")]),
                    style={
                        "width": "100%",
                        "height": "60px",
                        "lineHeight": "60px",
                        "borderWidth": "1px",
                        "borderStyle": "dashed",
                        "borderRadius": "5px",
                        "textAlign": "center",
                        "marginBottom": "5px",
                    },
                    multiple=False,
                ),
                html.Div(
                    id="upload-info",
                    style={
                        "fontSize": "12px",
                        "color": "#555",
                        "marginBottom": "10px",
                    },
                ),

                html.Div(
                    [
                        html.Label("File delimiter:"),
                        dcc.RadioItems(
                            id="upload-delimiter",
                            options=[
                                {"label": "Comma", "value": ","},
                                {"label": "Tab", "value": "tab"},
                            ],
                            value=",",
                            labelStyle={
                                "display": "inline-block",
                                "marginRight": "15px",
                            },
                            style={"marginTop": "4px"},
                        ),
                    ],
                    style={"marginTop": "10px", "marginBottom": "10px"},
                ),
                html.Div(
                    [
                        html.Button(
                            "Submit",
                            id="upload-submit",
                            n_clicks=0,
                            style={
                                "marginTop": "10px",
                                "padding": "6px 16px",
                                "backgroundColor": "#2c3e50",
                                "color": "white",
                                "border": "none",
                                "borderRadius": "4px",
                                "cursor": "pointer",
                            },
                        ),
                        html.Div(
                            id="upload-error",
                            style={
                                "fontSize": "12px",
                                "color": "#c0392b",
                                "marginTop": "8px",
                            },
                        ),
                    ],
                    style={"marginTop": "10px", "marginBottom": "10px"},
                ),
            ],
            style={"marginBottom": "25px", "display": "block"},
        ),
        # Viewer Page
        html.Div(
            id="viewer-page",
            children=[
                html.Button(
                    "← Back",
                    id="back-to-upload",
                    n_clicks=0,
                    style={
                        "marginBottom": "10px",
                        "padding": "4px 12px",
                        "backgroundColor": "#2c3e50",
                        "color": "#ffffff",
                        "border": "none",
                        "borderRadius": "4px",
                        "cursor": "pointer",
                    },
                ),
                # Filter Blocks
                html.Div(
                    [
                        # Source
                        html.Div(
                            [
                                html.Label("Source:", style=FILTER_LABEL_STYLE),
                                dcc.Checklist(
                                    id="source-filter",
                                    options=SOURCE_FILTER_OPTIONS,
                                    value=DEFAULT_SOURCE_CURATED_ONLY,
                                    labelStyle={
                                        "display": "inline-flex",
                                        "alignItems": "center",
                                        "marginRight": "10px",
                                    },
                                    inputStyle={"marginRight": "4px"},
                                    style=CHECKLIST_INLINE_STYLE,
                                ),
                                html.Div(
                                    "Toggle 'User data' to overlay or hide your uploaded variants.",
                                    style={
                                        "fontSize": "11px",
                                        "color": "#555",
                                        "marginTop": "4px",
                                    },
                                ),
                            ],
                            style={**FILTER_BOX_STYLE, "flex": "0 0 210px"},
                        ),
                        # Pathogenicity (ClinVar)
                        html.Div(
                            [
                                html.Label(
                                    "Pathogenicity (ClinVar):",
                                    style=FILTER_LABEL_STYLE,
                                ),
                                dcc.Checklist(
                                    id="clinvar-path-filter",
                                    options=CLINVAR_CLASS_FILTER_OPTIONS,
                                    value=[
                                        opt["value"]
                                        for opt in CLINVAR_CLASS_FILTER_OPTIONS
                                    ],
                                    labelStyle={
                                        "display": "inline-flex",
                                        "alignItems": "center",
                                        "marginRight": "10px",
                                    },
                                    inputStyle={"marginRight": "4px"},
                                    style=CHECKLIST_INLINE_STYLE,
                                ),
                                html.Div(
                                    "Note: Only ~30% of variants are classified within ClinVar.",
                                    style={
                                        "fontSize": "11px",
                                        "color": "#555",
                                        "marginTop": "4px",
                                    },
                                ),
                            ],
                            style={**FILTER_BOX_STYLE, "flex": "1 1 260px"},
                        ),
                        # Genetic Origin
                        html.Div(
                            [
                                html.Label(
                                    "Genetic Origin:", style=FILTER_LABEL_STYLE
                                ),
                                dcc.Checklist(
                                    id="genetic-origin-filter",
                                    options=GENETIC_ORIGIN_FILTER_OPTIONS,
                                    value=[
                                        opt["value"]
                                        for opt in GENETIC_ORIGIN_FILTER_OPTIONS
                                    ],
                                    labelStyle={
                                        "display": "inline-flex",
                                        "alignItems": "center",
                                        "marginRight": "10px",
                                    },
                                    inputStyle={"marginRight": "4px"},
                                    style=CHECKLIST_INLINE_STYLE,
                                ),
                            ],
                            style={**FILTER_BOX_STYLE, "flex": "0 0 180px"},
                        ),
                        # Laterality
                        html.Div(
                            [
                                html.Label(
                                    "Laterality:", style=FILTER_LABEL_STYLE
                                ),
                                dcc.Checklist(
                                    id="laterality-filter",
                                    options=LATERALITY_FILTER_OPTIONS,
                                    value=[
                                        opt["value"]
                                        for opt in LATERALITY_FILTER_OPTIONS
                                    ],
                                    labelStyle={
                                        "display": "inline-flex",
                                        "alignItems": "center",
                                        "marginRight": "10px",
                                    },
                                    inputStyle={"marginRight": "4px"},
                                    style=CHECKLIST_INLINE_STYLE,
                                ),
                            ],
                            style={**FILTER_BOX_STYLE, "flex": "0 0 180px"},
                        ),
                        # Molecular Consequence
                        html.Div(
                            [
                                html.Label(
                                    "Molecular Consequence:",
                                    style=FILTER_LABEL_STYLE,
                                ),
                                dcc.Checklist(
                                    id="mc-filter",
                                    options=MCONSEQ_FILTER_OPTIONS,
                                    value=[
                                        opt["value"]
                                        for opt in MCONSEQ_FILTER_OPTIONS
                                    ],
                                    labelStyle={
                                        "display": "inline-flex",
                                        "alignItems": "center",
                                        "marginRight": "10px",
                                    },
                                    inputStyle={"marginRight": "4px"},
                                    style=CHECKLIST_INLINE_STYLE,
                                ),
                            ],
                            style={**FILTER_BOX_STYLE, "flex": "1 1 260px"},
                        ),
                    ],
                    style=FILTER_PANEL_STYLE,
                ),
                # Search bar
                html.Div(
                    [
                        html.Div(
                            [
                                html.Label(
                                    "Search for a specific variant:",
                                    style={
                                        "fontWeight": "600",
                                        "fontSize": "18px",
                                        "marginRight": "8px",
                                    },
                                ),
                                dcc.Input(
                                    id="variant-search-input",
                                    type="text",
                                    placeholder="e.g., c.958C>T / p.Arg320* / chr13:48430295 C>T",
                                    style={
                                        "width": "380px",
                                        "marginRight": "8px",
                                        "height": "32px",
                                    },
                                ),
                                html.Button(
                                    "Search",
                                    id="variant-search-button",
                                    n_clicks=0,
                                    style={
                                        "padding": "4px 12px",
                                        "backgroundColor": "#2c3e50",
                                        "color": "#ffffff",
                                        "border": "none",
                                        "borderRadius": "4px",
                                        "cursor": "pointer",
                                        "marginRight": "6px",
                                    },
                                ),
                                html.Button(
                                    "Reset",
                                    id="variant-search-reset",
                                    n_clicks=0,
                                    style={
                                        "padding": "4px 12px",
                                        "backgroundColor": "#2c3e50",
                                        "color": "#ffffff",
                                        "border": "1px solid #d1d5db",
                                        "borderRadius": "4px",
                                        "cursor": "pointer",
                                    },
                                ),
                            ],
                            style={
                                "display": "flex",
                                "flexWrap": "wrap",
                                "alignItems": "center",
                                "gap": "6px",
                                "marginBottom": "4px",
                            },
                        ),
                        # Short status / error messages
                        html.Div(
                            id="variant-search-feedback",
                            style={
                                "fontSize": "11px",
                                "color": "#555",
                                "marginTop": "4px",
                            },
                        ),
                        # scrollable list of pinned variants
                        html.Div(
                            id="pinned-variants-list",
                            style={
                                "fontSize": "11px",
                                "color": "#111827",
                                "marginTop": "6px",
                                "maxHeight": "110px",
                                "overflowY": "auto",
                                "overflowX": "auto",
                                "borderTop": "1px solid #e5e7eb",
                                "paddingTop": "4px",
                            },
                        ),
                    ],
                    style={
                        "border": "1px solid #e1e5eb",
                        "borderRadius": "8px",
                        "padding": "8px 12px",
                        "backgroundColor": "#f9fafb",
                        "marginBottom": "10px",
                    },
                ),

                dcc.Loading(
                    id="mutation-plot-loading",
                    type="default",
                    fullscreen=True,
                    children=html.Div(
                        [

                            dcc.Store(id="base-figure-store"),

                            dcc.Graph(
                                id="mutation-plot",
                                config={
                                    "modeBarButtonsToRemove": ["select2d", "lasso2d"],
                                    "displaylogo": False,
                                },
                                clear_on_unhover=False,
                            ),
                        ],
                        style={"height": "100%"},
                    ),
                ),

                html.Button("Download CSV", id="btn_csv"),
                dcc.Download(id="download-dataframe-csv"),
            ],
            style={"display": "none"},
        ),
    ],
    style={
        "fontFamily": "Helvetica, Arial, sans-serif",
        "marginLeft": "40px",
        "marginRight": "40px",
    },
)


@app.callback(
    Output("page-state", "data"),
    Output("upload-error", "children"),
    Input("upload-submit", "n_clicks"),
    Input("view-lovd-only", "n_clicks"),
    Input("back-to-upload", "n_clicks"),
    State("page-state", "data"),
    State("upload-data", "contents"),
)
def update_page_state(
        n_submit, n_view_lovd, n_back, current_state, upload_contents
):
    ctx = dash.callback_context

    if not ctx.triggered:
        return (current_state or "upload"), ""

    trigger = ctx.triggered[0]["prop_id"].split(".")[0]

    if trigger == "upload-submit":
        if upload_contents is None:
            return "upload", "Please upload a file before clicking 'Submit'."
        else:
            return "viewer", ""

    elif trigger == "view-lovd-only":
        return "viewer", ""

    elif trigger == "back-to-upload":
        return "upload", ""

    return (current_state or "upload"), ""


@app.callback(
    Output("upload-page", "style"),
    Output("viewer-page", "style"),
    Input("page-state", "data"),
)
def apply_page_state(page_state):
    if page_state == "viewer":
        return (
            {"display": "none", "marginBottom": "25px"},
            {"display": "block"},
        )
    else:
        return (
            {"display": "block", "marginBottom": "25px"},
            {"display": "none"},
        )


@app.callback(
    Output("upload-info", "children"),
    Input("upload-data", "contents"),
    State("upload-data", "filename"),
    # State("upload-has-header", "value"),
    State("upload-delimiter", "value"),
)
def update_upload_info(
        contents, filename, delimiter_value
):
    if contents is None:
        return "No user file uploaded."
    name = filename or "uploaded file"
    # has_header = has_header_value == "yes"
    delim = delimiter_value or ","
    try:
        user_df = parse_user_variant_upload(
            contents,
            filename,
            has_header=True,
            delimiter=delim,
        )
    except Exception as e:
        return f"Error reading {name}: {e}"
    if user_df is None or user_df.empty:
        return f"{name}: no valid {CHROMOSOME} variants found."
    return f"{name}: loaded {len(user_df)} unique genomic position(s)."


@app.callback(
    Output("clinvar-path-filter", "value"),
    Output("clinvar-filter-store", "data"),
    Input("clinvar-path-filter", "value"),
    State("clinvar-filter-store", "data"),
)
def sync_clinvar_filter(current, previous):
    option_values = [opt["value"] for opt in CLINVAR_CLASS_FILTER_OPTIONS]
    updated = apply_all_logic(current, previous, option_values)
    return updated, updated


@app.callback(
    Output("mc-filter", "value"),
    Output("mc-filter-store", "data"),
    Input("mc-filter", "value"),
    State("mc-filter-store", "data"),
)
def sync_mc_filter(current, previous):
    option_values = [opt["value"] for opt in MCONSEQ_FILTER_OPTIONS]
    updated = apply_all_logic(current, previous, option_values)
    return updated, updated


@app.callback(
    Output("genetic-origin-filter", "value"),
    Output("genetic-origin-filter-store", "data"),
    Input("genetic-origin-filter", "value"),
    State("genetic-origin-filter-store", "data"),
)
def sync_go_filter(current, previous):
    option_values = [opt["value"] for opt in GENETIC_ORIGIN_FILTER_OPTIONS]
    updated = apply_all_logic(current, previous, option_values)
    return updated, updated


@app.callback(
    Output("laterality-filter", "value"),
    Output("laterality-filter-store", "data"),
    Input("laterality-filter", "value"),
    State("laterality-filter-store", "data"),
)
def sync_laterality_filter(current, previous):
    option_values = [opt["value"] for opt in LATERALITY_FILTER_OPTIONS]
    updated = apply_all_logic(current, previous, option_values)
    return updated, updated


@app.callback(
    Output("source-filter", "value"),
    Output("source-filter-store", "data"),
    Input("source-filter", "value"),
    State("source-filter-store", "data"),
)
def sync_source_filter(current, previous):
    option_values = [opt["value"] for opt in SOURCE_FILTER_OPTIONS]
    # "User data" should NOT be auto-selected/deselected by "__ALL__".
    updated = apply_all_logic(
        current,
        previous,
        option_values,
        special_exclude=["User data"],
    )
    return updated, updated


@app.callback(
    Output("base-figure-store", "data"),
    Output("curated-grouped-data", "data"),
    Input("clinvar-path-filter", "value"),
    Input("mc-filter", "value"),
    Input("genetic-origin-filter", "value"),
    Input("laterality-filter", "value"),
    Input("upload-data", "contents"),
    Input("source-filter", "value"),
    # Input("mutation-plot", "relayoutData"),
    State("upload-data", "filename"),
    # State("upload-has-header", "value"),
    State("upload-delimiter", "value"),
)
def update_plots(
        selected_clinvar_classes,
        selected_mc_list,
        selected_genetic_origins,
        selected_laterality,
        upload_contents,
        source_filter,
        # relayout_data,
        upload_filename,
        # upload_header_value,
        upload_delim_value,
):
    def normalize_all(vals):
        if not vals:
            return []
        if "__ALL__" in vals:
            return ["__ALL__"]
        return vals

    selected_clinvar_classes = normalize_all(selected_clinvar_classes)
    selected_mc_list = normalize_all(selected_mc_list)
    selected_genetic_origins = normalize_all(selected_genetic_origins)
    selected_laterality = normalize_all(selected_laterality)
    source_filter = source_filter or DEFAULT_SOURCE_CURATED_ONLY

    data = plot_df.copy()

    # ClinVar
    if (
            selected_clinvar_classes != ["__ALL__"]
            and "ClinVar_Category" in data.columns
    ):
        clinvar_set = set(selected_clinvar_classes)
        data = data[data["ClinVar_Category"].isin(clinvar_set)]

    # Molecular consequence
    if (
            selected_mc_list != ["__ALL__"]
            and "Molecular_Consequence" in data.columns
    ):
        mc_set = set(selected_mc_list)
        data = data[data["Molecular_Consequence"].isin(mc_set)]

    # Genetic origin
    if (
            selected_genetic_origins != ["__ALL__"]
            and "genetic_origin" in data.columns
    ):
        go_set = set(selected_genetic_origins)
        data = data[data["genetic_origin"].isin(go_set)]

    # Laterality
    if (
            selected_laterality != ["__ALL__"]
            and "laterality" in data.columns
    ):
        lat_set = set(selected_laterality)
        data = data[data["laterality"].isin(lat_set)]

    selected_sources = source_filter
    include_user = "User data" in selected_sources

    if SOURCE_COLUMN and SOURCE_COLUMN in data.columns:
        curated_source_choices = [
            s for s in selected_sources if s not in {"User data", "__ALL__"}
        ]

        if "__ALL__" not in selected_sources:
            if curated_source_choices:
                sel_set = set(curated_source_choices)
                data = data[data[SOURCE_COLUMN].isin(sel_set)]
            else:
                data = data.iloc[0:0]

    agg_dict_filtered = {
        "Frequency": list,
        "Allele_Freq": list,
        "Mutation_Label": list,
        "Region_Label": "first",
        "ClinVar_Category": lambda x: x.mode()[0]
        if not x.mode().empty
        else "Uncertain significance",
        "LOVD_Category": lambda x: x.mode()[0]
        if not x.mode().empty
        else "Unspecified",
        "Molecular_Consequence": lambda x: list(x),
    }
    if "Protein_Label" in data.columns:
        agg_dict_filtered["Protein_Label"] = list
    if "cDNA_Label_raw" in data.columns:
        agg_dict_filtered["cDNA_Label_raw"] = list
    if "gDNA_Label_raw" in data.columns:
        agg_dict_filtered["gDNA_Label_raw"] = list

    if "chrom" in data.columns:
        agg_dict_filtered["chrom"] = "first"
    if "hg38_genomic_start" in data.columns:
        agg_dict_filtered["hg38_genomic_start"] = "first"
    if "ref_38" in data.columns:
        agg_dict_filtered["ref_38"] = "first"
    if "alt_38" in data.columns:
        agg_dict_filtered["alt_38"] = "first"
    if "dna_hg38_combined" in data.columns:
        agg_dict_filtered["dna_hg38_combined"] = "first"

    if "genetic_origin" in data.columns:
        agg_dict_filtered["genetic_origin"] = lambda x: list(x)
    if "laterality" in data.columns:
        agg_dict_filtered["laterality"] = lambda x: list(x)
    if SOURCE_COLUMN and SOURCE_COLUMN in data.columns:
        agg_dict_filtered[SOURCE_COLUMN] = list
    if "gnomad_frequency" in data.columns:
        agg_dict_filtered["gnomad_frequency"] = lambda x: list(x)

    if "inheritance" in data.columns:
        agg_dict_filtered["inheritance"] = lambda x: list(x)
    if "methylation" in data.columns:
        agg_dict_filtered["methylation"] = lambda x: list(x)

    for _score_col in (
            "cadd_score",
            "spliceai_max",
            "alpha_missense_pathogenicity_score",
    ):
        if _score_col in data.columns:
            agg_dict_filtered[_score_col] = list

    if data.empty:
        grouped_filtered = pd.DataFrame(
            columns=[
                "GRCh38 coordinates",
                "Total_AC",
                "Hover_Text",
                "Jittered_Coord",
            ]
        )
    else:
        grouped_filtered = (
            data.groupby("GRCh38 coordinates")
            .agg(agg_dict_filtered)
            .reset_index()
        )

        if SOURCE_COLUMN and SOURCE_COLUMN in grouped_filtered.columns:

            def _clean_sources_atomic(src_list):

                if not isinstance(src_list, list):
                    src_list = [src_list]
                out = []
                for v in src_list:
                    try:
                        if pd.isna(v):
                            out.append(None)
                            continue
                    except Exception:
                        pass
                    s = str(v).strip()
                    if not s or s.lower() == "nan":
                        out.append(None)
                    else:
                        out.append(s)
                return out

            def _unique_sources(src_atomic):
                return sorted(
                    {
                        s
                        for s in (src_atomic or [])
                        if s is not None and str(s).strip()
                           and str(s).strip().lower() != "nan"
                    }
                )

            grouped_filtered["Source_Atomic_List"] = grouped_filtered[
                SOURCE_COLUMN
            ].apply(_clean_sources_atomic)

            grouped_filtered["Source_List"] = grouped_filtered[
                "Source_Atomic_List"
            ].apply(_unique_sources)

            def _main_source(src_list):
                if not src_list:
                    return "Unknown"
                if len(src_list) == 1:
                    return src_list[0]
                return "Multiple"

            grouped_filtered["Source_Category"] = grouped_filtered[
                "Source_List"
            ].apply(_main_source)
            grouped_filtered["Display_Source"] = grouped_filtered[
                "Source_Category"
            ]

            grouped_filtered.drop(columns=[SOURCE_COLUMN], inplace=True)

        else:
            grouped_filtered["Source_List"] = [[] for _ in range(len(grouped_filtered))]
            grouped_filtered["Source_Category"] = "Unknown"
            grouped_filtered["Display_Source"] = "Unknown"

        grouped_filtered["Hover_Text"] = grouped_filtered.apply(
            lambda row: build_hover_text(row, per_variant_map),
            axis=1,
        )

        grouped_filtered["Total_AC"] = grouped_filtered["Frequency"].apply(
            sum
        )

        grouped_filtered["Jittered_Coord"] = grouped_filtered[
            "GRCh38 coordinates"
        ].astype(float)

        if LEGEND_ORDER_MAP:
            grouped_filtered["legend_order"] = grouped_filtered[
                "Display_Source"
            ].map(LEGEND_ORDER_MAP).fillna(len(LEGEND_ORDER_MAP))
            grouped_filtered = grouped_filtered.sort_values(
                ["legend_order", "GRCh38 coordinates"]
            )

    data = grouped_filtered
    curated_store_data = data.to_dict("records")

    user_data = None
    user_has_haplotype = False
    if include_user and upload_contents is not None:
        # has_header = upload_header_value == "yes"
        delim = upload_delim_value or ","
        try:
            tmp_user = parse_user_variant_upload(
                upload_contents,
                upload_filename,
                has_header=True,
                delimiter=delim,
            )
            if tmp_user is not None and not tmp_user.empty:
                user_data = tmp_user

                if (
                        "Hyper_Count" in user_data.columns
                        and "Hypo_Count" in user_data.columns
                ):
                    nz_mask = (
                                      user_data["Hyper_Count"]
                                      .fillna(0)
                                      .astype(int)
                                      .abs()
                                      + user_data["Hypo_Count"]
                                      .fillna(0)
                                      .astype(int)
                                      .abs()
                              ) > 0
                    user_has_haplotype = bool(nz_mask.any())
        except Exception as e:
            print(f"[User upload] Failed to parse uploaded file: {e}")
            user_data = None
            user_has_haplotype = False

    total_rows = 7
    exon_row = 2
    maternal_row = 3
    paternal_row = 4
    protein_header_row = 5
    protein_lollipop_row = 6
    protein_boxes_row = 7

    # Used to anchor the Protein View y-axis so the baseline is exactly at 0
    protein_max_y = 1.0

    # Fractions of total height (sum ~ 1.0)
    row_heights = [0.41, 0.08, 0.12, 0.12, 0.04, 0.15, 0.08]
    fig = sp.make_subplots(
        rows=total_rows,
        cols=1,
        shared_xaxes=False,  # protein rows use AA coordinates; other rows use genomic
        row_heights=row_heights,
        vertical_spacing=0.06,
    )

    # Accumulate shapes and apply once at the end (much faster than many fig.add_shape calls)
    shapes = []

    # Main curated variants (row 1)
    seen_sources = set()
    if not data.empty:
        for _, row in data.iterrows():
            x = row["Jittered_Coord"]
            y = row["Total_AC"]
            src_cat = row.get("Display_Source", "Unknown")
            color = SOURCE_COLOR_MAP.get(src_cat, "#000000")

            showlegend = False
            if src_cat not in seen_sources:
                showlegend = True
                seen_sources.add(src_cat)

            # stick
            fig.add_trace(
                go.Scatter(
                    x=[x, x],
                    y=[0, y],
                    mode="lines",
                    line=dict(color=color, width=1),
                    hoverinfo="text",
                    text=row["Hover_Text"],
                    showlegend=False,
                ),
                row=1,
                col=1,
            )

            # head
            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[y],
                    mode="markers",
                    marker=dict(
                        size=9,
                        color=color,
                        line=dict(color="white", width=1),
                    ),
                    text=row["Hover_Text"],
                    hoverinfo="text",
                    name=src_cat,
                    showlegend=showlegend,
                ),
                row=1,
                col=1,
            )

    #  User variants on main panel
    user_legend_added = False
    if user_data is not None and not user_data.empty:
        for _, row in user_data.iterrows():
            ux = row["Jittered_Coord"]
            uy = row["Total_AC"]

            fig.add_trace(
                go.Scatter(
                    x=[ux, ux],
                    y=[0, uy],
                    mode="lines",
                    line=dict(color=USER_VARIANT_COLOR, width=3),
                    hoverinfo="text",
                    text=row["Hover_Text"],
                    showlegend=False,
                ),
                row=1,
                col=1,
            )

            fig.add_trace(
                go.Scatter(
                    x=[ux],
                    y=[uy],
                    mode="markers",
                    marker=dict(
                        size=9,
                        color=USER_VARIANT_COLOR,
                        line=dict(color="white", width=1),
                    ),
                    text=row["Hover_Text"],
                    hoverinfo="text",
                    name="User data",
                    showlegend=not user_legend_added,
                ),
                row=1,
                col=1,
            )
            user_legend_added = True

    max_y = 1
    if not data.empty:
        max_y = max(max_y, data["Total_AC"].max())
    if user_data is not None and not user_data.empty:
        max_y = max(max_y, user_data["Total_AC"].max())

    hyper_legend_added = False  # for maternal panel
    hypo_legend_added = False  # for paternal panel

    def _norm_meth(v):
        if v is None: return None
        s = str(v).strip().lower()
        if not s or s in {"nan", "none", "null"}: return None
        if "unmeth" in s or "hypo" in s: return "unmethylated"
        if "meth" in s or "hyper" in s: return "methylated"
        return None

    def _norm_inh(v):
        if v is None: return None
        s = str(v).strip().lower()
        if not s or s in {"nan", "none", "unknown"}: return None
        if "mater" in s: return "maternal"
        if "pater" in s: return "paternal"
        return None

    # Curated methylation points derived from the aggregated per-coordinate rows
    curated_points = []
    if not data.empty and ("methylation" in data.columns or "inheritance" in data.columns):
        for _, r in data.iterrows():
            pos = r.get("GRCh38 coordinates")
            x = r.get("Jittered_Coord")
            meth_list = r.get("methylation") if "methylation" in r.index else None
            inh_list = r.get("inheritance") if "inheritance" in r.index else None

            # Determine how many per-variant entries we have at this coordinate.
            _mut_list = r.get("Mutation_Label")
            if isinstance(_mut_list, list):
                _n_variants = len(_mut_list)
            elif _mut_list is None or (isinstance(_mut_list, float) and pd.isna(_mut_list)):
                _n_variants = 0
            else:
                _n_variants = 1

            if _n_variants <= 0:
                continue

            # Normalize per-variant methylation/inheritance lists so we can filter hover text.
            if isinstance(meth_list, list):
                meth_norm = [_norm_meth(vv) for vv in meth_list]
            else:
                meth_norm = [_norm_meth(meth_list)] * _n_variants
            if len(meth_norm) < _n_variants:
                meth_norm = (meth_norm + [None] * _n_variants)[:_n_variants]

            if isinstance(inh_list, list):
                inh_norm = [_norm_inh(vv) for vv in inh_list]
            else:
                inh_norm = [_norm_inh(inh_list)] * _n_variants
            if len(inh_norm) < _n_variants:
                inh_norm = (inh_norm + [None] * _n_variants)[:_n_variants]

            # Split into methylated/unmethylated subsets. If a site has BOTH (e.g., the same variant
            # appears as methylated in some rows and unmethylated in others), we draw lollipops in BOTH panels.
            idx_meth = [i for i, nv in enumerate(meth_norm) if nv == "methylated"]
            idx_unmeth = [i for i, nv in enumerate(meth_norm) if nv == "unmethylated"]

            # Family-history-only points (inheritance known but methylation missing/unknown at those rows)
            idx_fh_m = [i for i, iv in enumerate(inh_norm) if iv == "maternal" and meth_norm[i] is None]
            idx_fh_p = [i for i, iv in enumerate(inh_norm) if iv == "paternal" and meth_norm[i] is None]

            # Skip coordinates that have neither methylation nor inheritance information
            if (not idx_meth) and (not idx_unmeth) and (not idx_fh_m) and (not idx_fh_p):
                continue

            freq_list = r.get("Frequency")

            def _sum_counts(idxs):
                if not isinstance(freq_list, list) or not idxs:
                    return 1
                s = 0
                for i in idxs:
                    if 0 <= i < len(freq_list):
                        try:
                            s += int(freq_list[i])
                        except Exception:
                            s += 1
                return max(1, s)

            def _poo_from_idxs(idxs):
                # Prefer a parent-of-origin value that exists within this subset
                for i in idxs:
                    if 0 <= i < len(inh_norm) and inh_norm[i] is not None:
                        return inh_norm[i]
                # fallback: any non-null from the full list
                for v in inh_norm:
                    if v is not None:
                        return v
                return "unknown"

            def _build_hover(idxs, report_count, poo_txt, meth_txt):
                base_hover = None
                try:
                    base_hover = build_hover_text_subset(r, idxs)
                except Exception:
                    base_hover = None

                if not isinstance(base_hover, str) or not base_hover.strip():
                    base_hover = r.get("Hover_Text")

                if not isinstance(base_hover, str) or not base_hover.strip():
                    try:
                        pos_i = int(pos)
                    except Exception:
                        pos_i = pos
                    base_hover = "<b>RB1</b><br>" + f"<b>Genomic Position:</b> {pos_i}"

                return (
                        base_hover
                        + f"<br><b>Report Count:</b> {report_count}"
                        + f"<br><b>Parent of origin:</b> {poo_txt}"
                        + f"<br><b>Methylation:</b> {meth_txt}"
                )

            # --- Methylated (maternal panel) ---
            if idx_meth:
                y_count = _sum_counts(idx_meth)
                poo_txt = _poo_from_idxs(idx_meth)
                hover_text = _build_hover(idx_meth, y_count, poo_txt, "methylated")
                curated_points.append(("maternal", x, y_count, HYPER_METHYL_COLOR, hover_text))

            # --- Unmethylated (paternal panel) ---
            if idx_unmeth:
                y_count = _sum_counts(idx_unmeth)
                poo_txt = _poo_from_idxs(idx_unmeth)
                hover_text = _build_hover(idx_unmeth, y_count, poo_txt, "unmethylated")
                curated_points.append(("paternal", x, y_count, HYPO_METHYL_COLOR, hover_text))

            # --- Family history / inheritance-only (gray) ---
            # Draw these even if methylated/unmethylated also exist at the same coordinate,
            # but only for the specific rows where methylation is missing.
            if idx_fh_m:
                y_count = _sum_counts(idx_fh_m)
                hover_text = _build_hover(idx_fh_m, y_count, "maternal", "unknown")
                curated_points.append(("maternal", x, y_count, INHERITANCE_GRAY, hover_text))

            if idx_fh_p:
                y_count = _sum_counts(idx_fh_p)
                hover_text = _build_hover(idx_fh_p, y_count, "paternal", "unknown")
                curated_points.append(("paternal", x, y_count, INHERITANCE_GRAY, hover_text))

    hyper_max = 1
    hypo_max = 1
    for trk, _x, y, _c, _h in curated_points:
        if trk == "maternal":
            hyper_max = max(hyper_max, y)
        elif trk == "paternal":
            hypo_max = max(hypo_max, y)

    # Draw curated points first
    for trk, x, y, color, hover_text in curated_points:
        target_row = maternal_row if trk == "maternal" else paternal_row
        if target_row is None:
            continue
        # stick
        fig.add_trace(
            go.Scatter(
                x=[x, x], y=[0, y], mode="lines",
                line=dict(color=color, width=1),
                hoverinfo="text", text=hover_text, showlegend=False,
            ),
            row=target_row, col=1,
        )
        # head
        fig.add_trace(
            go.Scatter(
                x=[x], y=[y], mode="markers",
                marker=dict(size=9, color=color, line=dict(color='white', width=1)),
                hoverinfo="text", text=hover_text,
                name="Curated methylation",
                showlegend=False,
            ),
            row=target_row, col=1,
        )

    # Overlay user haplotype (black) if present
    if user_data is not None and not user_data.empty and (
            "Hyper_Count" in user_data.columns or "Hypo_Count" in user_data.columns):
        # update maxima with user counts
        for _, row_u in user_data.iterrows():
            hyper_max = max(hyper_max, int(row_u.get("Hyper_Count", 0) or 0))
            hypo_max = max(hypo_max, int(row_u.get("Hypo_Count", 0) or 0))

        for _, row_u in user_data.iterrows():
            x = row_u["Jittered_Coord"]
            pos = row_u["GRCh38 coordinates"]
            hyper = int(row_u.get("Hyper_Count", 0) or 0)
            hypo = int(row_u.get("Hypo_Count", 0) or 0)

            # Build methylation-specific hovers so each lollipop only shows the
            # variants that belong to that panel (methylated vs. unmethylated).
            def _build_user_meth_hover(which, count):
                try:
                    pos_i = int(pos)
                except Exception:
                    pos_i = pos

                lines = [
                    "<b>User variant(s)</b>",
                    f"<b>Genomic Position:</b> {pos_i}",
                ]

                region_u = row_u.get("Region_Label")
                if isinstance(region_u, str) and region_u.strip():
                    lines.append(f"<b>Region:</b> {region_u}")

                lines.append(f"<b>Report Count:</b> {int(count)}")

                labels = row_u.get("Labels")
                haplos = row_u.get("HaploList")
                target = "hyper" if which == "hyper" else "hypo"

                subset_labels = []
                if isinstance(labels, list) and isinstance(haplos, list):
                    for lbl, h in zip(labels, haplos):
                        if h == target:
                            subset_labels.append(lbl)

                label_counts = {}
                for lbl in subset_labels:
                    if lbl is None or (isinstance(lbl, float) and pd.isna(lbl)):
                        continue
                    s = str(lbl).strip()
                    if not s:
                        continue
                    label_counts[s] = label_counts.get(s, 0) + 1

                if label_counts:
                    lines.append("<b>Variant label(s):</b>")
                    for lbl, cnt in sorted(label_counts.items(), key=lambda kv: (-kv[1], kv[0])):
                        lines.append(f"{lbl}, Report Count: {cnt}")

                if which == "hyper":
                    lines.append(f"<b>Methylation:</b> Methylated: {int(count)}")
                else:
                    lines.append(f"<b>Methylation:</b> Unmethylated: {int(count)}")

                return "<br>".join(lines)

            # Fallback: use precomputed hover if something is missing.
            base_hover = row_u.get("Hover_Text")
            if not isinstance(base_hover, str) or not base_hover.strip():
                base_hover = "<b>User variant(s)</b><br>" + f"<b>Genomic Position:</b> {int(pos)}"

            if hyper > 0 and maternal_row is not None:
                try:
                    hover_text = _build_user_meth_hover("hyper", hyper)
                except Exception:
                    hover_text = base_hover + f"<br><b>Methylated report count:</b> {hyper}"
                fig.add_trace(
                    go.Scatter(x=[x, x], y=[0, hyper], mode="lines",
                               line=dict(color=USER_HAPLO_COLOR, width=3),
                               hoverinfo="text", text=hover_text, showlegend=False),
                    row=maternal_row, col=1,
                )
                fig.add_trace(
                    go.Scatter(x=[x], y=[hyper], mode="markers",
                               marker=dict(size=9, color=USER_HAPLO_COLOR, line=dict(color='white', width=1)),
                               hoverinfo="text", text=hover_text,
                               name="User methylation", showlegend=False),
                    row=maternal_row, col=1,
                )

            if hypo > 0 and paternal_row is not None:
                try:
                    hover_text = _build_user_meth_hover("hypo", hypo)
                except Exception:
                    hover_text = base_hover + f"<br><b>Unmethylated report count:</b> {hypo}"
                fig.add_trace(
                    go.Scatter(x=[x, x], y=[0, hypo], mode="lines",
                               line=dict(color=USER_HAPLO_COLOR, width=3),
                               hoverinfo="text", text=hover_text, showlegend=False),
                    row=paternal_row, col=1,
                )
                fig.add_trace(
                    go.Scatter(x=[x], y=[hypo], mode="markers",
                               marker=dict(size=9, color=USER_HAPLO_COLOR, line=dict(color='white', width=1)),
                               hoverinfo="text", text=hover_text,
                               name="User methylation", showlegend=False),
                    row=paternal_row, col=1,
                )

    # y-axes for haplotype panels (methylation tracks)
    # Show integer tick labels so "Report Count" is meaningful.
    def _dtick_for_counts(m):
        try:
            m = int(m)
        except Exception:
            m = 1
        if m <= 5:
            return 1
        if m <= 20:
            return 5
        if m <= 50:
            return 10
        if m <= 100:
            return 20
        if m <= 200:
            return 50
        return 100

    if maternal_row is not None:
        fig.update_yaxes(
            title_text="Methylated/<br>Maternal Allele",
            showgrid=False,
            zeroline=False,
            showticklabels=True,
            tick0=0,
            dtick=_dtick_for_counts(hyper_max),
            tickformat=",d",
            fixedrange=True,
            range=[0, max(hyper_max * 1.2, hyper_max + 1)],
            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="outside",
            tickcolor="#444",
            tickwidth=1,
            ticklen=4,
            row=maternal_row,
            col=1,
        )
    if paternal_row is not None:
        fig.update_yaxes(
            title_text="Unmethylated/<br>Paternal Allele",
            showgrid=False,
            zeroline=False,
            showticklabels=True,
            tick0=0,
            dtick=_dtick_for_counts(hypo_max),
            tickformat=",d",
            fixedrange=True,
            range=[0, max(hypo_max * 1.2, hypo_max + 1)],
            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="outside",
            tickcolor="#444",
            tickwidth=1,
            ticklen=4,
            row=paternal_row,
            col=1,
        )

    # Put panel labels at the *top* of each methylation track (inside the blue area)
    try:
        if maternal_row is not None:
            ykey = "yaxis" if maternal_row == 1 else f"yaxis{maternal_row}"
            dom = getattr(fig.layout, ykey).domain
            fig.add_annotation(
                xref="paper", yref="paper",
                x=0.0, y=dom[1],
                #text="<b>Methylated / Maternal</b>",
                text=" ",
                showarrow=False,
                xanchor="left", yanchor="top",
                xshift=8, yshift=-6,
                font=dict(size=12, color="#444"),
            )
        if paternal_row is not None:
            ykey = "yaxis" if paternal_row == 1 else f"yaxis{paternal_row}"
            dom = getattr(fig.layout, ykey).domain
            fig.add_annotation(
                xref="paper", yref="paper",
                x=0.0, y=dom[1],
                #text="<b>Unmethylated / Paternal</b>",
                text=" ",
                showarrow=False,
                xanchor="left", yanchor="top",
                xshift=8, yshift=-4,
                font=dict(size=12, color="#444"),
            )

        # Mini-legends for methylation tracks (right side): two sets
        if (maternal_row is not None) or (paternal_row is not None):
            try:
                # ensure enough right margin so annotations are not clipped
                try:
                    _r = fig.layout.margin.r if (fig.layout.margin and fig.layout.margin.r is not None) else 0
                    if _r < 170:
                        fig.update_layout(margin=dict(r=170))
                except Exception:
                    pass

                x_dot = 1.02
                x_txt = 1.04
                txt_col = "#2a3f5f"
                dot_size = 16  # visually closer to Plotly legend marker size

                # Methylated track legend: red = methylated, gray = maternal
                if maternal_row is not None:
                    ykey_m = "yaxis" if maternal_row == 1 else f"yaxis{maternal_row}"
                    dom_m = getattr(fig.layout, ykey_m).domain
                    y_top_m = dom_m[1]
                    span_m = dom_m[1] - dom_m[0]

                    legend_items_m = [
                        ("Methylated", HYPER_METHYL_COLOR),
                        ("Maternal", INHERITANCE_GRAY),
                    ]
                    for i, (lab, colr) in enumerate(legend_items_m):
                        y_leg = y_top_m - (0.34 + 0.23 * i) * span_m
                        fig.add_annotation(
                            xref="paper", yref="paper",
                            x=x_dot, y=y_leg,
                            text="●",
                            showarrow=False,
                            xanchor="left", yanchor="middle",
                            font=dict(size=dot_size, color=colr),
                        )
                        fig.add_annotation(
                            xref="paper", yref="paper",
                            x=x_txt, y=y_leg,
                            text=lab,
                            showarrow=False,
                            xanchor="left", yanchor="middle",
                            font=dict(size=12, color=txt_col),
                        )

                # Unmethylated track legend: blue = unmethylated, gray = paternal
                if paternal_row is not None:
                    ykey_p = "yaxis" if paternal_row == 1 else f"yaxis{paternal_row}"
                    dom_p = getattr(fig.layout, ykey_p).domain
                    y_top_p = dom_p[1]
                    span_p = dom_p[1] - dom_p[0]

                    legend_items_p = [
                        ("Unmethylated", HYPO_METHYL_COLOR),
                        ("Paternal", INHERITANCE_GRAY),
                    ]
                    span_ref = span_p
                    try:
                        if maternal_row is not None:
                            span_ref = span_m
                    except Exception:
                        span_ref = span_p
                    for i, (lab, colr) in enumerate(legend_items_p):
                        y_leg = y_top_p - (0.34 + 0.23 * i) * span_ref
                        fig.add_annotation(
                            xref="paper", yref="paper",
                            x=x_dot, y=y_leg,
                            text="●",
                            showarrow=False,
                            xanchor="left", yanchor="middle",
                            font=dict(size=dot_size, color=colr),
                        )
                        fig.add_annotation(
                            xref="paper", yref="paper",
                            x=x_txt, y=y_leg,
                            text=lab,
                            showarrow=False,
                            xanchor="left", yanchor="middle",
                            font=dict(size=12, color=txt_col),
                        )
            except Exception:
                pass

    except Exception:
        pass

    '''
    hyper_legend_added = False  # for maternal panel
    hypo_legend_added = False  # for paternal panel

    # Find max counts
    hyper_max = 0
    hypo_max = 0
    for _, row_u in user_data.iterrows():
        hyper_max = max(hyper_max, int(row_u.get("Hyper_Count", 0) or 0))
        hypo_max = max(hypo_max, int(row_u.get("Hypo_Count", 0) or 0))

    if hyper_max <= 0:
        hyper_max = 1
    if hypo_max <= 0:
        hypo_max = 1

    # Draw sticks in appropriate panel
    for _, row_u in user_data.iterrows():
        x = row_u["Jittered_Coord"]
        pos = row_u["GRCh38 coordinates"]
        hyper = int(row_u.get("Hyper_Count", 0) or 0)
        hypo = int(row_u.get("Hypo_Count", 0) or 0)

        base_hover = row_u.get("Hover_Text")
        if not isinstance(base_hover, str) or not base_hover.strip():
            base_hover = (
                "<b>User variant(s)</b><br>"
                f"<b>Genomic Position:</b> {int(pos)}"
            )

        # --- Methylated / maternal panel ---
        if hyper > 0 and maternal_row is not None:
            hover_text = (
                    base_hover
                    + f"<br><b>Methylated report count:</b> {hyper}"
            )

            # stick
            fig.add_trace(
                go.Scatter(
                    x=[x, x],
                    y=[0, hyper],
                    mode="lines",
                    line=dict(color=USER_HAPLO_COLOR, width=3),
                    hoverinfo="text",
                    text=hover_text,
                    showlegend=False,
                ),
                row=maternal_row,
                col=1,
            )

            # head
            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[hyper],
                    mode="markers",
                    marker=dict(size=8, color=USER_HAPLO_COLOR,
                                line=dict(color="white", width=1),
                                ),
                    hoverinfo="text",
                    text=hover_text,
                    name="User methylation",
                    showlegend=not hyper_legend_added,
                ),
                row=maternal_row,
                col=1,
            )
            hyper_legend_added = True

        # --- Unmethylated / paternal panel ---
        if hypo > 0 and paternal_row is not None:
            hover_text = (
                    base_hover
                    + f"<br><b>Hypomethylated report count:</b> {hypo}"
            )

            # stick
            fig.add_trace(
                go.Scatter(
                    x=[x, x],
                    y=[0, hypo],
                    mode="lines",
                    line=dict(color=USER_HAPLO_COLOR, width=3),
                    hoverinfo="text",
                    text=hover_text,
                    showlegend=False,
                ),
                row=paternal_row,
                col=1,
            )

            # head
            fig.add_trace(
                go.Scatter(
                    x=[x],
                    y=[hypo],
                    mode="markers",
                    marker=dict(size=8, color=USER_HAPLO_COLOR,
                                line=dict(color="white", width=1),
                                ),
                    hoverinfo="text",
                    text=hover_text,
                    name="User methylation",
                    showlegend=not hypo_legend_added,
                ),
                row=paternal_row,
                col=1,
            )
            hypo_legend_added = True

    # Set y-axes for the two haplotype panels
    if maternal_row is not None:
        fig.update_yaxes(
            title_text="Methylated",
            showgrid=True, zeroline=False,
            showticklabels=False,
            fixedrange=True,
            range=[0, max(hyper_max * 1.2, hyper_max + 1)],
            row=maternal_row,
            col=1,
        )
    if paternal_row is not None:
        fig.update_yaxes(
            title_text="Unmethylated",
            showgrid=True, zeroline=False,
            showticklabels=False,
            fixedrange=True,
            range=[0, max(hypo_max * 1.2, hypo_max + 1)],
            row=paternal_row,
            col=1,
        )
    '''

    # Exon/intron panel
    if rb1_exons:
        min_exon_start = min(s for s, _ in rb1_exons)
        max_exon_end = max(e for _, e in rb1_exons)

        exon_xref = f"x{exon_row}"
        exon_yref = f"y{exon_row}"

        if PROMOTER_START is not None and PROMOTER_END is not None:
            p_start = min(PROMOTER_START, PROMOTER_END)
            p_end = max(PROMOTER_START, PROMOTER_END)

            shapes.append(dict(
                type="rect",
                xref=exon_xref,
                yref=exon_yref,
                x0=p_start,
                x1=p_end,
                y0=0,
                y1=1,
                fillcolor=PROMOTER_COLOR,
                line=dict(width=0),
                layer="above",
            ))

            promoter_label_html = (
                f"<span style='color:{PROMOTER_COLOR}; font-weight:bold;'>"
                f"Promoter</span>"
            )

            fig.add_trace(
                go.Scatter(
                    x=[p_start, p_end, p_end, p_start, p_start],
                    y=[0, 0, 1, 1, 0],
                    mode="lines",
                    line=dict(width=0),
                    fill="toself",
                    hoveron="fills",
                    hoverinfo="text",
                    text=(
                            promoter_label_html
                            + f"<br><b>Genomic:</b> {p_start}–{p_end}"
                    ),
                    opacity=0,  # invisible but still catches hover
                    showlegend=False,
                ),
                row=exon_row,
                col=1,
            )

        # Exons + introns (drawn first, so DMR can sit on top)
        for i, (start, end) in enumerate(rb1_exons):
            exon_num = i + 1

            shapes.append(dict(
                type="rect",
                xref=exon_xref,
                yref=exon_yref,
                x0=start,
                x1=end,
                y0=0,
                y1=1,
                fillcolor=EXON_COLOR,
                line=dict(width=0),
                layer="above",
            ))

            exon_hover_text = EXON_HOVER_MAP.get(
                exon_num,
                (
                    f"<span style='color:{EXON_COLOR}; font-weight:bold;'>"
                    f"Exon {exon_num}</span>"
                    f"<br><b>Genomic:</b> {start}–{end}"
                ),
            )

            fig.add_trace(
                go.Scatter(
                    x=[start, end, end, start, start],
                    y=[0, 0, 1, 1, 0],
                    mode="lines",
                    line=dict(width=0),
                    fill="toself",
                    hoveron="fills",
                    hoverinfo="text",
                    text=exon_hover_text,
                    opacity=0,
                    showlegend=False,
                ),
                row=exon_row,
                col=1,
            )

            # Intron between this exon and the next one
            if i < len(rb1_exons) - 1:
                intron_start = end
                intron_end = rb1_exons[i + 1][0]

                shapes.append(dict(
                    type="rect",
                    xref=exon_xref,
                    yref=exon_yref,
                    x0=intron_start,
                    x1=intron_end,
                    y0=0,
                    y1=1,
                    fillcolor=INTRON_COLOR,
                    line=dict(width=0),
                    layer="above",
                ))

                intron_label_html = (
                    f"<span style='color:{INTRON_COLOR}; font-weight:bold;'>"
                    f"Intron {i + 1}</span>"
                )

                intron_hover_text = (
                        intron_label_html
                        + f"<br><b>Genomic:</b> {intron_start}–{intron_end}"
                )

                fig.add_trace(
                    go.Scatter(
                        x=[intron_start, intron_end, intron_end, intron_start, intron_start],
                        y=[0, 0, 1, 1, 0],
                        mode="lines",
                        line=dict(width=0),
                        fill="toself",
                        hoveron="fills",
                        hoverinfo="text",
                        text=intron_hover_text,
                        opacity=0,
                        showlegend=False,
                    ),
                    row=exon_row,
                    col=1,
                )

        # DMR region box OVERLAY
        if DMR_START is not None and DMR_END is not None:
            d_start = min(DMR_START, DMR_END)
            d_end = max(DMR_START, DMR_END)

            # Purple rectangle over exons/introns
            shapes.append(dict(
                type="rect",
                xref=exon_xref,
                yref=exon_yref,
                x0=d_start,
                x1=d_end,
                y0=0,
                y1=1,
                fillcolor=DMR_COLOR,
                line=dict(width=1, color=DMR_COLOR),
                opacity=0.4,
                layer="above",
            ))

            # Hoverable polygon for rich tooltip
            dmr_hover_text = (
                f"<span style='color:{DMR_COLOR}; font-weight:bold;'>"
                f"RB1 intronic DMR</span>"
                f"<br><b>Genomic (GRCh38):</b> {d_start}–{d_end}"
                f"<br><b>Location:</b> intron 2"
            )

            fig.add_trace(
                go.Scatter(
                    x=[d_start, d_end, d_end, d_start, d_start],
                    y=[0, 0, 1, 1, 0],
                    mode="lines",
                    line=dict(width=0),
                    fill="toself",
                    hoveron="fills",
                    hoverinfo="text",
                    text=dmr_hover_text,
                    opacity=0,  # invisible but catches hover
                    showlegend=False,
                    meta={"pinned_overlay": False},
                ),
                row=exon_row,
                col=1,
            )

            # Visible text inside the purple box (CpG85 in white)
            '''
            fig.add_annotation(
                xref=exon_xref,
                yref=exon_yref,
                x=(d_start + d_end) / 2.0,
                y=0.5,
                text=f"<b>{DMR_LABEL}</b>",
                showarrow=False,
                font=dict(size=11, color="white"),
                align="center",
                bgcolor="rgba(0,0,0,0)",
                bordercolor="rgba(0,0,0,0)",
            )
            '''

    if protein_boxes_row is not None and protein_lollipop_row is not None and CDS_EXON_AA_RANGES:
        PROT_EXON_BG = "#FFFFFF"
        PROT_EXON_LINE = "#111827"
        PROT_DOMAIN_OPACITY = 0.85

        PROT_BAR_Y0, PROT_BAR_Y1 = 0.10, 0.90
        PROT_STICK_Y0 = 1.06
        PROT_MIN_HEAD_Y = 1.14
        PROT_MAX_HEAD_Y = 1.70  #

        # gnomAD-like exon tiles
        EXON_PAD = 0.0
        CODON_JITTER = 0.00

        # Build exon boxes in true protein AA coordinates (widths reflect real CDS exon lengths)
        ex_vis = []
        for i, ex in enumerate(CDS_EXON_AA_RANGES):
            aa0 = int(ex.get("aa_start", 0))
            aa1 = int(ex.get("aa_end", 0))
            if aa1 < aa0:
                continue

            # Use half-AA boundaries so adjacent exons touch with no visual gaps
            vx0 = aa0 - 0.5
            vx1 = aa1 + 0.5

            ex_vis.append({
                "aa0": aa0, "aa1": aa1,
                "vx0": vx0, "vx1": vx1,
                "idx": i,
                "exon_display": ex.get("exon_display", i + 1),
                "ex": ex
            })

            # Draw exon fill (no border)
            shapes.append(dict(
                type="rect",
                x0=vx0, x1=vx1,
                y0=PROT_BAR_Y0, y1=PROT_BAR_Y1,
                xref=f"x{protein_boxes_row}", yref=f"y{protein_boxes_row}",
                fillcolor=PROT_EXON_BG,
                opacity=1.0,
                line=dict(width=0),
                layer="below",
            ))

        def aa_to_vx(aa: int):
            if aa is None:
                return None
            return float(aa)

        if uniprot_dom_df is not None and len(uniprot_dom_df) > 0 and ex_vis:
            for _, dom in uniprot_dom_df.iterrows():
                try:
                    ds = int(dom.get("start", 0))
                    de = int(dom.get("end", 0))
                except Exception:
                    continue
                if de <= 0:
                    continue

                name = str(dom.get("description", "")).strip()
                col = DOMAIN_COLOR_MAP.get(name, "#FDE68A")  # fallback soft yellow

                for ex in ex_vis:
                    seg0 = max(ds, ex["aa0"])
                    seg1 = min(de, ex["aa1"])
                    if seg1 < seg0:
                        continue

                    x0 = aa_to_vx(seg0) - 0.5
                    x1 = aa_to_vx(seg1) + 0.5
                    if x0 is None or x1 is None:
                        continue
                    if x1 < x0:
                        x0, x1 = x1, x0

                    shapes.append(dict(
                        type="rect",
                        x0=x0, x1=x1,
                        y0=PROT_BAR_Y0, y1=PROT_BAR_Y1,
                        xref=f"x{protein_boxes_row}", yref=f"y{protein_boxes_row}",
                        fillcolor=col,
                        opacity=PROT_DOMAIN_OPACITY,
                        line=dict(width=0),
                        layer="below",
                    ))

        if ex_vis:
            full_x0 = ex_vis[0]["vx0"]
            full_x1 = ex_vis[-1]["vx1"]

            # Outer border
            shapes.append(dict(
                type="rect",
                x0=full_x0, x1=full_x1,
                y0=PROT_BAR_Y0, y1=PROT_BAR_Y1,
                xref=f"x{protein_boxes_row}", yref=f"y{protein_boxes_row}",
                fillcolor="rgba(0,0,0,0)",
                line=dict(color=PROT_EXON_LINE, width=1),
                layer="above",
            ))

            for i in range(len(ex_vis) - 1):
                boundary_x = ex_vis[i]["vx1"]
                shapes.append(dict(
                    type="line",
                    x0=boundary_x, x1=boundary_x,
                    y0=PROT_BAR_Y0, y1=PROT_BAR_Y1,
                    xref=f"x{protein_boxes_row}", yref=f"y{protein_boxes_row}",
                    line=dict(color=PROT_EXON_LINE, width=1),
                    layer="above",
                ))

        exon_stats = {ex["exon_display"]: {"count": 0, "sum_ac": 0.0, "examples": []} for ex in ex_vis}
        aa_groups = {}
        gpos_to_gstr = {}

        def _clean_label(v):
            if v is None:
                return None
            if isinstance(v, float) and math.isnan(v):
                return None
            s = str(v).strip()
            if not s:
                return None
            if s.lower() in {"nan", "none", "n/a", "na"}:
                return None
            return s

        def _first_clean(v):
            if isinstance(v, list):
                for it in v:
                    s = _clean_label(it)
                    if s:
                        return s
                return None
            return _clean_label(v)

        _cols_meta = [c for c in (
            '_row_id', 'GRCh38 coordinates', 'Display_Source',
            'cDNA_Label_raw', 'Protein_Label', 'Mutation_Label', 'gDNA_Label_raw',
            'Total_AC', 'Frequency'
        ) if c in data.columns]
        _records = data[_cols_meta].to_dict('records') if _cols_meta else data.to_dict('records')

        for rec in _records:
            gpos = rec.get("GRCh38 coordinates", None)
            if gpos is None:
                continue
            try:
                gpos_i = int(gpos)
            except Exception:
                continue

            if gpos_i not in gpos_to_gstr:
                try:
                    _gs = format_genomic_change_for_row(rec)
                    if _gs:
                        gpos_to_gstr[gpos_i] = _gs
                except Exception:
                    pass

            src = rec.get("Display_Source", "Curated")

            fallback_aa = None
            mapped = POS_TO_AA_OFFSET.get(gpos_i)
            if mapped:
                try:
                    fallback_aa = int(mapped[0])
                except Exception:
                    fallback_aa = None

            variant_meta = _get_variant_meta_cached(rec)

            if not variant_meta:
                variant_meta = [
                    {
                        "genomic": gpos_to_gstr.get(gpos_i)
                                   or (format_genomic_change_for_row(rec) or f"{CHROMOSOME}:{gpos_i}"),
                        "coding": _first_clean(rec.get("cDNA_Label_raw", None)),
                        "protein": _first_clean(rec.get("Protein_Label", None)),
                        "total_count": rec.get("Total_AC", rec.get("Frequency", 1)),
                    }
                ]

            for vm in variant_meta:
                # Per-variant count
                ac = vm.get("total_count", vm.get("count", 1))
                try:
                    ac_f = float(ac)
                    if not math.isfinite(ac_f):
                        ac_f = 1.0
                except Exception:
                    ac_f = 1.0
                ac_f = max(ac_f, 1.0)

                cdna = _clean_label(vm.get("coding")) or _first_clean(rec.get("cDNA_Label_raw", None))
                g_lbl = _clean_label(vm.get("genomic"))

                prot_raw = vm.get("protein", None)
                prot_pairs = extract_protein_aa_pairs(prot_raw)

                if not prot_pairs:
                    if fallback_aa is None:
                        continue
                    prot_pairs = [(fallback_aa, prot_raw)]

                # Deduplicate (AA, label) pairs
                seen = set()
                prot_pairs_uniq = []
                for aa_pos_i, prot_lbl in prot_pairs:
                    try:
                        aa_pos_i = int(aa_pos_i)
                    except Exception:
                        continue
                    key = (aa_pos_i, str(prot_lbl) if prot_lbl is not None else "")
                    if key in seen:
                        continue
                    seen.add(key)
                    prot_pairs_uniq.append((aa_pos_i, prot_lbl))

                if not prot_pairs_uniq:
                    continue

                for aa_pos_i, prot_lbl in prot_pairs_uniq:
                    vx = aa_to_vx(aa_pos_i)
                    if vx is None:
                        continue

                    # Protein label for hover/table
                    prot = _clean_label(prot_lbl)

                    # Exon display: infer from AA ranges (fast)
                    exon_disp = None
                    for ex in ex_vis:
                        if ex["aa0"] <= aa_pos_i <= ex["aa1"]:
                            exon_disp = ex["exon_display"]
                            break

                    # Domains at this AA position
                    doms = _domains_for_aa(aa_pos_i)

                    # update exon stats (for exon hover)
                    if exon_disp is not None and exon_disp in exon_stats:
                        exon_stats[exon_disp]["count"] += 1
                        exon_stats[exon_disp]["sum_ac"] += ac_f
                        if len(exon_stats[exon_disp]["examples"]) < 10:
                            exon_stats[exon_disp]["examples"].append(
                                (ac_f, gpos_i, aa_pos_i, cdna, prot)
                            )

                    # group by AA
                    entry = aa_groups.get(aa_pos_i)
                    if entry is None:
                        entry = {
                            "x": vx,
                            "sum_ac": 0.0,
                            "src_ac": {},
                            "gmin": gpos_i,
                            "gmax": gpos_i,
                            "exon": exon_disp,
                            "domains": doms,
                            "examples": [],
                        }
                        aa_groups[aa_pos_i] = entry

                    entry["sum_ac"] += ac_f
                    entry["src_ac"][src] = entry["src_ac"].get(src, 0.0) + ac_f
                    entry["gmin"] = min(entry["gmin"], gpos_i)
                    entry["gmax"] = max(entry["gmax"], gpos_i)

                    if len(entry["examples"]) < 12:
                        # Include a per-variant genomic label so AA hovers don't "borrow" the wrong variant
                        entry["examples"].append((ac_f, gpos_i, g_lbl, cdna, prot, src))

        # Add hoverable exon tiles
        for ex in ex_vis:
            exon_disp = ex["exon_display"]
            g0 = ex["ex"]["g_start"]
            g1 = ex["ex"]["g_end"]
            aa0 = ex["aa0"]
            aa1 = ex["aa1"]

            # domains overlapping this exon (by AA overlap)
            exon_doms = _domains_for_range(aa0, aa1)

            hover_lines = [
                f"<b>Exon</b> {exon_disp}",
                f"<b>AA range:</b> {aa0}–{aa1}",
                # f"<b>Genomic range (GRCh38):</b> {int(g0):,}–{int(g1):,}",
            ]

            if exon_doms:
                show = exon_doms[:6]
                dom_txt = ", ".join(show)
                if len(exon_doms) > 6:
                    dom_txt += f", … (+{len(exon_doms) - 6} more)"
                hover_lines.append(f"<b>Protein domains:</b> {dom_txt}")

            st = exon_stats.get(exon_disp, None)
            if st and st["count"] > 0:
                # hover_lines.append(
                # f"<b>Variants (current filters):</b> {st['count']}  |  sum count ≈ {st['sum_ac']:.0f}")
                '''
                ex_rows = sorted(st["examples"], key=lambda t: t[0], reverse=True)[:6]

                headers = ["Genomic", "AA", "cDNA", "Protein", "Count"]
                rows = []
                for ac_f, gpos_i, aa_pos_i, cdna, prot in ex_rows:
                    # Use the same "Genomic" string as the variant hover
                    g_str = None
                    # Prefer the same genomic string used in variant hovers
                    try:
                        gpos_key = int(gpos_i)
                    except Exception:
                        gpos_key = gpos_i
                    if isinstance(per_variant_map, dict):
                        _pv = per_variant_map.get(gpos_key) or per_variant_map.get(gpos_i)
                        if isinstance(_pv, list) and _pv:
                            _cdna_norm = str(cdna or "").strip()
                            _prot_norm = str(prot or "").strip()
                            _best = None
                            for _item in _pv:
                                if not isinstance(_item, dict):
                                    continue
                                _cd = str(_item.get("cdna") or "").strip()
                                _pr = str(_item.get("prot") or "").strip()
                                if _cdna_norm and _cd == _cdna_norm and _prot_norm and _pr == _prot_norm:
                                    _best = _item
                                    break
                                if _cdna_norm and _cd == _cdna_norm:
                                    _best = _item
                            if _best is None:
                                for _item in _pv:
                                    if isinstance(_item, dict):
                                        _best = _item
                                        break
                            if _best:
                                g_str = _best.get("g_str") or _best.get("genomic") or _best.get("g")
                        elif isinstance(_pv, dict):
                            g_str = _pv.get("g_str") or _pv.get("genomic") or _pv.get("g")
                    if not g_str:
                        g_str = gpos_to_gstr.get(gpos_key) or gpos_to_gstr.get(gpos_i)
                    if not g_str:
                        g_str = f"{CHROMOSOME}:{gpos_key}"

                    prot = prot or "-"

                    rows.append((str(g_str), str(aa_pos_i), str(cdna), str(prot), str(int(ac_f))))

                # Compute per-column widths for alignment
                widths = []
                for i, h in enumerate(headers):
                    max_len = len(h)
                    for r in rows:
                        if len(r[i]) > max_len:
                            max_len = len(r[i])
                    widths.append(max_len)

                fmt = (
                    f"{{:<{widths[0]}}}  "
                    f"{{:>{widths[1]}}}  "
                    f"{{:<{widths[2]}}}  "
                    f"{{:<{widths[3]}}}  "
                    f"{{:>{widths[4]}}}"
                )

                table_lines = [fmt.format(*headers), fmt.format(*["─" * w for w in widths])]
                for r in rows:
                    table_lines.append(fmt.format(*r))

                hover_lines.append(
                    "<span style='font-family:monospace; font-size:11px; white-space:pre; line-height:1.1;'>"
                    + "<br>".join(html_lib.escape(s) for s in table_lines)
                    + "</span>"
                )

                if st["count"] > len(ex_rows):
                    hover_lines.append(
                        f"<span style='font-size:11px;'>… and {st['count'] - len(ex_rows)} more in this exon</span>")
                '''
            else:
                hover_lines.append("<b>Variants (current filters):</b> 0")

            fig.add_trace(
                go.Scatter(
                    x=[ex["vx0"], ex["vx1"], ex["vx1"], ex["vx0"], ex["vx0"]],
                    y=[PROT_BAR_Y0, PROT_BAR_Y0, PROT_BAR_Y1, PROT_BAR_Y1, PROT_BAR_Y0],
                    fill="toself",
                    mode="lines",
                    line=dict(width=0),
                    fillcolor="rgba(0,0,0,0)",
                    hoveron="fills",
                    text="<br>".join(hover_lines),
                    hoverinfo="text",
                    showlegend=False,
                ),
                row=protein_boxes_row, col=1
            )

        pts = []
        for aa_pos_i, entry in aa_groups.items():
            # dominant source (for color)
            if entry["src_ac"]:
                best_src = max(entry["src_ac"].items(), key=lambda kv: kv[1])[0]
            else:
                best_src = "Curated"
            col = FIXED_SOURCE_COLORS.get(best_src, "rgba(17,24,39,0.45)")

            hover_lines = [f"<b>Protein position (AA):</b> {aa_pos_i}"]
            if entry.get("exon") is not None:
                hover_lines.append(f"<b>Exon:</b> {entry['exon']}")
            if entry.get("domains"):
                show = entry["domains"][:6]
                dom_txt = ", ".join(show)
                if len(entry["domains"]) > 6:
                    dom_txt += f", … (+{len(entry['domains']) - 6} more)"
                hover_lines.append(f"<b>Protein domain:</b> {dom_txt}")

            # hover_lines.append(
            # f"<b>Genomic positions (GRCh38):</b> {int(entry['gmin']):,}–{int(entry['gmax']):,}"
            # )
            hover_lines.append(f"<b>Report Count :</b> {entry['sum_ac']:.0f}")

            # Examples table (aligned)
            ex_rows = sorted(entry["examples"], key=lambda t: t[0], reverse=True)[:8]
            if ex_rows:
                table_rows = []
                for ac_f, gpos_i, g_lbl, cdna, prot, src in ex_rows:
                    cdna = cdna or "-"
                    prot = prot or "-"
                    src = src or "-"
                    g_str = g_lbl

                    # If we don't already have a per-variant genomic label, try to match it from per_variant_map
                    if not g_str:
                        try:
                            gpos_key = int(gpos_i)
                        except Exception:
                            gpos_key = gpos_i

                        if isinstance(per_variant_map, dict):
                            _pv = per_variant_map.get(gpos_key) or per_variant_map.get(gpos_i)
                            if isinstance(_pv, list) and _pv:
                                _cdna_norm = str(cdna or "").strip()
                                _prot_norm = str(prot or "").strip()
                                _best = None
                                for _item in _pv:
                                    if not isinstance(_item, dict):
                                        continue
                                    _cd = str(_item.get("cdna") or "").strip()
                                    _pr = str(_item.get("prot") or "").strip()
                                    if _cdna_norm and _cd == _cdna_norm and _prot_norm and _pr == _prot_norm:
                                        _best = _item
                                        break
                                    if _cdna_norm and _cd == _cdna_norm:
                                        _best = _item
                                if _best is None:
                                    for _item in _pv:
                                        if isinstance(_item, dict):
                                            _best = _item
                                            break
                                if _best:
                                    g_str = _best.get("g_str") or _best.get("genomic") or _best.get("g")

                    if not g_str:
                        try:
                            gpos_key = int(gpos_i)
                        except Exception:
                            gpos_key = gpos_i
                        g_str = gpos_to_gstr.get(gpos_key) or gpos_to_gstr.get(gpos_i) or f"{CHROMOSOME}:{gpos_key}"

                    table_rows.append((g_str, cdna, prot, src, int(ac_f)))

                # Hover-safe monospace "table"
                col_headers = ["Genomic", "cDNA", "Protein", "Source", "Count"]

                def _trim(_s, _n):
                    _s = str(_s)
                    if _n is None or _n <= 0:
                        return _s
                    return (_s[:_n - 1] + "…") if (len(_s) > _n and _n >= 2) else _s

                # Trim long fields to keep hover boxes compact
                _rows = [(_trim(g, 34), _trim(c, 20), _trim(p, 20), _trim(s, 18), str(cnt)) for (g, c, p, s, cnt) in
                         table_rows]

                w = [len(h) for h in col_headers]
                for r in _rows:
                    for i, v in enumerate(r):
                        w[i] = max(w[i], len(str(v)))

                hdr = "{:<{}}  {:<{}}  {:<{}}  {:<{}}  {:>{}}".format(col_headers[0], w[0], col_headers[1], w[1],
                                                                      col_headers[2], w[2], col_headers[3], w[3],
                                                                      col_headers[4], w[4])
                sep = "─" * len(hdr)
                lines = [hdr, sep]
                for g, c, pr, s, cnt in _rows:
                    lines.append(
                        "{:<{}}  {:<{}}  {:<{}}  {:<{}}  {:>{}}".format(g, w[0], c, w[1], pr, w[2], s, w[3], cnt, w[4]))

                hover_lines.append(
                    "<span style='font-family:monospace; font-size:11px; white-space:pre; line-height:1.15;'>"
                    + "<br>".join(html_lib.escape(x) for x in lines)
                    + "</span>"
                )

            pts.append({"x": entry["x"], "ac": float(entry["sum_ac"]), "color": col, "hover": "<br>".join(hover_lines)})

        # Max y for protein lollipop (used to anchor y=0 at the x-axis baseline)
        protein_max_y = 1.0
        if pts:
            try:
                protein_max_y = max(protein_max_y, max(float(p.get("ac", 0.0) or 0.0) for p in pts))
            except Exception:
                protein_max_y = 1.0

        MAX_PROT_AA_POINTS = 5000
        if len(pts) > MAX_PROT_AA_POINTS:
            pts.sort(key=lambda d: d["ac"], reverse=True)
            pts = pts[:MAX_PROT_AA_POINTS]

        x_line, y_line = [], []
        x_head, y_head, c_head, h_head = [], [], [], []

        for p in pts:
            y = max(float(p["ac"]), 0.0)
            x = p["x"]
            x_line.extend([x, x, None])
            y_line.extend([0.0, y, None])
            x_head.append(x)
            y_head.append(y)
            c_head.append(p["color"])
            h_head.append(p["hover"])

        # sticks (uniform line color for readability)
        fig.add_trace(
            go.Scatter(
                x=x_line,
                y=y_line,
                mode="lines",
                line=dict(color="rgba(17,24,39,0.35)", width=1),
                hoverinfo="skip",
                showlegend=False,
            ),
            row=protein_lollipop_row, col=1
        )

        # heads (colored by dominant source)
        fig.add_trace(
            go.Scatter(
                x=x_head,
                y=y_head,
                mode="markers",
                marker=dict(size=7, color=c_head, line=dict(width=0)),
                hoverinfo="text",
                text=h_head,
                showlegend=False,
            ),
            row=protein_lollipop_row, col=1
        )

        # x-axis: show AA ticks (true AA coordinate) so the axis reads like protein position
        if ex_vis:
            aa_max = max(ex["aa1"] for ex in ex_vis)
            # Store exon-to-AA mapping in layout.meta so clientside can densify AA ticks on zoom
            _meta = dict(fig.layout.meta) if getattr(fig.layout, 'meta', None) else {}
            protein_lollipop_yaxis_key = "yaxis" if protein_lollipop_row == 1 else f"yaxis{protein_lollipop_row}"
            protein_boxes_yaxis_key = "yaxis" if protein_boxes_row == 1 else f"yaxis{protein_boxes_row}"
            protein_lollipop_xaxis_key = "xaxis" if protein_lollipop_row == 1 else f"xaxis{protein_lollipop_row}"
            protein_boxes_xaxis_key = "xaxis" if protein_boxes_row == 1 else f"xaxis{protein_boxes_row}"
            _meta.update({
                "protein_ex_vis": ex_vis,
                "protein_aa_max": int(aa_max),
                "protein_lollipop_yaxis_key": protein_lollipop_yaxis_key,
                "protein_boxes_yaxis_key": protein_boxes_yaxis_key,
                "protein_lollipop_xaxis_key": protein_lollipop_xaxis_key,
                "protein_boxes_xaxis_key": protein_boxes_xaxis_key
            })
            fig.update_layout(meta=_meta)
            aa_ticks = [1, 200, 400, 600, 800, aa_max]
            aa_ticks = [a for a in aa_ticks if 1 <= a <= aa_max]
            tickvals, ticktext = [], []
            for a in aa_ticks:
                vx = aa_to_vx(int(a))
                if vx is None:
                    continue
                tickvals.append(vx)
                ticktext.append(str(int(a)))
            xmax = ex_vis[-1]["vx1"]
            fig.update_xaxes(
                row=protein_boxes_row, col=1,
                range=[-EXON_PAD, xmax + EXON_PAD],
                tickmode="array",
                tickvals=tickvals,
                ticktext=ticktext,
                title_text="Protein position (AA)",
                title_standoff=10,
                automargin=True,
                zeroline=False,
                showgrid=False,
                showline=True,
                linecolor="#444",
                linewidth=1,
                ticks="outside",
                tickcolor="#444",
                tickwidth=1,
                ticklen=4,
            )
            # Keep the protein lollipop row perfectly aligned with the protein boxes row
            fig.update_xaxes(
                row=protein_lollipop_row, col=1,
                matches=f"x{protein_boxes_row}",
                showticklabels=False,
                zeroline=False,
                showgrid=False,
                showline=True,
                linecolor="#444",
                linewidth=1,
                ticks="outside",
                tickcolor="#444",
                tickwidth=1,
                ticklen=4,
            )
    fig.update_yaxes(
        title_text="RB1 Gene<br>(NM_000321.3)",
        showticklabels=False,
        showgrid=False,
        fixedrange=True,
        zeroline=False,
        range=[0, 1],
        row=exon_row,
        col=1,
    )

    # Protein lollipop row (AA coordinates): report counts on Y
    if protein_lollipop_row is not None:
        fig.update_yaxes(
            title_text="Report<br>Count",
            showgrid=True,
            fixedrange=True,
            zeroline=False,
            range=[0, protein_max_y * 1.1],

            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="",
            ticklen=0,
            ticklabelstandoff=0,
            row=protein_lollipop_row,
            col=1,
        )

    # Protein boxes row (AA coordinates)
    if protein_boxes_row is not None:
        fig.update_yaxes(
            title_text="prB protein<br>(NP_000321.2)",
            showticklabels=False,
            showgrid=False,
            fixedrange=True,
            range=[0.0, 1.0],
            row=protein_boxes_row,
            col=1,
        )

    fig.update_yaxes(
        title_text="Report<br>Count",
        showgrid=True,
        zeroline=False,
        fixedrange=True,
        range=[0, max_y * 1.1],
        row=1,
        col=1,
    )
    # Main variants row
    fig.update_xaxes(
        row=1,
        col=1,
        range=[GENOMIC_MIN, GENOMIC_MAX],
        showticklabels=False,
        tickformat=GENOMIC_TICKFORMAT,
        separatethousands=True,
        exponentformat="none",
        showgrid=False,
        showline=True,
        linecolor="#444",
        linewidth=1,
        ticks="outside",
        tickcolor="#444",
        tickwidth=1,
        ticklen=4,
    )
    '''
    # Exon / intron row: share x-axis but no ticks or title
    if exon_row is not None:
        fig.update_xaxes(
            row=exon_row,
            col=1,
            matches="x",
            showticklabels=False,
            tickformat=GENOMIC_TICKFORMAT,
            separatethousands=True,
            exponentformat="none",
            showgrid=False,
            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="outside",
            tickcolor="#444",
            tickwidth=1,
            ticklen=4,
        )
    '''

    if exon_row is not None:
        fig.update_xaxes(
            row=exon_row,
            col=1,
            matches="x",
            showticklabels=True,
            # title_text="Genomic position (GRCh38)",
            title_standoff=2,
            automargin=True,
            tickformat=GENOMIC_TICKFORMAT,
            separatethousands=True,
            exponentformat="none",
            showgrid=False,
            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="outside",
            tickcolor="#444",
            tickwidth=1,
            ticklen=4,
        )

    # Protein rows: AA axis (independent of genomic panels)
    if protein_boxes_row is not None:
        fig.update_xaxes(
            row=protein_boxes_row,
            col=1,
            showticklabels=True,
            tickformat="d",
            showgrid=False,
            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="outside",
            tickcolor="#444",
            tickwidth=1,
            ticklen=4,
        )
    if protein_lollipop_row is not None and protein_boxes_row is not None:
        fig.update_xaxes(
            row=protein_lollipop_row,
            col=1,
            matches=f"x{protein_boxes_row}",
            showticklabels=False,
            showgrid=False,
            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="outside",
            tickcolor="#444",
            tickwidth=1,
            ticklen=4,
        )

    # Maternal methylation row
    if maternal_row is not None:
        fig.update_xaxes(
            row=maternal_row,
            col=1,
            matches="x",
            showticklabels=False,
            tickformat=GENOMIC_TICKFORMAT,
            separatethousands=True,
            exponentformat="none",
            showgrid=False,
            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="outside",
            tickcolor="#444",
            tickwidth=1,
            ticklen=4,
        )

    # Paternal (bottom) methylation row
    if paternal_row is not None:
        fig.update_xaxes(
            row=paternal_row,
            col=1,
            matches="x",
            showticklabels=True,
            title_text="Genomic position (GRCh38)",
            tickformat=GENOMIC_TICKFORMAT,
            separatethousands=True,
            exponentformat="none",
            showgrid=False,
            showline=True,
            linecolor="#444",
            linewidth=1,
            ticks="outside",
            tickcolor="#444",
            tickwidth=1,
            ticklen=4,
        )

    nav_text = (
        "To navigate the viewer, click and drag to zoom in. "
        "Double click to zoom out. Hold Shift and drag to move along the gene. "
        "Hover over the lollipop to view variant details."
    )

    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=0.5,
        y=-0.18,
        text=nav_text,
        showarrow=False,
        font=dict(size=12, color="#444"),
        align="center",
    )

    layout_kwargs = dict(
        height=1050,
        dragmode="zoom",
        xaxis=dict(fixedrange=False),
        yaxis=dict(fixedrange=True),
        legend_title_text="Source",
        uirevision="rb1_combined",
        margin=dict(t=80, r=40, b=120, l=150),
    )
    layout_kwargs.update(
        xaxis2=dict(fixedrange=False),
        yaxis2=dict(fixedrange=True),
        xaxis3=dict(fixedrange=False),
        yaxis3=dict(fixedrange=True),
        xaxis4=dict(fixedrange=False),
        yaxis4=dict(fixedrange=True),
        xaxis5=dict(fixedrange=False),
        yaxis5=dict(fixedrange=True),
        xaxis6=dict(fixedrange=False),
        yaxis6=dict(fixedrange=True),
        xaxis7=dict(fixedrange=False),
        yaxis7=dict(fixedrange=True),
    )

    fig.update_layout(**layout_kwargs)
    fig["layout"]["yaxis"]["fixedrange"] = True

    try:
        # Default gap between rows (matches the make_subplots vertical_spacing)
        _default_gap = 0.06
        _gap_overrides = {(1, 2): 0.03, (3, 4): 0.03}
        _gaps = [
            _gap_overrides.get((i, i + 1), _default_gap)
            for i in range(1, total_rows)
        ]
        _available_h = 1.0 - sum(_gaps)

        _rh_sum = float(sum(row_heights)) if row_heights else 1.0
        _rh_norm = [(float(h) / _rh_sum) for h in row_heights]

        # Compute domains from top to bottom
        _domains = []
        _y_top = 1.0
        for i in range(total_rows):
            _h = _available_h * _rh_norm[i]
            _y_bottom = _y_top - _h
            _domains.append([_y_bottom, _y_top])
            if i < total_rows - 1:
                _y_top = _y_bottom - _gaps[i]

        for r in range(1, total_rows + 1):
            _key = "yaxis" if r == 1 else f"yaxis{r}"
            if _key in fig.layout:
                fig.layout[_key].domain = _domains[r - 1]
    except Exception as _e:
        print(f"[Row spacing] Could not apply custom domains: {_e}")

    try:
        import re
        for ax_key in list(fig.layout):
            if not str(ax_key).startswith("yaxis"):
                continue
            yax = fig.layout[ax_key]
            title_txt = None
            try:
                title_txt = yax.title.text
            except Exception:
                title_txt = None
            if title_txt:
                dom = getattr(yax, "domain", None)
                if dom and len(dom) == 2:
                    y_center = (dom[0] + dom[1]) / 2
                    # Clear the original axis title
                    try:
                        yax.title.text = ""
                    except Exception:
                        pass

                    xshift = -22

                    norm = str(title_txt)
                    norm = re.sub(r"<br\s*/?>", " ", norm, flags=re.IGNORECASE)
                    norm = re.sub(r"\s+", " ", norm).strip().lower()
                    if norm == "report count":
                        xshift = -50

                    fig.add_annotation(
                        xref="paper",
                        yref="paper",
                        x=0,
                        y=y_center,
                        xanchor="right",
                        yanchor="middle",
                        xshift=xshift,
                        text=str(title_txt),
                        showarrow=False,
                        align="center",
                        font=dict(size=12, color="#444"),
                    )
    except Exception as e:
        print(f"[Y-axis titles] Failed to make horizontal: {e}")

    # Section header above the main lollipop plot
    try:

        dom = fig.layout.yaxis.domain
        y_header = dom[1] + 0.001

        fig.add_annotation(
            xref="paper",
            yref="paper",
            x=0.0,
            y=y_header,
            xanchor="left",
            yanchor="bottom",
            text="<b>Gene view</b>",
            showarrow=False,
            font=dict(size=14),
        )
    except Exception as e:
        print(f"[Gene header] Could not add section header: {e}")

    #  Parent of Origin header
    try:
        if maternal_row is not None:
            _yaxis_key = "yaxis" if maternal_row == 1 else f"yaxis{maternal_row}"
            _dom = fig.layout[_yaxis_key].domain
            _y_header = min(0.995, _dom[1] + 0.001)

            fig.add_annotation(
                xref="paper",
                yref="paper",
                x=0.0,
                y=_y_header,
                xanchor="left",
                yanchor="bottom",
                text="<b>Parent of Origin</b>",
                showarrow=False,
                font=dict(size=14),
            )
    except Exception as e:
        print(f"[Parent of Origin header] Could not add section header: {e}")

    try:
        # Keep the spacer row axis hidden
        fig.update_xaxes(visible=False, showgrid=False, zeroline=False, row=protein_header_row, col=1)
        fig.update_yaxes(visible=False, showgrid=False, zeroline=False, row=protein_header_row, col=1)

        yaxis_key = "yaxis" if protein_lollipop_row == 1 else f"yaxis{protein_lollipop_row}"
        dom = fig.layout[yaxis_key].domain  # [bottom, top] in paper coordinates
        y_header = min(0.995, dom[1] + 0.01)  # slightly above the protein lollipop row

        fig.add_annotation(
            xref="paper",
            yref="paper",
            x=0.0,
            y=y_header,
            xanchor="left",
            yanchor="bottom",
            text="<b>Protein View</b>",
            showarrow=False,
            font=dict(size=14),
        )
    except Exception as e:
        print(f"[Protein header] Could not add section header: {e}")

    # Apply shapes in one shot (avoids O(n^2) overhead inside add_shape)
    if shapes:
        fig.update_layout(shapes=shapes)

    return fig.to_dict(), curated_store_data


app.clientside_callback(
    r"""
    function(baseFig, pinnedData, curatedData, relayoutData) {
        // If we don't have a base figure yet, just return empty.
        if (!baseFig) {
            return {};
        }

        // Deep copy so we don't mutate the stored base figure.
        var fig = JSON.parse(JSON.stringify(baseFig));

        // --- Dynamic AA ticks on protein axis (densify when you zoom) ---
        try {
            if (fig.layout && fig.layout.meta && fig.layout.meta.protein_ex_vis) {
                var exVis = fig.layout.meta.protein_ex_vis;
                var aaMaxAll = fig.layout.meta.protein_aa_max || null;

                // Find the protein x-axis (the one titled "Protein position (AA)")
                var protAxisKey = null;
                Object.keys(fig.layout).forEach(function(k) {
                    if (k.indexOf("xaxis") === 0) {
                        var ax = fig.layout[k];
                        if (ax && ax.title && ax.title.text && ax.title.text.indexOf("Protein position") !== -1) {
                            protAxisKey = k;
                        }
                    }
                });

                                function vxToAa(vx) {
                    if (vx === null || vx === undefined) return null;
                    // Protein track x-axis is in true AA coordinates; round to nearest residue index.
                    return Math.round(vx);
                }

                                function aaToVx(aa) {
                    if (aa === null || aa === undefined) return null;
                    // Protein track x-axis is in true AA coordinates.
                    return aa;
                }

                function niceStep(span, targetTicks) {
                    if (!isFinite(span) || span <= 0) return 1;
                    var rough = span / targetTicks;
                    var pow10 = Math.pow(10, Math.floor(Math.log10(rough)));
                    var base = rough / pow10;
                    var niceBase = (base <= 1) ? 1 : (base <= 2) ? 2 : (base <= 5) ? 5 : 10;
                    return niceBase * pow10;
                }

                if (protAxisKey) {
                    // Prefer the current zoom window from relayoutData if it refers to the protein axis
                    var relBase = (protAxisKey === "xaxis") ? "xaxis" : protAxisKey;
                    var vx0 = null, vx1 = null;

                    if (relayoutData) {
                        var k0 = relBase + ".range[0]";
                        var k1 = relBase + ".range[1]";
                        if (relayoutData.hasOwnProperty(k0) && relayoutData.hasOwnProperty(k1)) {
                            vx0 = relayoutData[k0];
                            vx1 = relayoutData[k1];
                        }
                    }

                    // Fallback to axis range if relayoutData is for a different panel
                    if (vx0 === null || vx1 === null) {
                        var axRange = fig.layout[protAxisKey] && fig.layout[protAxisKey].range;
                        if (axRange && axRange.length === 2) {
                            vx0 = axRange[0];
                            vx1 = axRange[1];
                        }
                    }

                    if (vx0 !== null && vx1 !== null) {
                        var vmin = Math.min(vx0, vx1);
                        var vmax = Math.max(vx0, vx1);

                        var aaMin = vxToAa(vmin);
                        var aaMax = vxToAa(vmax);
                        if (aaMin !== null && aaMax !== null) {
                            aaMin = Math.max(1, Math.floor(aaMin));
                            aaMax = Math.ceil(aaMax);
                            if (aaMaxAll) aaMax = Math.min(aaMaxAll, aaMax);

                            var span = Math.max(1, aaMax - aaMin);
                            var step = Math.max(1, niceStep(span, 7));

                            var start = Math.ceil(aaMin / step) * step;
                            var tickvals = [];
                            var ticktext = [];
                            for (var a = start; a <= aaMax; a += step) {
                                var vx = aaToVx(a);
                                if (vx === null) continue;
                                tickvals.push(vx);
                                ticktext.push(String(a));
                                if (tickvals.length > 30) break; // safety
                            }

                            if (tickvals.length >= 2) {
                                fig.layout[protAxisKey].tickmode = "array";
                                fig.layout[protAxisKey].tickvals = tickvals;
                                fig.layout[protAxisKey].ticktext = ticktext;
                            }
                        }
                    }
                }
            }
        } catch (e) {
            // If anything goes wrong, just keep the default ticks.
        }

        // --- Keep protein Y-axis ticks tight (Plotly can re-enable tick marks / padding after reset) ---
        function _tightenAxisLayout(axisKey) {
            if (!axisKey || !fig.layout || !fig.layout[axisKey]) return;
            fig.layout[axisKey].ticks = "";
            fig.layout[axisKey].ticklen = 0;
            fig.layout[axisKey].ticklabelstandoff = 0;
            fig.layout[axisKey].showline = true;
        }

        try {
            if (fig.layout && fig.layout.meta) {
                _tightenAxisLayout(fig.layout.meta.protein_lollipop_yaxis_key);
                _tightenAxisLayout(fig.layout.meta.protein_boxes_yaxis_key);
            }
        } catch (e) {
            // ignore
        }

        // Enforce on the live Plotly graph too (uirevision can preserve a "bad" axis state after zoom/reset)
        try {
            var host = document.getElementById("mutation-plot");
            var gd = host && (host.querySelector(".js-plotly-plot") || host);
            if (gd && window.Plotly && gd.layout) {
                function _needsTight(axisKey) {
                    if (!axisKey || !gd.layout[axisKey]) return false;
                    var ax = gd.layout[axisKey];
                    var t = (ax.ticks === undefined) ? null : ax.ticks;
                    var tl = (ax.ticklen === undefined) ? null : ax.ticklen;
                    var ts = (ax.ticklabelstandoff === undefined) ? null : ax.ticklabelstandoff;
                    return (t !== "" || tl !== 0 || ts !== 0);
                }

                function _needsXAxisFix(axisKey) {
                    if (!axisKey || !gd.layout[axisKey]) return false;
                    var ax = gd.layout[axisKey];
                    var z = (ax.zeroline === undefined) ? null : ax.zeroline;
                    var sg = (ax.showgrid === undefined) ? null : ax.showgrid;
                    var ar = (ax.autorange === undefined) ? null : ax.autorange;
                    var r = ax.range;
                    // "bad" state we see after zoom-reset: autorange or range dipping < 0, which makes a zeroline/gridline at x=0 appear.
                    var badRange = (r && r.length === 2 && r[0] < 0);
                    return (z !== false || sg !== false || ar === true || badRange);
                }

                var keys = [];
                if (fig.layout && fig.layout.meta) {
                    if (fig.layout.meta.protein_lollipop_yaxis_key) keys.push(fig.layout.meta.protein_lollipop_yaxis_key);
                    if (fig.layout.meta.protein_boxes_yaxis_key) keys.push(fig.layout.meta.protein_boxes_yaxis_key);
                }

                var xkeys = [];
                if (fig.layout && fig.layout.meta) {
                    if (fig.layout.meta.protein_lollipop_xaxis_key) xkeys.push(fig.layout.meta.protein_lollipop_xaxis_key);
                    if (fig.layout.meta.protein_boxes_xaxis_key) xkeys.push(fig.layout.meta.protein_boxes_xaxis_key);
                }

                var upd = {};
                var doIt = false;
                keys.forEach(function(k) {
                    if (_needsTight(k)) {
                        upd[k + ".ticks"] = "";
                        upd[k + ".ticklen"] = 0;
                        upd[k + ".ticklabelstandoff"] = 0;
                        upd[k + ".showline"] = true;
                        doIt = true;
                    }
                });


                xkeys.forEach(function(k) {
                    if (!_needsXAxisFix(k)) return;
                    upd[k + ".zeroline"] = false;
                    upd[k + ".showgrid"] = false;

                    var ax = gd.layout[k] || {};
                    // If Plotly reset to autorange (double click), snap back to the figure's intended range.
                    if (ax.autorange === true || !ax.range) {
                        if (fig.layout && fig.layout[k] && fig.layout[k].range && fig.layout[k].range.length === 2) {
                            upd[k + ".range[0]"] = fig.layout[k].range[0];
                            upd[k + ".range[1]"] = fig.layout[k].range[1];
                            upd[k + ".autorange"] = false;
                        }
                    } else if (ax.range && ax.range.length === 2 && ax.range[0] < 0) {
                        upd[k + ".range[0]"] = 0;
                        upd[k + ".range[1]"] = ax.range[1];
                        upd[k + ".autorange"] = false;
                    }
                    doIt = true;
                });
if (doIt) {
                    // Relayout can fire another relayout event; this converges quickly.
                    window.Plotly.relayout(gd, upd);
                }
            }
        } catch (e) {
            // ignore
        }

        if (!fig.data) fig.data = [];
        if (!fig.layout) fig.layout = {};
        if (!fig.layout.annotations) fig.layout.annotations = [];
        if (!fig.layout.shapes) fig.layout.shapes = [];

        // Remove previous pinned overlay traces/annotations/shapes.
        fig.data = fig.data.filter(function(tr) {
            return !(tr.meta && tr.meta.pinned_overlay);
        });
        fig.layout.annotations = fig.layout.annotations.filter(function(ann) {
            return !(ann.meta && ann.meta.pinned_overlay);
        });
        fig.layout.shapes = fig.layout.shapes.filter(function(shape) {
            return !(shape.meta && shape.meta.pinned_overlay);
        });

        // If nothing is pinned, just return the base figure.
        if (!pinnedData || !Array.isArray(pinnedData) || pinnedData.length === 0) {
            return fig;
        }
        if (!curatedData || !Array.isArray(curatedData) || curatedData.length === 0) {
            return fig;
        }

        // Map genomic coord -> record
        var coordToRecord = {};
        curatedData.forEach(function(rec) {
            var c = rec["GRCh38 coordinates"];
            if (c === undefined || c === null) return;
            c = parseInt(c);
            if (!isNaN(c)) {
                coordToRecord[c] = rec;
            }
        });

        // Color map matching Python's FIXED_SOURCE_COLORS
        var sourceColorMap = {
            "LOVD": "#3D6C88",
            "Publication": "#7B9041",
            "COSMIC": "#E4BC33",
            "Unknown": "#7f7f7f",
            "Multiple": "#141414"
        };

        // Helper function to extract AA positions from protein labels
        function extractAARange(proteinLabels) {
            if (!proteinLabels) return null;

            var labels = Array.isArray(proteinLabels) ? proteinLabels : [proteinLabels];
            var positions = [];

            // Regex to extract AA position from protein notation (e.g., p.Leu569Leu -> 569)
            var aaRegex = /p\.\(?[A-Za-z]{1,3}(\d+)/gi;

            labels.forEach(function(label) {
                if (!label) return;
                var match;
                while ((match = aaRegex.exec(String(label))) !== null) {
                    var pos = parseInt(match[1]);
                    if (!isNaN(pos)) {
                        positions.push(pos);
                    }
                }
            });

            if (positions.length === 0) return null;

            return {
                min: Math.min.apply(null, positions),
                max: Math.max.apply(null, positions)
            };
        }

        // Helper function to find exon from CDS_EXON_AA_RANGES
        function findExonFromAARange(aaRange, exVis) {
            if (!aaRange || !exVis) return null;

            for (var i = 0; i < exVis.length; i++) {
                var ex = exVis[i];
                // Check if AA range overlaps with this exon
                if (aaRange.min <= ex.aa1 && aaRange.max >= ex.aa0) {
                    return ex;
                }
            }
            return null;
        }

        // Get exon visualization data from metadata
        var exVis = (fig.layout && fig.layout.meta && fig.layout.meta.protein_ex_vis) || null;

        // Helper to convert AA to visual x coordinate
                function aaToVx(aa) {
            if (aa === null || aa === undefined) return null;
            // Protein track x-axis is in true AA coordinates.
            return aa;
        }

        // Build pinned point list with protein info
        var pinnedPoints = [];
        var exonsToHighlight = new Set();
        var aaRangesToHighlight = [];

        pinnedData.forEach(function(pin) {
            if (!pin || pin.coord === undefined || pin.coord === null) return;
            var c = parseInt(pin.coord);
            if (isNaN(c)) return;

            var rec = coordToRecord[c];
            if (!rec) return;

            var x = rec["Jittered_Coord"];
            if (x === undefined || x === null) x = rec["GRCh38 coordinates"];
            var y = rec["Total_AC"] || 1;
            var hover = rec["Hover_Text"] || ("Position " + c);

            var src = rec["Display_Source"] || rec["Source_Category"] || "Unknown";
            var color = sourceColorMap[src] || "#111827";

            pinnedPoints.push({coord: c, x: x, y: y, hover: hover, color: color, rec: rec});

            // Extract protein AA range and find corresponding exon
            var proteinLabels = rec["Protein_Label"];
            var aaRange = extractAARange(proteinLabels);

            if (aaRange && exVis) {
                var exon = findExonFromAARange(aaRange, exVis);
                if (exon) {
                    exonsToHighlight.add(exon.exon_display);
                }
                aaRangesToHighlight.push(aaRange);
            }
        });

        if (pinnedPoints.length === 0) {
            return fig;
        }

        // Highlight exons in protein boxes row (row 7)
        var proteinBoxesRow = 7;
        exonsToHighlight.forEach(function(exonDisplay) {
            var exon = exVis.find(function(ex) { return ex.exon_display === exonDisplay; });
            if (!exon) return;

            // Add orange background highlight for the exon
            fig.layout.shapes.push({
                type: "rect",
                xref: "x" + proteinBoxesRow,
                yref: "y" + proteinBoxesRow,
                x0: exon.vx0,
                x1: exon.vx1,
                y0: 0.10,
                y1: 0.90,
                fillcolor: "rgba(255, 165, 0, 0.4)",  // Orange with transparency
                line: {
                    color: "#FF8C00",  // Dark orange border
                    width: 4
                },
                layer: "above",
                meta: { pinned_overlay: true }
            });

            // Add diagonal hatching pattern (///) to make it more visible
            var hatchSpacing = 6;  // Space between diagonal lines
            var boxHeight = 0.90 - 0.10;  // Height of the exon box (0.80)
            var boxWidth = exon.vx1 - exon.vx0;

            // Calculate how many diagonal lines we need
            var numLines = Math.ceil(boxWidth / hatchSpacing) + Math.ceil(boxHeight * 10);

            for (var i = 0; i < numLines; i++) {
                // Start position for diagonal line
                var startX = exon.vx0 + (i * hatchSpacing);
                var startY = 0.10;

                // End position (diagonal at ~45 degrees, scaled to box proportions)
                var endX = startX + boxHeight;
                var endY = 0.90;

                // Clip the line to stay within exon boundaries
                if (startX < exon.vx0) {
                    // Adjust start point if it's before exon start
                    var diff = exon.vx0 - startX;
                    startX = exon.vx0;
                    startY = 0.10 + diff;
                }

                if (endX > exon.vx1) {
                    // Clip end point if it extends beyond exon end
                    var diff = endX - exon.vx1;
                    endX = exon.vx1;
                    endY = 0.90 - diff;
                }

                // Only draw if line is within bounds
                if (startX <= exon.vx1 && endX >= exon.vx0 && startY <= 0.90 && endY >= 0.10) {
                    fig.layout.shapes.push({
                        type: "line",
                        xref: "x" + proteinBoxesRow,
                        yref: "y" + proteinBoxesRow,
                        x0: startX,
                        y0: startY,
                        x1: endX,
                        y1: endY,
                        line: {
                            color: "#FF8C00",
                            width: 2
                        },
                        layer: "above",
                        meta: { pinned_overlay: true }
                    });
                }
            }
        });

        // Highlight AA positions in protein lollipop row (row 6)
        var proteinLollipopRow = 6;

        // Helper function to convert visual x back to approximate AA
                function vxToAa(vx) {
            if (vx === null || vx === undefined) return null;
            // Protein track x-axis is in true AA coordinates; round to nearest residue index.
            return Math.round(vx);
        }

        // Find all lollipop points that fall within highlighted AA ranges
        var lollipopPointsToHighlight = [];

        // Only proceed if we have AA ranges to highlight and exon data
        if (aaRangesToHighlight.length > 0 && exVis) {
            // Expected axis names for row 6
            var expectedXAxis = "x" + proteinLollipopRow;
            var expectedYAxis = "y" + proteinLollipopRow;

            // Look through existing traces to find lollipop head markers
            fig.data.forEach(function(trace) {
                // Skip if not a marker trace or already marked as overlay
                if (!trace.mode || trace.mode.indexOf("markers") === -1) return;
                if (trace.meta && trace.meta.pinned_overlay) return;

                // Check if this trace is assigned to row 6's axes
                // trace.xaxis and trace.yaxis are strings like "x6", "y6"
                var traceXAxis = trace.xaxis;
                var traceYAxis = trace.yaxis;

                // Skip if not in the protein lollipop row
                if (traceXAxis !== expectedXAxis || traceYAxis !== expectedYAxis) return;

                if (!trace.x || !trace.y) return;

                // Check each marker point
                for (var i = 0; i < trace.x.length; i++) {
                    var vx = trace.x[i];
                    var vy = trace.y[i];

                    // Skip null/undefined values
                    if (vx === null || vx === undefined || vy === null || vy === undefined) continue;

                    // Convert visual x back to AA position
                    var aa = vxToAa(vx);
                    if (aa === null) continue;

                    // Check if this AA position is within any highlighted range
                    // Add small tolerance (±0.5) for rounding issues
                    var inRange = false;
                    for (var j = 0; j < aaRangesToHighlight.length; j++) {
                        var range = aaRangesToHighlight[j];
                        if (aa >= (range.min - 0.5) && aa <= (range.max + 0.5)) {
                            inRange = true;
                            break;
                        }
                    }

                    if (inRange) {
                        lollipopPointsToHighlight.push({
                            x: vx,
                            y: vy,
                            hover: trace.text ? (Array.isArray(trace.text) ? trace.text[i] : trace.text) : null
                        });
                    }
                }
            });
        }

        // Add orange circles around the highlighted lollipop points
        if (lollipopPointsToHighlight.length > 0) {
            var x_coords = lollipopPointsToHighlight.map(function(p) { return p.x; });
            var y_coords = lollipopPointsToHighlight.map(function(p) { return p.y; });
            var hovers = lollipopPointsToHighlight.map(function(p) { return p.hover || ""; });

            // Add larger orange circles to highlight
            fig.data.push({
                x: x_coords,
                y: y_coords,
                mode: "markers",
                marker: {
                    size: 14,
                    color: "rgba(255, 165, 0, 0.6)",  // Semi-transparent orange fill
                    line: {
                        color: "#FF8C00",  // Dark orange border
                        width: 3
                    }
                },
                hoverinfo: "text",
                text: hovers,
                showlegend: false,
                xaxis: expectedXAxis,
                yaxis: expectedYAxis,
                meta: { pinned_overlay: true }
            });
        }

        // Sort by x for clustering
        pinnedPoints.sort(function(a, b) { return a.x - b.x; });

        // Cluster points that are close in x (so boxes don't totally overlap)
        var xMin = null, xMax = null;
        if (fig.layout && fig.layout.xaxis && Array.isArray(fig.layout.xaxis.range)) {
            xMin = fig.layout.xaxis.range[0];
            xMax = fig.layout.xaxis.range[1];
        }
        var clusterThreshold = 20;
        if (xMin !== null && xMax !== null) {
            var span = xMax - xMin;
            var thr = span * 0.002;
            if (thr > clusterThreshold) clusterThreshold = thr;
        }

        var clusters = [];
        var currentCluster = [pinnedPoints[0]];
        var lastX = pinnedPoints[0].x;

        for (var i = 1; i < pinnedPoints.length; i++) {
            var pt = pinnedPoints[i];
            if (Math.abs(pt.x - lastX) <= clusterThreshold) {
                currentCluster.push(pt);
            } else {
                clusters.push(currentCluster);
                currentCluster = [pt];
            }
            lastX = pt.x;
        }
        clusters.push(currentCluster);

        // Max y for small vertical offsets
        var maxY = 1;
        curatedData.forEach(function(rec) {
            var t = rec["Total_AC"];
            if (t && t > maxY) maxY = t;
        });

        var pinnedLegendAdded = false;

        clusters.forEach(function(cluster) {
            cluster.forEach(function(pt, j) {
                // Circle highlight around pinned lollipop (same color as the lollipop)
                fig.data.push({
                    x: [pt.x],
                    y: [pt.y],
                    mode: "markers",
                    marker: {
                        size: 14,
                        color: "rgba(0,0,0,0)",
                        line: { color: pt.color, width: 2 },
                        symbol: "circle-open"
                    },
                    hoverinfo: "text",
                    text: pt.hover,
                    name: pinnedLegendAdded ? null : "Pinned / search match",
                    showlegend: !pinnedLegendAdded,
                    xaxis: "x",
                    yaxis: "y",
                    meta: { pinned_overlay: true }
                });
                pinnedLegendAdded = true;

                // Hover-style annotation above the lollipop, not on top of it
                var baseShift = -40;
                var perClusterShift = -18;
                var yshift = baseShift + j * perClusterShift;

                fig.layout.annotations.push({
                    x: pt.x,
                    y: pt.y,
                    xref: "x",
                    yref: "y",
                    xanchor: "left",
                    yanchor: "bottom",
                    yshift: yshift,
                    text: pt.hover,
                    showarrow: false,
                    bordercolor: pt.color,
                    borderwidth: 1,
                    bgcolor: "rgba(255,255,255,0.97)",
                    opacity: 0.98,
                    align: "left",
                    font: { size: 10 },
                    meta: { pinned_overlay: true }
                });
            });
        });

        return fig;
    }
    """,
    Output("mutation-plot", "figure"),
    Input("base-figure-store", "data"),
    Input("pinned-variants", "data"),
    Input("curated-grouped-data", "data"),
    Input("mutation-plot", "relayoutData"),
)

'''
app.clientside_callback(
    r"""
    function(baseFig, pinnedData, curatedData, relayoutData) {
        // If we don't have a base figure yet, just return empty.
        if (!baseFig) {
            return {};
        }

        // Deep copy so we don't mutate the stored base figure.
        var fig = JSON.parse(JSON.stringify(baseFig));

        // --- Dynamic AA ticks on protein axis (densify when you zoom) ---
        try {
            if (fig.layout && fig.layout.meta && fig.layout.meta.protein_ex_vis) {
                var exVis = fig.layout.meta.protein_ex_vis;
                var aaMaxAll = fig.layout.meta.protein_aa_max || null;

                // Find the protein x-axis (the one titled "Protein position (AA)")
                var protAxisKey = null;
                Object.keys(fig.layout).forEach(function(k) {
                    if (k.indexOf("xaxis") === 0) {
                        var ax = fig.layout[k];
                        if (ax && ax.title && ax.title.text && ax.title.text.indexOf("Protein position") !== -1) {
                            protAxisKey = k;
                        }
                    }
                });

                                function vxToAa(vx) {
                    if (vx === null || vx === undefined) return null;
                    // Protein track x-axis is in true AA coordinates; round to nearest residue index.
                    return Math.round(vx);
                }

                                function aaToVx(aa) {
                    if (aa === null || aa === undefined) return null;
                    // Protein track x-axis is in true AA coordinates.
                    return aa;
                }

                function niceStep(span, targetTicks) {
                    if (!isFinite(span) || span <= 0) return 1;
                    var rough = span / targetTicks;
                    var pow10 = Math.pow(10, Math.floor(Math.log10(rough)));
                    var base = rough / pow10;
                    var niceBase = (base <= 1) ? 1 : (base <= 2) ? 2 : (base <= 5) ? 5 : 10;
                    return niceBase * pow10;
                }

                if (protAxisKey) {
                    // Prefer the current zoom window from relayoutData if it refers to the protein axis
                    var relBase = (protAxisKey === "xaxis") ? "xaxis" : protAxisKey;
                    var vx0 = null, vx1 = null;

                    if (relayoutData) {
                        var k0 = relBase + ".range[0]";
                        var k1 = relBase + ".range[1]";
                        if (relayoutData.hasOwnProperty(k0) && relayoutData.hasOwnProperty(k1)) {
                            vx0 = relayoutData[k0];
                            vx1 = relayoutData[k1];
                        }
                    }

                    // Fallback to axis range if relayoutData is for a different panel
                    if (vx0 === null || vx1 === null) {
                        var axRange = fig.layout[protAxisKey] && fig.layout[protAxisKey].range;
                        if (axRange && axRange.length === 2) {
                            vx0 = axRange[0];
                            vx1 = axRange[1];
                        }
                    }

                    if (vx0 !== null && vx1 !== null) {
                        var vmin = Math.min(vx0, vx1);
                        var vmax = Math.max(vx0, vx1);

                        var aaMin = vxToAa(vmin);
                        var aaMax = vxToAa(vmax);
                        if (aaMin !== null && aaMax !== null) {
                            aaMin = Math.max(1, Math.floor(aaMin));
                            aaMax = Math.ceil(aaMax);
                            if (aaMaxAll) aaMax = Math.min(aaMaxAll, aaMax);

                            var span = Math.max(1, aaMax - aaMin);
                            var step = Math.max(1, niceStep(span, 7));

                            var start = Math.ceil(aaMin / step) * step;
                            var tickvals = [];
                            var ticktext = [];
                            for (var a = start; a <= aaMax; a += step) {
                                var vx = aaToVx(a);
                                if (vx === null) continue;
                                tickvals.push(vx);
                                ticktext.push(String(a));
                                if (tickvals.length > 30) break; // safety
                            }

                            if (tickvals.length >= 2) {
                                fig.layout[protAxisKey].tickmode = "array";
                                fig.layout[protAxisKey].tickvals = tickvals;
                                fig.layout[protAxisKey].ticktext = ticktext;
                            }
                        }
                    }
                }
            }
        } catch (e) {
            // If anything goes wrong, just keep the default ticks.
        }


        if (!fig.data) fig.data = [];
        if (!fig.layout) fig.layout = {};
        if (!fig.layout.annotations) fig.layout.annotations = [];

        // Remove previous pinned overlay traces/annotations.
        fig.data = fig.data.filter(function(tr) {
            return !(tr.meta && tr.meta.pinned_overlay);
        });
        fig.layout.annotations = fig.layout.annotations.filter(function(ann) {
            return !(ann.meta && ann.meta.pinned_overlay);
        });

        // If nothing is pinned, just return the base figure.
        if (!pinnedData || !Array.isArray(pinnedData) || pinnedData.length === 0) {
            return fig;
        }
        if (!curatedData || !Array.isArray(curatedData) || curatedData.length === 0) {
            return fig;
        }

        // Map genomic coord -> record
        var coordToRecord = {};
        curatedData.forEach(function(rec) {
            var c = rec["GRCh38 coordinates"];
            if (c === undefined || c === null) return;
            c = parseInt(c);
            if (!isNaN(c)) {
                coordToRecord[c] = rec;
            }
        });

        // Color map matching Python's FIXED_SOURCE_COLORS
        var sourceColorMap = {
            "LOVD": "#3D6C88",
            "Publication": "#7B9041",
            "COSMIC": "#E4BC33",
            "Unknown": "#7f7f7f",
            "Multiple": "#141414"
        };

        // Build pinned point list
        var pinnedPoints = [];
        pinnedData.forEach(function(pin) {
            if (!pin || pin.coord === undefined || pin.coord === null) return;
            var c = parseInt(pin.coord);
            if (isNaN(c)) return;

            var rec = coordToRecord[c];
            if (!rec) return;

            var x = rec["Jittered_Coord"];
            if (x === undefined || x === null) x = rec["GRCh38 coordinates"];
            var y = rec["Total_AC"] || 1;
            var hover = rec["Hover_Text"] || ("Position " + c);

            var src = rec["Display_Source"] || rec["Source_Category"] || "Unknown";
            var color = sourceColorMap[src] || "#111827";

            pinnedPoints.push({coord: c, x: x, y: y, hover: hover, color: color});
        });

        if (pinnedPoints.length === 0) {
            return fig;
        }

        // Sort by x for clustering
        pinnedPoints.sort(function(a, b) { return a.x - b.x; });

        // Cluster points that are close in x (so boxes don't totally overlap)
        var xMin = null, xMax = null;
        if (fig.layout && fig.layout.xaxis && Array.isArray(fig.layout.xaxis.range)) {
            xMin = fig.layout.xaxis.range[0];
            xMax = fig.layout.xaxis.range[1];
        }
        var clusterThreshold = 20;
        if (xMin !== null && xMax !== null) {
            var span = xMax - xMin;
            var thr = span * 0.002;
            if (thr > clusterThreshold) clusterThreshold = thr;
        }

        var clusters = [];
        var currentCluster = [pinnedPoints[0]];
        var lastX = pinnedPoints[0].x;

        for (var i = 1; i < pinnedPoints.length; i++) {
            var pt = pinnedPoints[i];
            if (Math.abs(pt.x - lastX) <= clusterThreshold) {
                currentCluster.push(pt);
            } else {
                clusters.push(currentCluster);
                currentCluster = [pt];
            }
            lastX = pt.x;
        }
        clusters.push(currentCluster);

        // Max y for small vertical offsets
        var maxY = 1;
        curatedData.forEach(function(rec) {
            var t = rec["Total_AC"];
            if (t && t > maxY) maxY = t;
        });

        var pinnedLegendAdded = false;

                clusters.forEach(function(cluster) {
            cluster.forEach(function(pt, j) {
                // Circle highlight around pinned lollipop (same color as the lollipop)
                fig.data.push({
                    x: [pt.x],
                    y: [pt.y],
                    mode: "markers",
                    marker: {
                        size: 14,
                        color: "rgba(0,0,0,0)",
                        line: { color: pt.color, width: 2 },
                        symbol: "circle-open"
                    },
                    hoverinfo: "text",
                    text: pt.hover,
                    name: pinnedLegendAdded ? null : "Pinned / search match",
                    showlegend: !pinnedLegendAdded,
                    xaxis: "x",
                    yaxis: "y",
                    meta: { pinned_overlay: true }
                });
                pinnedLegendAdded = true;

                // Hover-style annotation above the lollipop, not on top of it
                // Negative yshift moves the box upward in *pixels*,
                // and yanchor="bottom" makes the whole box sit above the point.
                var baseShift = -40;        // how far above the lollipop (in pixels)
                var perClusterShift = -18;  // extra upwards spacing for stacked boxes
                var yshift = baseShift + j * perClusterShift;

                fig.layout.annotations.push({
                    x: pt.x,
                    y: pt.y,
                    xref: "x",
                    yref: "y",
                    xanchor: "left",
                    yanchor: "bottom",
                    yshift: yshift,
                    text: pt.hover,
                    showarrow: false,
                    bordercolor: pt.color,
                    borderwidth: 1,
                    bgcolor: "rgba(255,255,255,0.97)",
                    opacity: 0.98,
                    align: "left",
                    font: { size: 10 },
                    meta: { pinned_overlay: true }
                });
            });
        });

        return fig;
    }
    """,
    Output("mutation-plot", "figure"),
    Input("base-figure-store", "data"),
    Input("pinned-variants", "data"),
    Input("curated-grouped-data", "data"),
    Input("mutation-plot", "relayoutData"),
)
'''


@app.callback(
    Output("pinned-variants", "data"),
    Output("variant-search-feedback", "children"),
    Output("variant-search-input", "value"),
    Input("variant-search-button", "n_clicks"),
    Input("variant-search-reset", "n_clicks"),
    Input("mutation-plot", "clickData"),
    Input("variant-search-input", "n_submit"),
    State("variant-search-input", "value"),
    State("pinned-variants", "data"),
    State("curated-grouped-data", "data"),
    prevent_initial_call=True,
)
def handle_variant_search(
        n_search,
        n_reset,
        click_data,
        n_submit,
        query,
        pinned_data,
        curated_data,
):
    ctx = dash.callback_context
    pinned_data = pinned_data or []
    records = curated_data or []

    if not ctx.triggered:
        return pinned_data, "", query

    trigger = ctx.triggered[0]["prop_id"].split(".")[0]

    if trigger == "variant-search-reset":
        return [], "", ""

    if trigger == "mutation-plot":
        if not click_data or not click_data.get("points"):
            return pinned_data, "", query

        point = click_data["points"][0]
        x_clicked = point.get("x")
        if x_clicked is None or not records:
            return pinned_data, "", query

        best_rec = None
        best_delta = None
        x_clicked = float(x_clicked)
        for rec in records:
            j = rec.get("Jittered_Coord")
            if j is None:
                continue
            delta = abs(float(j) - x_clicked)
            if best_delta is None or delta < best_delta:
                best_delta = delta
                best_rec = rec

        if best_rec is None or best_delta is None or best_delta > 10:
            return pinned_data, "", query

        try:
            coord = int(best_rec.get("GRCh38 coordinates"))
        except (TypeError, ValueError):
            return pinned_data, "", query

        existing_coords = {
            int(item.get("coord"))
            for item in pinned_data
            if item and "coord" in item
        }

        if coord in existing_coords:
            new_pins = [
                item
                for item in pinned_data
                if not (
                        item
                        and "coord" in item
                        and int(item["coord"]) == coord
                )
            ]
            msg = f"Unpinned variant at {coord}. {len(new_pins)} pinned."
            return new_pins, msg, query

        new_entry = {"coord": coord, "query": "click"}
        new_pins = pinned_data + [new_entry]
        msg = f"Pinned variant at {coord} by click. {len(new_pins)} total pinned."
        return new_pins, msg, query

    if trigger in ("variant-search-button", "variant-search-input"):
        q = (query or "").strip()

        if not q:
            return (
                pinned_data,
                "Please enter a cDNA (e.g., c.1333C>T), protein (p.Thr489Ile), "
                "or genomic variant (chr13:48303995 C>G).",
                "",
            )

        if not records:
            return (
                pinned_data,
                "No curated variants are currently displayed under the selected filters, "
                "so the search could not find any matches.",
                "",
            )

        matches = search_in_curated_records(q, records)

        if not matches:
            return pinned_data, f"No variant in the current view matched '{q}'.", ""

        existing_coords = {
            int(item.get("coord"))
            for item in pinned_data
            if item and "coord" in item
        }

        new_entries = []
        for rec in matches:
            try:
                coord = int(rec.get("GRCh38 coordinates"))
            except (TypeError, ValueError):
                continue
            if coord in existing_coords:
                continue
            new_entries.append({"coord": coord, "query": q})
            existing_coords.add(coord)

        if not new_entries:
            return (
                pinned_data,
                f"All variants matching '{q}' are already pinned.",
                "",
            )

        updated_pins = pinned_data + new_entries
        msg = (
            f"Pinned {len(new_entries)} variant(s) for '{q}'. "
            f"{len(updated_pins)} total pinned."
        )
        return updated_pins, msg, ""

    return pinned_data, "", query


@app.callback(
    Output("download-dataframe-csv", "data"),
    Input("btn_csv", "n_clicks"),
    State("clinvar-path-filter", "value"),
    State("mc-filter", "value"),
    prevent_initial_call=True,
)
def download_csv(n_clicks, selected_clinvar_classes, selected_mc_list):
    full_df = pd.read_excel(str(EXCEL_FILE), sheet_name=INPUT_SHEET)
    full_df.columns = full_df.columns.str.strip()

    if "clinical_classification_combined" in full_df.columns:
        base_cls = full_df["clinical_classification_combined"]
    elif "clinvar_variant_classification" in full_df.columns:
        base_cls = full_df["clinvar_variant_classification"]
    else:
        base_cls = None

    if base_cls is not None:
        full_df["ClinVar_Category"] = (
            base_cls.astype(str)
            .str.strip()
            .replace({"nan": "Unspecified", "": "Unspecified"})
        )
    else:
        full_df["ClinVar_Category"] = "Unspecified"

    if "lovd_classification" in full_df.columns:
        lovd_base = full_df["lovd_classification"]
        full_df["LOVD_Category"] = (
            lovd_base.astype(str)
            .str.strip()
            .replace({"nan": "Unspecified", "": "Unspecified"})
        )
    else:
        full_df["LOVD_Category"] = "Unspecified"

    def normalize_all(vals):
        if not vals:
            return []
        if "__ALL__" in vals:
            return ["__ALL__"]
        return vals

    selected_clinvar_classes = normalize_all(selected_clinvar_classes)
    selected_mc_list = normalize_all(selected_mc_list)

    if (
            selected_clinvar_classes != ["__ALL__"]
            and "ClinVar_Category" in full_df.columns
    ):
        full_df = full_df[
            full_df["ClinVar_Category"].isin(selected_clinvar_classes)
        ]

    if (
            selected_mc_list != ["__ALL__"]
            and "molecular_consequence" in full_df.columns
    ):
        mc_set = set(selected_mc_list)
        full_df = full_df[
            full_df["molecular_consequence"]
            .astype(str)
            .str.strip()
            .isin(mc_set)
        ]

    return dcc.send_data_frame(
        full_df.to_csv,
        f"{GENE_NAME}_filtered_variants.csv",
        index=False,
    )


def _format_count_dict(counts):
    if not counts:
        return "-"
    parts = []
    for name, val in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0] or "")):
        label = name if (name is not None and name != "") else "unknown"
        parts.append(f"{label} ({int(val)})")
    return ", ".join(parts)


def _format_score_dict(counts):
    if not counts:
        return "-"
    keys = [
        k
        for k in counts.keys()
        if k is not None and str(k).strip().lower() not in {"", "nan", "none", "n/a"}
    ]
    if not keys:
        return "-"
    if len(keys) == 1:
        return str(keys[0])
    return _format_count_dict(counts)


@app.callback(
    Output("pinned-variants-list", "children"),
    Input("pinned-variants", "data"),
    Input("curated-grouped-data", "data"),
)
def render_pinned_variant_list(pinned_data, curated_data):
    pinned_data = pinned_data or []
    curated_data = curated_data or []

    if not pinned_data:
        return "No variants pinned yet. Click on a lollipop or use the search box to pin variants."

    if not curated_data:
        return "Pinned variants will appear here once variants are loaded in the current view."

    coord_to_rec = {}
    for rec in curated_data:
        coord = rec.get("GRCh38 coordinates")
        if coord is None:
            continue
        try:
            c_int = int(coord)
        except Exception:
            continue
        coord_to_rec[c_int] = rec

    # Styles (same as your current ones)
    table_style = {
        "width": "100%",
        "borderCollapse": "collapse",
        "fontSize": "10px",
        "fontFamily": "inherit",
        "marginTop": "2px",
    }
    th_style = {
        "borderBottom": "1px solid #9ca3af",
        "textAlign": "left",
        "padding": "2px 4px",
        "fontWeight": "600",
        "verticalAlign": "bottom",
        "whiteSpace": "nowrap",
        # Keep table headers visible while scrolling the pinned table
        "position": "sticky",
        "top": "0px",
        "backgroundColor": "#f9fafb",
        "zIndex": 5,
    }
    td_style = {
        "borderBottom": "1px solid #e5e7eb",
        "padding": "2px 4px",
        "verticalAlign": "top",
        "whiteSpace": "normal",
        "wordBreak": "break-word",
    }
    group_header_style = {
        "padding": "4px 4px 2px 4px",
        "backgroundColor": "#f9fafb",
        "borderBottom": "1px solid #d1d5db",
        "fontSize": "10px",
        "fontWeight": "600",
        "color": "#111827",
    }

    other_variants_style = {
        "padding": "2px 4px",
        "backgroundColor": "#ffffff",
        "borderBottom": "1px solid #e5e7eb",
        "fontSize": "10px",
        "fontStyle": "italic",
        "color": "#374151",
    }

    def _norm_text(s):
        s = '' if s is None else str(s)
        s = s.strip().lower()
        # remove all whitespace for robust matching (e.g., 'c.861G>A' vs 'c.861G>A ')
        s = re.sub(r'\s+', '', s)
        return s

    def _combo_matches_query(combo, q_norm):
        if not q_norm:
            return False
        g_norm = _norm_text(combo.get('genomic', ''))
        c_norm = _norm_text(combo.get('coding', ''))
        p_norm = _norm_text(combo.get('protein', ''))
        return (
                (q_norm == g_norm) or (q_norm == c_norm) or (q_norm == p_norm)
                or (q_norm in g_norm) or (q_norm in c_norm) or (q_norm in p_norm)
        )

    body_rows = []

    for idx, pin in enumerate(pinned_data, start=1):
        coord = pin.get("coord")
        if coord is None:
            continue
        try:
            c_int = int(coord)
        except Exception:
            continue

        rec = coord_to_rec.get(c_int)
        if rec is None:
            # Pinned variant might be filtered out in the current view
            continue

        row = pd.Series(rec)

        # Use existing helper to get combos
        combo_items = get_combo_counts_for_row(
            row,
            per_variant_map=per_variant_map,
            use_full_genomic=True,
            return_meta=True,
        )
        if not combo_items:
            continue

        region_label = rec.get("Region_Label", "-")

        # Group header row (one per pinned variant)
        header_text = f"{idx}. {c_int}"
        if region_label:
            header_text += f" ({region_label})"

        body_rows.append(
            html.Tr(
                [
                    html.Td(
                        header_text,
                        colSpan=15,  # now 11 columns total
                        style=group_header_style,
                    )
                ]
            )
        )

        # Data rows for each genomic/cDNA/protein combo
        pin_query = pin.get('query')
        is_search_pin = bool(pin_query) and str(pin_query) != 'click'
        q_norm = _norm_text(pin_query) if is_search_pin else ''

        def _append_combo_row(combo):
            g_str = combo.get("genomic", "-")
            c_str = combo.get("coding", "-")
            p_str = combo.get("protein") or "-"
            mc_txt = _format_count_dict(combo.get("molconseq_counts"))
            # Try a couple of possible keys for Nomad/gnomAD frequencyF
            nomad_txt = _format_count_dict(
                combo.get("nomad_counts") or combo.get("af_counts")
            )
            path_txt = _format_count_dict(combo.get("path_class_counts"))

            cadd_txt = _format_score_dict(combo.get("cadd_counts"))
            spliceai_txt = _format_score_dict(combo.get("spliceai_counts"))
            alphamissense_txt = _format_score_dict(combo.get("alphamissense_counts"))
            origin_txt = _format_count_dict(combo.get("origin_counts"))
            inher_txt = _format_count_dict(
                combo.get("inheritance_counts") or combo.get("inheritance_counts")
            )
            later_txt = _format_count_dict(combo.get("laterality_counts"))
            source_txt = _format_count_dict(combo.get("source_counts"))
            count_val = int(combo.get("total_count", 0) or 0)

            body_rows.append(
                html.Tr(
                    [
                        # 1–3: variant labels
                        html.Td(g_str, style=td_style),  # Genomic
                        html.Td(c_str, style=td_style),  # coding
                        html.Td(p_str, style=td_style),  # protein
                        # 4: region
                        html.Td(region_label or "-", style=td_style),  # Region
                        # 5–7: consequence, gnomAD, ClinVar
                        html.Td(mc_txt, style=td_style),  # Molecular consequence
                        html.Td(nomad_txt, style=td_style),  # Nomad Frequency
                        html.Td(path_txt, style=td_style),  # Clinvar Pathogenicity
                        html.Td(cadd_txt, style=td_style),  # CADD
                        html.Td(spliceai_txt, style=td_style),  # SpliceAI
                        html.Td(alphamissense_txt, style=td_style),  # AlphaMissense
                        html.Td(str(count_val), style=td_style),  # Count
                        # 9–11: origin / inheritance / laterality
                        html.Td(origin_txt, style=td_style),  # Genetic Origin
                        html.Td(inher_txt, style=td_style),  # Inheritance
                        html.Td(later_txt, style=td_style),  # Laterality
                        # 12: source
                        html.Td(source_txt, style=td_style),  # Source
                    ]
                )
            )

        if is_search_pin and q_norm:
            primary = []
            others = []
            for combo in combo_items:
                if _combo_matches_query(combo, q_norm):
                    primary.append(combo)
                else:
                    others.append(combo)

            if primary:
                for combo in primary:
                    _append_combo_row(combo)
                if others:
                    body_rows.append(
                        html.Tr(
                            [
                                html.Td(
                                    'Other variants in this position',
                                    colSpan=15,
                                    style=other_variants_style,
                                )
                            ]
                        )
                    )
                    for combo in others:
                        _append_combo_row(combo)
            else:
                # If we can't identify a specific match, fall back to the current behavior
                for combo in combo_items:
                    _append_combo_row(combo)
        else:
            # Click pins (or unknown pins) keep the current behavior: show everything at that position
            for combo in combo_items:
                _append_combo_row(combo)

    if not body_rows:
        return "Pinned variants will appear here once variants match the current filters."

    # Header row in the same order
    header_row = html.Tr(
        [
            html.Th("Genomic", style=th_style),
            html.Th("Coding", style=th_style),
            html.Th("Protein", style=th_style),
            html.Th("Region", style=th_style),
            html.Th("Molecular consequence", style=th_style),
            html.Th("gnomAD Frequency", style=th_style),
            html.Th("Clinvar Pathogenicity", style=th_style),
            html.Th("CADD", style=th_style),
            html.Th("SpliceAI", style=th_style),
            html.Th("AlphaMissense", style=th_style),
            html.Th("Count", style=th_style),
            html.Th("Genetic Origin", style=th_style),
            html.Th("Inheritance", style=th_style),
            html.Th("Laterality", style=th_style),
            html.Th("Source", style=th_style),
        ]
    )

    table = html.Table(
        [html.Thead(header_row), html.Tbody(body_rows)],
        style=table_style,
    )

    return table


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", "8051")))
