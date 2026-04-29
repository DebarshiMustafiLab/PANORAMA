library(tidyverse)
library(janitor)
library(readxl)
library(broom)
library(gt)
library(writexl)
library(glue)
library(scales)
library(rlang)
library(pwr)
library(ggrepel)
library(ggpattern)
library(colorspace)
library(htmltools)
library(patchwork)
library(ggtext)
library(grid)


###LOAD DATA
fpath <- "/Users/rooks/Desktop/Coding/VS\ Code/panorama_rb1/retinoblastoma_data_v1.xlsx"
rb <- read_excel(fpath, sheet = 1) %>% clean_names()

rb_clean_germline <- rb %>% filter(genetic_origin == "germline")

rb_clean_unique_germline <- rb_clean_germline %>%
  distinct(cdna, .keep_all = TRUE)


###CHI SQUARED AND FISHERS ANALYSES

#CHI SQUARED FOR GENETIC ORIGIN AND LATERALITY
origin_laterality <- rb %>%
  filter(genetic_origin %in% c("germline","somatic"),
         laterality %in% c("bilateral","unilateral")) %>%
  count(genetic_origin, laterality) %>%
  pivot_wider(names_from = laterality, values_from = n, values_fill = 0)

chisq_table <- origin_laterality %>%
  column_to_rownames("genetic_origin") %>%
  as.matrix()

chisq_table

chisq_test <- chisq.test(chisq_table, correct = FALSE)
chisq_test


a <- chisq_table["germline", "bilateral"]
b <- chisq_table["germline", "unilateral"]
c <- chisq_table["somatic",  "bilateral"]
d <- chisq_table["somatic",  "unilateral"]

or  <- (a * d) / (b * c)
log_or <- log(or)
se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
z <- qnorm(0.975)

ci_lower <- exp(log_or - z * se_log_or)
ci_upper <- exp(log_or + z * se_log_or)

data.frame(
  OR = or,
  CI_lower = ci_lower,
  CI_upper = ci_upper,
  chisq = as.numeric(chisq_test$statistic),
  df = chisq_test$parameter,
  p = chisq_test$p.value
)


#CHI SQUARED FOR FAMILY HISTORY AND LATERALITY

rb_subset_laterality_family <- rb_clean_germline %>%
  filter(
    laterality %in% c("unilateral", "bilateral"),
    family_history %in% c("isolated", "familial")
  )

contingency_laterality_family <- table(
  rb_subset_laterality_family$laterality,
  rb_subset_laterality_family$family_history
)

contingency_laterality_family

chisq_laterality_family <- chisq.test(contingency_laterality_family)

chisq_laterality_family

#FISHERS FOR PARENT OF ORIGIN AND LATERALITY
parent_laterality <- rb_clean_germline %>%
  filter(!is.na(inheritance),
         inheritance %in% c("maternal", "paternal"),
         laterality %in% c("bilateral", "unilateral"))

table_parent_laterality <- table(parent_laterality$inheritance,
                                 parent_laterality$laterality)

fisher_parent_laterality <- fisher.test(table_parent_laterality)

fisher_parent_laterality

table_parent_laterality


#MOST COMMON GERMLINE VARIANTS
top_germline <- rb_clean_germline %>%
  group_by(cdna) %>%
  summarise(
    count = n(),
    molecular_consequence = first(molecular_consequence),
    exon_intron = first(exon_intron_combined),
    clinvar_variant_classification = first(clinvar_variant_classification),
    mean_cadd = round(mean(cadd_score, na.rm = TRUE), 2),
    gnomad_frequency = formatC(mean(gnomad_frequency, na.rm = TRUE), format = "e", digits = 2),
    .groups = "drop"
  ) %>%
  arrange(desc(count)) %>%
  slice_head(n = 15)

list(top_germline = top_germline)



###TABLE 1: COHORT STUDIES INCLUDED

tab <- read_excel(
  "/Users/rooks/Desktop/Coding/VS Code/panorama_rb1/retinoblastoma_data_v1.xlsx",
  sheet = 3,
  range = "A1:D15"
)

tab %>%
  gt() %>%
  tab_header(
    title = md("**Table 1. Cohort studies included**")
  ) %>%
  tab_options(
    heading.align = "left",
    table.font.size = 12,
    heading.title.font.size = 12,
    heading.title.font.weight = "bold",
    table.border.top.color = "grey",
    table.border.bottom.color = "grey",
    data_row.padding = px(4),
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  ) %>%
  opt_table_font(font = "Times New Roman") %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

# ###TABLE 2. CHARACTERISTICS OF RB1 VARIANTS
# 
# mk_block <- function(data, group_col, levels_vec, label_map, section_label) {
#   
#   wide <- data %>%
#     mutate(
#       Category_raw = factor(.data[[group_col]], levels = levels_vec)
#     ) %>%
#     count(Category_raw, genetic_origin, name = "n") %>%
#     tidyr::pivot_wider(
#       names_from = genetic_origin,
#       values_from = n,
#       values_fill = 0
#     ) %>%
#     mutate(
#       germline = coalesce(germline, 0L),
#       somatic  = coalesce(somatic, 0L),
#       unknown  = coalesce(unknown, 0L),
#       Total    = germline + somatic + unknown,
#       Section  = section_label,
#       Category = label_map[as.character(Category_raw)]
#     ) %>%
#     arrange(Category_raw) %>%
#     select(-Category_raw)
#   
#   den_g     <- sum(wide$germline, na.rm = TRUE)
#   den_s     <- sum(wide$somatic,  na.rm = TRUE)
#   den_u     <- sum(wide$unknown,  na.rm = TRUE)
#   den_total <- sum(wide$Total,    na.rm = TRUE)
#   
#   fmt_inline <- function(x, den) {
#     ifelse(
#       den > 0 & x > 0,
#       sprintf("%d (%.1f)", x, 100 * x / den),
#       sprintf("%d", x)
#     )
#   }
#   
#   wide %>%
#     transmute(
#       Section,
#       Category,
#       `Germline, n (%)` = fmt_inline(germline, den_g),
#       `Somatic, n (%)`  = fmt_inline(somatic,  den_s),
#       `Unknown, n (%)`  = fmt_inline(unknown,  den_u),
#       `Total, n (%)`    = fmt_inline(Total,    den_total)
#     )
# }
# 
# tbl_sources <- rb %>%
#   mutate(source_clean = case_when(
#     str_detect(str_to_lower(str_trim(as.character(source))), "lovd")        ~ "Lovd",
#     str_detect(str_to_lower(str_trim(as.character(source))), "publication") ~ "Publication",
#     str_detect(str_to_lower(str_trim(as.character(source))), "cosmic")      ~ "Cosmic",
#     TRUE ~ NA_character_
#   )) %>%
#   filter(!is.na(source_clean)) %>%
#   mk_block(
#     group_col   = "source_clean",
#     levels_vec  = c("Lovd", "Publication", "Cosmic"),
#     label_map   = c("Lovd" = "LOVD",
#                     "Publication" = "Publications",
#                     "Cosmic" = "COSMIC"),
#     section_label = "Sources"
#   )
# 
# tbl_laterality <- rb %>%
#   mk_block(
#     group_col   = "laterality",
#     levels_vec  = c("bilateral", "unilateral", "unknown"),
#     label_map   = c("bilateral" = "Bilateral",
#                     "unilateral" = "Unilateral",
#                     "unknown" = "Unknown"),
#     section_label = "Laterality"
#   )
# 
# tbl_famhx <- rb %>%
#   mk_block(
#     group_col   = "family_history",
#     levels_vec  = c("familial", "isolated", "unknown"),
#     label_map   = c("familial" = "Yes",
#                     "isolated" = "No",
#                     "unknown"  = "Unknown"),
#     section_label = "Family History"
#   )
# 
# dat_all <- rb %>%
#   mutate(
#     molecular_consequence = case_when(
#       molecular_consequence %in% c("stop_gained", "nonsense") ~ "Nonsense",
#       molecular_consequence %in% c("frameshift_variant", "frameshift") ~ "Frameshift",
#       molecular_consequence %in% c("splice_donor", "splice_acceptor") ~ "Splice Site",
#       molecular_consequence %in% c("missense_variant", "missense") ~ "Missense",
#       molecular_consequence %in% c("intron_variant", "intron") ~ "Intronic",
#       molecular_consequence %in% c("inframe_deletion", "inframe_insertion", "inframe_indel") ~ "Inframe",
#       molecular_consequence %in% c("synonymous_variant", "synonymous") ~ "Synonymous",
#       molecular_consequence %in% c("5_prime_UTR_variant", "5_prime_UTR") ~ "5 prime UTR",
#       molecular_consequence %in% c("3_prime_UTR_variant", "3_prime_UTR") ~ "3 prime UTR",
#       molecular_consequence %in% c("start_lost", "start_lost_variant") ~ "Start Lost",
#       molecular_consequence %in% c("deletion", "large_deletion", "structural_deletion") ~ "Large Deletion",
#       TRUE ~ molecular_consequence
#     ) %>%
#       str_replace_all("_", " ") %>%
#       str_squish()
#   )
# 
# mc_levels <- dat_all %>%
#   count(molecular_consequence, sort = TRUE) %>%
#   pull(molecular_consequence)
# 
# tbl_consequence <- dat_all %>%
#   mk_block(
#     group_col   = "molecular_consequence",
#     levels_vec  = mc_levels,
#     label_map   = setNames(mc_levels, mc_levels),
#     section_label = "Molecular consequence"
#   )
# 
# combined <- bind_rows(tbl_sources, tbl_laterality, tbl_famhx, tbl_consequence)
# 
# grand_total <- nrow(rb)
# total_den_g <- sum(rb$genetic_origin == "germline", na.rm = TRUE)
# total_den_s <- sum(rb$genetic_origin == "somatic",  na.rm = TRUE)
# total_den_u <- sum(rb$genetic_origin == "unknown",  na.rm = TRUE)
# 
# fmt_inline_total <- function(x) {
#   ifelse(
#     grand_total > 0 & x > 0,
#     sprintf("%d (%.1f)", x, 100 * x / grand_total),
#     sprintf("%d", x)
#   )
# }
# 
# total_row <- tibble(
#   Section = "",
#   Category = "Total",
#   `Germline, n (%)` = fmt_inline_total(total_den_g),
#   `Somatic, n (%)`  = fmt_inline_total(total_den_s),
#   `Unknown, n (%)`  = fmt_inline_total(total_den_u),
#   `Total, n (%)`    = sprintf("%d", grand_total)
# )
# 
# final_tbl <- bind_rows(combined, total_row)
# 
# final_tbl %>%
#   gt(groupname_col = "Section") %>%
#   tab_header(title = md("**Table 2. Characteristics of RB1 Variants**")) %>%
#   tab_options(heading.align = "left") %>%
#   cols_align(align = "left", columns = everything()) %>%
#   opt_table_font(font = "Times New Roman") %>%
#   cols_width(
#     Category ~ px(190)
#   ) %>%
#   tab_options(
#     table.width = px(700),
#     table.font.size = 12,
#     heading.title.font.size = 12,
#     heading.title.font.weight = "bold",
#     table.border.top.color = "grey",
#     table.border.bottom.color = "grey",
#     data_row.padding = px(4)
#   ) %>%
#   # bold top header row
#   tab_style(
#     style = cell_text(weight = "bold"),
#     locations = cells_column_labels(everything())
#   ) %>%
#   # bold section labels (Sources, Laterality, etc)
#   tab_style(
#     style = cell_text(weight = "bold"),
#     locations = cells_row_groups()
#   ) %>%
#   # bold Total row
#   tab_style(
#     style = cell_text(weight = "bold"),
#     locations = cells_body(rows = Category == "Total")
#   ) %>%
#   # indent non-Total category labels
#   tab_style(
#     style = cell_text(indent = px(16)),
#     locations = cells_body(
#       columns = "Category",
#       rows = Category != "Total"
#     )
#   )



#TABLE 2 UPDATED WITH GERMLINE MOSAIC

mk_block <- function(data, group_col, levels_vec, label_map, section_label) {
  
  data <- data %>%
    mutate(
      origin_group = case_when(
        genetic_origin_source == "germline (mosaic)" ~ "germline mosaic",
        genetic_origin == "germline" ~ "germline",
        genetic_origin == "somatic" ~ "somatic",
        genetic_origin == "unknown" ~ "unknown",
        TRUE ~ as.character(genetic_origin)
      )
    )
  
  wide <- data %>%
    mutate(
      Category_raw = factor(.data[[group_col]], levels = levels_vec)
    ) %>%
    count(Category_raw, origin_group, name = "n") %>%
    tidyr::pivot_wider(
      names_from = origin_group,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      germline = coalesce(germline, 0L),
      `germline mosaic` = coalesce(`germline mosaic`, 0L),
      somatic  = coalesce(somatic, 0L),
      unknown  = coalesce(unknown, 0L),
      Total    = germline + `germline mosaic` + somatic + unknown,
      Section  = section_label,
      Category = label_map[as.character(Category_raw)]
    ) %>%
    arrange(Category_raw) %>%
    select(-Category_raw)
  
  den_g     <- sum(wide$germline, na.rm = TRUE)
  den_gm    <- sum(wide$`germline mosaic`, na.rm = TRUE)
  den_s     <- sum(wide$somatic,  na.rm = TRUE)
  den_u     <- sum(wide$unknown,  na.rm = TRUE)
  den_total <- sum(wide$Total,    na.rm = TRUE)
  
  fmt_inline <- function(x, den) {
    ifelse(
      den > 0 & x > 0,
      sprintf("%d (%.1f)", x, 100 * x / den),
      sprintf("%d", x)
    )
  }
  
  wide %>%
    transmute(
      Section,
      Category,
      `Germline, n (%)` = fmt_inline(germline, den_g),
      `Germline mosaic, n (%)` = fmt_inline(`germline mosaic`, den_gm),
      `Somatic, n (%)`  = fmt_inline(somatic,  den_s),
      `Unknown, n (%)`  = fmt_inline(unknown,  den_u),
      `Total, n (%)`    = fmt_inline(Total,    den_total)
    )
}

tbl_sources <- rb %>%
  mutate(source_clean = case_when(
    str_detect(str_to_lower(str_trim(as.character(source))), "lovd")        ~ "Lovd",
    str_detect(str_to_lower(str_trim(as.character(source))), "publication") ~ "Publication",
    str_detect(str_to_lower(str_trim(as.character(source))), "cosmic")      ~ "Cosmic",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(source_clean)) %>%
  mk_block(
    group_col   = "source_clean",
    levels_vec  = c("Lovd", "Publication", "Cosmic"),
    label_map   = c("Lovd" = "LOVD",
                    "Publication" = "Publications",
                    "Cosmic" = "COSMIC"),
    section_label = "Sources"
  )

tbl_laterality <- rb %>%
  mk_block(
    group_col   = "laterality",
    levels_vec  = c("bilateral", "unilateral", "unknown"),
    label_map   = c("bilateral" = "Bilateral",
                    "unilateral" = "Unilateral",
                    "unknown" = "Unknown"),
    section_label = "Laterality"
  )

tbl_famhx <- rb %>%
  mk_block(
    group_col   = "family_history",
    levels_vec  = c("familial", "isolated", "unknown"),
    label_map   = c("familial" = "Yes",
                    "isolated" = "No",
                    "unknown"  = "Unknown"),
    section_label = "Family History"
  )

dat_all <- rb %>%
  mutate(
    molecular_consequence = case_when(
      molecular_consequence %in% c("stop_gained", "nonsense") ~ "Nonsense",
      molecular_consequence %in% c("frameshift_variant", "frameshift") ~ "Frameshift",
      molecular_consequence %in% c("splice_donor", "splice_acceptor") ~ "Splice Site",
      molecular_consequence %in% c("missense_variant", "missense") ~ "Missense",
      molecular_consequence %in% c("intron_variant", "intron") ~ "Intronic",
      molecular_consequence %in% c("inframe_deletion", "inframe_insertion", "inframe_indel") ~ "Inframe",
      molecular_consequence %in% c("synonymous_variant", "synonymous") ~ "Synonymous",
      molecular_consequence %in% c("5_prime_UTR_variant", "5_prime_UTR") ~ "5 prime UTR",
      molecular_consequence %in% c("3_prime_UTR_variant", "3_prime_UTR") ~ "3 prime UTR",
      molecular_consequence %in% c("start_lost", "start_lost_variant") ~ "Start Lost",
      molecular_consequence %in% c("deletion", "large_deletion", "structural_deletion") ~ "Large Deletion",
      TRUE ~ molecular_consequence
    ) %>%
      str_replace_all("_", " ") %>%
      str_squish()
  )

mc_levels <- dat_all %>%
  count(molecular_consequence, sort = TRUE) %>%
  pull(molecular_consequence)

tbl_consequence <- dat_all %>%
  mk_block(
    group_col   = "molecular_consequence",
    levels_vec  = mc_levels,
    label_map   = setNames(mc_levels, mc_levels),
    section_label = "Molecular consequence"
  )

combined <- bind_rows(tbl_sources, tbl_laterality, tbl_famhx, tbl_consequence)

rb_totals <- rb %>%
  mutate(
    origin_group = case_when(
      genetic_origin_source == "germline (mosaic)" ~ "germline mosaic",
      genetic_origin == "germline" ~ "germline",
      genetic_origin == "somatic" ~ "somatic",
      genetic_origin == "unknown" ~ "unknown",
      TRUE ~ as.character(genetic_origin)
    )
  )

grand_total <- nrow(rb_totals)
total_den_g  <- sum(rb_totals$origin_group == "germline", na.rm = TRUE)
total_den_gm <- sum(rb_totals$origin_group == "germline mosaic", na.rm = TRUE)
total_den_s  <- sum(rb_totals$origin_group == "somatic", na.rm = TRUE)
total_den_u  <- sum(rb_totals$origin_group == "unknown", na.rm = TRUE)

fmt_inline_total <- function(x) {
  ifelse(
    grand_total > 0 & x > 0,
    sprintf("%d (%.1f)", x, 100 * x / grand_total),
    sprintf("%d", x)
  )
}

total_row <- tibble(
  Section = "",
  Category = "Total",
  `Germline, n (%)` = fmt_inline_total(total_den_g),
  `Germline mosaic, n (%)` = fmt_inline_total(total_den_gm),
  `Somatic, n (%)`  = fmt_inline_total(total_den_s),
  `Unknown, n (%)`  = fmt_inline_total(total_den_u),
  `Total, n (%)`    = sprintf("%d", grand_total)
)

final_tbl <- bind_rows(combined, total_row)

final_tbl %>%
  gt(groupname_col = "Section") %>%
  tab_header(title = md("**Table 2. Characteristics of RB1 Variants**")) %>%
  tab_options(heading.align = "left") %>%
  cols_align(align = "left", columns = everything()) %>%
  opt_table_font(font = "Times New Roman") %>%
  cols_width(
    Category ~ px(190)
  ) %>%
  tab_options(
    table.width = px(850),
    table.font.size = 12,
    heading.title.font.size = 12,
    heading.title.font.weight = "bold",
    table.border.top.color = "grey",
    table.border.bottom.color = "grey",
    data_row.padding = px(4)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = Category == "Total")
  ) %>%
  tab_style(
    style = cell_text(indent = px(16)),
    locations = cells_body(
      columns = "Category",
      rows = Category != "Total"
    )
  )



#========= Variant consequence by laterality: counts and within-laterality percentages ========

rb_lat <- rb_clean_germline %>%
  dplyr::filter(
    laterality %in% c("unilateral", "bilateral"),
    !is.na(molecular_consequence),
    molecular_consequence != ""
  ) %>%
  dplyr::mutate(
    laterality = factor(laterality, levels = c("unilateral", "bilateral")),
    consequence_raw = stringr::str_trim(as.character(molecular_consequence)),
    consequence_binned = dplyr::case_when(
      consequence_raw %in% c("splice_donor", "splice_acceptor", "splice_site") ~ "splice_site",
      consequence_raw %in% c("inframe_deletion", "inframe_insertion", "inframe_indel", "inframe_variant") ~ "inframe_variant",
      TRUE ~ consequence_raw
    )
  )

by_lat <- rb_lat %>%
  dplyr::count(laterality, consequence_binned, name = "N") %>%
  dplyr::group_by(laterality) %>%
  dplyr::mutate(Total_lat = sum(N), Percent = 100 * N / Total_lat) %>%
  dplyr::ungroup()

wide <- by_lat %>%
  dplyr::mutate(n_pct_lbl = sprintf("%d (%.1f)", N, Percent)) %>%
  dplyr::select(consequence_binned, laterality, n_pct_lbl) %>%
  tidyr::pivot_wider(
    names_from  = laterality,
    values_from = n_pct_lbl
  ) %>%
  dplyr::rename(
    `Unilateral, n (%)` = unilateral,
    `Bilateral, n (%)`  = bilateral
  )

tot_by_conseq <- rb_lat %>%
  dplyr::count(consequence_binned, name = "Total (n)")

tbl <- wide %>%
  dplyr::left_join(tot_by_conseq, by = "consequence_binned") %>%
  dplyr::arrange(dplyr::desc(`Total (n)`)) %>%
  dplyr::rename(Variant_Consequence = consequence_binned)

# Totals row
tot_unilat <- by_lat %>% dplyr::filter(laterality == "unilateral") %>% dplyr::summarise(N = sum(N)) %>% dplyr::pull(N)
tot_bilat  <- by_lat %>% dplyr::filter(laterality == "bilateral")  %>% dplyr::summarise(N = sum(N)) %>% dplyr::pull(N)

tot_row <- tibble::tibble(
  Variant_Consequence = "Total",
  `Unilateral, n (%)` = if (!is.na(tot_unilat) && tot_unilat > 0) sprintf("%d (100.0)", tot_unilat) else "—",
  `Bilateral, n (%)`  = if (!is.na(tot_bilat)  && tot_bilat  > 0) sprintf("%d (100.0)",  tot_bilat)  else "—",
  `Total (n)`         = nrow(rb_lat)
)

tbl_out <- dplyr::bind_rows(tbl, tot_row)

gt::gt(tbl_out) %>%
  gt::tab_header(title = md("**Table 7. Variant Consequence versus Tumor Laterality**")) %>%
  gt::cols_label(
    Variant_Consequence = "Variant Consequence",
    `Unilateral, n (%)` = "Unilateral, n (%)",
    `Bilateral, n (%)`  = "Bilateral, n (%)",
    `Total (n)`         = "Total (n)"
  ) %>%
  gt::cols_align(align = "left", columns = everything()) %>%
  gt::opt_table_font(font = "Times New Roman") %>%
  gt::fmt_missing(columns = everything(), missing_text = "—") %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = list(
      gt::cells_body(rows = Variant_Consequence == "Total")
    )
  ) %>%
  gt::tab_options(
    heading.align = "left",
    table.font.size = 12,
    heading.title.font.size = 12,
    heading.title.font.weight = "bold",
    table.border.top.color = "grey",
    table.border.bottom.color = "grey",
    data_row.padding = gt::px(4), 
    table.width = gt::pct(70)   # <- make table scale with window/container
  )



#========= avg cadd score vs laterality (wilcoxon and table)
cadd_laterality <- rb_clean_unique_germline %>%
  filter(!is.na(cadd_score),
         laterality %in% c("bilateral", "unilateral"))

cadd_summary <- cadd_laterality %>%
  group_by(laterality) %>%
  summarise(
    n = n(),
    mean_cadd = mean(cadd_score, na.rm = TRUE),
    median_cadd = median(cadd_score, na.rm = TRUE),
    sd_cadd = sd(cadd_score, na.rm = TRUE),
    .groups = "drop"
  )

wilcox_test <- wilcox.test(cadd_score ~ laterality, data = cadd_laterality)
wilcox_text <- glue(
  "Wilcoxon rank-sum test: W = {round(wilcox_test$statistic,1)}, p = {signif(wilcox_test$p.value,3)}"
)

gt_tbl <- cadd_summary %>%
  mutate(
    mean_cadd = sprintf("%.2f", mean_cadd),
    median_cadd = sprintf("%.2f", median_cadd),
    sd_cadd = sprintf("%.2f", sd_cadd)
  ) %>%
  rename(
    Laterality = laterality,
    `n (variants)` = n,
    `Mean CADD` = mean_cadd,
    `Median CADD` = median_cadd,
    `SD` = sd_cadd
  ) %>%
  gt() %>%
  tab_header(title = "CADD Scores by Laterality") %>%
  fmt_number(columns = c(`Mean CADD`, `Median CADD`, `SD`), decimals = 2) %>%
  opt_table_font(font = "Times New Roman") %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.title.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    data_row.padding = px(4)
  ) %>%
  tab_footnote(
    footnote = wilcox_text,
    locations = cells_title(groups = "title")
  )

gt_tbl


#========= Mean cadd scores by latearialty looking at variant consequence -----


rb_lat_cadd <- rb_clean_germline %>%
  filter(
    !is.na(cadd_score),
    laterality %in% c("unilateral", "bilateral"),
    !is.na(molecular_consequence)
  ) %>%
  mutate(
    laterality = factor(laterality, levels = c("unilateral", "bilateral")),
    consequence_raw = str_trim(as.character(molecular_consequence)),
    consequence_binned = case_when(
      consequence_raw %in% c("splice_donor", "splice_acceptor", "splice_site") ~ "splice_site",
      consequence_raw %in% c("inframe_deletion", "inframe_insertion", "inframe_indel", "inframe_variant") ~ "inframe_variant",
      TRUE ~ consequence_raw
    )
  )

summarise_ci <- function(df) {
  n  <- nrow(df)
  mu <- mean(df$cadd_score)
  sd <- sd(df$cadd_score)
  se <- sd / sqrt(n)
  ci_low  <- mu - 1.96 * se
  ci_high <- mu + 1.96 * se
  tibble(N = n, Mean = mu, SD = sd, SE = se, CI_low = ci_low, CI_high = ci_high)
}

by_lat <- rb_lat_cadd %>%
  group_by(consequence_binned, laterality) %>%
  group_modify(~ summarise_ci(.x)) %>%
  ungroup() %>%
  mutate(
    Mean_lbl = ifelse(N > 0, sprintf("%.2f", Mean), "—"),
    CI_lbl   = ifelse(N > 1, sprintf("%.2f–%.2f", CI_low, CI_high), "NA–NA")
  ) %>%
  select(consequence_binned, laterality, N, Mean_lbl, CI_lbl) %>%
  pivot_wider(
    names_from = laterality,
    values_from = c(N, Mean_lbl, CI_lbl),
    names_sep = "_"
  ) %>%
  rename(
    N_unilat       = N_unilateral,
    Mean_unilat    = Mean_lbl_unilateral,
    CI_unilat      = CI_lbl_unilateral,
    N_bilat        = N_bilateral,
    Mean_bilat     = Mean_lbl_bilateral,
    CI_bilat       = CI_lbl_bilateral
  )

by_total <- rb_lat_cadd %>%
  group_by(consequence_binned) %>%
  group_modify(~ summarise_ci(.x)) %>%
  ungroup() %>%
  mutate(
    Total_N = N,
    Total_Mean_lbl = sprintf("%.2f", Mean),
    Total_CI_lbl   = ifelse(Total_N > 1, sprintf("%.2f–%.2f", CI_low, CI_high), "NA–NA")
  ) %>%
  select(consequence_binned, Total_N, Total_Mean_lbl, Total_CI_lbl)

tbl <- by_lat %>%
  left_join(by_total, by = "consequence_binned") %>%
  arrange(desc(Total_N)) %>%
  select(
    Variant_Consequence = consequence_binned,
    N_unilat, Mean_unilat, CI_unilat,
    N_bilat,  Mean_bilat,  CI_bilat,
    `Total (n)` = Total_N,
    `Total Mean CADD` = Total_Mean_lbl,
    `Total 95% CI` = Total_CI_lbl
  )

tot_unilat <- rb_lat_cadd %>% filter(laterality == "unilateral")
tot_bilat  <- rb_lat_cadd %>% filter(laterality == "bilateral")

tot_row <- tibble(
  Variant_Consequence = "Total",
  N_unilat = nrow(tot_unilat),
  Mean_unilat = ifelse(nrow(tot_unilat) > 0, sprintf("%.2f", mean(tot_unilat$cadd_score)), "—"),
  CI_unilat = ifelse(
    nrow(tot_unilat) > 1,
    {
      sd <- sd(tot_unilat$cadd_score); se <- sd / sqrt(nrow(tot_unilat))
      sprintf("%.2f–%.2f", mean(tot_unilat$cadd_score) - 1.96 * se,
              mean(tot_unilat$cadd_score) + 1.96 * se)
    },
    "NA–NA"
  ),
  N_bilat = nrow(tot_bilat),
  Mean_bilat = ifelse(nrow(tot_bilat) > 0, sprintf("%.2f", mean(tot_bilat$cadd_score)), "—"),
  CI_bilat = ifelse(
    nrow(tot_bilat) > 1,
    {
      sd <- sd(tot_bilat$cadd_score); se <- sd / sqrt(nrow(tot_bilat))
      sprintf("%.2f–%.2f", mean(tot_bilat$cadd_score) - 1.96 * se,
              mean(tot_bilat$cadd_score) + 1.96 * se)
    },
    "NA–NA"
  ),
  `Total (n)` = nrow(rb_lat_cadd),
  `Total Mean CADD` = sprintf("%.2f", mean(rb_lat_cadd$cadd_score)),
  `Total 95% CI` = {
    sd <- sd(rb_lat_cadd$cadd_score); se <- sd / sqrt(nrow(rb_lat_cadd))
    sprintf("%.2f–%.2f", mean(rb_lat_cadd$cadd_score) - 1.96 * se,
            mean(rb_lat_cadd$cadd_score) + 1.96 * se)
  }
)

tbl_out <- bind_rows(tbl, tot_row)

gt(tbl_out) %>%
  tab_header(title = md("**Table 4. Variant Consequence versus Laterality with Mean CADD Scores**")) %>%
  tab_spanner(label = "Unilateral Cases", columns = c(N_unilat, Mean_unilat, CI_unilat)) %>%
  tab_spanner(label = "Bilateral Cases",  columns = c(N_bilat,  Mean_bilat,  CI_bilat)) %>%
  cols_label(
    Variant_Consequence = "Variant Consequence",
    N_unilat = "Count",
    Mean_unilat = "Mean CADD",
    CI_unilat = "95% CI",
    N_bilat  = "Count",
    Mean_bilat = "Mean CADD",
    CI_bilat = "95% CI"
  ) %>%
  cols_align(align = "left", columns = everything()) %>%
  opt_table_font(font = "Times New Roman") %>%
  fmt_missing(everything(), missing_text = "—") %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_body(columns = c(N_unilat, N_bilat, `Total (n)`)),
      cells_body(rows = Variant_Consequence == "Total")
    )
  ) %>%
  tab_options(
    heading.align = "left",
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.title.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    data_row.padding = px(4)
  ) 





#========= lovd and clinvar variant classifications confusion matrix 
five <- c("Benign","Likely benign","Uncertain significance","Likely pathogenic","Pathogenic")

cm <- rb_clean_unique_germline %>%
  dplyr::filter(!is.na(clinical_classification_combined),
                !is.na(clinvar_variant_classification)) %>%
  dplyr::mutate(
    clinvar_classification = dplyr::case_when(
      clinvar_variant_classification %in% c("Likely benign") ~ "Likely benign",
      clinvar_variant_classification %in% c("Likely pathogenic") ~ "Likely pathogenic",
      clinvar_variant_classification %in% c("Benign","Benign/Likely benign") ~ "Benign",
      clinvar_variant_classification %in% c("Pathogenic","Pathogenic/Likely pathogenic") ~ "Pathogenic",
      clinvar_variant_classification %in% c("Uncertain significance","Conflicting classifications of pathogenicity") ~ "Uncertain significance",
      TRUE ~ clinvar_variant_classification
    ),
    lovd_classification = dplyr::case_when(
      clinical_classification_combined %in% c("likely benign") ~ "Likely benign",
      clinical_classification_combined %in% c("likely pathogenic") ~ "Likely pathogenic",
      clinical_classification_combined %in% c("benign") ~ "Benign",
      clinical_classification_combined %in% c("pathogenic") ~ "Pathogenic",
      clinical_classification_combined %in% c("VUS") ~ "Uncertain significance",
      TRUE ~ clinical_classification_combined
    )
  ) %>%
  dplyr::count(lovd_classification, clinvar_classification, name = "n") %>%
  tidyr::pivot_wider(names_from = clinvar_classification, values_from = n, values_fill = 0) %>%
  dplyr::mutate(.row_order = match(lovd_classification, five)) %>%
  dplyr::arrange(.row_order, lovd_classification) %>%
  dplyr::select(-.row_order) %>%
  dplyr::select(lovd_classification, dplyr::any_of(five), dplyr::everything())

total_n <- sum(as.matrix(cm[, five]), na.rm = TRUE)

cm_fmt <- cm
for (col in five) {
  cm_fmt[[col]] <- ifelse(
    is.na(cm[[col]]) | cm[[col]] == 0,
    "0",
    paste0(cm[[col]], " (", round(100 * cm[[col]] / total_n, 1), ")")
  )
}

row_totals <- rowSums(cm[, five, drop = FALSE], na.rm = TRUE)
cm_fmt$Total <- ifelse(
  row_totals == 0,
  "0",
  paste0(row_totals, " (", round(100 * row_totals / total_n, 1), ")")
)

col_totals <- colSums(cm[, five, drop = FALSE], na.rm = TRUE)
total_col_total <- sum(row_totals, na.rm = TRUE)
total_row <- data.frame(
  lovd_classification = "Total",
  t(setNames(
    ifelse(col_totals == 0,
           "0",
           paste0(col_totals, " (", round(100 * col_totals / total_n, 1), ")")),
    five
  )),
  Total = paste0(total_col_total, " (100%)"),
  check.names = FALSE
)

tbl <- dplyr::bind_rows(cm_fmt, total_row)

gt::gt(tbl) %>%
  gt::tab_header(
    title = gt::md("**Table 4. Comparison of LOVD and ClinVar Variant Classifications**")
  ) %>%
  gt::tab_options(
    heading.align = "left",
    table.font.size = 12,
    data_row.padding = gt::px(4),
    column_labels.font.weight = "normal"
  ) %>%
  gt::cols_label(
    lovd_classification = gt::md("**LOVD Classification**"),
    Benign = "Benign",
    `Likely benign` = "Likely Benign",
    `Uncertain significance` = "Uncertain Significance",
    `Likely pathogenic` = "Likely Pathogenic",
    Pathogenic = "Pathogenic",
    Total = "Total"
  ) %>%
  gt::tab_spanner(
    label = gt::md("**ClinVar Classification**"),
    columns = c(Benign, `Likely benign`, `Uncertain significance`, `Likely pathogenic`, Pathogenic)
  ) %>%
  gt::opt_table_font(font = "Times New Roman") %>%
  gt::tab_source_note(
    gt::md("ClinVar classifications Pathogenic/Likely Pathogenic were grouped with Pathogenic and Benign/Likely Benign were grouped with Benign.")
  )








#========= ClinVar classification ----
cadd_by_clinvar <- rb_clean_unique_germline %>%
  filter(!is.na(cadd_score), !is.na(clinvar_variant_classification)) %>%
  mutate(
    clinvar_classification = case_when(
      clinvar_variant_classification %in% c("Benign", "Benign/Likely benign") ~ "Benign",
      clinvar_variant_classification %in% c("Likely benign") ~ "Likely benign",
      clinvar_variant_classification %in% c("Uncertain significance", "Conflicting classifications of pathogenicity") ~ "Uncertain significance",
      clinvar_variant_classification %in% c("Likely pathogenic") ~ "Likely pathogenic",
      clinvar_variant_classification %in% c("Pathogenic", "Pathogenic/Likely pathogenic") ~ "Pathogenic",
      TRUE ~ clinvar_variant_classification
    ),
    clinvar_classification = factor(
      clinvar_classification,
      levels = c("Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic")
    )
  ) %>%
  group_by(clinvar_classification) %>%
  summarise(
    Mean_CADD = round(mean(cadd_score, na.rm = TRUE), 2),
    SD_CADD   = round(sd(cadd_score,   na.rm = TRUE), 2),
    N = n(),
    SE = SD_CADD / sqrt(N),
    CI_low  = round(Mean_CADD - 1.96 * SE, 2),
    CI_high = round(Mean_CADD + 1.96 * SE, 2),
    .groups = "drop"
  ) %>%
  mutate(CI = paste0("[", CI_low, "–", CI_high, "]")) %>%
  dplyr::select(clinvar_classification, N, Mean_CADD, SD_CADD, CI)

total_n <- sum(cadd_by_clinvar$N)
cadd_by_clinvar <- cadd_by_clinvar %>%
  mutate(
    Percent = round(100 * N / total_n, 1),
    Count_label = paste0(N, " (", Percent, "%)")
  )

total_row <- rb_clean_unique_germline %>%
  filter(!is.na(cadd_score), !is.na(clinvar_variant_classification)) %>%
  summarise(
    clinvar_classification = "Total",
    Mean_CADD = round(mean(cadd_score, na.rm = TRUE), 2),
    SD_CADD   = round(sd(cadd_score,   na.rm = TRUE), 2),
    N = n(),
    SE = SD_CADD / sqrt(N),
    CI_low  = round(Mean_CADD - 1.96 * SE, 2),
    CI_high = round(Mean_CADD + 1.96 * SE, 2)
  ) %>%
  mutate(
    CI = paste0("[", CI_low, "–", CI_high, "]"),
    Count_label = paste0(N, " (100%)")
  ) %>%
  dplyr::select(clinvar_classification, Count_label, Mean_CADD, SD_CADD, CI)

cadd_by_clinvar_total <- cadd_by_clinvar %>%
  dplyr::select(clinvar_classification, Count_label, Mean_CADD, SD_CADD, CI) %>%
  bind_rows(total_row)

cadd_by_clinvar_total %>%
  gt() %>%
  tab_header(
    title = md("**Table 5. ClinVar Classification and CADD Scores**")
  ) %>%
  tab_options(
    heading.align = "left",
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.title.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    data_row.padding = px(4)
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  ) %>%
  cols_label(
    clinvar_classification = md("**ClinVar Classification**"),
    Count_label            = md("**Count (n, %)**"),
    Mean_CADD              = md("**Mean CADD**"),
    SD_CADD                = md("**SD**"),
    CI                     = md("**95% CI**")
  ) %>%
  fmt_number(columns = c(Mean_CADD, SD_CADD), decimals = 2) %>%
  opt_table_font(font = "Times New Roman") %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = clinvar_classification == "Total")
  ) %>%
  sub_missing(everything(), missing_text = "—") %>%
  tab_source_note(
    md("Mean CADD per ClinVar classification. One variant (c.979_1033dup, ClinVar Pathogenic) lacked a CADD score and was excluded. Pathogenic/Likely Pathogenic grouped with Pathogenic; Benign/Likely Benign grouped with Benign.")
  )



#========= benign likely benign clinvar classifications --> cadd

benign_low_cadd <- rb_clean_unique_germline %>%
  filter(
    str_to_lower(clinvar_variant_classification) %in%
      c("benign", "likely benign", "benign/likely benign"),
    cadd_score < 20)


#========= alpha missense breakdown by cadd score
missense_summary <- rb_clean_unique_germline %>%
  # Keep only missense variants
  filter(molecular_consequence == "missense") %>%
  # Convert and bin CADD scores
  mutate(
    cadd_score = suppressWarnings(as.numeric(cadd_score)),
    cadd_bin = case_when(
      cadd_score > 20 ~ "CADD > 20 (n)",
      cadd_score <= 20 ~ "CADD ≤ 20 (n)",
      TRUE ~ NA_character_
    ),
    alpha_missense_pathogenicity_class = str_to_title(
      str_replace_all(alpha_missense_pathogenicity_class, "_", " ")
    )
  ) %>%
  filter(!is.na(cadd_bin)) %>%
  # Count variants by AlphaMissense class and CADD bin
  group_by(alpha_missense_pathogenicity_class, cadd_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  # Reshape to wide format
  pivot_wider(
    names_from = cadd_bin,
    values_from = n,
    values_fill = 0
  ) %>%
  # Compute total column
  mutate(`Total (n)` = `CADD > 20 (n)` + `CADD ≤ 20 (n)`) %>%
  # Reorder AlphaMissense classes
  mutate(alpha_missense_pathogenicity_class = factor(
    alpha_missense_pathogenicity_class,
    levels = c("Likely Pathogenic", "Ambiguous", "Likely Benign")
  )) %>%
  arrange(alpha_missense_pathogenicity_class)

missense_summary %>%
  gt() %>%
  tab_header(
    title = md("**Table 6. Missense Variants by AlphaMissense Classification and CADD Score**")
  ) %>%
  cols_label(
    alpha_missense_pathogenicity_class = "AlphaMissense Classification",
    `CADD > 20 (n)` = html("CADD &gt; 20 (n)"),
    `CADD ≤ 20 (n)` = html("CADD ≤ 20 (n)"),
    `Total (n)` = "Total (n)"
  ) %>%
  fmt_number(columns = where(is.numeric), decimals = 0) %>%
  tab_options(
    table.font.size = 14,
    table.font.name = "Times New Roman",
    heading.title.font.size = 16,
    heading.align = "left",
    data_row.padding = px(6)
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  ) %>%
  tab_source_note(
    md("Distribution of unique germline missense variants according to AlphaMissense classification and CADD score. Variants were grouped by AlphaMissense-predicted pathogenicity class and stratified by CADD score thresholds (>20 vs. ≤20).")
  )





#========= Likely benign missense variants with CADD < 20
likely_benign_lowCADD_df <- rb_clean_unique_germline %>%
  filter(
    molecular_consequence == "missense",
    alpha_missense_pathogenicity_class == "likely_benign",
    suppressWarnings(as.numeric(cadd_score)) < 20
  ) %>%
  mutate(
    # Convert numeric fields
    cadd_num  = suppressWarnings(as.numeric(cadd_score)),
    am_prob   = suppressWarnings(as.numeric(alpha_missense_pathogenicity_score)),
    acc_g     = suppressWarnings(as.numeric(spliceai_acc_gain)),
    acc_l     = suppressWarnings(as.numeric(spliceai_acc_loss)),
    don_g     = suppressWarnings(as.numeric(spliceai_don_gain)),
    don_l     = suppressWarnings(as.numeric(spliceai_don_loss)),
    
    # Compute single max SpliceAI score
    spliceai_max = {
      v <- pmax(acc_g, acc_l, don_g, don_l, na.rm = TRUE)
      ifelse(is.infinite(v), NA, v)
    },
    
    # Clean formatting
    cadd_score = ifelse(is.na(cadd_num), "-", sub("\\.?0+$", "", format(cadd_num, trim = TRUE))),
    spliceai_max = ifelse(is.na(spliceai_max), "-", format(signif(spliceai_max, 2), trim = TRUE)),
    alpha_missense_pathogenicity_score = ifelse(is.na(am_prob), "-", format(round(am_prob, 2), nsmall = 2, trim = TRUE))
  ) %>%
  dplyr::select(
    cdna,
    protein_combined,
    alpha_missense_pathogenicity_score,
    cadd_score,
    spliceai_max,
    clinvar_variant_classification,
    mmcif_structure,
    mmcif_plddt,
    uniprot_regions
  ) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "-", .))) %>%
  arrange(desc(as.numeric(alpha_missense_pathogenicity_score))) %>%
  mutate(Count = row_number()) %>%
  dplyr::select(
    Count,
    cdna,
    protein_combined,
    alpha_missense_pathogenicity_score,
    cadd_score,
    spliceai_max,
    clinvar_variant_classification,
    mmcif_structure,
    mmcif_plddt,
    uniprot_regions
  )

likely_benign_lowCADD_df %>%
  gt() %>%
  tab_header(title = "Likely Benign Missense Variants with CADD < 20") %>%
  cols_label(
    Count = "Count",
    cdna = "cDNA",
    protein_combined = "Protein",
    alpha_missense_pathogenicity_score = "AlphaMissense Probability",
    cadd_score = "CADD",
    spliceai_max = "SpliceAI (max)",
    clinvar_variant_classification = "ClinVar Classification",
    mmcif_structure = "Structure (mmCIF)",
    mmcif_plddt = "Local Confidence",
    uniprot_regions = "UniProt Region"
  ) %>%
  fmt_number(
    columns = c(alpha_missense_pathogenicity_score, spliceai_max),
    decimals = 2
  ) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 16,
    heading.align = "center"
  )




#========= synonymous variants 


synonymous_variants_df <- rb_clean_unique_germline %>%
  filter(
    !is.na(molecular_consequence),
    molecular_consequence %in% c("synonymous_variant", "synonymous")
  ) %>%
  dplyr::select(
    cdna,
    protein_combined,
    clinvar_variant_classification,
    cadd_score,
    synvep_score,
    synvep_effect,
    gnomad_frequency,
    spliceai_max,
    dist2splice
  ) %>%
  mutate(
    cadd_score        = suppressWarnings(as.numeric(cadd_score)),
    synvep_score      = suppressWarnings(as.numeric(synvep_score)),
    gnomad_frequency  = suppressWarnings(as.numeric(gnomad_frequency)),
    spliceai_max      = suppressWarnings(as.numeric(spliceai_max)),
    dist2splice       = suppressWarnings(as.numeric(dist2splice))
  ) %>%
  arrange(desc(cadd_score))

n_syn <- nrow(synonymous_variants_df)

legend_text <- paste(
  "Entries classified as Pathogenic or Likely Pathogenic in ClinVar are in bold.",
  "CADD scores > 20 and synVEP scores > 0.5 are also bolded.",
  sep = " "
)

synonymous_variants_df %>%
  gt() %>%
  tab_header(
    title = md(glue("**Table 7. Summary of {n_syn} Unique Germline Synonymous Variants**"))
  ) %>%
  cols_label(
    cdna = "cDNA",
    protein_combined = "Protein",
    clinvar_variant_classification = "ClinVar Classification",
    cadd_score = "CADD",
    synvep_score = "synVEP Score",
    synvep_effect = "synVEP Effect",
    gnomad_frequency = "gnomAD Frequency",
    spliceai_max = "SpliceAI Max",
    dist2splice = "Distance to Splice Site"
  ) %>%
  fmt_missing(everything(), missing_text = "-") %>%
  fmt_number(columns = c(cadd_score, synvep_score, spliceai_max), decimals = 2) %>%
  fmt_number(columns = gnomad_frequency, decimals = 5) %>%
  fmt_number(columns = dist2splice, decimals = 0) %>%
  # Bold ClinVar P/LP
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = clinvar_variant_classification,
      rows = clinvar_variant_classification %in% c(
        "Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"
      )
    )
  ) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = cadd_score, rows = cadd_score > 20)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = synvep_score, rows = synvep_score > 0.5)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  # Left-align all columns
  cols_align(align = "left", columns = everything()) %>%
  tab_source_note(source_note = md(legend_text)) %>%
  tab_options(
    table.width = gt::pct(100),
    table.font.size = 14,
    table.font.name = "Times New Roman",
    heading.title.font.size = 16,
    heading.align = "left",
    data_row.padding = px(6)
  )
















common_limits <- c(0.5, 27.5)

# Exon lengths
exon_lengths <- read_excel("/Users/rooks/Desktop/Coding/VS\ Code/panorama_rb1/retinoblastoma_data_v1.xlsx", sheet = 2) %>%
  clean_names() %>%
  mutate(exon_number = as.integer(str_extract(exon, "\\d+"))) %>%
  distinct(exon_number, .keep_all = TRUE) %>%
  select(exon_number, length) %>%
  filter(!is.na(exon_number), !is.na(length), length > 0)

# Counts using ALL germline records (not deduplicated)
exon_counts_all <- rb_clean_germline %>%
  filter(str_detect(
    exon_intron_combined,
    regex("^\\s*exon\\s*\\d+\\s*$", ignore_case = TRUE)
  )) %>%
  mutate(exon_number = as.integer(str_extract(exon_intron_combined, "\\d+"))) %>%
  filter(!is.na(exon_number)) %>%
  group_by(exon_number) %>%
  summarise(variant_count_all = n(), .groups = "drop")

# Counts using UNIQUE cDNA changes
exon_counts_unique <- rb_clean_germline %>%
  filter(str_detect(
    exon_intron_combined,
    regex("^\\s*exon\\s*\\d+\\s*$", ignore_case = TRUE)
  )) %>%
  mutate(exon_number = as.integer(str_extract(exon_intron_combined, "\\d+"))) %>%
  filter(!is.na(exon_number)) %>%
  group_by(exon_number) %>%
  summarise(variant_count_unique = n_distinct(cdna), .groups = "drop")

# Rates per kb for both
exon_rates_all <- exon_lengths %>%
  left_join(exon_counts_all, by = "exon_number") %>%
  mutate(
    variant_count_all   = coalesce(variant_count_all, 0L),
    variants_per_kb_all = variant_count_all / (length / 1000)
  ) %>%
  arrange(exon_number)

exon_rates_unique <- exon_lengths %>%
  left_join(exon_counts_unique, by = "exon_number") %>%
  mutate(
    variant_count_unique   = coalesce(variant_count_unique, 0L),
    variants_per_kb_unique = variant_count_unique / (length / 1000)
  ) %>%
  arrange(exon_number)

#----------------- Bottom panel: horizontal protein domain line ----------------
# Domains:
#  N-terminal: exons 1–11
#  Pocket A : exons 11–18
#  Pocket B : exons 19–22
#  C-term   : exons 22–27
# Outer endpoints 0.5 and 27.5; A/B boundary at 18.5

domain_line <- tibble::tibble(
  domain = c(
    "N-terminal domain",
    "Pocket domain A",
    "Pocket domain B",
    "C-terminal domain"
  ),
  xstart = c(0.5, 11, 18.5, 22),
  xend   = c(11,   18.5, 22,   27.5)
)

# Vertical boundaries at 0.5, 11, 18.5, 22, 27.5
domain_boundaries <- tibble::tibble(
  x = c(0.5, 11, 18.5, 22, 27.5)
)

# Label positions (midpoints of each segment) - placed just below the line
domain_labels <- domain_line %>%
  mutate(
    x = (xstart + xend) / 2,
    label = case_when(
      domain == "N-terminal domain" ~ "N-terminal\ndomain",
      domain == "Pocket domain A"   ~ "Pocket\ndomain A",
      domain == "Pocket domain B"   ~ "Pocket\ndomain B",
      domain == "C-terminal domain" ~ "C-terminal\ndomain",
      TRUE                          ~ domain
    ),
    y = -0.02
  )

p_domains <- ggplot() +
  geom_segment(
    data = domain_line,
    aes(x = xstart, xend = xend, y = 0.1, yend = 0.1),
    linewidth = 0.8
  ) +
  geom_segment(
    data = domain_boundaries,
    aes(x = x, xend = x, y = 0.15, yend = 0.05),
    linewidth = 0.6
  ) +
  geom_text(
    data = domain_labels,
    aes(x = x, y = y, label = label),
    family = "Times New Roman",
    size = 3.5,
    lineheight = 0.9
    
  ) +
  coord_cartesian(xlim = common_limits, ylim = c(-0.2, 0.1), clip = "off") +
  scale_x_continuous(
    limits = common_limits,
    breaks = 1:27,
    labels = 1:27,
    expand = c(0, 0)
  ) +
  theme_void() +
  theme(
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
  )

#----------------- Top and middle panels: density plots ------------------------

# Panel: all germline reports (tag a)
p2a <- ggplot(exon_rates_all,
              aes(x = exon_number, y = variants_per_kb_all)) +
  geom_col(width = 0.6, fill = "gray70") +
  scale_x_continuous(
    limits = common_limits,
    breaks = 1:27,
    labels = 1:27,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1850),
    breaks = seq(0, 1500, by = 500),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Exon",
    y = "Variants per kb",
    tag = "a"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 5)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# Panel: unique germline variants (tag b)
p2b <- ggplot(exon_rates_unique,
              aes(x = exon_number, y = variants_per_kb_unique)) +
  geom_col(width = 0.6, fill = "gray40") +
  scale_x_continuous(
    limits = common_limits,
    breaks = 1:27,
    labels = 1:27,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1850),
    breaks = seq(0, 1500, by = 500),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Exon",
    y = "Variants per kb",
    tag = "b"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 5)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#----------------- Caption and combined figure ---------------------------------

caption_text <- str_wrap(
  "Figure 3. Germline RB1 variant density by exon shown as variants per kilobase, with exons mapped to pRB domains a) Density of all germline variants by exon b) Variant density of unique germline variants by exon",
  width = 300
)

(p2a / p2b / p_domains) +
  plot_layout(heights = c(1, 1, 0.3)) +
  plot_annotation(
    caption = caption_text,
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 14,
        family = "Times New Roman"
      ),
      plot.caption = element_textbox_simple(
        size = 12,
        family = "Times New Roman",
        width = unit(9, "in"),
        margin = margin(t = 10),
        hjust = 0
      )
    )
  )




#========= laterality counts by variants


variant_laterality_counts <- rb_clean_germline %>%
  filter(laterality %in% c("bilateral", "unilateral")) %>%
  count(cdna, laterality, name = "n") %>%
  pivot_wider(
    names_from = laterality,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    both_present = bilateral > 0 & unilateral > 0
  ) %>%
  left_join(
    rb %>%
      dplyr::select(cdna, cadd_score, molecular_consequence, exon_intron_combined) %>%
      distinct(cdna, .keep_all = TRUE),
    by = "cdna"
  ) %>%
  arrange(desc(bilateral + unilateral))



###INDIVIDUAL VARIANT LATERALITY ANALYSIS (MULTIPLE COMPARISONS CORRECTED)

dat <- rb_clean_germline %>%
  filter(laterality %in% c("bilateral", "unilateral")) %>%
  select(cdna, laterality)

total_bilat  <- sum(dat$laterality == "bilateral")
total_unilat <- sum(dat$laterality == "unilateral")

res <- dat %>%
  group_by(cdna) %>%
  summarise(
    bilateral  = sum(laterality == "bilateral"),
    unilateral = sum(laterality == "unilateral"),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    ft = list(stats::fisher.test(
      matrix(c(bilateral, unilateral,
               total_bilat - bilateral,
               total_unilat - unilateral),
             nrow = 2, byrow = TRUE)
    )),
    fisher_p   = ft$p.value,
    odds_ratio = unname(ft$estimate),
    conf_low   = ft$conf.int[1],
    conf_high  = ft$conf.int[2]
  ) %>%
  ungroup() %>%
  select(-ft) %>%
  mutate(
    p_bonf = p.adjust(fisher_p, method = "bonferroni")
  ) %>%
  arrange(p_bonf, fisher_p)

head(res, 20)


#Laterality vs fh chi squared

rb_subset_laterality_family <- rb_clean_germline %>%
  filter(
    laterality %in% c("unilateral", "bilateral"),
    family_history %in% c("isolated", "familial")
  )

contingency_laterality_family <- table(
  rb_subset_laterality_family$laterality,
  rb_subset_laterality_family$family_history
)

contingency_laterality_family

chisq_laterality_family <- chisq.test(contingency_laterality_family)

chisq_laterality_family



###MISSENSE IN CLINVAR

rb_missense_unique <- rb_clean_unique_germline %>%
  filter(str_detect(molecular_consequence, "missense"))

clinvar_breakdown_missense <- rb_missense_unique %>%
  mutate(
    clinvar_variant_classification = if_else(
      is.na(clinvar_variant_classification) | clinvar_variant_classification == "",
      "Unannotated",
      clinvar_variant_classification
    )
  ) %>%
  count(clinvar_variant_classification, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(prop = n / sum(n)) %>%
  bind_rows(
    tibble(
      clinvar_variant_classification = "Total",
      n = sum(.$n),
      prop = 1
    )
  )

clinvar_breakdown_missense





### FIGURES

# ----------------------------- PANEL A: exon density --------------------------
common_limits <- c(0.5, 27.5)

exon_lengths <- read_excel("/Users/rooks/Desktop/Coding/VS\ Code/panorama_rb1/retinoblastoma_data_v1.xlsx", sheet = 2) %>%
  clean_names() %>%
  mutate(exon_number = as.integer(str_extract(exon, "\\d+"))) %>%
  distinct(exon_number, .keep_all = TRUE) %>%
  select(exon_number, length) %>%
  filter(!is.na(exon_number), !is.na(length), length > 0)

exon_counts_all <- rb_clean_germline %>%
  filter(str_detect(exon_intron_combined, regex("^\\s*exon\\s*\\d+\\s*$", ignore_case = TRUE))) %>%
  mutate(exon_number = as.integer(str_extract(exon_intron_combined, "\\d+"))) %>%
  filter(!is.na(exon_number)) %>%
  count(exon_number, name = "variant_count_all")

exon_counts_unique <- rb_clean_germline %>%
  filter(str_detect(exon_intron_combined, regex("^\\s*exon\\s*\\d+\\s*$", ignore_case = TRUE))) %>%
  mutate(exon_number = as.integer(str_extract(exon_intron_combined, "\\d+"))) %>%
  filter(!is.na(exon_number)) %>%
  group_by(exon_number) %>%
  summarise(variant_count_unique = n_distinct(cdna), .groups = "drop")

exon_rates <- exon_lengths %>%
  left_join(exon_counts_all, by = "exon_number") %>%
  left_join(exon_counts_unique, by = "exon_number") %>%
  mutate(
    variant_count_all    = coalesce(variant_count_all, 0L),
    variant_count_unique = coalesce(variant_count_unique, 0L),
    variants_per_kb_all    = variant_count_all / (length / 1000),
    variants_per_kb_unique = variant_count_unique / (length / 1000)
  ) %>%
  arrange(exon_number)

plot_dat_a <- exon_rates %>%
  select(exon_number, variants_per_kb_all, variants_per_kb_unique) %>%
  pivot_longer(
    cols = c(variants_per_kb_all, variants_per_kb_unique),
    names_to = "series",
    values_to = "variants_per_kb"
  ) %>%
  mutate(
    series = recode(series,
                    variants_per_kb_all    = "All variants",
                    variants_per_kb_unique = "Unique variants"),
    series = factor(series, levels = c("All variants", "Unique variants")),
    # keep BOTH red + blue bars fully saturated for exons 4, 14, and 15
    alpha_val = if_else(exon_number %in% c(4, 14, 15), 1.00, 0.70)
  )

p_a <- ggplot(
  plot_dat_a,
  aes(x = exon_number, y = variants_per_kb, fill = series #, alpha = alpha_val
      )
) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7)) +
  scale_x_continuous(limits = common_limits, breaks = 1:27, labels = 1:27, expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 1850),
    breaks = seq(0, 1500, by = 500),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c("All variants" = "#D14625", "Unique variants" = "#3D6C88"),
    name = NULL
  ) +
  scale_alpha_identity(guide = "none") +
  labs(x = "Exon", y = "Variants per Kilobase", tag = "a") +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14)
  )

# ----------------------------- PANEL A bottom: domain bar ---------------------
domain_line <- tibble(
  domain = c("N-terminal domain", "Pocket domain A", "Pocket domain B", "C-terminal domain"),
  xstart = c(0.5, 11, 18.5, 22),
  xend   = c(11, 18.5, 22, 27.5)
)

domain_boundaries <- tibble(x = c(0.5, 11, 18.5, 22, 27.5))

domain_labels <- domain_line %>%
  mutate(
    x = (xstart + xend) / 2,
    label = case_when(
      domain == "N-terminal domain" ~ "N-terminal\ndomain",
      domain == "Pocket domain A"   ~ "Pocket\ndomain A",
      domain == "Pocket domain B"   ~ "Pocket\ndomain B",
      domain == "C-terminal domain" ~ "C-terminal\ndomain",
      TRUE ~ domain
    ),
    y = -0.02
  )

p_domains <- ggplot() +
  geom_segment(data = domain_line,
               aes(x = xstart, xend = xend, y = 0.1, yend = 0.1),
               linewidth = 0.8) +
  geom_segment(data = domain_boundaries,
               aes(x = x, xend = x, y = 0.15, yend = 0.05),
               linewidth = 0.6) +
  geom_text(data = domain_labels,
            aes(x = x, y = y, label = label),
            family = "Times New Roman", size = 3.5, lineheight = 0.9) +
  coord_cartesian(xlim = common_limits, ylim = c(-0.2, 0.1), clip = "off") +
  scale_x_continuous(limits = common_limits, breaks = 1:27, labels = 1:27, expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 5, b = 5, l = 5))

panel_a <- wrap_plots(p_a, p_domains, ncol = 1, heights = c(1, 0.28))



# ----------------------------- PANEL B: consequence density by domain ---------
domain_lengths <- tibble(
  uniprot_regions = c("N-terminal domain", "Pocket Domain A", "Pocket Domain B", "C-terminal domain"),
  domain_len_aa   = c(372, 207, 132, 157)
)

bar_dat <- rb_clean_unique_germline %>%
  filter(!is.na(uniprot_regions), str_squish(uniprot_regions) != "") %>%
  transmute(
    uniprot_regions = str_squish(uniprot_regions),
    molecular_consequence = str_squish(molecular_consequence),
    cdna = cdna
  ) %>%
  filter(molecular_consequence %in% c("stop_gained", "missense", "synonymous")) %>%
  distinct(uniprot_regions, molecular_consequence, cdna) %>%
  count(uniprot_regions, molecular_consequence, name = "n_unique") %>%
  left_join(domain_lengths, by = "uniprot_regions") %>%
  mutate(
    v_per_aa = n_unique / domain_len_aa,
    uniprot_regions = fct_relevel(
      uniprot_regions,
      "N-terminal domain", "Pocket Domain A", "Pocket Domain B", "C-terminal domain"
    ),
    consequence_lab = recode(
      molecular_consequence,
      "synonymous"  = "Synonymous",
      "missense"    = "Missense",
      "stop_gained" = "Stop gained"
    ),
    consequence_lab = factor(consequence_lab, levels = c("Synonymous", "Missense", "Stop gained")),
    label = sprintf("%.3f", v_per_aa)
  )

ymax_b <- max(bar_dat$v_per_aa, na.rm = TRUE) * 1.12

p_b <- ggplot(bar_dat, aes(x = uniprot_regions, y = v_per_aa, fill = consequence_lab)) +
  geom_col(
    width = 0.75,
    position = position_dodge(width = 0.8),
    color = NA
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.8),
    vjust = -0.35,
    family = "Times New Roman",
    size = 3.6
  ) +
  scale_fill_manual(
    name = "",
    values = c("Synonymous" = "#D8E1E8", "Missense" = "#8AA7B8", "Stop gained" = "#3D6C88")
  ) +
  scale_y_continuous(
    limits = c(0, ymax_b),
    expand = c(0, 0)
  ) +
  labs(
    x = "",
    y = "Variants per Amino Acid",
    tag = "b"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 12)
  )

caption_text <- str_wrap(
  "Figure 3. (a) Germline RB1 variant density across exon normalized by exon length; red bars represent all variants and blue bars represent unique variants. (b) Unique germline variants by each SNV consequence in each pRB domain, normalized by protein domain length.",
  width = 300
)

final_fig <- wrap_plots(panel_a, p_b, ncol = 1, heights = c(1.25, 1)) +
  plot_annotation(
    caption = caption_text,
    theme = theme(
      plot.caption = element_textbox_simple(
        size = 12,
        family = "Times New Roman",
        width = unit(9.8, "in"),
        margin = margin(t = 10),
        hjust = 0
      )
    )
  )

final_fig


## confidence intervals #========= 95% CI for bilateral / unilateral proportions by genetic origin -----

prop_ci <- function(x, n, conf.level = 0.95) {
  bt <- binom.test(x, n, conf.level = conf.level)
  tibble(
    count = x,
    total = n,
    percent = 100 * x / n,
    ci_low = 100 * bt$conf.int[1],
    ci_high = 100 * bt$conf.int[2]
  )
}

fmt_pct_ci <- function(p, lo, hi, digits = 1) {
  sprintf(
    paste0("%.", digits, "f%% (95%% CI, %.", digits, "f–%.", digits, "f%%)"),
    p, lo, hi
  )
}

laterality_counts <- rb %>%
  dplyr::filter(
    genetic_origin %in% c("germline", "somatic"),
    laterality %in% c("bilateral", "unilateral")
  ) %>%
  dplyr::count(genetic_origin, laterality, name = "n") %>%
  tidyr::pivot_wider(
    names_from = laterality,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    total = bilateral + unilateral
  )

laterality_ci <- purrr::map_dfr(
  seq_len(nrow(laterality_counts)),
  function(i) {
    row <- laterality_counts[i, ]
    
    bilat_ci <- prop_ci(row$bilateral, row$total)
    unilat_ci <- prop_ci(row$unilateral, row$total)
    
    tibble(
      genetic_origin = row$genetic_origin,
      bilateral_n = row$bilateral,
      unilateral_n = row$unilateral,
      total_n = row$total,
      bilateral_percent = round(bilat_ci$percent, 1),
      bilateral_ci_low = round(bilat_ci$ci_low, 1),
      bilateral_ci_high = round(bilat_ci$ci_high, 1),
      unilateral_percent = round(unilat_ci$percent, 1),
      unilateral_ci_low = round(unilat_ci$ci_low, 1),
      unilateral_ci_high = round(unilat_ci$ci_high, 1),
      bilateral_label = fmt_pct_ci(
        round(bilat_ci$percent, 1),
        round(bilat_ci$ci_low, 1),
        round(bilat_ci$ci_high, 1),
        digits = 1
      ),
      unilateral_label = fmt_pct_ci(
        round(unilat_ci$percent, 1),
        round(unilat_ci$ci_low, 1),
        round(unilat_ci$ci_high, 1),
        digits = 1
      )
    )
  }
)

laterality_ci


#### CODE TO ENUMERATE THE SNVS POSSIBLE: 

library(dplyr)
library(tidyr)
library(tibble)
library(Biostrings)

# ---------- INPUT: raw CDS string ----------
cds <- "ATGCCGCCCAAAACCCCCCGAAAAACGGCCGCCACCGCCGCCGCTGCCGCCGCGGAACCCCCGGCACCGCCGCCGCCGCCCCCTCCTGAGGAGGACCCAGAGCAGGACAGCGGCCCGGAGGACCTGCCTCTCGTCAG

GCT
TGAGTTTGAAGAAACAGAAGAACCTGATTTTACTGCATTATGTCAGAAATTAAAGATACCAGATCATGTC
AGAGAGAGAGCTTGGTTAACTTGGGAGAAAGTTTCATCTGTGGATGGAGTATTG

GGAGGTTATATTCAAAAGAAAAAGGAACTGTGGGGAATCTGTATCTTTATTGCAGCAGTTGACCTAGATGAGATGTCGTTCACTTT
TACTGAGCTACAGAAAAACATAGAAATCAGTGTCCATAAATTCTTTAACTTACTAAAAGAAATTGATACC
AGTACCAAAGTTGATAATGCTATGTCAAGACTGTTGAAGAAGTATGATGTATTGTTTGCACTCTTCAGCA
AATTGGAAAGGACATGTGAACTTATATATTTGACACAACCCAGCAGTTCGATATCTACTGAAATAAATTC
TGCATTGGTGCTAAAAGTTTCTTGGATCACATTTTTATTAGCTAAAGGGGAAGTATTACAAATGGAAGAT
GATCTGGTGATTTCATTTCAGTTAATGCTATGTGTCCTTGACTATTTTATTAAACTCTCACCTCCCATGT
TGCTCAAAGAACCATATAAAACAGCTGTTATACCCATTAATGGTTCACCTCGAACACCCAGGCGAGGTCA
GAACAGGAGTGCACGGATAGCAAAACAACTAGAAAATGATACAAGAATTATTGAAGTTCTCTGTAAAGAA
CATGAATGTAATATAGATGAGGTGAAAAATGTTTATTTCAAAAATTTTATACCTTTTATGAATTCTCTTG
GACTTGTAACATCTAATGGACTTCCAGAGGTTGAAAATCTTTCTAAACGATACGAAGAAATTTATCTTAA
AAATAAAGATCTAGATGCAAGATTATTTTTGGATCATGATAAAACTCTTCAGACTGATTCTATAGACAGT
TTTGAAACACAGAGAACACCACGAAAAAGTAACCTTGATGAAGAGGTGAATGTAATTCCTCCACACACTC
CAGTTAGGACTGTTATGAACACTATCCAACAATTAATGATGATTTTAAATTCAGCAAGTGATCAACCTTC
AGAAAATCTGATTTCCTATTTTAACAACTGCACAGTGAATCCAAAAGAAAGTATACTGAAAAGAGTGAAG
GATATAGGATACATCTTTAAAGAGAAATTTGCTAAAGCTGTGGGACAGGGTTGTGTCGAAATTGGATCAC
AGCGATACAAACTTGGAGTTCGCTTGTATTACCGAGTAATGGAATCCATGCTTAAATCAGAAGAAGAACG
ATTATCCATTCAAAATTTTAGCAAACTTCTGAATGACAACATTTTTCATATGTCTTTATTGGCGTGCGCT
CTTGAGGTTGTAATGGCCACATATAGCAGAAGTACATCTCAGAATCTTGATTCTGGAACAGATTTGTCTT
TCCCATGGATTCTGAATGTGCTTAATTTAAAAGCCTTTGATTTTTACAAAGTGATCGAAAGTTTTATCAA
AGCAGAAGGCAACTTGACAAGAGAAATGATAAAACATTTAGAACGATGTGAACATCGAATCATGGAATCC
CTTGCATGGCTCTCAGATTCACCTTTATTTGATCTTATTAAACAATCAAAGGACCGAGAAGGACCAACTG
ATCACCTTGAATCTGCTTGTCCTCTTAATCTTCCTCTCCAGAATAATCACACTGCAGCAGATATGTATCT
TTCTCCTGTAAGATCTCCAAAGAAAAAAGGTTCAACTACGCGTGTAAATTCTACTGCAAATGCAGAGACA
CAAGCAACCTCAGCCTTCCAGACCCAGAAGCCATTGAAATCTACCTCTCTTTCACTGTTTTATAAAAAAG
TGTATCGGCTAGCCTATCTCCGGCTAAATACACTTTGTGAACGCCTTCTGTCTGAGCACCCAGAATTAGA
ACATATCATCTGGACCCTTTTCCAGCACACCCTGCAGAATGAGTATGAACTCATGAGAGACAGGCATTTG
GACCAAATTATGATGTGTTCCATGTATGGCATATGCAAAGTGAAGAATATAGACCTTAAATTCAAAATCA
TTGTAACAGCATACAAGGATCTTCCTCATGCTGTTCAGGAGACATTCAAACGTGTTTTGATCAAAGAAGA
GGAGTATGATTCTATTATAGTATTCTATAACTCGGTCTTCATGCAGAGACTGAAAACAAATATTTTGCAG
TATGCTTCCACCAGGCCCCCTACCTTGTCACCAATACCTCACATTCCTCGAAGCCCTTACAAGTTTCCTA
GTTCACCCTTACGGATTCCTGGAGGGAACATCTATATTTCACCCCTGAAGAGTCCATATAAAATTTCAGA
AGGTCTGCCAACACCAACAAAAATGACTCCAAGATCAAGAATCTTAGTATCAATTGGTGAATCATTCGGG
ACTTCTGAGAAGTTCCAGAAAATAAATCAGATGGTATGTAACAGCGACCGTGTGCTCAAAAGAAGTGCTG
AAGGAAGCAACCCTCCTAAACCACTGAAAAAACTACGCTTTGATATTGAAGGATCAGATGAAGCAGATGG
AAGTAAACATCTCCCAGGAGAGTCCAAATTTCAGCAGAAACTGGCAGAAATGACTTCTACTCGAACACGA
ATGCAAAAGCAGAAAATGAATGATAGCATGGATACCTCAAACAAGGAAGAGAAA"

# clean
cds <- toupper(gsub("[^ACGT]", "", cds))

cat("CDS length:", nchar(cds), "\n")
cat("Divisible by 3:", nchar(cds) %% 3 == 0, "\n")

# ---------- SETUP ----------
bases <- c("A","C","G","T")
code  <- Biostrings::GENETIC_CODE

codons <- substring(cds,
                    seq(1, nchar(cds), by = 3),
                    seq(3, nchar(cds), by = 3))

aa_ref <- unname(code[codons])

# ---------- ENUMERATION ----------
per_codon <- lapply(seq_along(codons), function(i) {
  
  nts <- strsplit(codons[i], "")[[1]]
  aa0 <- aa_ref[i]
  
  do.call(rbind, lapply(1:3, function(pos) {
    
    alts <- setdiff(bases, nts[pos])
    
    do.call(rbind, lapply(alts, function(nb) {
      
      mut <- nts
      mut[pos] <- nb
      mut_codon <- paste0(mut, collapse = "")
      
      aa_mut <- unname(code[mut_codon])
      
      outcome <- if (aa_mut == "*") "nonsense"
      else if (aa_mut == aa0) "synonymous"
      else "missense"
      
      data.frame(
        codon_index = i,
        pos = pos,
        outcome = outcome,
        stringsAsFactors = FALSE
      )
    }))
  }))
})

all_variants <- do.call(rbind, per_codon)

summary_counts <- table(all_variants$outcome)
print(summary_counts)

cat("Number of codons:", length(codons), "\n")
cat("Total SNVs:", nrow(all_variants), "\n")


# ---------- Possible missense / nonsense / synonymous mutations by protein domain ----------

domain_map <- tibble(
  domain   = c("N-terminal domain", "Pocket Domain A", "Pocket Domain B", "C-terminal domain"),
  aa_start = c(1, 373, 640, 772),
  aa_stop  = c(372, 579, 771, 928)
)

domain_variant_table <- all_variants %>%
  mutate(aa_pos = codon_index) %>%
  rowwise() %>%
  mutate(
    domain = domain_map$domain[
      which(aa_pos >= domain_map$aa_start & aa_pos <= domain_map$aa_stop)[1]
    ]
  ) %>%
  ungroup() %>%
  filter(!is.na(domain)) %>%
  count(domain, outcome, name = "n") %>%
  pivot_wider(
    names_from = outcome,
    values_from = n,
    values_fill = 0
  ) %>%
  left_join(domain_map, by = "domain") %>%
  select(domain, aa_start, aa_stop, missense, nonsense, synonymous) %>%
  arrange(aa_start)

total_row <- domain_variant_table %>%
  summarise(
    domain = "Total",
    aa_start = NA_integer_,
    aa_stop = NA_integer_,
    missense = sum(missense, na.rm = TRUE),
    nonsense = sum(nonsense, na.rm = TRUE),
    synonymous = sum(synonymous, na.rm = TRUE)
  )

domain_variant_table <- bind_rows(domain_variant_table, total_row)

domain_variant_table_updated <- domain_variant_table %>%
  filter(domain != "Total") %>%
  mutate(
    domain_length = aa_stop - aa_start + 1L,
    Total = missense + nonsense + synonymous
  ) %>%
  mutate(
    pct_missense = missense / Total,
    pct_nonsense = nonsense / Total,
    pct_synonymous = synonymous / Total,
    pct_total = 1
  )

total_row <- domain_variant_table_updated %>%
  summarise(
    domain = "Total",
    domain_length = sum(domain_length, na.rm = TRUE),
    missense = sum(missense, na.rm = TRUE),
    nonsense = sum(nonsense, na.rm = TRUE),
    synonymous = sum(synonymous, na.rm = TRUE),
    Total = sum(Total, na.rm = TRUE)
  ) %>%
  mutate(
    pct_missense = missense / Total,
    pct_nonsense = nonsense / Total,
    pct_synonymous = synonymous / Total,
    pct_total = 1
  )

domain_variant_table_updated <- bind_rows(domain_variant_table_updated, total_row) %>%
  mutate(
    missense = sprintf("%d (%.1f)", missense, 100 * pct_missense),
    nonsense = sprintf("%d (%.1f)", nonsense, 100 * pct_nonsense),
    synonymous = sprintf("%d (%.1f)", synonymous, 100 * pct_synonymous),
    Total = sprintf("%d (%.1f)", Total, 100 * pct_total)
  ) %>%
  select(domain, domain_length, missense, nonsense, synonymous, Total)

domain_variant_table_updated

domain_variant_table_updated %>%
  mutate(
    line = sprintf(
      "%-20s  (n=%3d aa)  Missense: %-15s  Nonsense: %-15s  Synonymous: %-15s  Total: %s",
      domain, domain_length, missense, nonsense, synonymous, Total
    )
  ) %>%
  pull(line) %>%
  cat(sep = "\n")







#### SUPPLEMENTARY TABLES AND ANALYSES CODE
### code depends on enumeration above


library(tidyverse)
library(janitor)
library(readxl)
library(broom)
library(gt)
library(writexl)
library(glue)
library(scales)
library(rlang)
library(pwr)
library(ggrepel)
library(ggpattern)
library(colorspace)
library(htmltools)
library(patchwork)
library(ggtext)
library(grid)


###LOAD DATA
fpath <- "/Users/rooks/Desktop/Coding/VS\ Code/panorama_rb1/retinoblastoma_data_v1.xlsx"
rb <- read_excel(fpath, sheet = 1) %>% clean_names()

rb_clean_germline <- rb %>% filter(genetic_origin == "germline")

rb_clean_unique_germline <- rb_clean_germline %>%
  distinct(cdna, .keep_all = TRUE)


#========= SUPPLEMENTAL TABLE 1. ClinVar and LOVD classification by variant consequence

five <- c("Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic")

rb_class_by_consequence <- rb_clean_unique_germline %>%
  dplyr::filter(
    !is.na(molecular_consequence),
    molecular_consequence != ""
  ) %>%
  dplyr::mutate(
    consequence_raw = stringr::str_trim(as.character(molecular_consequence)),
    consequence_binned = dplyr::case_when(
      consequence_raw %in% c("splice_donor", "splice_acceptor", "splice_site") ~ "splice_site",
      consequence_raw %in% c("inframe_deletion", "inframe_insertion", "inframe_indel", "inframe_variant") ~ "inframe_variant",
      TRUE ~ consequence_raw
    )
  ) %>%
  dplyr::filter(
    consequence_binned %in% c(
      "stop_gained", "frameshift", "splice_site", "missense",
      "intron", "inframe_variant", "synonymous", "5_prime_UTR", "start_lost"
    )
  )

make_class_block <- function(data, class_var, section_label) {
  
  dat <- data %>%
    dplyr::mutate(
      classification = dplyr::case_when(
        !!rlang::sym(class_var) %in% c("Likely benign", "likely benign") ~ "Likely benign",
        !!rlang::sym(class_var) %in% c("Likely pathogenic", "likely pathogenic") ~ "Likely pathogenic",
        !!rlang::sym(class_var) %in% c("Benign", "Benign/Likely benign", "benign") ~ "Benign",
        !!rlang::sym(class_var) %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "pathogenic") ~ "Pathogenic",
        !!rlang::sym(class_var) %in% c(
          "Uncertain significance",
          "Conflicting classifications of pathogenicity",
          "VUS"
        ) ~ "Uncertain significance",
        TRUE ~ as.character(!!rlang::sym(class_var))
      )
    ) %>%
    dplyr::filter(
      !is.na(classification),
      classification %in% five
    )
  
  cm <- dat %>%
    dplyr::count(classification, consequence_binned, name = "n") %>%
    tidyr::pivot_wider(
      names_from = consequence_binned,
      values_from = n,
      values_fill = 0
    ) %>%
    dplyr::mutate(
      classification = factor(classification, levels = five)
    ) %>%
    dplyr::arrange(classification)
  
  consequence_cols <- c(
    "stop_gained", "frameshift", "splice_site", "missense",
    "intron", "inframe_variant", "synonymous", "5_prime_UTR", "start_lost"
  )
  
  for (col in consequence_cols) {
    if (!col %in% names(cm)) {
      cm[[col]] <- 0
    }
  }
  
  cm <- cm %>%
    dplyr::select(classification, dplyr::all_of(consequence_cols))
  
  row_totals <- rowSums(cm[, consequence_cols, drop = FALSE], na.rm = TRUE)
  col_totals <- colSums(cm[, consequence_cols, drop = FALSE], na.rm = TRUE)
  grand_total <- sum(row_totals, na.rm = TRUE)
  
  cm <- cm %>%
    dplyr::mutate(
      Section = section_label,
      Total = row_totals
    ) %>%
    dplyr::rename(Classification = classification) %>%
    dplyr::select(Section, Classification, dplyr::all_of(consequence_cols), Total)
  
  total_row <- tibble::tibble(
    Section = section_label,
    Classification = "Total",
    stop_gained = col_totals["stop_gained"],
    frameshift = col_totals["frameshift"],
    splice_site = col_totals["splice_site"],
    missense = col_totals["missense"],
    intron = col_totals["intron"],
    inframe_variant = col_totals["inframe_variant"],
    synonymous = col_totals["synonymous"],
    `5_prime_UTR` = col_totals["5_prime_UTR"],
    start_lost = col_totals["start_lost"],
    Total = grand_total
  )
  
  dplyr::bind_rows(cm, total_row)
}

tbl_clinvar <- make_class_block(
  data = rb_class_by_consequence %>%
    dplyr::filter(
      !is.na(clinvar_variant_classification),
      clinvar_variant_classification != ""
    ),
  class_var = "clinvar_variant_classification",
  section_label = "ClinVar Classification"
)

tbl_lovd <- make_class_block(
  data = rb_class_by_consequence %>%
    dplyr::filter(
      !is.na(clinical_classification_combined),
      clinical_classification_combined != ""
    ),
  class_var = "clinical_classification_combined",
  section_label = "LOVD Classification"
)

# variants unclassified by BOTH ClinVar and LOVD
unclassified_dat <- rb_class_by_consequence %>%
  dplyr::filter(
    (is.na(clinvar_variant_classification) | stringr::str_trim(clinvar_variant_classification) == "") &
      (is.na(clinical_classification_combined) | stringr::str_trim(clinical_classification_combined) == "")
  )

unclassified_row <- unclassified_dat %>%
  dplyr::count(consequence_binned, name = "n") %>%
  tidyr::pivot_wider(
    names_from = consequence_binned,
    values_from = n,
    values_fill = 0
  )

consequence_cols <- c(
  "stop_gained", "frameshift", "splice_site", "missense",
  "intron", "inframe_variant", "synonymous", "5_prime_UTR", "start_lost"
)

for (col in consequence_cols) {
  if (!col %in% names(unclassified_row)) {
    unclassified_row[[col]] <- 0
  }
}

unclassified_row <- unclassified_row %>%
  dplyr::select(dplyr::all_of(consequence_cols))

unclassified_row <- tibble::tibble(
  Section = "",
  Classification = "Unclassified",
  stop_gained = unclassified_row$stop_gained,
  frameshift = unclassified_row$frameshift,
  splice_site = unclassified_row$splice_site,
  missense = unclassified_row$missense,
  intron = unclassified_row$intron,
  inframe_variant = unclassified_row$inframe_variant,
  synonymous = unclassified_row$synonymous,
  `5_prime_UTR` = unclassified_row$`5_prime_UTR`,
  start_lost = unclassified_row$start_lost,
  Total = sum(unclassified_row[1, consequence_cols], na.rm = TRUE)
)

tbl_out <- dplyr::bind_rows(
  tbl_clinvar,
  tbl_lovd,
  unclassified_row
)

gt(tbl_out, groupname_col = "Section") %>%
  tab_header(
    title = md("**Supplemental Table 1. ClinVar and LOVD Classification by Variant Consequence**")
  ) %>%
  cols_label(
    Classification = "Classification",
    stop_gained = "Stop gained",
    frameshift = "Frameshift",
    splice_site = "Splice site",
    missense = "Missense",
    intron = "Intron",
    inframe_variant = "Inframe variant",
    synonymous = "Synonymous",
    `5_prime_UTR` = "5′ UTR",
    start_lost = "Start lost",
    Total = "Total"
  ) %>%
  tab_spanner(
    label = md("**Variant Consequence**"),
    columns = c(
      stop_gained, frameshift, splice_site, missense,
      intron, inframe_variant, synonymous, `5_prime_UTR`, start_lost
    )
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  ) %>%
  opt_table_font(font = "Times New Roman") %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = Classification %in% c("Total", "Unclassified"))
  ) %>%
  tab_options(
    heading.align = "left",
    table.font.size = 12,
    heading.title.font.size = 12,
    heading.title.font.weight = "bold",
    table.border.top.color = "grey",
    table.border.bottom.color = "grey",
    data_row.padding = px(4)
  ) %>%
  tab_source_note(
    md("Counts reflect unique germline variants. For ClinVar, Pathogenic/Likely pathogenic were grouped with Pathogenic, Benign/Likely benign were grouped with Benign, and conflicting classifications were grouped with Uncertain significance. For LOVD, VUS was grouped with Uncertain significance. The final Unclassified row includes variants lacking both a ClinVar classification and a LOVD classification.")
  )

#Remade Figure 3: 


# ----------------------------- PANEL A: exon density --------------------------
common_limits <- c(0.5, 27.5)

exon_lengths <- read_excel("/Users/rooks/Desktop/Coding/VS\ Code/panorama_rb1/retinoblastoma_data_v1.xlsx", sheet = 2) %>%
  clean_names() %>%
  mutate(exon_number = as.integer(str_extract(exon, "\\d+"))) %>%
  distinct(exon_number, .keep_all = TRUE) %>%
  select(exon_number, length) %>%
  filter(!is.na(exon_number), !is.na(length), length > 0)

exon_counts_all <- rb_clean_germline %>%
  filter(str_detect(exon_intron_combined, regex("^\\s*exon\\s*\\d+\\s*$", ignore_case = TRUE))) %>%
  mutate(exon_number = as.integer(str_extract(exon_intron_combined, "\\d+"))) %>%
  filter(!is.na(exon_number)) %>%
  count(exon_number, name = "variant_count_all")

exon_counts_unique <- rb_clean_germline %>%
  filter(str_detect(exon_intron_combined, regex("^\\s*exon\\s*\\d+\\s*$", ignore_case = TRUE))) %>%
  mutate(exon_number = as.integer(str_extract(exon_intron_combined, "\\d+"))) %>%
  filter(!is.na(exon_number)) %>%
  group_by(exon_number) %>%
  summarise(variant_count_unique = n_distinct(cdna), .groups = "drop")

exon_rates <- exon_lengths %>%
  left_join(exon_counts_all, by = "exon_number") %>%
  left_join(exon_counts_unique, by = "exon_number") %>%
  mutate(
    variant_count_all    = coalesce(variant_count_all, 0L),
    variant_count_unique = coalesce(variant_count_unique, 0L),
    variants_per_kb_all    = variant_count_all / (length / 1000),
    variants_per_kb_unique = variant_count_unique / (length / 1000)
  ) %>%
  arrange(exon_number)

plot_dat_a <- exon_rates %>%
  select(exon_number, variants_per_kb_all, variants_per_kb_unique) %>%
  pivot_longer(
    cols = c(variants_per_kb_all, variants_per_kb_unique),
    names_to = "series",
    values_to = "variants_per_kb"
  ) %>%
  mutate(
    series = recode(series,
                    variants_per_kb_all    = "All variants",
                    variants_per_kb_unique = "Unique variants"),
    series = factor(series, levels = c("All variants", "Unique variants")),
    # keep BOTH red + blue bars fully saturated for exons 4, 14, and 15
    alpha_val = if_else(exon_number %in% c(4, 14, 15), 1.00, 0.70)
  )

p_a <- ggplot(
  plot_dat_a,
  aes(x = exon_number, y = variants_per_kb, fill = series #, alpha = alpha_val
  )
) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7)) +
  scale_x_continuous(limits = common_limits, breaks = 1:27, labels = 1:27, expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 1850),
    breaks = seq(0, 1500, by = 500),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c("All variants" = "#D14625", "Unique variants" = "#3D6C88"),
    name = NULL
  ) +
  scale_alpha_identity(guide = "none") +
  labs(x = "Exon", y = "Variants per kilobase", tag = "a") +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14)
  )

# ----------------------------- PANEL A bottom: domain bar ---------------------
domain_line <- tibble(
  domain = c("N-terminal domain", "Pocket domain A", "Pocket domain B", "C-terminal domain"),
  xstart = c(0.5, 11, 18.5, 22),
  xend   = c(11, 18.5, 22, 27.5)
)

domain_boundaries <- tibble(x = c(0.5, 11, 18.5, 22, 27.5))

domain_labels <- domain_line %>%
  mutate(
    x = (xstart + xend) / 2,
    label = case_when(
      domain == "N-terminal domain" ~ "N-terminal\ndomain",
      domain == "Pocket domain A"   ~ "Pocket\ndomain A",
      domain == "Pocket domain B"   ~ "Pocket\ndomain B",
      domain == "C-terminal domain" ~ "C-terminal\ndomain",
      TRUE ~ domain
    ),
    y = -0.02
  )

p_domains <- ggplot() +
  geom_segment(data = domain_line,
               aes(x = xstart, xend = xend, y = 0.1, yend = 0.1),
               linewidth = 0.8) +
  geom_segment(data = domain_boundaries,
               aes(x = x, xend = x, y = 0.15, yend = 0.05),
               linewidth = 0.6) +
  geom_text(data = domain_labels,
            aes(x = x, y = y, label = label),
            family = "Times New Roman", size = 3.5, lineheight = 0.9) +
  coord_cartesian(xlim = common_limits, ylim = c(-0.2, 0.1), clip = "off") +
  scale_x_continuous(limits = common_limits, breaks = 1:27, labels = 1:27, expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 5, b = 5, l = 5))

panel_a <- wrap_plots(p_a, p_domains, ncol = 1, heights = c(1, 0.28))


# ----------------------------- PANEL B: consequence enrichment by domain ---------

# domain_variant_table is assumed to already exist with columns:
# domain, aa_start, aa_stop, missense, nonsense, synonymous

domain_possible_long <- domain_variant_table %>%
  filter(domain != "Total") %>%
  select(domain, missense, nonsense, synonymous) %>%
  pivot_longer(
    cols = c(missense, nonsense, synonymous),
    names_to = "molecular_consequence",
    values_to = "n_possible"
  )

bar_dat <- rb_clean_unique_germline %>%
  filter(!is.na(uniprot_regions), str_squish(uniprot_regions) != "") %>%
  transmute(
    domain = str_squish(uniprot_regions),
    molecular_consequence = str_squish(molecular_consequence),
    cdna = cdna
  ) %>%
  filter(molecular_consequence %in% c("stop_gained", "missense", "synonymous")) %>%
  distinct(domain, molecular_consequence, cdna) %>%
  count(domain, molecular_consequence, name = "n_unique") %>%
  mutate(
    molecular_consequence = recode(
      molecular_consequence,
      "stop_gained" = "nonsense",
      "missense" = "missense",
      "synonymous" = "synonymous"
    )
  ) %>%
  left_join(domain_possible_long, by = c("domain", "molecular_consequence")) %>%
  mutate(
    obs_over_possible = n_unique / n_possible,
    domain = factor(
      domain,
      levels = c("N-terminal domain", "Pocket Domain A", "Pocket Domain B", "C-terminal domain")
    ),
    consequence_lab = recode(
      molecular_consequence,
      "synonymous" = "Synonymous",
      "missense" = "Missense",
      "nonsense" = "Stop gained"
    ),
    consequence_lab = factor(
      consequence_lab,
      levels = c("Synonymous", "Missense", "Stop gained")
    ),
    label = sprintf("%.3f", obs_over_possible)
  )

ymax_b <- max(bar_dat$obs_over_possible, na.rm = TRUE) * 1.12

p_b <- ggplot(bar_dat, aes(x = domain, y = obs_over_possible, fill = consequence_lab)) +
  geom_col(
    width = 0.75,
    position = position_dodge(width = 0.8),
    color = NA
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.8),
    vjust = -0.35,
    family = "Times New Roman",
    size = 3.6
  ) +
  scale_fill_manual(
    name = "",
    values = c(
      "Synonymous" = "#D8E1E8",
      "Missense" = "#8AA7B8",
      "Stop gained" = "#3D6C88"
    )
  ) +
  scale_y_continuous(
    limits = c(0, ymax_b),
    expand = c(0, 0)
  ) +
  labs(
    x = "",
    y = "Observed / possible variants",
    tag = "b"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 12)
  )

caption_text <- str_wrap(
  "Figure 3. (a) Germline RB1 variant density across exons normalized by exon length; red bars represent all variants and blue bars represent unique variants. (b) Unique germline variants in each pRB domain, stratified by SNV consequence and normalized by the number of possible variants of that consequence within the domain.",
  width = 300
)

final_fig <- wrap_plots(panel_a, p_b, ncol = 1, heights = c(1.25, 1)) +
  plot_annotation(
    caption = caption_text,
    theme = theme(
      plot.caption = element_textbox_simple(
        size = 12,
        family = "Times New Roman",
        width = unit(9.8, "in"),
        margin = margin(t = 10),
        hjust = 0
      )
    )
  )

final_fig


#### WILCOXON to see if age at diagnosis is different between somatic and germline
dplyr::group_by(rb, genetic_origin) |>
  dplyr::summarise(
    n = sum(!is.na(age_at_diagnosis_months)),
    mean_age = mean(age_at_diagnosis_months, na.rm = TRUE),
    sd_age = sd(age_at_diagnosis_months, na.rm = TRUE),
    median_age = median(age_at_diagnosis_months, na.rm = TRUE),
    IQR_age = IQR(age_at_diagnosis_months, na.rm = TRUE)
  )

wilcox.test(
  age_at_diagnosis_months ~ genetic_origin,
  data = rb[rb$genetic_origin %in% c("germline", "somatic"), ],
  exact = FALSE
)

