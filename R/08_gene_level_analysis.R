# Sourcing functions
source("R/00_functions.R")

# Reading in a list of data frames
# dfs <- list(
#   AD = readr::read_rds("data/paired_diff_archs4/AD_paired_diff.rds") |> 
#     combine_lists_to_df(),
#   CVD = readr::read_rds("data/paired_diff_archs4/CVD_paired_diff.rds") |> 
#     combine_lists_to_df(),
#   COPD = readr::read_rds("data/paired_diff_archs4/COPD_paired_diff.rds") |> 
#     combine_lists_to_df(),
#   OB = readr::read_rds("data/paired_diff_archs4/OB_paired_diff.rds") |> 
#     combine_lists_to_df(),
#   RA = readr::read_rds("data/paired_diff_archs4/RA_paired_diff.rds") |> 
#     combine_lists_to_df(),
#   T2D = readr::read_rds("data/paired_diff_archs4/T2D_paired_diff.rds") |> 
#     combine_lists_to_df(),
#   TB = readr::read_rds("data/paired_diff_archs4/TB_paired_diff.rds") |> 
#     combine_lists_to_df(),
#   IBD = readr::read_rds("data/paired_diff_archs4/IBD_paired_diff.rds") |> 
#     combine_lists_to_df(),
#   CANCERS = readr::read_rds("data/paired_diff_archs4/Cancers_paired_diff_5MAY.rds") |> 
#     combine_lists_to_df() |> 
#     dplyr::mutate(gene = stringr::str_remove(gene, "\\.\\d+$"))
# )

new_dfs <- list(
  AD = readr::read_rds("data/pairedV2_diff_archs4/AD_paired_diff.rds") |> 
    combine_lists_to_df(),
  CVD = readr::read_rds("data/pairedV2_diff_archs4/CVD_paired_diff.rds") |> 
    combine_lists_to_df(),
  COPD = readr::read_rds("data/pairedV2_diff_archs4/COPD_paired_diff.rds") |> 
    combine_lists_to_df(),
  OB = readr::read_rds("data/pairedV2_diff_archs4/OB_paired_diff.rds") |> 
    combine_lists_to_df(),
  RA = readr::read_rds("data/pairedV2_diff_archs4/RA_paired_diff.rds") |> 
    combine_lists_to_df(),
  T2D = readr::read_rds("data/pairedV2_diff_archs4/T2D_paired_diff.rds") |> 
    combine_lists_to_df(),
  TB = readr::read_rds("data/pairedV2_diff_archs4/TB_paired_diff.rds") |> 
    combine_lists_to_df(),
  IBD = readr::read_rds("data/pairedV2_diff_archs4/IBD_paired_diff.rds") |> 
    combine_lists_to_df(),
  CANCERS = readr::read_rds("data/pairedV2_diff_archs4/Cancers_paired_diff.rds") |> 
    combine_lists_to_df() |> 
    dplyr::mutate(gene = stringr::str_remove(gene, "\\.\\d+$"))
)

# Preparing data with gene_symbol
mart <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

saveRDS(mart, file = "data/ensg_to_symbol/mart.rds")

# Collect all Ensembl IDs across all disease data frames
all_gene_ids <- new_dfs |>
  purrr::map(~ .x$gene) |>
  unlist() |>
  stringr::str_remove("\\.\\d+$") |>
  unique()

saveRDS(all_gene_ids, file = "data/ensg_to_symbol/all_gene_ids.rds")

# Get Ensembl ID -> gene symbol mapping
gene_map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = all_gene_ids,
  mart = mart
)

saveRDS(gene_map, file = "data/ensg_to_symbol/gene_mapping.rds")

# Replace ENSG IDs with gene symbols in each data frame
new_dfs <- purrr::map(new_dfs, function(df) {
  df |>
    dplyr::mutate(
      gene = stringr::str_remove(gene, "\\.\\d+$")
    ) |>
    dplyr::left_join(
      gene_map,
      by = c("gene" = "ensembl_gene_id"),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      gene = dplyr::if_else(
        is.na(hgnc_symbol) | hgnc_symbol == "",
        gene,
        hgnc_symbol
      )
    ) |>
    dplyr::select(-hgnc_symbol)
})

saveRDS(new_dfs, file = "data/ensg_to_symbol/new_dfs.rds")

# Expression plots
vulcano_exp <- plot_volcano_expression(
  new_dfs, 
  lfc_col = "lfc_expression", 
  pval_col = "pvalue_expression",
  label_top_n = 10
)
vulcano_exp$AD
vulcano_exp$COPD
vulcano_exp$CVD
vulcano_exp$IBD
vulcano_exp$OB
vulcano_exp$RA
vulcano_exp$TB
vulcano_exp$T2D
vulcano_exp$CANCERS

# Splicing plots
splicing_plots <- plot_splicing_scatter(
  new_dfs,
  gene_col = "gene",
  dif_col = "max_abs_dif_splicing",
  pval_col = "pvalue_splicing",
  pval_threshold = 0.05,
  dif_threshold = 0.1,
  n_labels = 37
)
splicing_plots$AD
splicing_plots$COPD
splicing_plots$CVD
splicing_plots$IBD
splicing_plots$OB
splicing_plots$RA
splicing_plots$TB
splicing_plots$T2D
splicing_plots$CANCERS
