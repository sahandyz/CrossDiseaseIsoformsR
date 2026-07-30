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

dfsV2 <- list(
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

# Replace ENSG IDs with gene symbols in each data frame
dfsV2 <- purrr::map(dfsV2, function(df) {
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

# Run fgsea - Cumulative
fgsea_cumu_d4_s4_genes <- run_fgsea_on_cutoff(dfs = dfsV2, 
                                              path_to_feature_cutoff = "data/cut_off_cumulative/genes_diseases4_studies4.rds")

fgsea_cumu_d4_s4_isoforms <- run_fgsea_on_cutoff(dfs = dfsV2, 
                                                 path_to_feature_cutoff = "data/cut_off_cumulative/isoforms_diseases4_studies4.rds")

fgsea_cumu_d9_s5_genes <- run_fgsea_on_cutoff(dfs = dfsV2,
                                              path_to_feature_cutoff = "data/cut_off_cumulative/genes_diseases9_studies5.rds")
saveRDS(object = fgsea_cumu_d9_s5_genes,
        file = "data/cut_off_fgsea_objects/fgsea_cumu_d9_s5_genes.rds")

fgsea_cumu_d9_s5_isoforms <- run_fgsea_on_cutoff(dfs = dfsV2, 
                                                 path_to_feature_cutoff = "data/cut_off_cumulative/isoforms_diseases9_studies5.rds")
saveRDS(object = fgsea_cumu_d9_s5_isoforms,
        file = "data/cut_off_fgsea_objects/fgsea_cumu_d9_s5_isoforms.rds")

fgsea_cumu_d9_s6_genes <- run_fgsea_on_cutoff(dfs = dfsV2,
                                              path_to_feature_cutoff = "data/cut_off_cumulative/genes_diseases9_studies6.rds")
saveRDS(object = fgsea_cumu_d9_s6_genes,
        file = "data/cut_off_fgsea_objects/fgsea_cumu_d9_s6_genes.rds")

fgsea_cumu_d9_s6_isoforms <- run_fgsea_on_cutoff(dfs = dfsV2, 
                                                 path_to_feature_cutoff = "data/cut_off_cumulative/isoforms_diseases9_studies6.rds")
saveRDS(object = fgsea_cumu_d9_s6_isoforms,
        file = "data/cut_off_fgsea_objects/fgsea_cumu_d9_s6_isoforms.rds")

## Continue with cumulative.##########################################
cumu_d9_s5_combined <- combine_fgsea_gene_isoform(gene_res = fgsea_cumu_d9_s5_genes,
                                                  isoform_res = fgsea_cumu_d9_s5_isoforms,
                                                  sigGenes_genes = readr::read_rds("data/cut_off_cumulative/genes_diseases9_studies5.rds"),
                                                  sigGenes_isoform = readr::read_rds("data/cut_off_cumulative/genes_diseases9_studies5.rds"),
                                                  save_list_path = "data/fgsea_ora_sig_objects/cumu_d9_s5_lists.rds",
                                                  save_enrichment_path = "data/fgsea_enrich_objects/enrich_cumu_d9_s5_combined.rds")
cumu_d9_s6_combined <- combine_fgsea_gene_isoform(gene_res = fgsea_cumu_d9_s6_genes,
                                                  isoform_res = fgsea_cumu_d9_s6_isoforms,
                                                  sigGenes_genes = readr::read_rds("data/cut_off_cumulative/genes_diseases9_studies6.rds"),
                                                  sigGenes_isoform = readr::read_rds("data/cut_off_cumulative/genes_diseases9_studies6.rds"),
                                                  save_list_path = "data/fgsea_ora_sig_objects/cumu_d9_s6_lists.rds",
                                                  save_enrichment_path = "data/fgsea_enrich_objects/enrich_cumu_d9_s6_combined.rds")
