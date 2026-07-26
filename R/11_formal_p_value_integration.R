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

# meta-p-vals 
meta_per_disease <- meta_pvalues_by_gene(new_dfs)

saveRDS(meta_per_disease, "data/formal_pvalue_integration_res/per_disease_fisher_pvalues_V2.rds")
#meta_per_disease_AD <- meta_per_disease$AD
meta_across_diseases <- meta_pvalues_across_diseases(meta_per_disease)
saveRDS(meta_across_diseases, "data/formal_pvalue_integration_res/across_diseases_pvalues_V2.rds")

plots_expr_per_disease <- plot_volcano_meta(
  meta_per_disease,
  lfc_col = "median_lfc_expression",
  padj_col = "expression_meta_fdr"
)
plots_expr_per_disease$AD
plots_expr_per_disease$COPD
plots_expr_per_disease$CVD
plots_expr_per_disease$IBD
plots_expr_per_disease$OB
plots_expr_per_disease$RA
plots_expr_per_disease$TB
plots_expr_per_disease$T2D
plots_expr_per_disease$CANCERS

plots_sp_per_disease <- plot_splicing_scatter(
  meta_per_disease,
  dif_col = "median_max_abs_dif_splicing",
  pval_col = "splicing_meta_fdr",
  dif_threshold = 0.1,
  pval_threshold = 0.05,
  gene_col = "gene"
)
plots_sp_per_disease$AD
plots_sp_per_disease$COPD
plots_sp_per_disease$CVD
plots_sp_per_disease$IBD
plots_sp_per_disease$OB
plots_sp_per_disease$RA
plots_sp_per_disease$TB
plots_sp_per_disease$T2D
plots_sp_per_disease$CANCERS

expr_fisher_ranks <- meta_across_diseases |>
  dplyr::mutate(
    expression_fisher_p = pmax(expression_fisher_p, .Machine$double.xmin),
    rank = -log10(expression_fisher_p)
  ) |>
  dplyr::select(gene, rank) |>
  tibble::deframe() |>
  sort(decreasing = TRUE)

sp_fisher_ranks <- meta_across_diseases |>
  dplyr::mutate(
    splicing_fisher_p = pmax(splicing_fisher_p, .Machine$double.xmin),
    rank = -log10(splicing_fisher_p)
  ) |>
  dplyr::select(gene, rank) |>
  tibble::deframe() |>
  sort(decreasing = TRUE)

expr_edgintons_ranks <- meta_across_diseases |>
  dplyr::mutate(
    expression_edgington_p = pmax(expression_edgington_p, .Machine$double.xmin),
    rank = -log10(expression_edgington_p)
  ) |>
  dplyr::select(gene, rank) |>
  tibble::deframe() |>
  sort(decreasing = TRUE)

sp_edgintons_ranks <- meta_across_diseases |>
  dplyr::mutate(
    splicing_edgington_p = pmax(splicing_edgington_p, .Machine$double.xmin),
    rank = -log10(splicing_edgington_p)
  ) |>
  dplyr::select(gene, rank) |>
  tibble::deframe() |>
  sort(decreasing = TRUE)

genesets <- readRDS("data/genesets/genesets.rds")

fgsea_expr_fisher <- run_fgsea_multilevel(pathways = genesets, counted_ranks = expr_fisher_ranks)
fgsea_expr_fisher <- fgsea_expr_fisher |> dplyr::filter(padj < 0.05)
fgsea_sp_fisher <- run_fgsea_multilevel(pathways = genesets, counted_ranks = sp_fisher_ranks)
fgsea_sp_fisher <- fgsea_sp_fisher |> dplyr::filter(padj < 0.05)

fgsea_expr_edgintons <- run_fgsea_multilevel(pathways = genesets, counted_ranks = expr_edgintons_ranks)
fgsea_expr_edgintons <- fgsea_expr_edgintons |> dplyr::filter(padj < 0.05)
fgsea_sp_edgintons <- run_fgsea_multilevel(pathways = genesets, counted_ranks = sp_edgintons_ranks)
fgsea_sp_edgintons <- fgsea_sp_edgintons |> dplyr::filter(padj < 0.05)
saveRDS(fgsea_expr_fisher, "data/formal_pvalue_fgseaMulti/fgsea_expr_fisher_V2.rds")
saveRDS(fgsea_sp_fisher, "data/formal_pvalue_fgseaMulti/fgsea_sp_fisher_V2.rds")
saveRDS(fgsea_expr_edgintons, "data/formal_pvalue_fgseaMulti/fgsea_expr_edgintons_V2.rds")
saveRDS(fgsea_sp_edgintons, "data/formal_pvalue_fgseaMulti/fgsea_sp_edgintons_V2.rds")


# GenesetSimplifyR object creation
fisher_enrichment <- dplyr::bind_rows(
  fgsea_expr_fisher |>
    dplyr::transmute(
      source = "Expression",
      pathway,
      enrichment_score = NES
    ),
  fgsea_sp_fisher |>
    dplyr::transmute(
      source = "Splicing",
      pathway,
      enrichment_score = NES
    )
)

fisher_gs_list <- geneSetSimplifyR::geneSetSimplifyR(
  geneSetsList = genesets,
  geneSetsDF = fisher_enrichment,
  verbose = FALSE,
  removeFirstWord = TRUE
)
saveRDS(fisher_gs_list, "data/formal_pvalue_integration_geneSetSimplifyR_objects/fisher_gs_list.rds")

edgington_enrichment <- dplyr::bind_rows(
  fgsea_expr_edgintons |>
    dplyr::transmute(
      source = "Expression",
      pathway,
      enrichment_score = NES
    ),
  fgsea_sp_edgintons |>
    dplyr::transmute(
      source = "Splicing",
      pathway,
      enrichment_score = NES
    )
)

edgington_gs_list <- geneSetSimplifyR::geneSetSimplifyR(
  geneSetsList = genesets,
  geneSetsDF = edgington_enrichment,
  verbose = FALSE,
  removeFirstWord = TRUE
)
saveRDS(edgington_gs_list, "data/formal_pvalue_integration_geneSetSimplifyR_objects/edgingtons_gs_list.rds")
