# Sourcing functions
source("R/00_functions.R")

# library
library(patchwork)

# Reading in a list of data frames
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

# meta-p-vals fisher
meta_per_disease <- meta_pvalues_by_gene(new_dfs)
saveRDS(meta_per_disease, "data/formal_pvalue_integration_res/per_disease_fisher_pvalues_V2.rds")

meta_across_diseases <- meta_pvalues_across_diseases(meta_per_disease)
saveRDS(meta_across_diseases, "data/formal_pvalue_integration_res/across_diseases_pvalues_V2.rds")

# Change ensg-id to hgnc mapping
new_dfs_fisher <- purrr::map(meta_per_disease, function(df) {
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

# Plotting fisher aggregated splicing scatter plot
plots_expr_per_disease <- plot_volcano_meta(
  new_dfs_fisher,
  gene_col = "gene",
  x_metric = "median_lfc_expression",
  padj_col = "expression_meta_p",
  padj_threshold = 0.05,
  log2fc_threshold = 1,
  
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

# Saving cancer plot
ggplot2::ggsave(filename = "plots/gene_analysis_plots/Cancer_Fisher_Volcano.png",
                plot = plots_expr_per_disease$CANCERS,
                scale = 1, 
                dpi = 300,
                bg = "white")

# Plotting fisher aggregated splicing scatter plot
plots_sp_per_disease <- plot_splicing_scatter(
  new_dfs_fisher,
  dif_col = "median_max_abs_dif_splicing",
  pval_col = "splicing_meta_p",
  dif_threshold = 0.1,
  pval_threshold = 0.05,
  gene_col = "gene",
  n_labels = 20
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

# Saving RA
ggplot2::ggsave(filename = "plots/gene_analysis_plots/RA_Fisher_scatter.png",
                plot = plots_sp_per_disease$RA,
                scale = 1, 
                dpi = 300,
                bg = "white")

# Saving combined plot
combined_plot <-
  (plots_expr_per_disease$CANCERS | plots_sp_per_disease$RA) +
  plot_annotation(tag_levels = "A")

ggplot2::ggsave(filename = "plots/gene_analysis_plots/combined_plot_cancer_RA.png",
                plot = combined_plot,
                scale = 1, 
                dpi = 300,
                bg = "white")



# Fisher Across diseases expression
expr_fisher_ranks <- meta_across_diseases |>
  dplyr::mutate(
    expression_fisher_p = pmax(expression_fisher_p, .Machine$double.xmin),
    rank = -log10(expression_fisher_p)
  ) |>
  dplyr::select(gene, rank) |>
  tibble::deframe() |>
  sort(decreasing = TRUE)

# Fisher Across diseases splicing
sp_fisher_ranks <- meta_across_diseases |>
  dplyr::mutate(
    splicing_fisher_p = pmax(splicing_fisher_p, .Machine$double.xmin),
    rank = -log10(splicing_fisher_p)
  ) |>
  dplyr::select(gene, rank) |>
  tibble::deframe() |>
  sort(decreasing = TRUE)

# Edgington Across diseases expression
expr_edgintons_ranks <- meta_across_diseases |>
  dplyr::mutate(
    expression_edgington_p = pmax(expression_edgington_p, .Machine$double.xmin),
    rank = -log10(expression_edgington_p)
  ) |>
  dplyr::select(gene, rank) |>
  tibble::deframe() |>
  sort(decreasing = TRUE)

# Edgington Across diseases splicing
sp_edgintons_ranks <- meta_across_diseases |>
  dplyr::mutate(
    splicing_edgington_p = pmax(splicing_edgington_p, .Machine$double.xmin),
    rank = -log10(splicing_edgington_p)
  ) |>
  dplyr::select(gene, rank) |>
  tibble::deframe() |>
  sort(decreasing = TRUE)

# Load genesets for fgsea
genesets <- readRDS("data/genesets/genesets.rds")

# fgsea performed for fisher across disease and filtering for padj afterwards
fgsea_expr_fisher <- run_fgsea_multilevel(pathways = genesets, counted_ranks = expr_fisher_ranks)
fgsea_expr_fisher <- fgsea_expr_fisher |> dplyr::filter(padj < 0.05)
fgsea_sp_fisher <- run_fgsea_multilevel(pathways = genesets, counted_ranks = sp_fisher_ranks)
fgsea_sp_fisher <- fgsea_sp_fisher |> dplyr::filter(padj < 0.05)

# fgsea performed for edgingtons across disease and filtering for padj afterwards
fgsea_expr_edgintons <- run_fgsea_multilevel(pathways = genesets, counted_ranks = expr_edgintons_ranks)
fgsea_expr_edgintons <- fgsea_expr_edgintons |> dplyr::filter(padj < 0.05)
fgsea_sp_edgintons <- run_fgsea_multilevel(pathways = genesets, counted_ranks = sp_edgintons_ranks)
fgsea_sp_edgintons <- fgsea_sp_edgintons |> dplyr::filter(padj < 0.05)
saveRDS(fgsea_expr_fisher, "data/formal_pvalue_fgseaMulti/fgsea_expr_fisher_V2.rds")
saveRDS(fgsea_sp_fisher, "data/formal_pvalue_fgseaMulti/fgsea_sp_fisher_V2.rds")
saveRDS(fgsea_expr_edgintons, "data/formal_pvalue_fgseaMulti/fgsea_expr_edgintons_V2.rds")
saveRDS(fgsea_sp_edgintons, "data/formal_pvalue_fgseaMulti/fgsea_sp_edgintons_V2.rds")


# Binding rows of both FCS analyses - fisher
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

# fisher_gs_list created for geneSetSimplifyR
fisher_gs_list <- geneSetSimplifyR::geneSetSimplifyR(
  geneSetsList = genesets,
  geneSetsDF = fisher_enrichment,
  verbose = FALSE,
  removeFirstWord = TRUE
)

# saved to path
saveRDS(fisher_gs_list, "data/formal_pvalue_integration_geneSetSimplifyR_objects/fisher_gs_list.rds")

# Binding rows of both FCS analyses - edgington
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

# edgington_gs_list created for geneSetSimplifyR
edgington_gs_list <- geneSetSimplifyR::geneSetSimplifyR(
  geneSetsList = genesets,
  geneSetsDF = edgington_enrichment,
  verbose = FALSE,
  removeFirstWord = TRUE
)

# saved to path
saveRDS(edgington_gs_list, "data/formal_pvalue_integration_geneSetSimplifyR_objects/edgingtons_gs_list.rds")