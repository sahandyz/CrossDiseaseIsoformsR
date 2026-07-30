# Sourcing functions
source("R/00_functions.R")

# Reading in a list of data frames
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

# Inspecting cut offs - exact needs to be used 
background_exact <- background_cutoff(dfs = dfsV2,
                                      studies_comparison = "==",
                                      datasets_comparison = "==")
p_exact <- background_exact$plot +
  ggplot2::ggtitle('Exact combinations of number disease and studies cut-off.') +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 28,
      face = "bold"
    )
  )

ggplot2::ggsave(filename = "plots/exact_cut_off_heatmap.png",
                plot = p_exact,
                width = 14,
                height = 10,
                units = "in",
                dpi = 300,
                bg = "white"
                )

# Inspecting cut offs - for visualization only
background_cumu <- background_cutoff(dfs = dfsV2,
                                     studies_comparison = ">=",
                                     datasets_comparison = ">=")

saveRDS(background_cumu, "data/cut_off/background.rds")

p_cumu <- background_cumu$plot +
  ggplot2::ggtitle('Cumulative number (">=") disease and studies cut-off.') +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 28,
      face = "bold"
    )
  )

ggplot2::ggsave(filename = "plots/cumu_cut_off_heatmap.png",
                plot = p_cumu,
                width = 14,
                height = 10,
                units = "in",
                dpi = 300,
                bg = "white"
                )


# Sampled cut-off utilizing background and the full pool
sampled <- sampled_cutoff(dfs = dfsV2, feature_col = "gene")
saveRDS(sampled, "data/cut_off/sampled.rds")
sampled$plot

p_sampled <- sampled$plot +
  ggplot2::ggtitle('Sampled cut-off.') +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 28,
      face = "bold"
    )
  )

ggplot2::ggsave(filename = "plots/sampled_cut_off_heatmap.png",
                plot = p_sampled,
                width = 14,
                height = 10,
                units = "in",
                dpi = 300,
                bg = "white"
)

threshold_res <- background_exact$threshold_results
sampled_res <- sampled$threshold_results_median

# Enrichment of background and samples data
enrich_res <- enrich_cutoff(observed_results = threshold_res, 
                            sampled_summary = sampled_res) |> 
  dplyr::mutate(
    log2_enrichment = dplyr::case_when(
      is.infinite(log2_enrichment) | is.nan(log2_enrichment) ~ 0,
      TRUE ~ log2_enrichment
    )
  )

# heatmap of log2 enrichment score
p_enrich <- ggplot2::ggplot(
  enrich_res,
  ggplot2::aes(
    x = min_studies,
    y = min_datasets,
    fill = log2_enrichment
  )
  ) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = round(log2_enrichment, 2)
    ),
    color = "white",
    angle = -45
  ) +
  ggplot2::facet_wrap(~ signal_type) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(enrich_res$min_studies))
  ) +
  ggplot2::scale_y_continuous(
    breaks = sort(unique(enrich_res$min_datasets))
  ) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_minimal() +
  ggplot2::ggtitle('Enrichment cut-off.') +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 28,
      face = "bold"
    )
  )

ggplot2::ggsave(filename = "plots/enrich_cut_off_heatmap.png",
                plot = p_enrich,
                width = 14,
                height = 10,
                units = "in",
                dpi = 300,
                bg = "white"
)

# heatmap of z-score
ggplot2::ggplot(
  enrich_res,
  ggplot2::aes(
    x = min_studies,
    y = min_datasets,
    fill = z_score
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = round(z_score, 1)
    ),
    color = "white"
  ) +
  ggplot2::facet_wrap(~ signal_type) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_minimal()


### Cumulative extraction ###
# min_studies 6 and min_datasets (disease) 9
# Genes
extract_background_genes(dfs = dfsV2,
                         padj_col = "padj_expression",
                         min_studies = 6,
                         min_datasets = 9,
                         feature_col = "gene",
                         cumulative = TRUE,
                         new_filename = "data/cut_off_cumulative/genes_diseases9_studies6.rds")
extract_background_genes(dfs = dfsV2,
                         padj_col = "padj_splicing",
                         min_studies = 6,
                         min_datasets = 9,
                         feature_col = "gene",
                         cumulative = TRUE,
                         new_filename = "data/cut_off_cumulative/isoforms_diseases9_studies6.rds")

# min_studies 5 and min_datasets (disease) 9
# Genes
extract_background_genes(dfs = dfsV2,
                         padj_col = "padj_expression",
                         min_studies = 5,
                         min_datasets = 9,
                         feature_col = "gene",
                         cumulative = TRUE,
                         new_filename = "data/cut_off_cumulative/genes_diseases9_studies5.rds")
extract_background_genes(dfs = dfsV2,
                         padj_col = "padj_splicing",
                         min_studies = 5,
                         min_datasets = 9,
                         feature_col = "gene",
                         cumulative = TRUE,
                         new_filename = "data/cut_off_cumulative/isoforms_diseases9_studies5.rds")




