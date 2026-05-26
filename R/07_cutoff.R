# Sourcing functions
source("R/00_functions.R")

# Reading in a list of data frames
dfs <- list(
  AD = readr::read_rds("data/paired_diff_archs4/AD_paired_diff.rds") |> 
    combine_lists_to_df(),
  CVD = readr::read_rds("data/paired_diff_archs4/CVD_paired_diff.rds") |> 
    combine_lists_to_df(),
  COPD = readr::read_rds("data/paired_diff_archs4/COPD_paired_diff.rds") |> 
    combine_lists_to_df(),
  OB = readr::read_rds("data/paired_diff_archs4/OB_paired_diff.rds") |> 
    combine_lists_to_df(),
  RA = readr::read_rds("data/paired_diff_archs4/RA_paired_diff.rds") |> 
    combine_lists_to_df(),
  T2D = readr::read_rds("data/paired_diff_archs4/T2D_paired_diff.rds") |> 
    combine_lists_to_df(),
  TB = readr::read_rds("data/paired_diff_archs4/TB_paired_diff.rds") |> 
    combine_lists_to_df(),
  IBD = readr::read_rds("data/paired_diff_archs4/IBD_paired_diff.rds") |> 
    combine_lists_to_df(),
  CANCERS = readr::read_rds("data/paired_diff_archs4/Cancers_paired_diff_5MAY.rds") |> 
    combine_lists_to_df() |> 
    dplyr::mutate(gene = stringr::str_remove(gene, "\\.\\d+$"))
)

x <- background_cutoff(dfs = dfs)
x$plot

y <- sampled_cutoff(dfs = dfs, feature_col = "gene")
y$plot

enrich_res <- enrich_cutoff(observed_results = x$threshold_results, sampled_summary = y$threshold_results_median)

ggplot2::ggplot(
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
    color = "white"
  ) +
  ggplot2::facet_wrap(~ signal_type) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_minimal()

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







