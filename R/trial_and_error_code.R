# Sourcing functions
source("R/00_functions.R")

AD_paired_diff_df <- AD_paired_diff |> 
  combine_lists_to_df() |> 
  dplyr::filter(padj_splicing < 0.05) |> 
  dplyr::group_by(gene) |> 
  dplyr::mutate(n = dplyr::n()) |> 
  dplyr::arrange(n)


Cancers_paired_diff_4APR_df <- Cancers_paired_diff_4APR |> 
  combine_lists_to_df()

library(dplyr)
library(ggplot2)

Cancers_paired_diff_4APR_df |>
  filter(.id == "BRCA : 6 vs 1") |>
  filter(padj_splicing < 0.05) |> 
  ggplot(aes(x = lfc_expression, y = -log10(pvalue_expression))) +
  geom_point(alpha = 0.6) +
  theme_minimal()



cancer_subset <- Cancers_paired_diff_4APR_df |> 
  filter(.id == "BRCA : 6 vs 1")

# Define thresholds
padj_threshold <- 0.05  # Adjusted p-value threshold
log2fc_threshold <- 1   # Log2 fold change threshold

tidy_deseq_results <- cancer_subset |>
  mutate(
    regulation = case_when(
      padj_expression < padj_threshold & lfc_expression > log2fc_threshold ~ "Up-regulated",
      padj_expression < padj_threshold & lfc_expression < -log2fc_threshold ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

volcano <- ggplot(tidy_deseq_results, aes(x = lfc_expression, y = -log10(padj_expression), color = regulation)) +
  geom_point(alpha = 0.8, size = 2) +  
  scale_color_manual(values = c(
    "Up-regulated" = "#F8766D",
    "Down-regulated" = "#619CFF",
    "Not significant" = "gray"
  )) +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Regulation"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 12)
  )
volcano


library(dplyr)
library(purrr)
library(tidyr)
library(UpSetR)

# ------------------------------------------------------------
# 0. Put all your dataframes into a named list
# ------------------------------------------------------------



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
#saveRDS(object = dfs, file = "data/list_of_diseases/dfs.rds")
# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

padj_thresh <- 0.05
min_studies <- 3     # within a dataframe
min_datasets <- 2    # across dataframes

# ------------------------------------------------------------
# 2. Expression: filter + require >= min_studies per dataframe
# ------------------------------------------------------------

genes_expr_by_df <- dfs |>
  map(~ .x |>
        filter(!is.na(padj_expression),
               padj_expression < padj_thresh) |>
        distinct(.id, gene) |>
        count(gene, name = "n_studies") |>
        filter(n_studies >= min_studies) |>
        pull(gene)
  )

# ------------------------------------------------------------
# 3. Splicing: filter + require >= min_studies per dataframe
# ------------------------------------------------------------

genes_splice_by_df <- dfs |>
  map(~ .x |>
        filter(!is.na(padj_splicing),
               padj_splicing < padj_thresh) |>
        distinct(.id, gene) |>
        count(gene, name = "n_studies") |>
        filter(n_studies >= min_studies) |>
        pull(gene)
  )

# ------------------------------------------------------------
# 4. Genes present in >= min_datasets datasets
# ------------------------------------------------------------

expr_genes_summary <- genes_expr_by_df |>
  imap_dfr(~ tibble(gene = .x, dataset = .y)) |>
  distinct(gene, dataset) |>
  count(gene, name = "n_datasets") |>
  filter(n_datasets >= min_datasets) |>
  arrange(desc(n_datasets))

splice_genes_summary <- genes_splice_by_df |>
  imap_dfr(~ tibble(gene = .x, dataset = .y)) |>
  distinct(gene, dataset) |>
  count(gene, name = "n_datasets") |>
  filter(n_datasets >= min_datasets) |>
  arrange(desc(n_datasets))

expr_genes_summary
splice_genes_summary

# vectors of genes (if you need them)
expr_genes <- expr_genes_summary$gene
splice_genes <- splice_genes_summary$gene

# ------------------------------------------------------------
# 5. Optional: genes shared by ALL datasets (strict intersection)
# ------------------------------------------------------------

common_expr_all <- reduce(genes_expr_by_df, intersect)
common_splice_all <- reduce(genes_splice_by_df, intersect)

common_expr_all
common_splice_all

# ------------------------------------------------------------
# 6. Get filtered rows back from original dfs
# ------------------------------------------------------------

filtered_expr_dfs <- map(
  dfs,
  ~ .x |>
    filter(gene %in% expr_genes,
           !is.na(padj_expression),
           padj_expression < padj_thresh)
)

filtered_splice_dfs <- map(
  dfs,
  ~ .x |>
    filter(gene %in% splice_genes,
           !is.na(padj_splicing),
           padj_splicing < padj_thresh)
)

# ------------------------------------------------------------
# 7. Presence/absence tables (nice for inspection)
# ------------------------------------------------------------

expr_presence <- genes_expr_by_df |>
  imap_dfr(~ tibble(gene = .x, dataset = .y)) |>
  distinct() |>
  mutate(present = TRUE) |>
  pivot_wider(names_from = dataset,
              values_from = present,
              values_fill = FALSE)

splice_presence <- genes_splice_by_df |>
  imap_dfr(~ tibble(gene = .x, dataset = .y)) |>
  distinct() |>
  mutate(present = TRUE) |>
  pivot_wider(names_from = dataset,
              values_from = present,
              values_fill = FALSE)

expr_presence
splice_presence

# ------------------------------------------------------------
# 8. UpSet plots
# ------------------------------------------------------------

upset(fromList(genes_expr_by_df),
      order.by = "freq",
      nsets = length(genes_expr_by_df),
      mainbar.y.label = "Expression intersections")

upset(fromList(genes_splice_by_df),
      order.by = "freq",
      nsets = length(genes_splice_by_df),
      mainbar.y.label = "Splicing intersections")

# ------------------------------------------------------------
# 9. Quick diagnostics (VERY useful if things look empty)
# ------------------------------------------------------------

map_int(genes_expr_by_df, length)
map_int(genes_splice_by_df, length)

library(dplyr)

padj_thresh <- 0.05

precompute_gene_counts <- function(padj_col) {
  dfs |>
    purrr::imap_dfr(~ .x |> 
                      dplyr::filter(
                        !is.na(.data[[padj_col]]),
                        .data[[padj_col]] < padj_thresh
                      ) |>
                      dplyr::distinct(.id, gene) |> 
                      dplyr::count(gene, name = "n_studies") |> 
                      dplyr::mutate(dataset = .y)
    ) |> 
    dplyr::group_by(gene) |> 
    tidyr::nest(dataset_counts = c(dataset, n_studies)) |> 
    dplyr::ungroup()
}

expression_counts <- precompute_gene_counts("padj_expression")
splicing_counts   <- precompute_gene_counts("padj_splicing")

genes_for_threshold <- function(counts, min_studies, min_datasets) {
  counts |>
    dplyr::mutate(
      n_datasets = purrr::map_int(
        dataset_counts,
        ~ sum(.x$n_studies == min_studies)
      )
    ) |>
    dplyr::filter(n_datasets == min_datasets) |>
    dplyr::pull(gene)
}

threshold_results <- tidyr::expand_grid(
  min_studies = 1:6,
  min_datasets = 1:length(dfs)
) |>
  dplyr::mutate(
    expression_genes = purrr::map2(
      min_studies, min_datasets,
      ~ genes_for_threshold(expression_counts, .x, .y)
    ),
    splicing_genes = purrr::map2(
      min_studies, min_datasets,
      ~ genes_for_threshold(splicing_counts, .x, .y)
    ),
    n_expression_genes = purrr::map_int(expression_genes, length),
    n_splicing_genes = purrr::map_int(splicing_genes, length),
    n_overlap_genes = purrr::map2_int(
      expression_genes,
      splicing_genes,
      ~ length(intersect(.x, .y))
    )
  ) |>
  dplyr::select(-expression_genes, -splicing_genes)

threshold_results

threshold_results_long <- threshold_results |>
  tidyr::pivot_longer(
    cols = c(n_expression_genes, n_splicing_genes, n_overlap_genes),
    names_to = "signal_type",
    values_to = "n_genes"
  ) |>
  dplyr::mutate(
    signal_type = dplyr::recode(
      signal_type,
      n_expression_genes = "Expression",
      n_splicing_genes = "Splicing",
      n_overlap_genes = "Overlap"
    ),
    signal_type = factor(
      signal_type,
      levels = c("Expression", "Splicing", "Overlap")
    )
  )

ggplot2::ggplot(
  threshold_results_long,
  ggplot2::aes(
    x = min_studies,
    y = min_datasets,
    fill = n_genes
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = n_genes),
    color = "white",
    size = 5
  ) +
  ggplot2::facet_wrap(~ signal_type) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(threshold_results_long$min_studies))
  ) +
  ggplot2::scale_y_continuous(
    breaks = sort(unique(threshold_results_long$min_datasets))
  ) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 12),
    strip.text = ggplot2::element_text(size = 14)
  ) +
  ggplot2::labs(
    x = "Minimum studies",
    y = "Minimum diseases",
    fill = "Genes retained"
  )

ggplot2::ggplot(
  threshold_results_long,
  ggplot2::aes(
    x = min_studies,
    y = n_genes,
    color = factor(min_datasets)
  )
) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~ signal_type) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(threshold_results_long$min_studies))
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    x = "Minimum studies",
    y = "Genes retained",
    color = "Min datasets"
  )


ggplot2::ggplot(
  threshold_results,
  ggplot2::aes(
    x = n_splicing_genes,
    y = n_expression_genes
  )
) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0("S", min_studies, "_D", min_datasets)),
    vjust = -0.7,
    size = 3
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    x = "Number of splicing genes",
    y = "Number of expression genes",
    title = "Expression vs Splicing gene counts across thresholds"
  )

gene_study_counts <- dfs |>
  imap_dfr(~ .x |>
             select(.id, gene, padj_expression, padj_splicing) |>
             mutate(dataset = .y)
  )

expr_counts <- gene_study_counts |>
  filter(!is.na(padj_expression),
         padj_expression < padj_thresh) |>
  distinct(dataset, .id, gene) |>
  count(dataset, gene, name = "n_studies")

splice_counts <- gene_study_counts |>
  filter(!is.na(padj_splicing),
         padj_splicing < padj_thresh) |>
  distinct(dataset, .id, gene) |>
  count(dataset, gene, name = "n_studies")


target_min_studies <- 6
target_min_datasets <- 9

splicing_genes_target <- splice_counts |>
  filter(n_studies >= target_min_studies) |>
  distinct(dataset, gene) |>
  count(gene, name = "n_datasets") |>
  filter(n_datasets >= target_min_datasets) |>
  pull(gene)

splicing_genes_target



library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

set.seed(123)

padj_thresh <- 0.05
n_repeats <- 99

test_thresholds <- tidyr::expand_grid(
  min_studies = 1:6,
  min_datasets = 1:length(dfs)
)

get_sampling_x <- function(dfs, padj_col, padj_thresh = 0.05) {
  dfs |>
    purrr::imap_dfr(~ .x |>
                      dplyr::filter(
                        !is.na(.data[[padj_col]]),
                        .data[[padj_col]] < padj_thresh
                      ) |>
                      dplyr::distinct(gene) |>
                      dplyr::summarise(
                        dataset = .y,
                        n_sig_genes = dplyr::n(),
                        .groups = "drop"
                      )
    ) |>
    dplyr::summarise(x = floor(min(n_sig_genes, na.rm = TRUE))) |>
    dplyr::pull(x)
}

x_expression <- get_sampling_x(dfs, "padj_expression", padj_thresh)
x_splicing   <- get_sampling_x(dfs, "padj_splicing", padj_thresh)

sample_dfs_by_gene <- function(dfs, padj_col, x, padj_thresh = 0.05) {
  dfs |>
    purrr::map(~ {
      df_filtered <- .x |>
        dplyr::filter(
          !is.na(.data[[padj_col]]),
          .data[[padj_col]] < padj_thresh
        )
      
      unique_genes <- df_filtered |>
        dplyr::distinct(gene)
      
      sampled_genes <- unique_genes |>
        dplyr::slice_sample(n = min(x, nrow(unique_genes))) |>
        dplyr::pull(gene)
      
      df_filtered |>
        dplyr::filter(gene %in% sampled_genes)
    })
}

summarise_sampled_dfs <- function(sampled_dfs) {
  sampled_dfs |>
    purrr::imap_dfr(~ .x |>
                      dplyr::distinct(.id, gene) |>
                      dplyr::count(gene, name = "n_studies") |>
                      dplyr::mutate(dataset = .y)
    )
}

count_from_summary <- function(sample_summary, min_studies, min_datasets) {
  sample_summary |>
    dplyr::filter(n_studies == min_studies) |>
    dplyr::count(gene, name = "n_datasets") |>
    dplyr::filter(n_datasets == min_datasets) |>
    nrow()
}

run_one_repeat_fast <- function(dfs,
                                padj_col,
                                x,
                                repeat_id,
                                test_thresholds,
                                padj_thresh = 0.05) {
  
  sampled_dfs <- sample_dfs_by_gene(
    dfs = dfs,
    padj_col = padj_col,
    x = x,
    padj_thresh = padj_thresh
  )
  
  sample_summary <- summarise_sampled_dfs(sampled_dfs)
  
  test_thresholds |>
    dplyr::mutate(
      n_genes = purrr::map2_int(
        min_studies,
        min_datasets,
        ~ count_from_summary(sample_summary, .x, .y)
      ),
      repeat_id = repeat_id
    )
}

expr_repeat_results <- purrr::map_dfr(
  seq_len(n_repeats),
  ~ run_one_repeat_fast(
    dfs = dfs,
    padj_col = "padj_expression",
    x = x_expression,
    repeat_id = .x,
    test_thresholds = test_thresholds,
    padj_thresh = padj_thresh
  )
) |>
  dplyr::mutate(signal_type = "expression")

splice_repeat_results <- purrr::map_dfr(
  seq_len(n_repeats),
  ~ run_one_repeat_fast(
    dfs = dfs,
    padj_col = "padj_splicing",
    x = x_splicing,
    repeat_id = .x,
    test_thresholds = test_thresholds,
    padj_thresh = padj_thresh
  )
) |>
  dplyr::mutate(signal_type = "splicing")

repeat_results <- dplyr::bind_rows(
  expr_repeat_results,
  splice_repeat_results
)

threshold_results_median <- repeat_results |>
  dplyr::group_by(signal_type, min_studies, min_datasets) |>
  dplyr::summarise(
    median_n_genes = stats::median(n_genes),
    mean_n_genes = mean(n_genes),
    sd_n_genes = stats::sd(n_genes),
    .groups = "drop"
  )

threshold_results_median

ggplot2::ggplot(
  threshold_results_median,
  ggplot2::aes(
    x = min_studies,
    y = min_datasets,
    fill = median_n_genes
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = median_n_genes),
    color = "white",
    size = 3
  ) +
  ggplot2::facet_wrap(~ signal_type) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(threshold_results_median$min_studies))
  ) +
  ggplot2::scale_y_continuous(
    breaks = sort(unique(threshold_results_median$min_datasets))
  ) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    x = "Minimum studies",
    y = "Minimum diseases",
    fill = "Median genes retained"
  )


library(dplyr)
# Counting analysis
# Expression objects
metaGenesExp_list <- purrr::imap(dfs, ~ {
  obj <- .x |>
    dplyr::filter(padj_expression < 0.05) |>
    dplyr::count(gene, name = "n") |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::slice_head(n = 5)
  
  assign(
    paste0("metaGenesExp_", .y),
    obj,
    envir = .GlobalEnv
  )
  
  obj
})

# Splicing objects
metaGenesSp_list <- purrr::imap(dfs, ~ {
  obj <- .x |>
    dplyr::filter(padj_splicing < 0.05) |>
    dplyr::count(gene, name = "n") |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::slice_head(n = 5)
  
  assign(
    paste0("metaGenesSp_", .y),
    obj,
    envir = .GlobalEnv
  )
  
  obj
})

get_disease_patterns <- function(counts,
                                 min_studies = 6,
                                 min_datasets = 9) {
  
  counts |>
    tidyr::unnest(dataset_counts) |>
    dplyr::group_by(gene) |>
    dplyr::filter(sum(n_studies >= min_studies) >= min_datasets) |>
    dplyr::arrange(dataset, .by_group = TRUE) |>
    dplyr::summarise(
      diseases = paste(dataset, collapse = " | "),
      study_counts = paste(n_studies, collapse = " | "),
      combo = paste0(dataset, ":", n_studies, collapse = " | "),
      .groups = "drop"
    )
}
expression_combos <- get_disease_patterns(
  expression_counts,
  min_studies = 1,
  min_datasets = 9
)

splicing_combos <- get_disease_patterns(
  splicing_counts,
  min_studies = 1,
  min_datasets = 9
)

expression_combos
splicing_combos


library(dplyr)
library(tidyr)
library(ggplot2)

count_retained_by_threshold <- function(counts,
                                        min_datasets = 9,
                                        thresholds = 1:20,
                                        label = "expression") {
  
  counts_long <- counts |>
    unnest(dataset_counts)
  
  purrr::map_dfr(thresholds, function(thr) {
    counts_long |>
      group_by(gene) |>
      summarise(
        n_datasets_passing = sum(n_studies >= thr),
        .groups = "drop"
      ) |>
      filter(n_datasets_passing >= min_datasets) |>
      summarise(
        n_retained = n(),
        threshold = thr,
        type = label
      )
  })
}
plot_df <- bind_rows(
  count_retained_by_threshold(
    expression_counts,
    min_datasets = 9,
    thresholds = 1:20,
    label = "Expression"
  ),
  count_retained_by_threshold(
    splicing_counts,
    min_datasets = 9,
    thresholds = 1:20,
    label = "Splicing"
  )
)

ggplot(plot_df, aes(x = threshold, y = n_retained, color = type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 10)) +
  labs(
    x = "Minimum studies threshold",
    y = "Retained genes / isoforms"
  ) +
  theme_minimal(base_size = 14)
