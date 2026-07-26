# Sourcing functions
source("R/00_functions.R")

combine_paired_diff_lists(
  input_dir = "data/pairedV2_diff_tcga",
  output_dir = "data/pairedV2_diff_archs4",
  output_file = "Cancers_paired_diff.rds"
)

# preliminary feature plots
plot_summary_feature_boxplots(input_dir = "data/pairedV2_diff_archs4")

# preliminary gene plots
geneset_results <- plot_summary_geneset_boxplots(
  input_dir = "data/ora_archs4/",
  pattern = "_ora_raw\\.rds$",
  output_dir = "plots_genesets"
)

# preliminary ggridges plot
make_ridge_from_ora(path = "data/ora_archs4/")


