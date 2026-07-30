# Sourcing functions
source("R/00_functions.R")

combine_paired_diff_lists(
  input_dir = "data/pairedV2_diff_tcga",
  output_dir = "data/pairedV2_diff_archs4",
  output_file = "Cancers_paired_diff.rds"
)

# preliminary feature plots
plot_summary_feature_boxplots(input_dir = "data/pairedV2_diff_archs4")




