# Sourcing functions
source("R/00_functions.R")

# Combining all paired differential analysis cancer results and putting it next all ARCHS4 data
combine_paired_diff_lists(
  input_dir = "data/pairedV2_diff_tcga",
  output_dir = "data/pairedV2_diff_archs4",
  output_file = "Cancers_paired_diff.rds"
)

# preliminary feature boxplots
plot_summary_feature_boxplots(input_dir = "data/pairedV2_diff_archs4")




