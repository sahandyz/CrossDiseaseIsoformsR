# Sourcing functions
source("R/00_functions.R")


# preliminary feature plots
plot_summary_feature_boxplots(input_dir = "data/paired_diff_archs4")

# preliminary gene plots
geneset_results <- plot_summary_geneset_boxplots(
  input_dir = "data/ora_archs4/",
  pattern = "_ora_raw\\.rds$",
  output_dir = "plots_genesets"
)

# preliminary ggridges plot
make_ridge_from_ora(path = "data/ora_archs4/")


