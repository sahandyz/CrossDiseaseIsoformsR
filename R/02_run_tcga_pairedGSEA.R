source(file = "R/00_functions.R")

# Run functions
run_all_tcga(
  tcga_dir = "data/curated_tcga/",
  n_cores = 12,
  out_dir = "data/paired_diff_tcga"
  )