source("R/00_functions.R")


# Run all tcga subtypes
run_all_tcga(
  tcga_dir = "data/curated_tcga_immune/",
  out_dir = "data/pairedV2_diff_tcga/"
  )
