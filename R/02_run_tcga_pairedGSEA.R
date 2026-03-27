source(file = "/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/R/00_functions.R")

# Run functions
# run_all_tcga(
#   tcga_dir = "data/curated_tcga/",
#   n_cores = 12,
#   out_dir = "data/paired_diff_tcga"
#   )
ctrl_vs_BRCA_LumA <- run_paired_diff(
  curated_tcga_file_path = "/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/data/curated_tcga/BRCA_LumA.xlsx",
  archs4 = FALSE
)

saveRDS(ctrl_vs_BRCA_LumA, "/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/data/paired_diff_tcga/BRCA_LumA_paired_diff.rds")