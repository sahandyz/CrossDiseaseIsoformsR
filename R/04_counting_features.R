source("R/00_functions.R")

# Counting features for expression and saving them in data/counting_features_archs4
counting_features(
  path = "data/paired_diff_archs4/",
  padj_score = "padj_expression", 
  new_path = "data/counting_features_archs4",
  n_cores = 8
)

# Counting features for splicing and saving them in data/counting_features_archs4
counting_features(
  path = "data/paired_diff_archs4/",
  padj_score = "padj_splicing", 
  new_path = "data/counting_features_archs4",
  n_cores = 8
)

# Counting features for expression and saving them in data/counting_features_archs4
counting_features(
  path = "data/pairedV2_diff_archs4/",
  padj_score = "padj_expression", 
  new_path = "data/counting_features_archs4_V2",
  n_cores = 8
)
AD_counts <- AD_feature_counts_padj_expression |> 
  dplyr::inner_join(gene_map, by = c("gene" = "ensembl_gene_id"))
# Counting features for splicing and saving them in data/counting_features_archs4
counting_features(
  path = "data/pairedV2_diff_archs4/",
  padj_score = "padj_splicing", 
  new_path = "data/counting_features_archs4_V2",
  n_cores = 8
)
