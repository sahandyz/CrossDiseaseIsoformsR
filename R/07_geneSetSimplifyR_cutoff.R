# Sourcing functions
source("R/00_functions.R")

# packages needed for geneSetSimplifyR
library(geneSetSimplifyR)
library(tidyseurat)
library(ggraph)   

# Create gsList objects for each cutoff
cumu_d9_s5_gsList <- geneSetSimplifyR(geneSetsList = readr::read_rds("data/fgsea_ora_sig_objects/cumu_d9_s5_lists.rds"),
                                      geneSetsDF = readr::read_rds("data/fgsea_enrich_objects/enrich_cumu_d9_s5_combined.rds"),
                                      verbose = FALSE,
                                      removeFirstWord = TRUE)
#saveRDS(cumu_d9_s5_gsList, "data/geneSetSimplifyR_objects/cumu_d9_s5_gsList.rds")

cumu_d9_s6_gsList <- geneSetSimplifyR(geneSetsList = readr::read_rds("data/fgsea_ora_sig_objects/cumu_d9_s6_lists.rds"),
                                      geneSetsDF = readr::read_rds("data/fgsea_enrich_objects/enrich_cumu_d9_s6_combined.rds"),
                                      verbose = FALSE,
                                      removeFirstWord = TRUE)
#saveRDS(cumu_d9_s6_gsList, "data/geneSetSimplifyR_objects/cumu_d9_s6_gsList.rds")


# Find appropriate split and thereby clustering resolution
cumu_d9_s5_Clustree <- plotClustree(updated_cumu_d9_s5_gsList) +
  ggtitle("Gene-set Clustree of cut-off with 9 diseases and 5 studies")

ggplot2::ggsave(filename = "plots/cut_off_geneSetSimplifyR/cumu_d9_s5_Clustree.png",
                plot = cumu_d9_s5_Clustree,
                scale = 1,
                dpi = 300,
                bg = "white"
)

plotUMAP(
  geneSetsList = cumu_d9_s5_gsList,
  resolution = 1.1
  ) +
  ggtitle("9 disease and 5 studies")

plotProportions(
  geneSetsList = cumu_d9_s5_gsList,
  resolution = 1.1,
  removeClusterId = TRUE
) +
  ggtitle("9 disease and 5 studies") 

plotClusterLabels(
  geneSetsList = cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 17,
  maxChar = 35 
) +
  theme_bw(base_size = 12)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 3,
  newClusterName = 'Fatty Acid Metabolic Process'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 4,
  newClusterName = 'T Cell Mediated Immunity'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 5,
  newClusterName = 'Actin Filament Regulation'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 7,
  newClusterName = 'Membrane Protein Transport'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 8,
  newClusterName = 'Intrinsic Apoptotic Signaling'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 9,
  newClusterName = 'Ion Transmembrane Transport'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 10,
  newClusterName = 'Developmental Cell Growth'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 11,
  newClusterName = 'Mitotic Cell Cycle'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 15,
  newClusterName = 'Large Ribosomal Subunit'
)

updated_cumu_d9_s5_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  clusterNumber = 16,
  newClusterName = 'Pattern Recognition Receptor'
)


saveRDS(updated_cumu_d9_s5_gsList, "data/geneSetSimplifyR_objects/updated_cumu_d9_s5_gsList.rds")

# Plots for cumu_d4_s4_gsList - best resolution 1.1
cumu_d9_s5_UMAP <- plotUMAP(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  removeClusterId = TRUE,
  highlightClusterNo = c(4, 6, 8),
  labelSize = 4
) +
  ggtitle("Gene-set cluster UMAP of cut-off with 9 diseases and 5 studies") +
  theme(
    text = element_text(size = 16)
  )

ggplot2::ggsave(filename = "plots/cut_off_geneSetSimplifyR/cumu_d9_s5_UMAP.png",
                plot = cumu_d9_s5_UMAP,
                scale = 1,
                dpi = 300,
                bg = "white"
)

cumu_d9_s5_Proportions <- plotProportions(
  geneSetsList = updated_cumu_d9_s5_gsList,
  resolution = 1.1,
  removeClusterId = TRUE,
  highlightClusterNo = c(4, 6, 8)
) +
  ggtitle("Gene-set cluster Proportions of cut-off 9 diseases and 5 studies") +
  theme(
    text = element_text(size = 14)
  )

ggplot2::ggsave(filename = "plots/cut_off_geneSetSimplifyR/cumu_d9_s5_Proportions.png",
                plot = cumu_d9_s5_Proportions,
                scale = 1,
                dpi = 300,
                bg = "white"
)

# Strict cut-off
# Find appropriate split and thereby clustering resolution
cumu_d9_s6_Clustree <- plotClustree(updated_cumu_d9_s6_gsList) +
  ggtitle("Gene-set Clustree of cut-off with 9 diseases and 6 studies")

ggplot2::ggsave(filename = "plots/cut_off_geneSetSimplifyR/cumu_d9_s6_Clustree.png",
                plot = cumu_d9_s6_Clustree,
                scale = 1,
                dpi = 300,
                bg = "white"
)

plotUMAP(
  geneSetsList = cumu_d9_s6_gsList,
  resolution = 0.6
) +
  ggtitle("9 disease and 6 studies")

plotProportions(
  geneSetsList = cumu_d9_s6_gsList,
  resolution = 0.6,
  removeClusterId = TRUE
) +
  ggtitle("9 disease and 6 studies") 

plotClusterLabels(
  geneSetsList = cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = 15,
  maxChar = 35
) +
  theme_bw(base_size = 12)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = 1,
  newClusterName = 'Intrinsic Apoptotic Signaling'
)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = "2",
  newClusterName = 'Leukocyte Mediated Immune Response'
)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = "3",
  newClusterName = 'Cancer Targets'
)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = "4",
  newClusterName = 'Lipid Metabolic Process'
)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = "5",
  newClusterName = 'Actin Filament Regulation'
)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = "8",
  newClusterName = 'PD L1 Mediated Degradation'
)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = "9",
  newClusterName = 'Amino Acid Transport'
)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = "13",
  newClusterName = 'Translation Initiation'
)

updated_cumu_d9_s6_gsList <- updateClusterLabels(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  clusterNumber = "14",
  newClusterName = 'Large Ribosomal Subunit'
)

saveRDS(updated_cumu_d9_s6_gsList, "data/geneSetSimplifyR_objects/updated_cumu_d9_s6_gsList.rds")

# Plots for cumu_d9_s6_gsList - best resolution 1.6
cumu_d9_s6_UMAP <- plotUMAP(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  removeClusterId = TRUE,
  highlightClusterNo = c(2, 4, 5),
  labelSize = 4
) +
  ggtitle("Gene-set cluster UMAP of cut-off with 9 diseases and 6 studies") +
  theme(
    text = element_text(size = 16)
  )
  

ggplot2::ggsave(filename = "plots/cut_off_geneSetSimplifyR/cumu_d9_s6_UMAP.png",
                plot = cumu_d9_s6_UMAP,
                scale = 1,
                dpi = 300,
                bg = "white"
)

cumu_d9_s6_Proportions <- plotProportions(
  geneSetsList = updated_cumu_d9_s6_gsList,
  resolution = 0.6,
  removeClusterId = TRUE,
  highlightClusterNo = c(2, 4, 5)
) +
  ggtitle("Gene-set cluster Proportions of cut-off 9 diseases and 6 studies") +
  theme(
    text = element_text(size = 14)
  )

ggplot2::ggsave(filename = "plots/cut_off_geneSetSimplifyR/cumu_d9_s6_Proportions.png",
                plot = cumu_d9_s6_Proportions,
                scale = 1,
                dpi = 300,
                bg = "white"
)
