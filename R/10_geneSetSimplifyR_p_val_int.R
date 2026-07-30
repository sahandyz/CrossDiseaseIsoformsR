# Sourcing functions
source("R/00_functions.R")

# packages needed for geneSetSimplifyR
library(geneSetSimplifyR)
library(tidyseurat)
library(ggraph)   

### --- Fisher --- ### Res 0.6
plotClustree(fisher_gs_list) + 
  ggtitle("Double Fisher p-value Integration")

plotUMAP(fisher_gs_list,
         resolution = 0.6) + 
  ggtitle("Fisher p-value integration")

plotProportions(fisher_gs_list,
                resolution = 0.6) + 
  ggtitle("Fisher p-value integration")
  
plotClusterLabels(
  geneSetsList = fisher_gs_list,
  resolution = 0.6,
  clusterNumber = 17,
  maxChar = 35 
) +
  theme_bw(base_size = 12)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "0",
  newClusterName = 'Breast Cancer Targets'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "1",
  newClusterName = 'Cancer Targets'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "2",
  newClusterName = 'Membrane Protein Localization'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "4",
  newClusterName = 'Biosynthetic Metabolic Process'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "5",
  newClusterName = 'T Cell Immune Response'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "6",
  newClusterName = 'Cell Differentiation'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "7",
  newClusterName = 'PD L1 Mediated Degradation'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "8",
  newClusterName = 'Cell Maturation Targets'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "9",
  newClusterName = 'Mitotic Cell Cycle'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "10",
  newClusterName = 'Gtpase Cycle Regulation'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "11",
  newClusterName = 'RNA Polymerase II'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "12",
  newClusterName = 'MRNA Catabolic Process'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "13",
  newClusterName = 'Ion Transmembrane Transport'
)

updated_fisher_gs_list <- updateClusterLabels(
  geneSetsList = updated_fisher_gs_list,
  resolution = 0.6,
  clusterNumber = "15",
  newClusterName = 'Protein Catabolic Process'
)

plotUMAP(
  updated_fisher_gs_list,
  resolution = 0.6,
  removeClusterId = TRUE
) +
  ggtitle("UMAP of Double Fisher p-value Integration")

plotProportions(
  updated_fisher_gs_list,
  resolution = 0.6,
  removeClusterId = TRUE
) +
  ggtitle("Proportions of Double Fisher p-value Integration")


### --- Edginton --- ### Res 1.6
# Using basic vignette workflow of geneSetSimplifyR
plotClustree(edgingtons_gs_list) + 
  ggtitle("Gene-set Clustree of Edgingtons p-value aggregation method")

plotClusterLabels(
  geneSetsList = edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = 22,
  maxChar = 35 
) +
  theme_bw(base_size = 12)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "1",
  newClusterName = 'Cancer Targets'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "3",
  newClusterName = 'Bronchial Epithelial Cells'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "5",
  newClusterName = 'Muscle Cell Differentiation'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "6",
  newClusterName = 'Perturbation Targets'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "7",
  newClusterName = 'Vesicle Mediated Transport'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "10",
  newClusterName = 'Gtpase Cycle Regulation'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "11",
  newClusterName = 'Defense Immune Response'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "12",
  newClusterName = 'Small Ribosomal Subunit'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "14",
  newClusterName = 'Translation Initiation Factor'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "15",
  newClusterName = 'Transmembran Ion Transport'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "16",
  newClusterName = 'T Cell Mediated Immunity'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "17",
  newClusterName = 'Protein Catabolic Process'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "18",
  newClusterName = 'Mitotic Cell Cycle'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "19",
  newClusterName = 'Wound Healing Regulation'
)

updated_edgingtons_gs_list <- updateClusterLabels(
  geneSetsList = updated_edgingtons_gs_list,
  resolution = 1.6,
  clusterNumber = "20",
  newClusterName = 'Adaptive Immune System'
)

# Save relevant updated plots
edgingtons_clustree <- plotClustree(updated_edgingtons_gs_list) + 
  ggtitle("Gene-set Clustree of Edgingtons p-value aggregation method")

ggplot2::ggsave(filename = "plots/formal_pval_int_geneSetSimplifyR/edgingtons_Clustree.png",
                plot = edgingtons_clustree,
                scale = 1,
                dpi = 300,
                bg = "white"
)

edgingtons_UMAP <- plotUMAP(
  updated_edgingtons_gs_list,
  resolution = 1.6,
  removeClusterId = TRUE,
  highlightClusterNo = c(10, 14, 21),
  labelSize = 4
) +
  ggtitle("Gene-set cluster UMAP with Edgingtons p-value aggregation method") +
  theme(
    text = element_text(size = 16)
  )

ggplot2::ggsave(filename = "plots/formal_pval_int_geneSetSimplifyR/edgingtons_UMAP.png",
                plot = edgingtons_UMAP,
                scale = 1,
                dpi = 300,
                bg = "white"
)

edgingtons_proportions <- plotProportions(
  updated_edgingtons_gs_list,
  resolution = 1.6,
  removeClusterId = TRUE,
  highlightClusterNo = c(10, 14, 21)
) +
  ggtitle("Gene-set cluster Proportions with Edgingtons p-value aggregation method") +
  theme(
    text = element_text(size = 14)
  )

ggplot2::ggsave(filename = "plots/formal_pval_int_geneSetSimplifyR/edgingtons_Proportions.png",
                plot = edgingtons_proportions,
                scale = 1,
                dpi = 300,
                bg = "white"
)
