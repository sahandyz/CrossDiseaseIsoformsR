# Collects a dataframe of wanted geneset categories
my_go <- msigdbr::msigdbr(species = 'Homo sapiens') |> 
  
  dplyr::filter(
    
    gs_collection %in% c(
      
      'C5', # Gene ontology
      'H'   # hallmarks
      
    ) | gs_subcollection %in% c(
      
      'CP:KEGG', # KEGG DB
      'CP:REACTOME', # REACTOME DB
      'CP:WIKIPATHWAYS', # WIKIPATHWAYS
      'CGP'  # Chemical and genetic perterbutation
      
    )
  )

# Removal of all genesets that are phenotype (HP_geneset_name)
exclude_phenotypes <- my_go |> 
  dplyr::filter(!stringr::str_starts(gs_name, "HP_"))

# Create a geneset + gene list
m_sig_db_list <- split(
  exclude_phenotypes$ensembl_gene,
  exclude_phenotypes$gs_name
)

# Save mSigDB as rds object
saveRDS(
  object = m_sig_db_list,
  file = "/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/data/genesets/genesets.rds"
)
