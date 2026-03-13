# Load libraries and functions
source("R/00_functions.R")

# Run data_mine_disease - TB
data_mine_disease(metadata_path = "data/metadata_archs4/01_archs4_filtered_metadata.rds",
                  keywords = c("tuberculosis",
                               "Granuloma",
                               "Interferon gamma",
                               "Latent TB infection",
                               "Acid-fast bacilli",
                               "Caseous necrosis",
                               "Caseation",
                               "Sputum smear",
                               "Pulmonary cavitation",
                               "Lung cavitary lesions",
                               "Rifampicin",
                               "Isoniazid"
                               ),
                  disease_name = "TB")

# Run data_mine_disease - RA
data_mine_disease(metadata_path = "data/metadata_archs4/01_archs4_filtered_metadata.rds",
                  keywords = c("Synovial inflammation",
                               "synovitis",
                               "Autoantibodies",
                               "Anti-CCP",
                               "Rheumatoid",
                               "Joint erosion",
                               "Bone erosion",
                               "Tumor necrosis factor alpha",
                               "TNF-alpha",
                               "Interleukin-6",
                               "Pannus formation",
                               "Invasive synovial tissue",
                               "Morning stiffness",
                               "DMARDs"
                               ), 
                  disease_name = "RA")

# Run data_mine_disease - COPD
data_mine_disease(metadata_path = "data/metadata_archs4/01_archs4_filtered_metadata.rds",
                  keywords = c("Emphysema",
                               "Alveolar destruction",
                               "Chronic bronchitis",
                               "Productive chronic cough",
                               "Airflow limitation",
                               "Obstructive ventilatory defect",
                               "Forced expiratory volume in 1 second",
                               "Hyperinflation",
                               "Air trapping",
                               "Neutrophilic inflammation",
                               "hypersecretion",
                               "Goblet cell hyperplasia",
                               "Exacerbations",
                               "Acute COPD flare-ups",
                               "Hypoxemia",
                               "Low arterial oxygen",
                               "COPD"
                               ),
                  disease_name = "COPD_no_dist",
                  max.distance = 0.05)

# Run data_mine_disease - IBD
data_mine_disease(metadata_path = "data/metadata_archs4/01_archs4_filtered_metadata.rds",
                  keywords = c("Crohn’s disease",
                               "Regional enteritis",
                               "Ulcerative colitis",
                               "Intestinal inflammation",
                               "Mucosal ulceration",
                               "Tumor necrosis factor alpha",
                               "TNF-alpha",
                               "Interleukin-23",
                               "Dysbiosis",
                               "Microbial imbalance",
                               "Intestinal permeability",
                               "Leaky gut",
                               "Hematochezia",
                               "Biologic therapy"
                               ), 
                  disease_name = "IBD")
