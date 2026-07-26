# Sourcing functions
source("R/00_functions.R")

# Read ChEMBL
chembl_db <- readr::read_tsv("/home/databases/ChEMBL/02_ChEMBL_drug_gene_interactions.tsv")

# Load saved mapping of ensg to symbols
gene_map <- readRDS("data/ensg_to_symbol/gene_mapping.rds")

# Intersection with Edgingtons p-val aggregation results
multiLevel_res <- run_chembl_targets(
  data = readRDS("data/formal_pvalue_fgseaMulti/fgsea_sp_edgintons_V2.rds"),
  method = "fgseaMultilevel",
  chembl = df,
  gene_map = gene_map,
  padj_cutoff = 0.05
)

comp_scores <- multiLevel_res$compound_scores

top_50 <- multiLevel_res$compound_scores |> 
  dplyr::slice_head(n = 50)

multiLevel_res$hits |> 
  dplyr::slice_head(n = 10)

multiLevel_res$targets |> 
  dplyr::arrange(desc(NES)) |> 
  dplyr::slice_head(n = 10)

multiLevel_res$targets |> dplyr::arrange(NES)


library(readr)
library(dplyr)
library(tidyr)

compoundinfo <- read_tsv("/home/databases/CMap_links_2020/compoundinfo_beta.txt", show_col_types = FALSE)


library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)

# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------

get_chembl_json <- function(url, query = NULL) {
  
  response <- GET(
    url,
    query = query,
    user_agent("R ChEMBL API script")
  )
  
  stop_for_status(response)
  
  fromJSON(
    content(
      response,
      as = "text",
      encoding = "UTF-8"
    ),
    simplifyVector = FALSE
  )
}


# ------------------------------------------------------------
# 1. Download ALL human targets
# ------------------------------------------------------------

target_url <- paste0(
  "https://www.ebi.ac.uk/chembl/api/data/",
  "target.json"
)

human_targets <- list()

limit <- 1000
offset <- 0

repeat {
  
  message(
    "Downloading human targets — offset: ",
    offset
  )
  
  data <- get_chembl_json(
    target_url,
    query = list(
      organism__exact = "Homo sapiens",
      limit = limit,
      offset = offset
    )
  )
  
  human_targets <- c(
    human_targets,
    data$targets
  )
  
  total_count <- data$page_meta$total_count
  
  message(
    "Downloaded human targets: ",
    length(human_targets),
    " / ",
    total_count
  )
  
  offset <- offset + limit
  
  if (offset >= total_count) {
    break
  }
}


# ------------------------------------------------------------
# 2. Create human target lookup table
# ------------------------------------------------------------

human_target_df <- tibble(
  target_chembl_id = map_chr(
    human_targets,
    ~ .x$target_chembl_id %||% NA_character_
  ),
  
  hgnc_symbol = map_chr(
    human_targets,
    function(target) {
      
      components <- target$target_components
      
      if (
        is.null(components) ||
        length(components) == 0
      ) {
        return(NA_character_)
      }
      
      xrefs <- map(
        components,
        function(component) {
          component$target_component_xrefs %||% list()
        }
      ) |>
        flatten()
      
      hgnc <- keep(
        xrefs,
        function(x) {
          
          identical(
            toupper(
              x$xref_src_db %||% ""
            ),
            "HGNC"
          )
        }
      )
      
      if (length(hgnc) == 0) {
        return(NA_character_)
      }
      
      as.character(
        hgnc[[1]]$xref_name %||%
          NA_character_
      )
    }
  )
)


# ------------------------------------------------------------
# 3. Download ALL mechanism records
# ------------------------------------------------------------

mechanism_url <- paste0(
  "https://www.ebi.ac.uk/chembl/api/data/",
  "mechanism.json"
)

mechanisms <- list()

limit <- 1000
offset <- 0

repeat {
  
  message(
    "Downloading mechanisms — offset: ",
    offset
  )
  
  data <- get_chembl_json(
    mechanism_url,
    query = list(
      limit = limit,
      offset = offset
    )
  )
  
  mechanisms <- c(
    mechanisms,
    data$mechanisms
  )
  
  total_count <- data$page_meta$total_count
  
  message(
    "Downloaded mechanisms: ",
    length(mechanisms),
    " / ",
    total_count
  )
  
  offset <- offset + limit
  
  if (offset >= total_count) {
    break
  }
}


# ------------------------------------------------------------
# 4. Create mechanism dataframe
# ------------------------------------------------------------

mechanism_df <- tibble(
  
  compound_chembl_id = map_chr(
    mechanisms,
    ~ .x$molecule_chembl_id %||%
      NA_character_
  ),
  
  target_chembl_id = map_chr(
    mechanisms,
    ~ .x$target_chembl_id %||%
      NA_character_
  ),
  
  phase = map_dbl(
    mechanisms,
    ~ as.numeric(
      .x$max_phase %||%
        NA_real_
    )
  )
)


# ------------------------------------------------------------
# 5. Keep only mechanisms linked to human targets
# ------------------------------------------------------------

df <- mechanism_df |>
  inner_join(
    human_target_df,
    by = "target_chembl_id"
  ) |>
  select(
    compound_chembl_id,
    target_chembl_id,
    hgnc_symbol,
    phase
  )


# ------------------------------------------------------------
# 6. Inspect results
# ------------------------------------------------------------

head(df)

nrow(df)

summary(df$phase)

sum(is.na(df$hgnc_symbol))


# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------

write.csv(
  df,
  "chembl_human_mechanisms.csv",
  row.names = FALSE
)

message("Done!")