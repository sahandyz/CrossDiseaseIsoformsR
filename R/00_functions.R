# Function to extract metadata from ARCHS4 database h5-file - inputs described below
## *metadata_path* is the meta information such as title and characteristics columns from ARCHS4
## *keywords* are words chosen to be associated with the disease in question
## *max.distance* is the fraction specified as distance that each word that can differ from the original word (insertions, deletions, and substitutions).

data_mine_diseases <- function(metadata_path, keywords, max.distance = 0.2) {
  metadata <- read_rds(metadata_path) |> 
    filter(
      rowSums(sapply(keywords, function(keyword) {
        agrepl(pattern = keyword, x = title, ignore.case = TRUE, max.distance = max.distance) |
          agrepl(pattern = keyword, x = characteristics, ignore.case = TRUE, max.distance = max.distance)
      })) > 0
    )
}

