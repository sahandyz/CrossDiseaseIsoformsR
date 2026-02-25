# nice /home/ctools/opt/R-4.4.1/bin/R
    
    library(tidyverse)
    library(rhdf5)
    
    ### Define paths to h5 files
    archs4gene <- "/home/databases/archs4/v2.latest/human_gene_v2.latest.h5"
    archs4TxCount <- "/home/databases/archs4/v2.latest/human_transcript_v2.latest.h5"
    file.exists(archs4gene)
    file.exists(archs4TxCount)
    
    
    ### Extract hq samples/studies
    if(TRUE) {
        ### Extract meta data
        arch4meta <- 
            data.frame(
                study             = h5read(archs4gene, "/meta/samples/series_id"),
                id                = h5read(archs4gene, "/meta/samples/geo_accession"),
                organism          = h5read(archs4gene, "/meta/samples/organism_ch1"),
                source            = h5read(archs4gene, "/meta/samples/source_name_ch1"),
                title             = h5read(archs4gene, "/meta/samples/title"),
                characteristics   = h5read(archs4gene, "/meta/samples/characteristics_ch1"),
                singleCellProb    = h5read(archs4gene, "/meta/samples/singlecellprobability"),
                type              = h5read(archs4gene, "/meta/samples/type"),
                library_selection = h5read(archs4gene, "/meta/samples/library_selection"),
                channel_count     = h5read(archs4gene, "/meta/samples/channel_count"),
                readsAlligned     = h5read(archs4gene, "/meta/samples/alignedreads"), # v2 beta
                submissionDate    = h5read(archs4gene, "/meta/samples/submission_date"),
                stringsAsFactors  = FALSE
            )
        
        ### Filter data
        arch4metaFilt <- 
          arch4meta |> 
          as_tibble() |> 
          group_by(study) |> 
          mutate(
            medianScProb = median(singleCellProb),
            n_samples = n()
          ) |> 
          filter(
            organism == 'Homo sapiens',
            type == 'SRA',
            library_selection == 'cDNA',
            medianScProb <= 1/3,
            readsAlligned >= 10e6
          ) |> 
          ungroup() |> 
          arrange(study)
        
        nrow(arch4metaFilt)
        # 143927
        
        n_distinct(arch4metaFilt$study)
        
    }
    # arch4metaFilt
    
saveRDS(object = arch4metaFilt, file = "data/01_archs4_filtered_metadata.rds")
