## Code to prepare Ravel_2011_bv_taxonomy_table.tsv

library(magrittr)
library(taxize)

fname <- system.file(
    ## The code to prepare this file is an article within the vignettes directory.
    ## Look for the save abundance matrix code chunk in the
    ## Ravel_2011_16S_nugent_score.Rmd file.
    "extdata/Ravel_2011_bv_abundance_matrix.tsv",
    package = "MicrobiomeBenchmarkDataAnalyses"
)

abMat <- read.table(fname, sep = "\t", row.names = 1)

taxa_names <- rownames(abMat) %>%
    sub("(Incertae_sedis_5_[1-2])$", "\\1dontdelete", .) %>%
    sub("_[0-9]+$", "", .) %>% 
    sub("_Incertae_Sedis", "", .) %>% 
    sub("_j$", "", .) %>% 
    sub("_c$", "", .) %>% 
    sub("_genera_incertae_sedis", "", .) %>% 
    sub("^TM7$", "Candidatus Saccharibacteria", .)

## Some taxa are Bacteria incertae sedis, so we have to add that manually
taxa_names[grepl("dontdelete", taxa_names)] <- "Bacteria incertae sedis"
    
taxonomy <- taxize::classification(taxa_names, db = "ncbi")
taxonomy_table <- taxize_classification_to_taxonomy_table(taxonomy)
taxonomy_table <- taxonomy_table[,-1]
rownames(taxonomy_table) <- rownames(abMat)
write.table(
    taxonomy_table,
    file = "inst/extdata/Ravel_2011_bv_taxonomy_table.tsv",
    sep = "\t",
)


