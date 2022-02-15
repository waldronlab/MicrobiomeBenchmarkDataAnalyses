## code to prepare `abundance_matrix_bacterial_vaginosis_dataset_1_ravel_2011.tsv` dataset goes here

library(magrittr)
library(taxize)
library(MicrobiomeBenchmarkDataAnalyses)

fname <- system.file(
    "extdata/abundance_matrix_bacterial_vaginosis_dataset_1_ravel_2011.tsv",
    package = "MicrobiomeBenchmarkDataAnalyses")

abMat <- read.table(fname, sep = "\t", row.names = 1)

taxa_names <- rownames(abMat) %>%
    sub("(Incertae_sedis_5_[1-2])$", "\\1dontdelete", .) %>%
    sub("_[0-9]+$", "", .) %>% 
    sub("_Incertae_Sedis", "", .) %>% 
    sub("_j$", "", .) %>% 
    sub("_c$", "", .) %>% 
    sub("_genera_incertae_sedis", "", .)

## Some taxa are Bacteria incertae sedis, so we have to add that manually
taxa_names[grepl("dontdelete", taxa_names)] <- "Bacteria incertae sedis"
    
taxonomy <- taxize::classification(taxa_names, db = "ncbi")
taxonomy_table <- taxize_classification_to_taxonomy_table(taxonomy)
taxonomy_table <- taxonomy_table[,-1]
rownames(taxonomy_table) <- rownames(abMat)
write.table(
    taxonomy_table,
    file = "inst/extdata/abundance_matrix_bacterial_vaginosis_dataset_1_ravel_2011.tsv",
    sep = "\t",
    
)


