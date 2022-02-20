library(MicrobiomeBenchmarkDataAnalyses)
library(taxize)

## This table was created in of the articles under the
## vignettes directory
abundance_matrix <- read.table(
    file = "inst/extdata/Ravel_2013_bv_abundance_matrix.tsv",
    sep = "\t", header = TRUE, row.names = 1
)

taxa_names_original <- rownames(abundance_matrix)

## TM7 is in the NCBI taxonomy browser as Candidatus Saccharibacteria
tm7_pos <- grep("TM7", taxa_names)
taxa_names <- sub("^(\\w+)\\s.*.*", "\\1", taxa_names)
taxa_names[tm7_pos] <- "Candidatus Saccharibacteria"
taxonomy <- classification(taxa_names, db = "ncbi")
taxonomy_table <- taxize_classification_to_taxonomy_table(taxonomy)
taxonomy_table$species <- NULL
rownames(taxonomy_table) <- taxa_names_original
taxonomy_table$query <- NULL

write.table(
    x = taxonomy_table,
    file = "inst/extdata/Ravel_2013_bv_taxonomy_table.tsv",
    sep = "\t", col.names = TRUE, row.names = TRUE
)


