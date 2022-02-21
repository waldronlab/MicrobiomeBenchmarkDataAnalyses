library(bugphyzz)
library(purrr)
library(dplyr)

## Aerophilicity data from bugphyzz
aer <- as_tibble(physiologies("aerophilicity")[[1]])
aer <- aer %>% 
    filter(
        Rank == "genus",
        Attribute %in% c("aerobic", "anaerobic", "facultatively anaerobic")
    ) %>% 
    select(Taxon_name, Attribute) %>% 
    distinct()

readr::write_tsv(
    aer, "inst/extdata/bugphyzz_aerophilicity.tsv"
)

## Gram stain data from bugphyzz
gs <- as_tibble(physiologies("gram stain")[[1]])
gs <- filter(gs, Rank == "genus") %>% 
    select(Taxon_name, Attribute) %>% 
    distinct()

duplicated_taxa_in_gs <- unique(gs[duplicated(gs$Taxon_name), ]$Taxon_name)
duplicated_gs <- gs %>% 
    filter(Taxon_name %in% duplicated_taxa_in_gs)
gs_filtered <- gs %>% 
    filter(
        !Taxon_name %in% duplicated_gs,
        grepl("(negative|positive)$", Attribute)
    )

readr::write_tsv(
    gs_filtered, "inst/extdata/bugphyzz_gram_stain.tsv"
)

## Aerophilicty data from nychanes biosis
biosis_url <- "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/biosis.tsv"
biosis <- readr::read_tsv(biosis_url) %>% 
    magrittr::set_colnames(c("genera", "biosis"))

readr::write_tsv(
    biosis, "inst/extdata/nychanes_biosis.tsv"
)
