## code to prepare `method_classification` dataset goes here

url <- 'https://docs.google.com/spreadsheets/d/e/2PACX-1vTrwtp0Br2q79Iq_Ar1rQeq3cll4aab0l6lkRBqcX0I7EbWnunRkQVkIYUG58cKlYfoSFY4UCPYUdef/pub?gid=0&single=true&output=tsv'

fname <- system.file(
    'extdata/method_classification.tsv',
    package = 'MicrobiomeBenchmarkDataAnalyses'
)

methods_classification <- read.table(
    file = fname, header = TRUE, sep = '\t', row.names = NULL
)


usethis::use_data(methods_classification, overwrite = TRUE, internal = TRUE)

