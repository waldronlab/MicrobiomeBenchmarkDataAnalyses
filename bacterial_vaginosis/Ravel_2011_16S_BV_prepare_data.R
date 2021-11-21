library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(purrr)

supp_tbl_s3_url <- "https://www.pnas.org/highwire/filestream/603190/field_highwire_adjunct_files/4/st04.xlsx"
supp_tbl_s3_file <- tempfile()
download.file(url = supp_tbl_s3_url, destfile = supp_tbl_s3_file)

## Import sample metadata
sample_metadata <- readxl::read_xlsx(supp_tbl_s3_file, range = "A3:G397")

## Import count data
otu_table <- readxl::read_xlsx(supp_tbl_s3_file, range = "I3:IU397")
otu_table <- bind_cols(sample_metadata["Subject ID"], otu_table)

## relative abundance matrix
relative_abundance_matrix <- otu_table %>% 
    column_to_rownames(var = "Subject ID") %>% 
    as.data.frame() %>% 
    as.matrix() %>% 
    t()

## Count matrix
count_matrix <- 
    t(t(relative_abundance_matrix)*sample_metadata$`Total number readsd`/100)

## col_data for SummarizedExperiment
col_data <- sample_metadata %>% 
    column_to_rownames(var = "Subject ID") %>% 
    as.data.frame() %>% 
    S4Vectors::DataFrame()

se <- SummarizedExperiment(
    assays = S4Vectors::SimpleList(
        counts = count_matrix, relative_abundance = relative_abundance_matrix
    ),
    colData = col_data
)

# taxa_names <- rownames(se) %>% 
#     sub("^L. ", "Lactobacillus ", .) %>% 
#     sub("_.+$", "", .)
# 
# 
# taxonomies <- vector("list", length(taxa_names))
# names(taxonomies) <- taxa_names
# for (i in seq_along(taxonomies)) {
#     taxonomies[[i]] <- 
#         tryCatch(
#             error = function(e) e,
#             taxizedb::classification(taxa_names[i], db = "ncbi")[[1]]
#         )
#         
# }
# 
# positions_with_error <- which(map_lgl(taxonomies, ~ "error" %in% class(.x)))
# 
# unique_names <- unique(names(positions_with_error))
# 
# unique_name1 <- taxize::classification(unique_names[1], db = "ncbi")
# unique_name2 <- taxize::classification(unique_names[2], db = "ncbi")
# unique_name3 <- taxize::classification(unique_names[3], db = "ncbi")
# unique_name4 <- taxize::classification(unique_names[4], db = "ncbi")
# unique_name5 <- taxize::classification(unique_names[5], db = "ncbi")
# unique_name6 <- taxize::classification(unique_names[6], db = "ncbi")
# 
# unique_taxonomies <- list(
#     unique_name1[[1]],
#     unique_name2[[1]],
#     unique_name3[[1]],
#     unique_name4[[1]],
#     unique_name5[[1]],
#     unique_name6[[1]]
#     
# )
# names(unique_taxonomies) <- unique_names
# 
# for (i in seq_along(taxonomies)) {
#     for (j in seq_along(unique_taxonomies)) {
#         name1 <- names(taxonomies)[i]
#         name2 <- names(unique_taxonomies)[j]
#         if (name1 == name2) {
#             taxonomies[[i]] <- unique_taxonomies[[j]]
#         }
#     }
# 
# }
# 
# 
# 
# 
# taxonomy_table <- map(taxonomies, ~ {
#     
#     valid_ranks <- c(
#         "superkingdom", "phylum", "class", "order", "family", "genus"
#     )
#     
#     if (is.logical(.x)) {
#         n <- length(valid_ranks)
#         mat <- matrix(rep(NA, n), nrow = 1)
#         colnames(mat) <- valid_ranks
#         df <- tibble::as_tibble(as.data.frame(mat))
#     } else if (is.data.frame(.x)) {
#         x <- .x[.x$rank %in% valid_ranks, c("rank", "name")]
#         df <- as.data.frame(matrix(x$name, ncol = length(x$name)))
#         colnames(df) <- x$rank
#         df <- tibble::as_tibble(df)
#         df
#     }
# }) %>% 
#     bind_rows(.id = "taxa_name")
#     
rownames(se) <- sub("^L. ", "Lactobacillus ", rownames(se))
# taxonomy_table$taxa_name <- rownames(se)
# 
readr::write_tsv(
    taxonomy_table, "bacterial_vaginosis/Ravel_2011_16S_BV_taxonomy_table.tsv"
)


taxonomy_table_tsv <- readr::read_tsv(
    file = "bacterial_vaginosis/Ravel_2011_16S_BV_taxonomy_table.tsv"
)

row_data <- taxonomy_table_tsv %>% 
    column_to_rownames(var = "taxa_name") %>% 
    as.data.frame() %>% 
    S4Vectors::DataFrame()

rowData(se) <- row_data
saveRDS(se, file = "bacterial_vaginosis/Ravel_2011_16S_BV_se.rds")
