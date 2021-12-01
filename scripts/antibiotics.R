## In this script I test several datasets of antibiotics

library(dplyr)
library(tibble)
library(magrittr)
library(S4Vectors)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(bugphyzz)
library(ggplot2)
library(purrr)
source("R/functions_da.R")

file_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-020-00991-x/MediaObjects/40168_2020_991_MOESM5_ESM.xlsx"
temp_file <- tempfile()
download.file(file_url, temp_file)

count_matrix <- readxl::read_xlsx(temp_file, sheet = "OUT_table") %>% 
    column_to_rownames(var = "OTU") %>%
    as.data.frame() %>% 
    as.matrix()

col_data <- readxl::read_xlsx(temp_file, sheet = "Sample_data") %>% 
    column_to_rownames(var = "Sample_ID") %>% 
    as.data.frame() %>% 
    DataFrame()

row_data <- readxl::read_xlsx(temp_file, sheet = "Taxa", col_names = FALSE) %>% 
    set_colnames(c("taxa", "taxonomy")) %>% 
    column_to_rownames(var = "taxa") %>% 
    as.data.frame() %>% 
    DataFrame()

tse <- TreeSummarizedExperiment(
    assays = SimpleList(counts = count_matrix),
    colData = col_data,
    rowData = row_data
)

## FMT stands for fecal microbiota transplantation
table(colData(tse))


# Get taxonomy ranks ------------------------------------------------------
## Get taxonomy ranks and NCBI IDs
tax = rowData(tse)$taxonomy
rank_names <- map(tax, ~ {
    str <- stringr::str_split(.x, ";")[[1]]
    str <- tail(str, 1)
    str
    }
) %>% 
    flatten_chr() %>% 
    set_names(rownames(tse))
unclassified_terms <- c("Bacteria", "Unclassified")
sum(rank_names %in% unclassified_terms) / length(rank_names) * 100

tax_summary = taxizedb::name2taxid(rank_names, db = "ncbi", out_type = "summary")

append_taxa <- as.data.frame(rank_names) %>%
    rownames_to_column(var = "taxa") %>% 
    left_join(tax_summary, by = c("rank_names" = "name"))

ranks <- map_chr(append_taxa$id, ~{
    tryCatch(
        error = function(e) NA,
        taxizedb::taxid2rank(.x, db = "ncbi")
    )
})
append_taxa$rank <- ranks

new_append_taxa <- append_taxa %>% 
    group_by(taxa) %>% 
    arrange(rank_names) %>% 
    slice_head()
row_data <- rowData(tse) %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "otu")

new_row_data <- left_join(row_data, new_append_taxa, by = c("otu" = "taxa")) %>% 
    column_to_rownames(var = "otu") %>% 
    as.data.frame() %>% 
    DataFrame()
rowData(tse) <- new_row_data


# Analysis ----------------------------------------------------------------


## Comparisons
## DSS_FMT (35) vs DSS_FMT_Metronidazole (6)
## DSS_FMT (35) vs DSS_FMT_Streptomycin (6)
## DSS_FMT (35) vs DSS_FMT_Vancomycin (3)

## Select
van_sam_1 <- colnames(tse[,colData(tse)$Treatment == "DSS_FMT_Vancomycin"])
van_sam_2 <- colnames(tse[,sample(which(colData(tse)$Treatment == "DSS_FMT"), 3)])

tse_van <- tse[,c(van_sam_1, van_sam_2)]
grp <- "Treatment"
ref <- "DSS_FMT"

desesq_results <- deseq2Poscounts(tse_van, grp, ref)
edger_results <- edgerTmm(tse_van, grp, ref)
limma_results <- limmaVoomTmm(tse_van, grp, ref)
ancombc_results <- ancom_bc(tse_van, grp, ref) # there are a few warnings
#corncob_results <- corncob(tse_van, grp, ref) # there are a few errors
# mgs_results <- metagenomeSeq(tse_van, grp, ref) # there are a few errors
aldex2_results <- aldex2(tse_van, grp, ref)

p1 <- volcano_plot(desesq_results) + ggtitle("deseq2")
p2 <- volcano_plot(edger_results) + ggtitle("edger")
p3 <- volcano_plot(limma_results) + ggtitle("limma")
p4 <- volcano_plot(ancombc_results) + ggtitle("ancom-bc")
# p5 <- volcano_plot(corncob_results) + ggtitle("corncob")
p6 <- volcano_plot(aldex2_results) + ggtitle("aldex2")


ggpubr::ggarrange(p1, p2, p3, p4, p6)


dres <- desesq_results %>% 
    left_join(new_append_taxa, by = c("TAXA" = "taxa")) %>% 
    dplyr::mutate(
        sig = ifelse(
            .data$ADJPVAL > 0.1 | is.na(.data$ADJPVAL), "unchanged", "DA"
        )
    ) 
    
dres_subset <- dres %>% 
    dplyr::filter(sig == "DA", FC > 0)

dres %>% 
    ggplot2::ggplot(ggplot2::aes(.data$FC, -log10(.data$PVAL))) +
    ggplot2::geom_point(
        ggplot2::aes(color = .data$sig), shape = 1, size = 2, stroke = 0.8
    ) +
    ggplot2::labs(y = "-log10(P-value)", x = "log2(Fold Change)") +
    ggrepel::geom_text_repel(
        data = dres_subset, 
        mapping = aes(FC, -log10(PVAL), label = rank_names)
        ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
