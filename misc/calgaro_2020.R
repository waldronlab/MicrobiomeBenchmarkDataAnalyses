## Script for analysis for the paper of calgaro

packages <- c(
    "MicrobiomeBenchmarkData",
    "ggplot2",
    "dplyr",
    "ggpubr"
)

for (pkg in packages) suppressMessages(library(pkg, character.only = TRUE))
source("R/functions_da.R")

## Load the calgaro dataset
tse <- getDataset("HMP_2012_16S_gingival_V35_subset", FALSE)[[1]]
grp <-  "hmp_body_subsite"
ref <- "subgingival_plaque"
gg_title <- function(method) {
    ggplot2::ggtitle(
        paste0("Sugingival (reference) vs Supragingival - ", method)
    )
}


# DA methods --------------------------------------------------------------

## 1. DESEQ2
dds <- deseq2Poscounts(tse, grp, ref)
dds_p <- dds %>% 
    volcano_plot() +
    gg_title("DESeq2-poscounts") 

## 2. edgeR
dge <- edgerTmm(tse, grp, ref) 
dge_p <- dge %>% 
    volcano_plot() +
    gg_title("edgeR-TMM")

## 3. Limma voom
lv <- limmaVoomTmm(tse, grp, ref)
lv_p <- lv %>% 
    volcano_plot() +
    gg_title("Limma-Voom-TMM")

## 4. ALDEx2
al2 <- aldex2(tse, grp, ref)
al2_p <- al2 %>% 
    volcano_plot() + 
    gg_title("ALDEx2")

## 5. MetagenomoeSeq
mgs <- metagenomeSeq(tse, grp, ref)
mgs_p <- mgs %>% 
    volcano_plot() +
    gg_title("MetagenomeSeq")

## 6. Corncob
cc <- corncob(tse, grp, ref, fdr_cutoff = 0.1)
cc_p <- cc %>% 
    volcano_plot() +
    gg_title("corncob")

## 7. ANCOM-BC
ambc <- ancom_bc(tse, grp, ref)
ambc_p <- ambc %>% 
    volcano_plot() +
    gg_title("ANCOM-BC")

ggarrange(
    dds_p, dge_p, lv_p, al2_p, mgs_p, cc_p, ambc_p,
    common.legend = TRUE
)


# Combine datasets --------------------------------------------------------

row_data <- rowData(tse) %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "TAXA") %>% 
    select(TAXA, GENUS, BIOSIS)

results <- list(deseq2 = dds, edger = dge, limma_voom = lv,
    aldex2 = al2, metagenomeseq = mgs, cornconb = cc, ancombc = ambc) %>% 
    bind_rows(.id = "method") %>% 
    mutate(
        SIG = case_when(
            ADJPVAL <= 0.1 & FC > 0 ~ "over",
            ADJPVAL <= 0.1 & FC < 0 ~ "under",
            TRUE ~ "unchanged"
            )
    ) %>% 
    left_join(row_data, by = "TAXA")

results %>% 
    filter(SIG != "unchanged") %>% 
    ggplot(aes(method, FC)) +
    geom_boxplot() +
    geom_point(aes(color = SIG), position = "jitter") +
    labs(x = "DA method", y = "log2(Fold Change)") +
    theme_bw()


results %>% 
    filter(SIG != "unchanged") %>% 
    count(method, SIG) %>% 
    ggplot(aes(method, n)) +
    geom_col(aes(fill = SIG), position = position_dodge2(0.9, preserve = "single")) +
    coord_flip() +
    theme_bw()

results %>% 
    filter(SIG != "unchanged") %>% 
    count(method, SIG, BIOSIS) %>% 
    mutate(n = ifelse(SIG == "under", -n, n )) %>% 
    ggplot(aes(method, n)) +
    geom_col(aes(fill = BIOSIS), position = position_dodge(0.9, preserve = "single")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


# positives and negatives -------------------------------------------------

positives_1 <- results %>% 
    mutate(method = as.factor(method)) %>%
    filter(BIOSIS == "Aerobic", SIG == "over") %>% 
    count(method, .drop = FALSE)
positives_2 <- results %>% 
    mutate(method = as.factor(method)) %>%
    filter(BIOSIS == "Anaerobic", SIG == "under") %>% 
    count(method, .drop = FALSE)
negatives_1 <- results %>% 
    mutate(method = as.factor(method)) %>%
    filter(BIOSIS == "Aerobic", SIG == "under") %>% 
    count(method, .drop = FALSE)
negatives_2 <- results %>% 
    mutate(method = as.factor(method)) %>%
    filter(BIOSIS == "Anaerobic", SIG == "over") %>% 
    count(method, .drop = FALSE)

ranks <- (positives_1$n + positives_2$n) - (negatives_1$n + negatives_2$n)
names(ranks) <- levels(positives_1$method)
ranks <- ranks %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "method") %>% 
    magrittr::set_colnames(c("method", "rank"))

ranks %>% 
    ggplot(aes(reorder(method, rank), rank)) +
    geom_col(fill = "firebrick") +
    labs(x = "DA method", y = "Putative TP - Putative FP",
        title = "Putative Positives - Putative Negatives",
        subtitle = "Aerobic vs Anaerobic - Adjusted P-value 0.1") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
