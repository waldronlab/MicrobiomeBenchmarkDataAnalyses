
library(MicrobiomeBenchmarkData)
library(phyloseq)
library(mia)
library(benchdamic)
library(purrr)

# Import data -------------------------------------------------------------

tse <- getDataset("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)[[1]]

grp <- "hmp_body_subsite"
conditions <- c("subgingival_plaque", "supragingival_plaque")

# prior knowledge (biological info) ---------------------------------------

row_data <- as.data.frame(rowData(tse))
taxaNames <- rownames(row_data)
newNames <- paste0(taxaNames, "|", row_data$GENUS) 
priorInfo <- data.frame(
    taxaNames = taxaNames,
    genus = row_data$GENUS,
    newNames = newNames,
    type = row_data$BIOSIS 
)

rownames(priorInfo) <- taxaNames

# Convert to phyloseq -----------------------------------------------------

ps <- makePhyloseqFromTreeSummarizedExperiment(tse)

## This step places the first condition as reference
sample_data(ps)[[grp]] <- factor(
    sample_data(ps)[[grp]], levels = conditions
)

# Normalization -----------------------------------------------------------
norm_par <- tibble::tribble(
    ~ fun, ~ method,
    "norm_edgeR", "TMM",
    "norm_edgeR", "none",
    "norm_DESeq2", "poscounts",
    "norm_TSS", "TSS",
    "norm_CSS", "median"
)

norms <- setNormalizations(fun = norm_par$fun, method = norm_par$method)

ps <- 
    runNormalizations(normalization_list = norms, object = ps, verbose = FALSE)


# Estimate weights --------------------------------------------------------

## This section is only for DESeq2, edgeR, and limma-voom
## Calculate zinbwave weights
zinbWeights <- weights_ZINB(object = ps, design = grp)


# DESEq2 ------------------------------------------------------------------

output_DESeq2 <- DA_DESeq2(
    object = ps, pseudo_count = FALSE, 
    design = as.formula(paste0("~", grp)),
    contrast = c(grp, conditions[2], conditions[1]), 
    norm = "poscounts", 
    verbose = FALSE
)

output_DESeq2_zinbweights <- DA_DESeq2(
    object = ps, pseudo_count = FALSE,
    design = as.formula(paste0("~", grp)),
    contrast = c(grp, conditions[2], conditions[1]),
    norm = "poscounts", weights = zinbWeights,
    verbose = FALSE
)

# edgeR -------------------------------------------------------------------

output_edgeR <- DA_edgeR(
    object = ps, pseudo_count = FALSE, group_name = grp,
    design = as.formula(paste0("~", grp)), norm = "TMM",
    verbose = FALSE
)

output_edgeR_zinbweihgts <- DA_edgeR(
    object = ps, pseudo_count = FALSE, group_name = grp,
    design = as.formula(paste0("~", grp)), norm = "TMM",
    weights = zinbWeights,
    verbose = FALSE
)

# limma benchdamic ---------------------------------------------------------

output_limma <- DA_limma(
    object = ps, pseudo_count = FALSE, design = grp, norm = "TMM",
    verbose = FALSE,
)

output_limma_zinbweights <- DA_limma(
   object = ps, pseudo_count = FALSE, design = grp, norm = "TMM" ,
   weights = zinbWeights, verbose = FALSE 
)


# Combine results in single object ----------------------------------------

DA_results <- list(
    output_DESeq2,
    output_DESeq2_zinbweights,
    output_edgeR,
    output_edgeR_zinbweihgts,
    output_limma,
    output_limma_zinbweights
)
names(DA_results) <- map_chr(DA_results, ~ .x$name)



# Direction ---------------------------------------------------------------

direction <- c(
    DESeq2.poscounts = "log2FoldChange",
    DESeq2.poscounts.weighted = "log2FoldChange",
    edgeR.TMM = "logFC",
    edgeR.TMM.weighted = "logFC",
    limma.TMM = "logFC",
    limma.TMM.weighted = "logFC"
)

# Enrichment --------------------------------------------------------------

enrichment <- createEnrichment(
    object = DA_results,
    priorKnowledge = priorInfo,
    enrichmentCol = "type",
    namesCol = "newNames",
    slot = "pValMat", colName = "adjP", type = "pvalue",
    direction = direction,
    threshold_pvalue = 0.1,
    threshold_logfc = 0,
    top = NULL,
    alternative = "greater",
    verbose = FALSE 
)


plotEnrichment(
    enrichment = enrichment, 
    enrichmentCol = "type",
    levels_to_plot = c("Aerobic", "Anaerobic", "F Anaerobic")
)





