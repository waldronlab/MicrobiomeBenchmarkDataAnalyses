## This is a file for testing some features of the benchdamic package
## for performing benchmarking of differential abundance methods

library(benchdamic)
library(phyloseq)

ps <- ps_plaque_16S

sample_data(ps)$HMP_BODY_SUBSITE <- factor(sample_data(ps)$HMP_BODY_SUBSITE)
sample_data(ps)$HMP_BODY_SUBSITE <- relevel(
    x = sample_data(ps)$HMP_BODY_SUBSITE, 
    ref = "Subgingival Plaque"
)

## Prepare normalization methods
norm_pars <- tibble::tribble(
    ~fun, ~method,
    "norm_edgeR", "none",
    "norm_edgeR", "TMM",
    "norm_DESeq2", "poscounts",
    "norm_CSS", "median"
)

## Set normalization 
my_norm <- setNormalizations(fun = norm_pars$fun, method = norm_pars$method)

## Run normalization
ps <- runNormalizations(normalization_list = my_norm, object = ps, verbose = T)

## Zero-inflated negative binomial weights
zinbweights <- weights_ZINB(
    object = ps,
    K = 0,
    design = "~ 1"
)

## Prepare methods

# Some variables
grp <- "HMP_BODY_SUBSITE"
contrast <- c("Subgingival Plaque", "Supragingival Plaque")

# edgeR
my_edger <- set_edgeR(
    group_name = grp,
    design = as.formula(paste0("~", grp)),
    norm = "TMM",
    coef = 2
)

# DESeq2
my_deseq2 <- set_DESeq2(
    contrast = c(grp, contrast),
    design = as.formula(paste0("~", grp)),
    norm = "poscounts"
)

# limma 
my_limma <- set_limma( # I get a warning
    design = as.formula(paste0("~", grp)),
    norm = c("TMM", "CSSmedian"),
    coef = 2
)

# metagenomeSeq
my_metagenomeseq <- set_metagenomeSeq(
    design = as.formula(paste0("~", grp)),
    norm = "CSSmedian",
    coef = 2
)

# ALDEx2
my_aldex2 <- set_ALDEx2(
    conditions = grp,
    test = "t",
    norm = "none"
)

# corncob
my_corncob <- set_corncob(
    formula = as.formula(paste0("~", grp)),
    phi.formula = as.formula(paste0("~", grp)),
    formula_null = ~ 1,
    phi.formula_null = as.formula(paste0("~", grp)),
    test = "Wald",
    coefficient = paste0(grp, contrast[2]),
    norm = "none"
)

# MAST
my_mast <- set_MAST(
    rescale = "median",
    design = as.formula(paste0("~", grp)),
    coefficient = paste0(grp, contrast[2]),
    norm = "none"
)

# Seurat
my_seurat <- set_Seurat(
    test.use = "wilcox",
    contrast = c(grp, contrast),
    norm = "none"
)

my_methods <- c(
    my_edger, my_deseq2, my_limma, my_metagenomeseq, my_aldex2,
    my_corncob, my_mast, my_seurat
)


DA <- runDA(method_list = my_methods, object = ps, weights = zinbweights)


## A prior knowledge
mb <- microbial_metabolism
genera <- tax_table(ps)[,"GENUS"]
rownames(mb) <- mb$Genus
priorInfo <- data.frame(genera, Type = mb[genera, "Type"])
unknown_metabolism <- is.na(priorInfo$Type)
priorInfo[unknown_metabolism, "Type"] <- "Unknown"
priorInfo$Type <- factor(priorInfo$Type, levels = c("Aerobic","Anaerobic","F Anaerobic","Unknown"), labels = c("Aerobic","Anaerobic","F_Anaerobic","Unknown"))
priorInfo[, "newNames"] <- paste0(rownames(priorInfo), "|", 
    priorInfo[, "GENUS"])

direction <- c(
    edgeR.TMM = "logFC", 
    DESeq2.poscounts = "log2FoldChange",
    limma.CSSmedian = "logFC",
    limma.TMM = "logFC",
    metgenomeSeq.CSSmedian = "HMP_BODY_SUBSITESupragingival Plaque",
    ALDEx2.none = "effect",
    corncob.none = "Estimate",
    MAST.none = "logFC",
    Seurat.none = "avg_log2FC"
)
names(direction) <- NULL
## Visualization
enrichment <- createEnrichment(
    object = DA,
    priorKnowledge = priorInfo,
    enrichmentCol = "Type",
    namesCol = "newNames",
    slot = "pValMat",
    colName = "adjP",
    type = "pvalue",
    direction = direction,
    threshold_pvalue = 0.1,
    threshold_logfc = 0,
    top = NULL,
    alternative = "greater",
    verbose = TRUE
)

plotContingency(enrichment = enrichment, levels_to_plot = c("Aerobic", "Anaerobic"), method = "metagenomeSeq.CSSmedian")
