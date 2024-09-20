library(benchdamic)
library(ggpubr)
data("ps_plaque_16S")
data("microbial_metabolism")

# Extract genera from the phyloseq tax_table slot
genera <- phyloseq::tax_table(ps_plaque_16S)[, "GENUS"]
# Genera as rownames of microbial_metabolism data.frame
rownames(microbial_metabolism) <- microbial_metabolism$Genus
# Match OTUs to their metabolism
priorInfo <- data.frame(genera,
                        "Type" =  microbial_metabolism[genera, "Type"])
# Unmatched genera becomes "Unknown"
unknown_metabolism <- is.na(priorInfo$Type)
priorInfo[unknown_metabolism, "Type"] <- "Unknown"
priorInfo$Type <- factor(priorInfo$Type)
# Add a more informative names column
priorInfo[, "newNames"] <- paste0(rownames(priorInfo), priorInfo[, "GENUS"])

# Add some normalization/scaling factors to the phyloseq object
my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
                             method = c("TMM", "CSS"))
ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
                                   object = ps_plaque_16S)

# Initialize some limma based methods
my_limma <- set_limma(design = ~ 1 + RSID + HMP_BODY_SUBSITE, 
                      coef = "HMP_BODY_SUBSITESupragingival Plaque",
                      norm = c("TMM", "CSS"))

# Make sure the subject ID variable is a factor
phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
    phyloseq::sample_data(ps_plaque_16S)[["RSID"]])

# Perform DA analysis
Plaque_16S_DA <- runDA(method_list = my_limma, object = ps_plaque_16S)



## Fitlering by adjusted P-value of 0.1 and 0 logFC
enrichment1 <- createEnrichment(
    object = Plaque_16S_DA,
    priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
    slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
    threshold_pvalue = 0.1, threshold_logfc = 0, top = NULL, verbose = TRUE
)
p1 <- plotEnrichment(enrichment1, enrichmentCol = "Type")

## Filtering by raw P-value of 0.1 and 0 logFC
enrichment2 <- createEnrichment(
    object = Plaque_16S_DA,
    priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
    slot = "pValMat", colName = "rawP", type = "pvalue", direction = "logFC",
    threshold_pvalue = 0.1, threshold_logfc = 0, top = NULL, verbose = TRUE
)
p2 <- plotEnrichment(enrichment2, enrichmentCol = "Type")


## Filtering by threshold of p-value of 0? (ignored?) and 0 logFC
enrichment3 <- createEnrichment(
    object = Plaque_16S_DA,
    priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
    slot = "statInfo", colName = "logFC", type = "logfc", direction = "logFC",
    threshold_pvalue = 0, threshold_logfc = 0, top = NULL, verbose = TRUE
)
p3 <- plotEnrichment(enrichment3, enrichmentCol = "Type")

pL <- list(p1, p2, p3)
ggarrange(plotlist = pL, nrow = 1)





