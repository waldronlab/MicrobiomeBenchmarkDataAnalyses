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

# Count TPs and FPs, from the top 1 to the top 20 features.
# As direction is supplied, features are ordered by "logFC" absolute values.


## is ranking based on raw pvalue?
positives1 <- createPositives(
    object = Plaque_16S_DA,
    priorKnowledge = priorInfo, enrichmentCol = "Type", 
    namesCol = "newNames", 
    slot = "pValMat", colName = "rawP", 
    type = "pvalue", direction = "logFC",
    threshold_pvalue = 1, 
    threshold_logfc = 0, 
    top = 1:20, alternative = "greater", 
    verbose = FALSE,
    TP = list(c("DOWN Abundant", "Anaerobic"), c("UP Abundant", "Aerobic")),
    FP = list(c("DOWN Abundant", "Aerobic"), c("UP Abundant", "Anaerobic"))
)
p1 <- plotPositives(positives = positives1)


## Is ranking based on logFC?
positives2 <- createPositives(
    object = Plaque_16S_DA,
    priorKnowledge = priorInfo, enrichmentCol = "Type", 
    namesCol = "newNames", 
    slot = "statInfo", colName = "logFC", 
    type = "logfc", direction = "logFC",
    threshold_pvalue = 1, 
    threshold_logfc = 0, 
    top = 1:20, alternative = "greater", 
    verbose = FALSE,
    TP = list(c("DOWN Abundant", "Anaerobic"), c("UP Abundant", "Aerobic")),
    FP = list(c("DOWN Abundant", "Aerobic"), c("UP Abundant", "Anaerobic"))
)
p2 <- plotPositives(positives = positives2)


pL <- list(p1, p2)
ggarrange(plotlist = pL, nrow = 1)

sessioninfo::session_info()