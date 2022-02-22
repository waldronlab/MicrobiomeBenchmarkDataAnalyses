
library(MicrobiomeBenchmarkData)
library(MicrobiomeBenchmarkDataAnalyses)
library(mia)
library(benchdamic)
library(phyloseq)
library(magrittr)

se <- getDataset("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)[[1]]

## Get prior knowledge
prior_knowledge <- as.data.frame(rowData(se))
prior_knowledge$newNames <- 
    paste0(rownames(prior_knowledge), "|", prior_knowledge$GENUS)
prior_knowledge <- subset(prior_knowledge, select = c(newNames, BIOSIS))
prior_knowledge$BIOSIS[is.na(prior_knowledge$BIOSIS)] <- "Unknown"

ps <- makePhyloseqFromTreeSummarizedExperiment(se)

## Set normalization
ps_normalized <- runNormalizations(set_normalizations, object = ps)


grp <- "hmp_body_subsite"
contrast <- c("subgingival_plaque", "supragingival_plaque")


## set methods
my_edger <- set_edgeR(
    group_name = grp,
    design = as.formula(paste0("~", grp)),
    norm = "TMM",
    coef = 2
)

# DESeq2
my_deseq2 <- set_DESeq2(
    contrast = c(grp, "supragingival_plaque", "subgingival_plaque"),
    design = as.formula(paste0("~", grp)),
    norm = "poscounts"
)

my_zinq <- list(
    zinq.1 = list(method = "zinq", grp = "hmp_body_subsite", ref = "supragingival_plaque")
)


my_methods <- c(my_edger, my_deseq2, my_zinq)

output <- runDA(my_methods, ps_normalized)
      
#



x <- c(my_edger, my_deseq2)
y <- list(my_edger, my_deseq2)
z <- as.list(my_edger, my_deseq2)

