
library(MicrobiomeBenchmarkData)
library(mia)
library(benchdamic)

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

# limma benchdamic ---------------------------------------------------------

ps <- makePhyloseqFromTreeSummarizedExperiment(tse)

## Add normalization factors
ps <- norm_edgeR(object = ps, method = "TMM", verbose = TRUE)

## Calculate zinbwave weights
zinbweights <- weights_ZINB(
    object = ps,
    K = 0,
    design = "~ 1"
)

## Run without zinbwave weights
output_limma <- DA_limma(
    object = ps, pseudo_count = FALSE, design = grp, norm = "TMM",
    verbose = TRUE,
)

## Run with zinwave weights
output_limma_zinbweights <- DA_limma(
   object = ps, pseudo_count = FALSE, design = grp, norm = "TMM" ,
   weights = zinbweights, verbose = TRUE
)
