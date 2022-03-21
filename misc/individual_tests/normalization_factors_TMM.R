
## BiocManager::install("waldronlab/MicrobiomeBenchmarkData")
suppressMessages({
    library(MicrobiomeBenchmarkData)
    library(benchdamic)
    library(mia)
    library(phyloseq)
})

tse <- getDataset("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)[[1]]
ps <- makePhyloseqFromTreeSummarizedExperiment(tse)

## Add normalization factor to sample metadata (NF.TSS column)
ps <- norm_edgeR(ps, "TMM", verbose = FALSE)

counts <-as(otu_table(ps), "matrix")
NF <- sample_data(ps)[["NF.TMM"]]

## Perform matrix normalization
## What I think:
# norm_counts <- t(apply(t(counts), 2, function(x) round(x * NF * 1e6 )))
# summary(colSums(norm_counts))

## What I find in the functions:
NF2 <- NF * colSums(counts)
NF2 <- NF2/exp(mean(log(NF2)))
norm_counts_2 <- round(counts %*% diag(1/NF2), digits = 0)
colnames(norm_counts_2) <- colnames(counts)
summary(colSums(norm_counts_2))

## Which are the same values as no normalization
summary(colSums(counts))
all(counts == norm_counts_2)


