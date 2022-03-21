library(MicrobiomeBenchmarkData)
library(mia)
library(phyloseq)

tse <- getDataset('HMP_2012_16S_gingival_V35_subset', dryrun = FALSE)[[1]]
ps <- makePhyloseqFromTreeSummarizedExperiment(tse)

grp <- "hmp_body_subsite"
conditions <- c(condA = 'subgingival_plaque', condB = 'supragingival_plaque')


output_wilcox_none <- DA_wilcox(
    object = ps, norm = 'none', design = grp, denom = conditions[2]
)

output_wilcox_clr <- DA_wilcox(
    object = ps, norm = 'CLR', design = grp, denom = conditions[2],
    verbose = TRUE
)

output_wilcox_tss <- DA_wilcox(
    object = ps, norm = 'TSS', design = grp, denom = conditions[2]
)
