
## BiocManager::install("waldronlab/MicrobiomeBenchmarkData")

suppressMessages({
    library(MicrobiomeBenchmarkData)
    library(benchdamic)
    library(mia)
    library(phyloseq)
})

## Get dataset from MicrobiomeBenchmarkData
tse <- getDataset("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)[[1]]

grp <- "hmp_body_subsite"
conditions <- c("subgingival_plaque", "supragingival_plaque")

## Convert to phyloseq
ps <- makePhyloseqFromTreeSummarizedExperiment(tse)

sample_data(ps)[[grp]] <- factor(
    sample_data(ps)[[grp]], levels = conditions
)

## Normalization
norm_par <- tibble::tribble(
    ~ fun, ~ method,
    "norm_edgeR", "none",
    "norm_TSS", "TSS",
    "norm_CSS", "median"
)

norms <- setNormalizations(fun = norm_par$fun, method = norm_par$method)

ps <- 
    runNormalizations(normalization_list = norms, object = ps, verbose = FALSE)

## Run ALDEx2

set.seed(123)
output_aldex2_t_tss <- DA_ALDEx2(
    object = ps, pseudo_count = FALSE, conditions = grp,
    test = "t", denom = "iqlr", norm = "TSS",
    verbose = FALSE
)

set.seed(123)
output_aldex2_t_none <- DA_ALDEx2(
    object = ps, pseudo_count = FALSE, conditions = grp,
    test = "t", denom = "iqlr", norm = "none",
    verbose = FALSE
)

set.seed(123)
output_aldex2_t_css <- DA_ALDEx2(
    object = ps, pseudo_count = FALSE, conditions = grp,
    test = "t", denom = "iqlr", norm = "CSSmedian",
    verbose = FALSE
)

## Comparing main outputs with and without TSS normalization

all(output_aldex2_t_tss$pValMat$rawP == output_aldex2_t_none$pValMat$rawP)
all(output_aldex2_t_tss$pValMat$adjP == output_aldex2_t_none$pValMat$adjP)
all(output_aldex2_t_tss$statInfo$effect == output_aldex2_t_none$statInfo$effect)

