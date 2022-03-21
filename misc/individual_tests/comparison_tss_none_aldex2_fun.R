
## BiocManager::install("waldronlab/MicrobiomeBenchmarkData")

suppressMessages({
    library(MicrobiomeBenchmarkData)
    library(benchdamic)
    library(mia)
    library(phyloseq)
})

## modified function
my_fun <- function(object, pseudo_count = FALSE, conditions = NULL,
                      mc.samples = 128, test = c("t","wilcox"), denom = "iqlr", norm = c("TMM",
                                                                                         "TMMwsp", "RLE", "upperquartile", "posupperquartile", "none", "ratio",
                                                                                         "poscounts", "iterate", "TSS", "CSSmedian", "CSSdefault"), verbose = TRUE){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "ALDEx2"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential",
             " abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop("Can't find the ", NF.col," column in your object.",
             " Make sure to add the normalization factors column in your",
             " object first.")}
    name <- paste(name, ".", norm, sep = "")
    NFs = unlist(metadata[, NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    # if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          # "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        # NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    if(is.null(conditions))
        stop("Please supply the name of the variable of interest or the",
             " entire character vector.")
    else if(length(conditions) == 1)
        conditions = unlist(metadata[, conditions])
    name <- paste(name, ".", denom, sep = "")
    if(!is.element(test, c("t","wilcox")) | length(test) != 1)
        stop("Please choose between p-values produced by Welch t-test (t) or",
             " by the Wilcoxon test (wilcox).")
    name <- paste(name, ".", test, sep = "")
    if(verbose){
        statInfo <- ALDEx2::aldex(reads = norm_counts, conditions = conditions,
                                  mc.samples = mc.samples, test = test, effect = TRUE,
                                  include.sample.summary = FALSE, denom = denom, verbose = verbose)
    } else {
        statInfo <- suppressMessages(ALDEx2::aldex(reads = norm_counts,
                                                   conditions = conditions, mc.samples = mc.samples, test = test,
                                                   effect = TRUE, include.sample.summary = FALSE, denom = denom,
                                                   verbose = verbose))
    }
    if(test == "t")
        pValMat <- data.frame(statInfo[, c("we.ep", "we.eBH")])
    else pValMat <- data.frame(statInfo[, c("wi.ep", "wi.eBH")])
    colnames(pValMat) <- c("rawP", "adjP")
    return(norm_counts)
    # return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_ALDEx2


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
    "norm_TSS", "TSS"
    
)

norms <- setNormalizations(fun = norm_par$fun, method = norm_par$method)

ps <- 
    runNormalizations(normalization_list = norms, object = ps, verbose = FALSE)


norm_counts_tss <- my_fun(
    object = ps, pseudo_count = FALSE, conditions = grp,
    test = "t", denom = "iqlr", norm = "TSS",
    verbose = FALSE
)


norm_counts_none <- my_fun(
    object = ps, pseudo_count = FALSE, conditions = grp,
    test = "t", denom = "iqlr", norm = "none",
    verbose = FALSE
)



