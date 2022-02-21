
#' Normalization methods
#' 
#' @importFrom tibble tribble
#' 
#' @export
normalization_methods <- tibble::tribble(
    ~fun, ~method,
    # "norm_CSS", "default",
    "norm_CSS", "median",
    "norm_edgeR", "TMM",
    # "norm_edgeR", "TMMwsp",
    # "norm_edgeR", "RLE", # This is the same as poscounts?
    # "norm_edgeR", "upperquartile",
    # "norm_edgeR", "posupperquartile",
    "norm_DESeq2", "poscounts", # RLE
    # "norm_DESeq2", "ratio",
    # "norm_DESeq2", "iterate",
    "norm_TSS", "TSS",
    "norm_edgeR", "none"
)

#' set normalization
#' 
#' @importFrom benchdamic setNormalizations
#' 
#' @export
set_normalizations <- benchdamic::setNormalizations(
    fun = normalization_methods$fun, 
    method = normalization_methods$method
)



