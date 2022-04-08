
#' Set list of normalization methods
#' 
#' \code{set_norm_list} is a wrapper of the 
#' \code{\link[benchdamic]{setNormalizations}} function with already predefined
#' normalization methods. The methods are: CSSdefault, CSSmedian, TMM,
#' poscounts, and none.
#' 
#' The output should be used as input of the
#' \code{\link[benchdamic]{runNormalizations}} function.
#' 
#' 
#' @return A list of normalization methods compatible with the benchdamic
#' framework.
#' 
#' @export
#' 
set_norm_list <- function() {
    
    normalization_methods <- tibble::tribble(
        ~fun, ~method,
        "norm_CSS", "default",
        "norm_CSS", "median",
        "norm_edgeR", "TMM",
        # "norm_edgeR", "TMMwsp",
        # "norm_edgeR", "RLE",
        # "norm_edgeR", "upperquartile",
        # "norm_edgeR", "posupperquartile",
        "norm_DESeq2", "poscounts", # RLE
        # "norm_DESeq2", "ratio",
        # "norm_DESeq2", "iterate",
        # "norm_TSS", "TSS",
        "norm_edgeR", "none"
    )
    
    benchdamic::setNormalizations(
        fun = normalization_methods$fun,
        method = normalization_methods$method
    )
}