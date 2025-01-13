
## Notes ------------------------------------------------------------------

## Note 1: Here, only the normalization methods implemented in benchdamic are
## used.

## Note 2: The TSS method will be implemented differently than in benchdamic,
## so it's not included in this script.

## Note 3: The CLR method will be added and implemented
## differently  than in benchdamic, so it's not included in this script.

# Code --------------------------------------------------------------------

#' Set list of normalization methods
#' 
#' \code{set_norm_list} is a wrapper of the 
#' \code{\link[benchdamic]{setNormalizations}} function.
#' 
#' The output should be used as input of the
#' \code{\link[benchdamic]{runNormalizations}} function.
#' 
#' @return A list of normalization methods compatible with the benchdamic
#' framework.
#' 
#' @export
#' 
set_norm_list <- function() {
    normalization_methods <- tibble::tribble(
        ~fun, ~method,
        # "norm_CSS", "default",
        # "norm_CSS", "median",
        "norm_CSS", "CSS",
        "norm_edgeR", "TMM",
        # "norm_edgeR", "TMMwsp",
        # "norm_edgeR", "RLE", # same as DESEq2's poscounts
        # "norm_edgeR", "upperquartile",
        # "norm_edgeR", "posupperquartile",
        "norm_DESeq2", "poscounts", # same as edgeR's RLE
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