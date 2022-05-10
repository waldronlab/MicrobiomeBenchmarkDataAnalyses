
#' Run differential abundance methods
#' 
#' \code{run_DA} runs a set of predefined differential abundance methods on a
#' phyloseq object. This function includes normalization, calculation of matrix
#' of zinbweights, and running of the DA methods.
#'
#' @param object A phyloseq object.
#' @param conditions_col A character string indicating the name of the column in
#' sample_metadata that contains the conditions to be compared.
#' @inheritParams set_DA_methods_list
#' @param verbose This argument is passed to the 
#' \code{\link[benchdamic]{runNormalizations}} and 
#' \code{\link[benchdamic]{runDA}} functions.
#'
#' @return A list containing the results of the DA methods. The output is 
#' compatible with the benchdamic framework.
#' @export
#' 
run_DA <- function(object, conditions_col, conditions, verbose = FALSE) {
    
    if (
        !length(names(conditions)) ||
        !all(names(conditions) == c('condB', 'condA'))
    ) {
        stop(paste0(
            'The `conditions` argument must be a named vector of length two', 
            ' and the names must be condB and condA. Example:',
            ' `c(condB = "control", condA = "Treatment")`'
        ),
        call. = FALSE
        )
    }
    
    phyloseq::sample_data(ps)[[conditions_col]] <- 
        factor(
            phyloseq::sample_data(ps)[[conditions_col]],
            levels = conditions 
        )
    
    ps <- benchdamic::runNormalizations(
        normalization_list = set_norm_list(), object = object, verbose = verbose 
    )
    
    zinbWeights <- benchdamic::weights_ZINB(object = ps, design = conditions_col)
    
    DA_methods <- set_DA_methods_list(conditions_col, conditions)
    
    benchdamic::runDA(
        method_list = DA_methods, object = ps, weights = zinbWeights,
        verbose = verbose 
    ) 
    
}