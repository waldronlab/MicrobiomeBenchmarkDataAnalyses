
#' ZINQ for benchdamic
#' 
#' \code{zinq} adapts the zinq method for integration with the \code{benchdamic}
#' package.
#' 
#' @section References:
#' Ling, W., Zhao, N., Plantinga, A.M. et al. Powerful and robust non-parametric
#' association testing for microbiome data via a zero-inflated quantile approach
#' (ZINQ). Microbiome 9, 181 (2021). https://doi.org/10.1186/s40168-021-01129-3
#'
#' @param object A (Tree)SummarizedExperiment or a phyloseq object.
#' @param pseudo_count Whether include or not a pseudo_count. Default = FALSE.
#' @param conditions_col The name of the grouping column, which is located in
#' the colData (SummarizedExperiment) or the sample_data (phyloseq).
#' @param conditions A character vector indicating the conditions. Example:
#' `c(condB = 'control', condA = 'treatment')`
#' @param norm A character string indicating the normalization method to be used
#' prior analysis. Default is 'none'.
#' @param pval_method A character string indicating the type of pvalue to use.
#' Options: Cauchy or MinP. Default: Cauchy.
#' @param y_CorD Indicate if data are counts 'D' or continuous 'C'.
#' @param verbose Whether include messages or not. Default is FALSE.
#'
#' @return
#' A list in the format used in the benhcdamic pipeline.
#' 
#' @export
#'
DA_ZINQ <- function(
    object, pseudo_count = FALSE, conditions_col, conditions, 
    norm = 'none', pval_method, y_CorD, verbose = FALSE 
) {
    
    if (class(object) != 'phyloseq')
        stop(
            'The `object` argument must be a phyloseq object.',
            call. = FALSE
        )
    
    name <- "ZINQ"
    
    abundances <- microbiome::abundances(object)
    taxa <- microbiome::taxa(object)
    sample_metadata <- microbiome::meta(object)
    
    ## Pseudo count or not
    if (pseudo_count) {
        if (verbose)
            message(
                'A pseudocount of 1 was added to the abundance matrix.',
                ' for ZINQ test.'
            )
        abundances <- abundances + 1
    }
    
    ## Normalize data 
    
    if (
        !length(names(conditions)) ||
        any(names(conditions) != c('condB', 'condA'))
    ) {
        stop(
            'Condtions must be a named character vector with "condB"',
            ' and "condA" as names. For example:',
            ' `c(condB = "control", condA = "Treatment")`'
        )
    }
    
    if (length(norm) != 1 || !is.character(norm))
        stop(
            '`norm` must be a single character string. Valid options:',
            ' none, TSS, or CLR.',
            call. = FALSE
        )
    
    if (norm == 'none') {
        if (verbose)
            message('No normalization applied for ZINQ.')
        name <- paste0(name, '.none')
    } else if (norm == 'CLR') {
        if (verbose)
            message('Applying CLR normalization for ZINQ.')
        name <- paste0(name, '.CLR')
        abundances <- norm_clr(abundances)
    } else if (norm == 'TSS') {
        if (verbose)
            message('Applying TSS normalization for ZINQ.')
        name <- paste0(name, '.TSS')
        abundances <- norm_tss(abundances)
    }
    
    ## Calculate log2 fold change
    condition_vector <- sample_metadata[[conditions_col]]
    denom <- conditions[['condB']]
    
    if (norm == 'CLR') {
        log2FoldChange <- 
            log2_fold_change(abundances, condition_vector, denom, log = FALSE)
        
    } else {
        log2FoldChange <- log2_fold_change(abundances, condition_vector, denom)
    }
    
    
    ## Sanity check from ZINQ. If there are any warnings they should be
    ## Captured and reported if verbose
    ## TODO ??
    
    sample_metadata[[conditions_col]] <- factor(
        ## ZINQ::ZINQ_check requires a factor with 0 and 1.
        sample_metadata[[conditions_col]], levels = conditions, labels = c(0,1)
    )

    abundances_t <- as.data.frame(t(abundances))
    
    rawP <- vector("double", ncol(abundances_t))
    names(rawP) <- colnames(abundances_t)
    
    form <- stats::as.formula('X ~ Y')
    
    for (i in seq_along(rawP)) {
        
        df <- data.frame(
            X = abundances_t[[i]], 
            Y = sample_metadata[[conditions_col]]
        )
        
        res <- tryCatch(
            error = function(e) NULL, {
                ZINQ::ZINQ_tests(
                    formula.logistic = form,
                    formula.quantile = form,
                    C = "Y", y_CorD = y_CorD, data = df
                )
            })
        
        rawP[[i]] <- ZINQ::ZINQ_combination(res, method = pval_method) 
        
    }
    
    name <- paste0(name, '.', pval_method)
    
    adjP <- stats::p.adjust(rawP, method = 'fdr')
    
    pValMat <- data.frame(rawP = rawP, adjP = adjP)
    rownames(pValMat) <- taxa
    
    statInfo <- data.frame(
        log2FoldChange = log2FoldChange,
        rawP = rawP,
        adjP = adjP
    )
    
    list(pValMat = pValMat, statInfo = statInfo, name = name)
}

#' Lefser method
#' 
#' \code{DA_lefser} is a modified version of the lefser package, which includes
#' the pvalues of the Kruskal test.
#'
#' @param object A phyloseq or (Tree)SummarizedExperiment object.
#' @param pseudo_count Add or not a pseudocount of 1. Default is FALSE.
#' @param conditions A named character vector of length two. Names must be
#' "condB" and "condA" in that order. "condB" must indicate
#' control/reference/denominator and "condA" must indicate 
#' treatment/target/numerator. Example:
#' `c(condB = 'control', condA = 'treatment')`
#' @param norm Normalization method. Options: none, CLR, and TSS.
#' @param groupCol Name of the column in colData with the conditions.
#' @param ... Parameters passed to \code{\link[lefser]{lefser}}. 
#' @param verbose If TRUE, messages are displayed on screen. Default is FALSE.
#'
#' @return An object ready to be included in the benchdamic framework.
#' @export
#'
DA_lefse <- function(
    object, pseudo_count = FALSE, conditions, norm = 'none', verbose = FALSE,
    groupCol, ...
) {
    
    if (class(object) != 'phyloseq') {
        stop(
            'The object argument must be a phyloseq object',
            call. = FALSE
        )
    }
    
    name <- 'lefse'
    
    se <- mia::makeTreeSummarizedExperimentFromPhyloseq(object)
    abundances <- SummarizedExperiment::assay(se)
   
    ## Pseudocount 
    if (pseudo_count) {
        if (verbose)
            message('Adding a pseudocount of 1.')
        abundances <- abundances + 1
    }
    
    ## Conditions
    if (
        !length(names(conditions)) ||
        !all(names(conditions) == c('condB', 'condA'))
    ) {
        stop(
            'the `conditions` argument must be a named vector, and the names',
            ' must be "condB" and "condA". "condB" must be',
            ' the reference/control/denominator. Example:',
            ' `c(condB = "control", condA = "treatment")`',
            call. = FALSE
        )
    }
    
    SummarizedExperiment::colData(se)[[groupCol]] <- 
        factor(
            SummarizedExperiment::colData(se)[[groupCol]], levels = conditions
        )
    
    ## Normalization
    if (
        length(norm) != 1 ||
        !norm %in% c('none', 'CLR', 'TSS')
    ) {
        stop(
            'Only one normalization method should be included.',
            ' Supported options are none, CLR, and TSS.',
            call. = FALSE
        )
    }
    
    if (norm == 'none') {
        if (verbose)
            message('No normalization applied for lefse')
        name <- paste0(name, '.none')
    } else if (norm == 'CLR') {
        if (verbose)
            message('Applying CLR normalization.')
        name <- paste0(name, '.CLR')
        abundances <- norm_clr(abundances)
    } else if (norm == 'TSS') {
        if (verbose)
            message('Applying TSS normalization.')
        name <- paste0(name, '.TSS')
        abundances <- norm_tss(abundances)
    }
    
    ## Analysis with lefser(2)
    
    SummarizedExperiment::assay(se) <- abundances
    
    statInfo <- lefser2(expr = se, groupCol = groupCol, ...) 
    statInfo$adj_pval <- stats::p.adjust(statInfo$kw_pvalues, method = "fdr")
    rownames(statInfo) <- statInfo[["Names"]]
    colnames(statInfo) <- c("Taxa", "LDA_scores", "rawP", "adjP")
    
   
    pValMat <- statInfo[, c("rawP", "adjP")]
    rownames(pValMat) <- statInfo[["Taxa"]]
   
    list(pValMat = pValMat, statInfo = statInfo, name = name) 
}


#' Wilcox test for differential abundance
#' 
#' \code{DA_wilcox} performs Wilcoxon test on a phyloseq object.
#'
#' @param object A phyloseq object.
#' @param pseudo_count Whether include a pseudocount or not. Default is FASLE.
#' @param norm String character. Choose between three normalization methods: 
#' 'none', 'CLR', or 'TSS'.
#' @param conditions_col String character.
#' The name of the column in sample data with the conditions.
#' @inheritParams set_DA_methods_list
#' @param verbose Print messages or not, Default is FALSE.
#'
#' @return A list with outputs compatible with the benchdamic framework.
#' @export
#'
DA_wilcox <- 
    function(
        object, pseudo_count = FALSE, norm = 'none', 
        conditions_col, conditions, 
        verbose = FALSE
    ) {
        if(class(object) != 'phyloseq')
            stop(
                'Object must be phyloseq.',
                call. = FALSE
            )
        
        name <- "wilcox"
        
        abundances <- microbiome::abundances(object)
        taxa <- microbiome::taxa(object)
        sample_metadata <- microbiome::meta(object)
        
        sample_metadata[[conditions_col]] <- 
            factor(sample_metadata[[conditions_col]], levels = conditions)
        
        if (pseudo_count) {
            if (verbose)
                message(
                    'A pseudocount of 1 was added to the abundance matrix.',
                    ' for wilcox test.'
                )
            abundances <- abundances + 1
        }
        
        ## Normalize data 
        
        if (
            !length(names(conditions)) ||
            any(names(conditions) != c('condB', 'condA'))
        ) {
            stop(
                'Condtions must be a named character vector with "condB"',
                ' and "condA" as names. For example:',
                ' `c(condB = "control", condA = "Treatment")`'
            )
        }
        
        if (length(norm) != 1 || !is.character(norm))
            stop(
                '`norm` must be a single character string. Valid options:',
                ' none, TSS, or CLR.',
                call. = FALSE
            )
        
        if (norm == 'none') {
            if (verbose)
                message('No normalization applied for wilcox test.')
            name <- paste0(name, '.none')
        } else if (norm == 'CLR') {
            if (verbose)
                message('Applying CLR normalization for wilcox test.')
            name <- paste0(name, '.CLR')
            abundances <- norm_clr(abundances, log = TRUE)
        } else if (norm == 'TSS') {
            if (verbose)
                message('Applying TSS normalization for wilcox test.')
            name <- paste0(name, '.TSS')
            abundances <- norm_tss(abundances)
        }
        
        ## Calculate log2 fold change
        condition_vector <- sample_metadata[[conditions_col]]
        denom <- conditions[['condB']]
        
        if (norm == 'CLR') {
            log2FoldChange <- 
                log2_fold_change(abundances, condition_vector, denom, log = TRUE)
            
        } else {
            log2FoldChange <- log2_fold_change(abundances, condition_vector, denom)
        }
        
        ## Perform Wilcox test 
        pvalues <- vector("double", length(taxa))
        
        for (i in seq_along(pvalues)) {
            df <- data.frame(
                condition = condition_vector, value = abundances[i,]
            )
            wi_res <- stats::wilcox.test(value ~ condition, data = df)
            pvalues[i] <- wi_res$p.value
        }
        
        adj_pvalues <- stats::p.adjust(pvalues, method = "fdr")
       
        ## Combine result and return output 
        statInfo <- data.frame(
            log2FoldChange = log2FoldChange, 
            rawP = pvalues, 
            adjP = adj_pvalues
        )
        rownames(statInfo) <- taxa
        
        pValMat <- statInfo[, c("rawP", "adjP")]
        rownames(pValMat) <- taxa
        
        list(pValMat = pValMat, statInfo = statInfo, name = name)

}

#' DA_ancombc
#' 
#' \code{DA_ancombc} is an adaptation of the \code{\link[ANCOMBC]{ancombc}}
#' function for integration in the benchdamic framework.
#'
#' @param object A phyloseq object.
#' @param pseudo_count Logical. Whether inlucde a pseudocount of 1 or not.
#' Default is FALSE.
#' @param norm Character string indicating the normalization to use.
#' Options: 'none' and 'TSS'. Default is 'none'.
#' @param conditions A named character vector of length 2 indicating the names
#' of the conditions to be compared. The names must be 'condB' for the reference
#' and 'condA' for the target, in that order. For example:
#' `c(condB = 'control', condA = 'treatment')`.
#' @param verbose Logical. If TRUE messages at each step are printed on screen.
#' Default is FALSE.
#' @param group Name of column with conditions/group information.
#' Same argument as in \code{\link[ANCOMBC]{ancombc}}. Check the original
#' function for more information.
#' @param formula Name of column with conditions/group information.
#' Same argument as in \code{\link[ANCOMBC]{ancombc}}. Check the original
#' function for more information.
#' @param ... Parameters passed to the \code{\link[ANCOMBC]{ancombc}} function.
#'
#' @return A list with results of the analysis ready to be integrated in the
#' benchadmic framework.
#' @export
#' 
#' @seealso 
#' \code{\link[ANCOMBC]{ancombc}}
#'
DA_ancombc <- function(
    object, pseudo_count = FALSE, norm = 'none', conditions, 
    verbose = TRUE, group, formula, ...
) {
    
    name <- 'ancombc'
    
    if (!phyloseq::taxa_are_rows(object)) {
        object <- t(object)
    }
    
    counts <- phyloseq::otu_table(object)
    
    ## Check and set conditions for 'control' and 'treatment'
    if (
        !length(names(conditions)) || any(names(conditions) != c('condB', 'condA'))
    ) {
        stop(
            'Condtions must be a named character vector with "condB"',
            ' and "condA" as names. For example:',
            ' `c(condB = "control", condA = "Treatment")`'
        )
    }
    
    conditions_col <- as.factor(phyloseq::sample_data(object)[[group]])
    
    if (length(levels(conditions_col)) != 2)
        stop(
            'Two and only two conditions must be present in the conditions',
            ' column of sample (meta)data.'
        )
    
    if (!any(conditions %in% as.character(conditions_col)))
        stop(
            'Conditions are not present in sample metadata. Check the right',
            ' column and conditions levels.'
        )
    
    if (verbose)
        message(
            paste0(conditions[['condB']], ' will be used as reference.')
        )
    
    phyloseq::sample_data(object)[[group]] <- 
        factor(conditions_col, levels = conditions)
    
    # phyloseq::sample_data(obsample_metadata[[group]] <- 
        # factor(sample_metadata[[group]], levels = conditions)
    
    ## Pseudocount
    if (any(counts == 0) && pseudo_count) {
        if (verbose)
            message('A pseudocount of 1 was added to the abundance matrix.')
        counts <- counts + 1
    }
    
    ## Normalization
    if (!norm %in% c('none', 'TSS'))
        stop('Normalization must be either `none` or `TSS`.')
    
    if (norm == 'none') {
        if (verbose)
            message('No normalization was applied.')
        name <- paste0(name, '.none')
        ## Nothing is done on the otu_table
        
    } else if (norm == 'TSS') {
        if (verbose)
            message('TSS normalization applied.')
        name <- paste0(name, '.TSS')
        counts <- norm_tss(counts)
    }
    
    ## Perform analysis with ancombc
    phyloseq::otu_table(object) <- counts ## replace otu table 
    res <- ANCOMBC::ancombc(
        phyloseq = object, group = group, formula = formula, ...
    )[['res']]
    
    # Create pValMat and statInfo
    
    # features_names <- rownames(res$p_val) # names are the same for all outputs
    features_names <- res$p_val[[1]] # change when code of the acnbombc package changed
    
    ## I had to change the column in p_val, q_val, etc from 1 to 3
    pValMat <- data.frame(rawP = res$p_val[[3]], adjP = res$q_val[[3]])
    pValMat <- as.data.frame(pValMat)
    rownames(pValMat) <- features_names
   
    statInfo <- do.call('cbind', lapply(res, function(x) x[[3]]))
    statInfo <- as.data.frame(statInfo)
    rownames(statInfo) <- features_names
    
    list(pValMat = pValMat, statInfo = statInfo, name = name)
}

