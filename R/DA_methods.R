
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
    norm = 'none', pval_method = "Cauchy", y_CorD, verbose = FALSE 
) {
    name <- "ZINQ"
    abundances <- microbiome::abundances(object)
    taxa <- microbiome::taxa(object)
    sample_metadata <- microbiome::meta(object)
    
    if (pseudo_count) {
        if (verbose)
            message(
                'A pseudocount of 1 was added to the abundance matrix.',
                ' for ZINQ test.'
            )
        abundances <- abundances + 1
    }
    
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
    
    condition_vector <- sample_metadata[[conditions_col]]
    denom <- conditions[['condB']]
    if (norm == 'CLR') {
        log2FoldChange <- log2_fold_change(
            abundances, condition_vector, denom, log = TRUE
        )
        
    } else {
        log2FoldChange <- log2_fold_change(abundances, condition_vector, denom)
    }
    
    sample_metadata[[conditions_col]] <- factor(
        sample_metadata[[conditions_col]], levels = conditions, labels = c(0,1)
    )
    
    abundances_t <- as.data.frame(t(abundances))
    
    rawP <- purrr::map_dbl(abundances_t, ~ {
        df <- data.frame(
            X = .x,
            Y = sample_metadata[[conditions_col]]
        )
        res <- tryCatch(
            error = function(e) NA, {
                ZINQ::ZINQ_tests(
                    formula.logistic =  stats::as.formula('X ~ Y'),
                    formula.quantile =  stats::as.formula('X ~ Y'),
                    C = "Y", y_CorD = y_CorD, data = df
                )
            })
        ZINQ::ZINQ_combination(res, method = pval_method)
    })
    
    # name <- paste0(name, '.', pval_method)
    
    adjP <- stats::p.adjust(rawP, method = 'fdr')
    
    pValMat <- data.frame(rawP = rawP, adjP = adjP)
    rownames(pValMat) <- taxa
    
    statInfo <- data.frame(
        log2FoldChange = log2FoldChange[names(rawP)],
        rawP = rawP,
        adjP = adjP
    )
    return(list(pValMat = pValMat, statInfo = statInfo, name = name))
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
    name <- 'LEfSe'
    
    # se <- mia::makeTreeSummarizedExperimentFromPhyloseq(object)
    se <- mia::convertFromPhyloseq(object)
    abundances <- SummarizedExperiment::assay(se)
    
    if (pseudo_count) {
        if (verbose)
            message('Adding a pseudocount of 1.')
        abundances <- abundances + 1
    }
    
    SummarizedExperiment::colData(se)[[groupCol]] <- 
        factor(
            SummarizedExperiment::colData(se)[[groupCol]], levels = conditions
        )
    
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
    
    SummarizedExperiment::assay(se) <- abundances
    
    statInfo <- lefser::lefser(se, classCol = groupCol, ...) 
    
    statInfo <- statInfo |> 
        dplyr::mutate(abs_score = abs(.data$scores)) |> 
        dplyr::arrange(abs_score)
    
    ## Add artificial p-values and adjusted p-values. This is for
    ## compatibility with the benchdamic workflow.
    ## I used these artificial values to order the lefse results by
    ## (adjusted) p-values, just like the other methods,
    ## That's why I'm using 'seq' below. 
    ## This is no longer necessary beacuse the results are now
    ## ordered according to LDA.
    ## However, I'm leaving the code here.
    
    statInfo$rawP <- seq(0.04, 0, length.out = nrow(statInfo))
    statInfo$adjP <- seq(0.09, 0, length.out = nrow(statInfo))
    rownames(statInfo) <- statInfo[["features"]]
    colnames(statInfo) <- c("Taxa", "LDA_scores", "abs_score", "rawP", "adjP")
    
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
DA_wilcox <- function(
        object, pseudo_count = FALSE, norm = 'none', 
        conditions_col, conditions, 
        verbose = FALSE
) {
    
    name <- "Wilcox"
    
    abundances <- microbiome::abundances(object)
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
        log2FoldChange <- log2_fold_change(
            abundances, condition_vector, denom, log = TRUE
        )
    } else {
        log2FoldChange <- log2_fold_change(
            abundances, condition_vector, denom
        )
    }
    
    abundances_t <- as.data.frame(t(abundances))
    pvalues <- purrr::map_dbl(abundances_t, ~ {
        df <- data.frame(
            condition = condition_vector,
            value = .x
        )
        wi_res <- stats::wilcox.test(value ~ condition, data = df)
        wi_res$p.value
    }) 
    
    adj_pvalues <- stats::p.adjust(pvalues, method = "fdr")
    names(adj_pvalues) <- names(pvalues)
    
    statInfo <- data.frame(
        log2FoldChange = log2FoldChange[names(pvalues)], 
        rawP = pvalues, 
        adjP = adj_pvalues
    )
    pValMat <- statInfo[, c("rawP", "adjP")]
    return(list(pValMat = pValMat, statInfo = statInfo, name = name))
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
    
    name <- 'ANCOM-BC'
    counts <- phyloseq::otu_table(object)
    
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
    
    if (any(counts == 0) && pseudo_count) {
        if (verbose)
            message('A pseudocount of 1 was added to the abundance matrix.')
        counts <- counts + 1
    }
    
    if (!norm %in% c('none', 'TSS'))
        stop('Normalization must be either `none` or `TSS`.')
    
    if (norm == 'none') {
        if (verbose)
            message('No normalization was applied.')
        name <- paste0(name, '.none')
        
    } else if (norm == 'TSS') {
        if (verbose)
            message('TSS normalization applied.')
        name <- paste0(name, '.TSS')
        counts <- norm_tss(counts)
    }
    
    phyloseq::otu_table(object) <- counts
    res <- ANCOMBC::ancombc(
        data = object, group = group, formula = formula, ...
    )[['res']]
    
    features_names <- res$p_val[[1]] # change when code of the acnbombc package changed
    
    pValMat <- data.frame(rawP = res$p_val[[3]], adjP = res$q_val[[3]])
    pValMat <- as.data.frame(pValMat)
    rownames(pValMat) <- features_names
   
    statInfo <- do.call('cbind', lapply(res, function(x) x[[3]]))
    statInfo <- as.data.frame(statInfo)
    rownames(statInfo) <- features_names
    
    return(list(pValMat = pValMat, statInfo = statInfo, name = name))
}

#' TSS normalization
#' 
#' \code{norm_TSS} Applies TSS normalization to a matrix of count data.
#' 
#' @param mat A numeric matrix of counts.
#' @param total_sum The tolal sum of the scaling, e.g. 100 or 1e.
#' 
#' @return A TSS-normalized matrix,
#' 
#' @export
#'
norm_tss <- function(mat, total_sum = 1e6) {
    apply(mat, 2, function(x) x / sum(x) * total_sum)
}

#' CLR normalization
#' 
#' \code{norm_clr} applies a centered-log ratio (CLR) transformation to a
#' matrix column-wise. Features (e.g. taxa, OTUs) must be in the rows and
#' samples in the columns.
#'
#' @param mat A count matrix.
#' @param pseudocount Pseudocount added to the matrix.
#' Default value is 0 (no pseudocount)..
#' @param log If TRUE, CLR will be logged. Default = TRUE. In most
#' cases this should be set to TRUE
#'
#' @return A CLR-transformed matrix.
#' @export
#'
norm_clr <- function(mat, pseudocount = 0, log = TRUE) {
    ## Centered log ratio transformation of a vector
    ## Sources: 
    ## + https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6755255/
    ## + https://www.r-bloggers.com/2021/08/calculate-geometric-mean-in-r/
    ## + https://doi.org/10.1038/s41467-022-28401-w
    if (any(is.na(mat) | any(mat < 0))) 
        stop("NA's or negative numbers are not allowed.", call. = FALSE)
    
    mat <- mat + pseudocount
    
    if (any(mat == 0)) {
        warning(
            "0s found in matrix. 1 was added as pseudocount.", call. = FALSE
        )
        mat <- mat + 1
    }
    
    if (log) {
        output <- apply(
            mat, 2, function(x) log(x / exp(mean(log(x))))
        ) # exp(mean(log(x))) is the geometric mean
    } else {
        output <- apply(mat, 2, function(x) x / exp(mean(log(x))))
    }
    return(output)
}

#' Calculate log2 fold change
#' 
#' \code{log2_fold_change} calculates the log2 fold change of a matrix with
#' features in the rows and samples in the columns.
#'
#' @param mat A matrix. Features in rows and samples in columns.
#' @param condition_vector A character vector or factor with the names of the 
#' conditions. The conditions must correspond to the samples, i.e. the exact
#' same order. Two and only two conditions (levels) are needed.
#' @param condB Condition used as reference. E.g. control condition.
#' @param log If log is TRUE, it's assumed that the matrix is already log
#' transformed.
#' @param pseudocount Numeric value indicating pseudocount to be added.
#'
#' @return
#' A named vector of log2 fold changes per feature.
log2_fold_change <- function(
        mat, condition_vector, condB = NULL, log = FALSE, pseudocount = 0
) {
    ## condB is control; condA is treated
    condition_vector <- as.factor(condition_vector)
    conditions <- levels(condition_vector)
    
    if (length(conditions) != 2)
        stop('Two and only two conditions are needed.', call. = FALSE)
    
    if (!is.null(condB)) {
        condA <- conditions[conditions != condB] # treated
        
    } else {
        condB <- conditions[1]
        condA <- conditions[2]
    }
    
    mat <- mat + pseudocount
    
    features <- rownames(mat)
    log2FoldChange <- vector("double", length(features))
    names(log2FoldChange) <- features 
    
    for (i in seq_along(features)) {
        
        mean_condB <- mean(mat[features[i], which(condition_vector == condB), drop = TRUE])
        mean_condA <- mean(mat[features[i], which(condition_vector == condA), drop = TRUE])
        
        if (log) { # CLR (already logged)
            log2FoldChange[i] <- mean_condA - mean_condB
            
        } else {
            if (mean_condA >= mean_condB) { # TSS - relative abundance
                log2FoldChange[i] <- log2(mean_condA / mean_condB)
                
            } else if (mean_condA < mean_condB) {
                log2FoldChange[i] <- -log2(mean_condB / mean_condA)
            }
        }
        
    }
    log2FoldChange
}