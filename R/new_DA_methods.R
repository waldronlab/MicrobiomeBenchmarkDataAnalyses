
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
#' @param grp The name of the grouping column, which is located in the colData
#' (SummarizedExperiment) or the sample_data (phyloseq).
#' @param ref A character string indicating the group that should be used as
#' reference.
#' @param pval_type A character string indicating the type of pvalue to use.
#' Options: Cauchy or MinP. Default: Cauchy.
#' @param ... Other parameters
#'
#' @return
#' A list in the format used in the benhcdamic pipeline.
#' 
#' @export
#'
zinq <- function(object, grp, ref = NULL, pval_type = "Cauchy", ...) {
    tse <- NULL
    if (class(object) == "phyloseq") {
        m <- microbiome::abundances(object)
        taxa <- rownames(m)
        metadata <- microbiome::meta(object)
        
    } else if (class(object) %in% c("SummarizedExperiment", "TreeSummarizedExperiment")) {
        m <- SummarizedExperiment::assay(object)
        taxa <- rownames(m)
        metadata <- as.data.frame(SummarizedExperiment::colData(object))
    }
    
    list_of_abundances <- vector("list", length(taxa))
    names(list_of_abundances) <- taxa
    
    for (i in seq_along(taxa)) {
        abundance <- data.frame(m[taxa[i],])
        covariates <- metadata[rownames(abundance), grp]
        df <- cbind(abundance = abundance, covariates)
        colnames(df) <- c("abundance", "covariate")
        
        if (!is.null(ref)) {
            df[["covariate"]] <- 
                stats::relevel(factor(df[["covariate"]]), ref = ref)
        } else {
            df[["covariate"]] <- factor(df[["covariate"]])
        }
        
        ## Calculate fold change
        num_lvl <- levels(df[["covariate"]])[1]
        denom_lvl <- levels(df[["covariate"]])[2]
        
        num <- mean(df$abundance[df$covariate == num_lvl])
        denom <- mean(df$abundance[df$covariate == denom_lvl])
    
        
        if (num >= denom) {
            log2FoldChange <- log2(num / denom)
        } else if (num < denom) {
            log2FoldChange <- -log2(denom / num)
        }
        
        output <- tryCatch(
            error = function(e) NULL, {
                
                ZINQ::ZINQ_tests(
                    formula.logistic = abundance ~ covariate, 
                    formula.quantile = abundance ~ covariate, 
                    C = "covariate", y_CorD = "D", data = df
                )
            })
        
            pvalues <- vector("double", 3)
            names(pvalues) <- c("log2FoldChange", "Cauchy", "MinP")
            
        if ( !is.null(output)) {
            pvalues[["log2FoldChange"]] <- log2FoldChange
            pvalues[["Cauchy"]] <- ZINQ::ZINQ_combination(output, method="Cauchy")
            pvalues[["MinP"]] <- ZINQ::ZINQ_combination(output, method="MinP")
            
        } else {
            pvalues[["log2FoldChange"]] <- log2FoldChange
            pvalues[["Cauchy"]] <- NA
            pvalues[["MinP"]] <- NA
        }
        
            list_of_abundances[[i]] <- pvalues
    }
    
    output <- as.data.frame(do.call("rbind", list_of_abundances))
    output[["adj_Cauchy"]] <- stats::p.adjust(output[["Cauchy"]], method = "fdr")
    output[["adj_MinP"]] <- stats::p.adjust(output[["MinP"]], method = "fdr")
    # return(output)
    statInfo <- output
    pValMat <- output[,c("Cauchy", "adj_Cauchy")]
    colnames(pValMat) <- c("rawP", "adjP")
    name <- "ZINQ"
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}

#' Lefser method
#' 
#' \code{DA_lefser} is a modified version of the lefser package, which includes
#' the pvalues of the Kruskal test.
#'
#' @param object A phyloseq or (Tree)SummarizedExperiment object.
#' @param grp A character string indicating the name of the column in colData
#' where the conditions are stored.
#' @param ref A character string indicating which condition should be used as
#' reference.
#' @param norm Normalization method. Options: CLR and norm.
#' @param ... Other parameters passed to benchdamic. 
#'
#' @return An object for benchdamic.
#' @export
#'
DA_lefse <- function(object, grp, ref = NULL, norm = 'none', ...) {
    
    if (class(object) == "phyloseq") {
        se <- mia::makeTreeSummarizedExperimentFromPhyloseq(object)
    } else if (any(grepl("SummarizedExperiment", class(object)))) {
        se <- object
    }
    
    condition_vector <- SummarizedExperiment::colData(se)[[grp]]
    
    if (!is.null(ref)) {
        condition_vector <- stats::relevel(factor(condition_vector), ref)
    } else {
        condition_vector <- factor(condition_vector)
    }
    
    SummarizedExperiment::colData(se)[[grp]] <- condition_vector
    
    
    if (norm == 'CLR') {
        name <- 'lefse.CLR'
        SummarizedExperiment::assay(se) <- 
            norm_clr(SummarizedExperiment::assay(se))
    } else {
        name <- 'lefse.none'
    }
    
    statInfo <- lefser::lefser2(expr = se, groupCol = grp, ...) 
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
#' @param object A phyloseq object or (Tree)SummarizedExperiment
#' @param norm String character. Choose between three normalization methods: 
#' 'none', 'CLR', or 'TSS'.
#' @param design String character.
#' The name of the column in sample data with the conditions.
#' @param denom The sample that will be used as denominator to calcluate 
#' log2FoldChange. Usually, the 'treated' sample is used here.
#' @param verbose Print messages or not, Default is FALSE.
#'
#' @return A list with outputs compatible with the benchdamic framework/package.
#' @export
#'
DA_wilcox <- 
    function(
        object, norm = 'none', design, denom = NULL, pseudocount = NULL,
        verbose = FALSE
    ) {
        
        ## Extract components 
        if (class(object) == "phyloseq") {
            m <- microbiome::abundances(object)
            taxa <- rownames(m)
            metadata <- microbiome::meta(object)
            
        } else if (class(object) %in% c("SummarizedExperiment", "TreeSummarizedExperiment")) {
            m <- SummarizedExperiment::assay(object)
            taxa <- rownames(m)
            metadata <- as.data.frame(SummarizedExperiment::colData(object))
        }
        
        ## Normalize data 
        if (norm == 'none') {
            if (verbose)
                message('No normalization applied.')
            m <- m
        } else if (norm == 'CLR') {
            if (verbose)
                message('Applying CLR normalization')
            m <- norm_clr(m)
        } else if (norm == 'TSS') {
            if (verbose)
                message('Applying TSS normalization')
            m <- norm_tss(m)
        }
        
        ## Calculate log2 fold change
        condition_vector <- metadata[[design]]
        
        if (norm == 'CLR') {
            log2FoldChange <- 
                log2_fold_change(m, condition_vector, denom, log = TRUE)
            
        } else {
            log2FoldChange <- log2_fold_change(m, condition_vector, denom)
        }
        
        ## Perfrom Wilcox test 
        taxa <- rownames(m)
        pvalues <- vector("double", length(taxa))
        
        for (i in seq_along(pvalues)) {
            df <- data.frame(condition = condition_vector, value = m[i,])
            wi_res <- stats::wilcox.test(value ~ condition, data = df)
            pvalues[i] <- wi_res$p.value
        }
        
        adj_pvalues <- stats::p.adjust(pvalues, method = "fdr")
       
        ## Combine all outputs 
        statInfo <- data.frame(
            log2FoldChange = log2FoldChange, 
            rawP = pvalues, 
            adjP = adj_pvalues
        )
        
        pValMat <- statInfo[, c("rawP", "adjP")]
        
        name <- paste("wilcox", norm, sep = ".")
        
        return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))

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
        if (vervose)
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
    
    features_names <- rownames(res$p_val) # names are the same for all outputs
    
    pValMat <- data.frame(rawP = res$p_val[[1]], adjP = res$q_val[[1]])
    rownames(pValMat) <- features_names
   
    statInfo <- do.call('cbind', lapply(res, function(x) x[[1]]))
    rownames(statInfo) <- features_names
    
    list(pValMat = pValMat, statInfo = statInfo, name = name)
}

