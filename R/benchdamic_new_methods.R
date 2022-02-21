
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
#' @param x A (Tree)SummarizedExperiment or a phyloseq object.
#' @param grp The name of the grouping column, which is located in the colData
#' (SummarizedExperiment) or the sample_data (phyloseq).
#' @param ref A character string indicating the group that should be used as
#' reference.
#' @param pval_type A character string indicating the type of pvalue to use.
#' Options: Cauchy or MinP. Default: Cauchy.
#'
#' @return
#' A list in the format used in the benhcdamic pipeline.
#' 
#' @export
#'
zinq <- function(x, grp, ref = NULL, pval_type = "Cauchy") {
    tse <- NULL
    if (class(x) == "phyloseq") {
        m <- microbiome::abundances(x)
        taxa <- rownames(m)
        metadata <- microbiome::meta(x)
        
    } else if (class(x) %in% c("SummarizedExperiment", "TreeSummarizedExperiment")) {
        m <- SummarizedExperiment::assay(tse)
        taxa <- rownames(m)
        metadata <- SummarizedExperiment::colData(tse) |>
            as.data.frame()
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
    
    output <- do.call("rbind", list_of_abundances) |>
        as.data.frame()
    output[["adj_Cauchy"]] <- stats::p.adjust(output[["Cauchy"]], method = "fdr")
    output[["adj_MinP"]] <- stats::p.adjust(output[["MinP"]], method = "fdr")
    # return(output)
    statInfo <- output
    pValMat <- output[,c("Cauchy", "adj_Cauchy")]
    colnames(pValMat) <- c("rawP", "adjP")
    name <- "ZINQ"
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}

#' ANCOMBC
#' 
#' \code{ancombc} performs ancombc
#'
#' @param ps A phyloseq object
#' @param formula A string.
#' @param group A string.
#'
#' @return ANCOMBC result.
#' @export
#'
ancombc <- function(ps, formula, group) {
    out <- ANCOMBC::ancombc(phyloseq = ps, formula = formula, 
                            p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
                            group = group, struc_zero = TRUE, neg_lb = TRUE, 
                            tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, 
                            global = TRUE)
    res <- out$res
    ### extract important statistics ###
    vector_of_pval <- res$p_val[[1]] # contains the p-values
    vector_of_adjusted_pval <- res$q_val[[1]] # contains the adjusted p-values
    name_of_your_features <- rownames(res$p_val) # contains the OTU, or ASV, or other feature 
    # names. Usually extracted from the rownames of 
    # the count data
    vector_of_logFC <- res$beta[[1]] # logos the logFCs
    vector_of_statistics <- res$beta[[1]] # contains other statistics
    
    ### prepare the output ###
    pValMat <- data.frame("rawP" = vector_of_pval,
                          "adjP" = vector_of_adjusted_pval)
    statInfo <- data.frame("logFC" = vector_of_logFC,
                           "statistics" = vector_of_statistics) 
    name <- "ANCOMBC"
    # Be sure that your method hasn't changed the order of the features. If it 
    # happens, you'll need to re-establish the original order.
    rownames(pValMat) <- rownames(statInfo) <- name_of_your_features 
    
    # Return the output as a list
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}




