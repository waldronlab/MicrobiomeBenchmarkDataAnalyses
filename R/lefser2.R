
## This is an adaptation of a couple of functions from the lefser package.
## This was done in order to include raw p-values from the kruskal test.

filterKruskal2 <- function(expr, group, p.value) {
    
    kw.res.pvalues <- apply(expr, 1L, function(x) {
        kruskal.test(x ~ group)[["p.value"]]
    })
    
    kw.res.pvalues <- kw.res.pvalues[kw.res.pvalues <= p.value] 
    kw.res.pvalues <- kw.res.pvalues[!is.na(kw.res.pvalues)]
    kw.res.pvalues <- as.data.frame(kw.res.pvalues)
    colnames(kw.res.pvalues) <- "kw.pval"
    
    taxa <- rownames(kw.res.pvalues)
    filtered_matrix <- expr[taxa, ]
    
    list(pvalues = kw.res.pvalues, submatrix = filtered_matrix)
    
}

#' R implementation of the LEfSe method 2
#'
#' \code{lefser2} is an adaptation from the
#' \code{\link[lefser]{lefser}} function. This adaptation allows the inclusion
#' of the raw p-values of the KW test on the output.
#'
#' @param expr A \code{\linkS4class{SummarizedExperiment}} with expression data.
#' @param kruskal.threshold numeric(1) The p-value for the Kruskal-Wallis Rank
#' Sum Test (default 0.05).
#' @param wilcox.threshold numeric(1) The p-value for the Wilcoxon Rank-Sum Test
#' when 'blockCol' is present (default 0.05).
#' @param lda.threshold numeric(1) The effect size threshold (default 2.0).
#' @param groupCol character(1) Column name in `colData(expr)` indicating
#' groups, usually a factor with two levels (e.g., `c("cases", "controls")`;
#' default "GROUP").
#' @param blockCol character(1) Optional column name in `colData(expr)`
#' indicating the blocks, usually a factor with two levels (e.g.,
#' `c("adult", "senior")`; default NULL).
#' @param assay The i-th assay matrix in the `SummarizedExperiment` ('expr';
#' default 1).
#' @param trim.names If `TRUE` extracts the most specific taxonomic rank of organism.
#' @return
#' The function returns a dataframe with two columns, which are
#' names of microorganisms and their LDA scores.
#'
#' @importFrom stats kruskal.test reorder rnorm
#' @importFrom coin pvalue statistic wilcox_test
#' @importFrom MASS lda
#' @importFrom methods as is
#' @importFrom stats setNames
#' @importFrom utils tail
#' @import SummarizedExperiment
#'
#' @export
#' 
lefser2 <-
    function(expr,
             kruskal.threshold = 0.05,
             wilcox.threshold = 0.05,
             lda.threshold = 2.0,
             groupCol = "GROUP",
             blockCol = NULL,
             assay = 1L,
             trim.names = FALSE, log = FALSE
             )
    {
        groupf <- colData(expr)[[groupCol]]
        if (is.null(groupf))
            stop("A valid group assignment 'groupCol' must be provided")
        groupf <- as.factor(groupf)
        groupsf <- levels(groupf)
        if (length(groupsf) != 2L)
            stop(
                "Group classification is not dichotomous:\n",
                "Found (", paste(groupsf, collapse = ", "), ")"
            )
        group <- .numeric01(groupf)
        groups <- 0:1
        expr_data <- assay(expr, i = assay)
        # expr_sub <- filterKruskal(expr_data, group, kruskal.threshold)
        
        ## Adding these lines
        kw_res <- filterKruskal2(expr_data, group, kruskal.threshold)
        kw_pval <- data.frame(
            Names = rownames(kw_res[[1]]), kw_pvalues = kw_res[[1]][[1]]
        )
        kw_pval <- .trunc(kw_pval, trim.names)
        expr_sub <- kw_res[[2]]
        ## Lines above were added
        
        if (!is.null(blockCol)) {
            block <- as.factor(colData(expr)[[blockCol]])
            expr_sub <- fillPmatZmat(groupf, block, expr_sub, wilcox.threshold)
        }
        
        # transposes matrix and add a "class" (i.e., group) column
        # matrix converted to dataframe
        
        if (log) {
            expr_sub_t <- t(exp(expr_sub)) ## added this line
        } else {
            expr_sub_t <- t(expr_sub)
        }
        expr_sub_t_df <- as.data.frame(expr_sub_t)
        expr_sub_t_df <- createUniqueValues(expr_sub_t_df, groupf)
        expr_sub_t_df <- cbind(expr_sub_t_df, class = group)
        
        # number of samples (i.e., subjects) in the dataframe
        lfk <- nrow(expr_sub_t_df)
        # rfk is the number of subject that will be used in linear discriminant analysis
        rfk <- floor(lfk * 2 / 3)
        # number of classes (two)
        ncl <- length(groups)
        # count samples in each class of the dataframe, select the number from the class with a smaller
        # count of samples and multiply that number by 2/*2/3*0.5
        min_cl <-
            as.integer(min(table(expr_sub_t_df$class)) * 2 / 3 * 2 / 3 *
                           0.5)
        # if min_cl is less than 1, then make it equal to 1
        min_cl <- max(min_cl, 1)
        
        # lda_fn repeated 30 times, producing a matrix of 30 scores per feature
        eff_size_mat <-
            replicate(30, suppressWarnings(ldaFunction(
                expr_sub_t_df, lfk, rfk, min_cl, ncl, groups
            )), simplify = TRUE)
        
        # mean of 30 scores per feature
        raw_lda_scores <- rowMeans(eff_size_mat)
        
        # processing of score
        processed_scores <-
            sign(raw_lda_scores) * log((1 + abs(raw_lda_scores)), 10)
        
        # sorting of scores
        processed_sorted_scores <- sort(processed_scores)
        scores_df <- data.frame(Names = names(processed_sorted_scores),
                                scores = as.vector(processed_sorted_scores),
                                stringsAsFactors = FALSE)
        
        scores_df <- .trunc(scores_df, trim.names)
        
        threshold_scores <- abs(scores_df$scores) >= lda.threshold
        
        ## Adding a these lines
        final_scores <- scores_df[threshold_scores, ]
        which_names <- which(kw_pval[["Names"]] %in% final_scores[["Names"]])
        kw_pval_subset <- kw_pval[which_names, ]
        merge(scores_df, kw_pval_subset, by = "Names")
    }