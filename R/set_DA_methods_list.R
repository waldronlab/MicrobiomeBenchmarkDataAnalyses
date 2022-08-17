
#' Set list of DA methods
#' 
#' \code{set_DA_methods_list} creates a predefined list of methods for the
#' benchdamic framework. The output of this function should input for the
#' \code{\link[benchdamic]{runDA}} function.
#'
#' @param conditions_col Character string indicating the name of the column
#' containing the conditions of the samples in sample_metadata/colData. 
#' @param conditions A named character vector. The names must be
#' "condB" and "condA". condB indicates the reference/numerator/control
#' and condA indicates the target/denominator/treatment. For example:
#' `c(condB = 'control', condA = 'treatment')`
#' @return A list of DA methods for benchdamic.
#' @export
#'
set_DA_methods_list <- function(conditions_col, conditions) {
    
    c(
        # edgeR
        my_edger <- benchdamic::set_edgeR(
            group_name = conditions_col,
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = FALSE
        ),
        
        # edgeR + weights
        
        my_edger <- benchdamic::set_edgeR(
            group_name = conditions_col,
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = TRUE
        ),
        
        
        # DESeq2
        my_deseq2 <- benchdamic::set_DESeq2(
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "poscounts",
            weights_logical = FALSE
        ),
        
        # DESeq2 + weights
        my_deseq2 <- benchdamic::set_DESeq2(
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "poscounts",
            weights_logical = TRUE
        ),
        
        # limma 
        my_limma <- benchdamic::set_limma( # I get a warning
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = FALSE
        ),
        
        # limma + weights 
        my_limma <- benchdamic::set_limma(
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "TMM", 
            coef = 2,
            weights_logical = TRUE
        ),
        
        ## metagenomeSeq
        my_metagenomeseq <- benchdamic::set_metagenomeSeq(
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = c("CSSmedian"),
            coef = 2
        ),

        ## ALDEx2 with t-test
        # my_aldex2 <- benchdamic::set_ALDEx2(
        #     conditions = conditions_col,
        #     test = "t",
        #     norm = "none"
        # ),

        ## ALDEx2 with wilcox
        my_aldex2 <- benchdamic::set_ALDEx2(
            conditions = conditions_col,
            test = "wilcox",
            norm = "none"
        ),

        ## corncob
        my_corncob <- benchdamic::set_corncob(
            formula = stats::as.formula(paste0("~", conditions_col)),
            phi.formula = stats::as.formula(paste0("~", conditions_col)),
            formula_null = ~ 1,
            phi.formula_null = stats::as.formula(paste0("~", conditions_col)),
            test = "Wald",
            coefficient = paste0(conditions_col, conditions['condA']),
            norm = "none"
        ),
        
        ## MAST
        my_mast <- benchdamic::set_MAST(
            rescale = "median",
            design = stats::as.formula(paste0("~", conditions_col)),
            coefficient = paste0(conditions_col, conditions['condA']),
            norm = "none"
        ),
        
        ## Seurat
        my_seurat <- benchdamic::set_Seurat(
            test.use = "wilcox",
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            norm = "none"
        ),
        
        new_methods <- list(
            ## ANCOMB-BC
            ancombc.1 = list(
               method = 'DA_ancombc',
               conditions = conditions,
               group = conditions_col,
               formula = conditions_col,
               norm = 'none',
               p_adj_method = 'fdr'
            ),
            ## Wilcox
            # wilcox.2 = list(
            #     method = 'DA_wilcox',
            #     norm = 'none',
            #     conditions_col = conditions_col,
            #     conditions = conditions
            # ),
            wilcox.3 = list(
                method = 'DA_wilcox',
                norm = 'TSS',
                conditions_col = conditions_col,
                conditions = conditions
            ),
            wilcox.4 = list(
                method = 'DA_wilcox',
                norm = 'CLR',
                conditions_col = conditions_col,
                conditions = conditions
            ),
            ## ZINQ
            # ZINQ.5 = list(
            #     method = 'DA_ZINQ',
            #     conditions_col = conditions_col,
            #     conditions = conditions, 
            #     norm = 'none', pval_method = 'MinP', y_CorD = 'D'
            # ),
            # ZINQ.6 = list(
            #     method = 'DA_ZINQ',
            #     conditions_col = conditions_col,
            #     conditions = conditions, 
            #     norm = 'TSS', pval_method = 'MinP', y_CorD = 'C'
            # ),
            # ZINQ.7 = list(
            #     method = 'DA_ZINQ',
            #     conditions_col = conditions_col,
            #     conditions = conditions, 
            #     norm = 'CLR', pval_method = 'MinP', y_CorD = 'C'
            # ),
            # ZINQ.8 = list(
            #     method = 'DA_ZINQ',
            #     conditions_col = conditions_col,
            #     conditions = conditions, 
            #     norm = 'none', pval_method = 'Cauchy', y_CorD = 'D'
            # ),
            ZINQ.9 = list(
                method = 'DA_ZINQ',
                conditions_col = conditions_col,
                conditions = conditions, 
                norm = 'TSS', pval_method = 'Cauchy', y_CorD = 'C'
            ),
            ZINQ.10 = list(
                method = 'DA_ZINQ',
                conditions_col = conditions_col,
                conditions = conditions, 
                norm = 'CLR', pval_method = 'Cauchy', y_CorD = 'C'
            ),
            # lefse.11 = list(
            #     method = 'DA_lefse',
            #     conditions = conditions,
            #     norm = 'none',
            #     groupCol = conditions_col,
            #     kruskal.threshold = 1,
            #     wilcox.threshold = 1,
            #     lda.threshold = 0
            # ),
            lefse.12 = list(
                method = 'DA_lefse',
                conditions = conditions,
                norm = 'CLR',
                groupCol = conditions_col,
                kruskal.threshold = 1,
                wilcox.threshold = 1,
                lda.threshold = 0
            ),
            lefse.13 = list(
                method = 'DA_lefse',
                conditions = conditions,
                norm = 'TSS',
                groupCol = conditions_col,
                kruskal.threshold = 1,
                wilcox.threshold = 1,
                lda.threshold = 0
            )
        )
        
    )

}