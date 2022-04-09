
#' Set list of DA methods
#' 
#' \code{set_DA_methods_list} creates a predefined list of methods for the
#' benchdamic workflow. This should be input to the
#' \code{\link[benchdamic]{runDA}} function.
#'
#' @param conditions_col Character string indicating the name of the column
#' containing the conditions of the samples in sample_metadata/colData. 
#' 
#' @param conditions A named character vector with names condB and condA. 
#' Reference is condB. For example, `c(condB = 'control', condA = 'treatment')`.
#'
#' @return A list of DA methods for benchdamic.
#' @export
#'
set_DA_methods_list <- function(conditions_col, conditions) {
    
    c(
        # edgeR
        my_edger <- benchdamic::set_edgeR(
            group_name = conditions_col,
            design = as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = FALSE
        ),
        
        # edgeR + weights
        
        my_edger <- benchdamic::set_edgeR(
            group_name = conditions_col,
            design = as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = TRUE
        ),
        
        
        # DESeq2
        my_deseq2 <- benchdamic::set_DESeq2(
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            design = as.formula(paste0("~", conditions_col)),
            norm = "poscounts",
            weights_logical = FALSE
        ),
        
        # DESeq2 + weights
        my_deseq2 <- benchdamic::set_DESeq2(
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            design = as.formula(paste0("~", conditions_col)),
            norm = "poscounts",
            weights_logical = TRUE
        ),
        
        # limma 
        my_limma <- benchdamic::set_limma( # I get a warning
            design = as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = FALSE
        ),
        
        # limma + weights 
        my_limma <- benchdamic::set_limma( # I get a warning
            design = as.formula(paste0("~", conditions_col)),
            norm = "TMM", 
            coef = 2,
            weights_logical = TRUE
        ),
        
        # metagenomeSeq
        my_metagenomeseq <- benchdamic::set_metagenomeSeq(
            design = as.formula(paste0("~", conditions_col)),
            norm = c("CSSmedian", "CSSdefault"),
            coef = 2
        ),
        
        # ALDEx2 with t-test
        my_aldex2 <- benchdamic::set_ALDEx2(
            conditions = conditions_col,
            test = "t",
            norm = "none"
        ),
        
        # ALDEx2 with wilcox
        my_aldex2 <- benchdamic::set_ALDEx2(
            conditions = conditions_col,
            test = "wilcox",
            norm = "none"
        ),
        
        # corncob
        my_corncob <- benchdamic::set_corncob(
            formula = as.formula(paste0("~", conditions_col)),
            phi.formula = as.formula(paste0("~", conditions_col)),
            formula_null = ~ 1,
            phi.formula_null = as.formula(paste0("~", conditions_col)),
            test = "Wald",
            coefficient = paste0(conditions_col, conditions['condA']),
            norm = "none"
        ),
        
        # MAST
        my_mast <- benchdamic::set_MAST(
            rescale = "median",
            design = as.formula(paste0("~", conditions_col)),
            coefficient = paste0(conditions_col, conditions['condA']),
            norm = "none"
        ),
        
        # Seurat
        my_seurat <- benchdamic::set_Seurat(
            test.use = "wilcox",
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            norm = "none"
        )
        
        # my_new_methods <- list(
        #     zinq.1 = list(method = "zinq", grp = grp, ref = contrast[2]),
        #     ancombc.2 = list(method = "ancombc", formula = grp, group = grp),
        #     wilcox_test.3 = list(
        #         method = "DA_wilcox_test", grp = grp, ref = contrast[2]
        #     ),
        #     wilcox_test_clr.4 = list(
        #         method = "DA_wilcox_test_clr", grp = grp, ref = contrast[2]
        #     ),
        #     kruskal_test.5 = list(
        #         method = "DA_kruskal_test", grp = grp, ref = contrast[2]
        #     ),
        #     kruskal_test_clr.6 = list(
        #         method = "DA_kruskal_test_clr", grp = grp, ref = contrast[2]
        #     ),
        #     lefse.7= list(
        #         method = "DA_lefse", grp = grp, ref = contrast[1]
        #     )
        #     
        # )
    )

}