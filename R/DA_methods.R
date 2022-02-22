
#' Set list of DA methods
#' 
#' \code{set_DA_methods_list} creates a list of methods for the benchdamic
#' workflow.
#'
#' @param grp Character string. Name of the column
#' @param contrast Character vector. Reference goes first.
#'
#' @return A list of methods for benchdamic
#' @export
#'
set_DA_methods_list <- function(grp, contrast) {
    
    # edgeR
    my_edger <- set_edgeR(
        group_name = grp,
        design = as.formula(paste0("~", grp)),
        norm = "TMM",
        coef = 2
    )
    
    # DESeq2
    my_deseq2 <- set_DESeq2(
        contrast = c(grp, contrast[2], contrast[1]),
        design = as.formula(paste0("~", grp)),
        norm = "poscounts"
    )
    
    # limma 
    my_limma <- set_limma( # I get a warning
        design = as.formula(paste0("~", grp)),
        norm = c("TMM", "CSSmedian"),
        coef = 2
    )
    
    # metagenomeSeq
    my_metagenomeseq <- set_metagenomeSeq(
        design = as.formula(paste0("~", grp)),
        norm = "CSSmedian",
        coef = 2
    )
    
    # ALDEx2
    my_aldex2 <- set_ALDEx2(
        conditions = grp,
        test = "t",
        norm = "none"
    )
    
    # corncob
    my_corncob <- set_corncob(
        formula = as.formula(paste0("~", grp)),
        phi.formula = as.formula(paste0("~", grp)),
        formula_null = ~ 1,
        phi.formula_null = as.formula(paste0("~", grp)),
        test = "Wald",
        coefficient = paste0(grp, contrast[2]),
        norm = "none"
    )
    
    # MAST
    my_mast <- set_MAST(
        rescale = "median",
        design = as.formula(paste0("~", grp)),
        coefficient = paste0(grp, contrast[2]),
        norm = "none"
    )
    
    # Seurat
    my_seurat <- set_Seurat(
        test.use = "wilcox",
        contrast = c(grp, contrast[2], contrast[1]),
        norm = "none"
    )
    
    my_new_methods <- list(
        zinq.1 = list(method = "zinq", grp = grp, ref = contrast[2]),
        ancombc.2 = list(method = "ancombc", formula = grp, group = grp),
        wilcox_test.3 = list(
            method = "wilcox_test", grp = grp, ref = contrast[2]
        ),
        wilcox_test_clr.4 = list(
            method = "wilcox_test_clr", grp = grp, ref = contrast[2]
        )
    )
    
    my_methods <- c(
        my_edger, my_deseq2, my_limma, my_metagenomeseq, my_aldex2,
        my_corncob, my_mast, my_seurat, my_new_methods
    )
    
    my_methods
}