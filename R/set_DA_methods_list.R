
#' Set list of DA methods
#' 
#' \code{set_DA_methods_list} creates a predefined list of methods for the
#' benchdamic framework. The output of this function should be the input
#' for the \code{\link[benchdamic]{runDA}} function.
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
        edgeR.TMM = benchdamic::set_edgeR(
            group_name = conditions_col,
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = FALSE
        ) |> 
            unname(),
        edgeR.TMM.w = benchdamic::set_edgeR(
            group_name = conditions_col,
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = TRUE
        ) |> 
            unname(),
        DESeq2.poscounts = benchdamic::set_DESeq2(
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "poscounts",
            weights_logical = FALSE
        ) |> 
            unname(),
        DESeq2.poscounts.w = benchdamic::set_DESeq2(
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "poscounts",
            weights_logical = TRUE
        ) |> 
            unname(),
        `Limma-Voom.TMM` = benchdamic::set_limma(
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "TMM",
            coef = 2,
            weights_logical = FALSE
        ) |> 
            unname(),
        `Limma-Voom.TMM.w` = benchdamic::set_limma(
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = "TMM", 
            coef = 2,
            weights_logical = TRUE
        ) |> 
            unname(),
        metagenomeSeq.CSS = benchdamic::set_metagenomeSeq(
            design = stats::as.formula(paste0("~", conditions_col)),
            norm = c("CSS"),
            coef = 2
        ) |> 
            unname(),
        `ALDEx2-Wilcox` = benchdamic::set_ALDEx2(
            contrast = c(conditions_col, conditions[[2]], conditions[[1]]),
            test = "wilcox",
            design = conditions_col
        ) |> 
            unname(),
        MAST = benchdamic::set_MAST(
            rescale = "median",
            design = stats::as.formula(paste0("~", conditions_col)),
            coefficient = paste0(conditions_col, conditions['condA'])
        ) |> 
            unname(),
        `Seurat-Wilcox` = benchdamic::set_Seurat(
            test = "wilcox",
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            # contrast = NULL,
            norm = "none"
        ) |> 
            unname(),
        `maaslin2` = benchdamic::set_Maaslin2(
            normalization = "TSS",
            transform = "LOG",
            analysis_method = "LM",
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            fixed_effects = conditions_col
        ) |> 
            unname(),
        `maaslin3` = benchdamic::set_maaslin3(
            normalization = "TSS",
            transform = "LOG",
            median_comparison_abundance =TRUE,
            subtract_median = TRUE,
            contrast = c(conditions_col, conditions['condA'], conditions['condB']),
            fixed_effects = conditions_col
        ) |> 
            unname(),
        `corncob` = benchdamic::set_corncob(
            formula = ~ conditions_col,
            phi.formula = ~ conditions_col,
            formula_null = ~1,
            phi.formula_null = ~ conditions_col,
            coefficient = conditions['condA'],
            test = "LRT",
            boot = TRUE
        ) |> 
            unname(),
        list(
            `ANCOM-BC` = list(
                method = 'DA_ancombc',
                conditions = conditions,
                group = conditions_col,
                formula = conditions_col,
                norm = 'none',
                p_adj_method = 'fdr'
            ),
            Wilcox.TSS = list(
                method = 'DA_wilcox',
                norm = 'TSS',
                conditions_col = conditions_col,
                conditions = conditions
            ),
            Wilcox.CLR = list(
                method = 'DA_wilcox',
                norm = 'CLR',
                conditions_col = conditions_col,
                conditions = conditions
            ),
            ZINQ.TSS = list(
                method = 'DA_ZINQ',
                conditions_col = conditions_col,
                conditions = conditions, 
                norm = 'TSS', pval_method = 'Cauchy', y_CorD = 'C'
            ),
            ZINQ.CLR = list(
                method = 'DA_ZINQ',
                conditions_col = conditions_col,
                conditions = conditions, 
                norm = 'CLR', pval_method = 'Cauchy', y_CorD = 'C'
            ),
            ## P-value thresholds are needed here since they are not provided
            ## in the output of the lefser package
            LEfSe.TSS = list(
                method = 'DA_lefse',
                conditions = conditions,
                norm = 'TSS',
                groupCol = conditions_col, ## classCol in lefser
                kruskal.threshold = 0.05,
                wilcox.threshold = 0.05,
                lda.threshold = 0
            ),
            LEfSe.CLR = list(
                method = 'DA_lefse',
                conditions = conditions,
                norm = 'CLR',
                groupCol = conditions_col, ## classCol in lefser
                kruskal.threshold = 0.05,
                wilcox.threshold = 0.05,
                lda.threshold = 0
            )
        )
    )
}
