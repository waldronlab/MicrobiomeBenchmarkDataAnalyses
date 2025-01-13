
#' Edit the taxa names by mia
#' 
#' \code{editMiaTaxaNames} corrects the taxonomic names of the 
#' `agglomerateByRank` function of the `mia` package. This is used in the
#' vignette of the Ravel dataset with bacterial vaginosis data.
#'
#' @param x A TreeSummarizedExperiment agglomerated by the `agglomerateByRank`
#' function of the mia package.
#'
#' @return A character vector with new names
#' @export
#'
editMiaTaxaNames <- function(x) {
    taxa_names <- rownames(x)
    taxa_ranks <- mia::taxonomyRanks(x)
    row_data <- as.data.frame(rowData(x))
    used_index <- 0
    counter <- length(taxa_ranks) + 1
    for (i in seq_along(taxa_ranks)) {
        x <- row_data[, taxa_ranks[counter - 1], drop = TRUE]
        index <- which(!is.na(x))
        index <- setdiff(index, used_index)
        taxa_names <-
            purrr::map_at(
                taxa_names, index, ~ paste0(taxa_ranks[counter-1], ':', .x)
            ) |> 
            as.character()
        used_index <- c(used_index, index)
        counter <- counter - 1
        
    }
    taxa_names |> 
        {\(y) gsub('_NA', '', y)}() |> 
        {\(y) sub('([a-z]:).*_([a-zA-Z]+)$', '\\1\\2', y)}()
}

#' Filter taxa based on abundance values per sample  
#' 
#' \code{filterTaxa} filters the number of taxa per sample based on a minimum
#' value of abundance. This functions works with both phyloseq and
#' TreeSummarizedExperiment objects.
#'
#' @param x A phyloseq or TreeSummarizedExperiment object with otu_table/assay
#' and sample_data/colData
#' @param min_ab The minimum value of abundance for taxon to be
#' considered as present in a sample. Default is 1. The default value of 1 
#' could be good for counts. Relative abundance or other data
#' transformations might require another threshold value.
#' @param min_per minimum percentage of samples in which each taxon must be
#' present in order to be kept in the data. Default is 0.2. Taxon presence is
#' dtermined by the `min_ab` argument (see above).
#'
#' @return The filtered phyloseq/TreeSummarized object
#' @export
#'
filterTaxa <- function(x, min_ab = 1, min_per = 0.2) { 
    
    if (is(x, 'TreeSummarizedExperiment') || is(x, 'SummarizedExperiment')) {
        m <- SummarizedExperiment::assay(x)
        n_samples <- ncol(m)
        min_n_samples <- round(n_samples * min_per)
        return(x[rowSums(m >= min_ab) >= min_n_samples,])
        
    } else if (is(x, 'phyloseq')) {
        m <- phyloseq::otu_table(x)
        n_samples <- ncol(m)
        min_n_samples <- round(n_samples * min_per)
        return(phyloseq::prune_taxa(rowSums(m >= min_ab) >= min_n_samples, x))
    }
}

#' Get direction columns
#' 
#' \code{get_direction_cols} gets the names of the columns with the
#' effect sizes for each DA method. The output is suitable for
#' enrichment analyses with benchdamic.
#'
#' @param x Output of \code{benchdamic::runDA}.
#' @inheritParams set_DA_methods_list
#'
#' @return A named vector ready to be used with the 
#' \code{\link[benchdamic]{createEnrichment}} and
#' \code{\link[benchdamic]{createPositives}} functions.
#' @export
#'
get_direction_cols <- function(x, conditions_col, conditions) {
    
    if (
        !length(names(conditions)) ||
        !any(names(conditions) == c('condB', 'condA'))
    ) {
        stop(
            'The `conditions` argument must be a named vector with names',
            ' "condB" and "condA, indicating reference/numerator and',
            ' target/denominator. For example:',
            '`c(condB = "control", condA = "condA"',
            call. = FALSE
        )
    }
    
    mgs <- paste0(conditions_col, conditions[['condA']])
    
    method_names <- names(x)
    index <-  match(method_names, methods_classification[['method']])
    effect_size_cols <- methods_classification[index,][['effect_size_col']]
    names(effect_size_cols) <- method_names
    
    for (i in seq_along(effect_size_cols)) {
        if (grepl('metagenomeSeq', names(effect_size_cols)[i])) {
            effect_size_cols[i] <- 'logFC'
        }
    }
    
    effect_size_cols
}

# Normalization -----------------------------------------------------------

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
        output <- apply(mat, 2, function(x) log(x / exp(mean(log(x))))) # exp(mean(log(x))) is the geometric mean
    } else {
        output <- apply(mat, 2, function(x) x / exp(mean(log(x))))
    }
    return(output)
}

# Calculations ------------------------------------------------------------

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


# Plotting ----------------------------------------------------------------

#' Get Method Classification
#' 
#' \code{getMethodClass} returns a tibble of the methods and their
#' classification. Helper function for plotting.
#'
#' @return A tibble.
#' @export
#'
getMethodClass <- function() {
    methods_classification |>
        dplyr::select(-.data$effect_size_col) |> 
        dplyr::relocate(.data$method_class, .data$base_method, .data$method) |> 
        dplyr::arrange(.data$method_class, .data$base_method, .data$method)
}
