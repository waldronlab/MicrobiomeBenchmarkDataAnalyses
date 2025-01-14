
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
    
    mgs <- paste0(conditions_col, conditions[['condA']])
    methods_classification <- getMethodClass()
    
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

#' Get Method Classification
#' 
#' \code{getMethodClass} returns a tibble of the methods and their
#' classification. Helper function for plotting.
#'
#' @return A tibble.
#' @export
#'
getMethodClass <- function() {
    system.file(
        "extdata", "method_classification.tsv",
        package = "MicrobiomeBenchmarkDataAnalyses", 
        mustWork = TRUE
    ) |> 
        readr::read_tsv(show_col_types = FALSE) |> 
        # dplyr::select(-.data$effect_size_col) |> 
        dplyr::relocate(.data$method_class, .data$base_method, .data$method) |> 
        dplyr::arrange(.data$method_class, .data$base_method, .data$method)
}
