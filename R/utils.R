
#' Quiet
#' 
#' \code{quiet} hides messages. This was taken from stackOverflow (reference
#' must be added).
#'
#' @param x Expression, command, etc.
#'
#' @return Nothing
#' 
#' @export
#'
quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
} 

#' Taxize classification to a taxonomy table
#'
#' \code{taxize_classification_to_taxonomy_table} converts the output of
#' \code{taxize::classification} to a data frame.
#' 
#' This fuction was taken from the misctoolsr package on GitHub. Repo:
#' sdgamboa/misctoolsr
#'
#' @param tax Output of `taxize::classification`.
#' @param id_type Either "name" or "id".
#'
#' @return A data frame
#' 
#' @importFrom dplyr bind_rows
#' 
#' @export
#'
taxize_classification_to_taxonomy_table <- function(tax, id_type = "name") {
    valid_ranks <- c(
        "superkingdom", "class", "order", "family", "genus", "species"
    )
    query_names <- as.data.frame(names(tax))
    colnames(query_names) <- "query"
    taxonomy_list <- lapply(tax, function(x) {
        if (any(is.na(x))) {
            data.frame(
                kingdom = NA, class = NA, order = NA, family = NA,
                genus = NA, species = NA
            )
        } else {
            df <- x[x$rank %in% valid_ranks, ]
            df <- df[,c("rank", id_type)]
            df <- as.data.frame(t(df))
            col_names <- as.character(df[1,])
            df <- as.data.frame(df[-1,])
            colnames(df) <- col_names
            colnames(df)[colnames(df) == "superkingdom"] <- "kingdom"
            rownames(df) <- NULL
            df
        }
    })
    taxonomy_table <- taxonomy_list %>%
        dplyr::bind_rows()
    cbind(query_names, taxonomy_table)
}

#' Fitler phyloseq object
#' 
#' \code{filter_phyloseq} filters both taxa and samples.
#'
#' @param ps A phyloseq object with otu_table and sample_data
#'
#' @return The filtered phyloseq object
#' @export
#'
filter_phyloseq <- function(ps) {
    m <- phyloseq::otu_table(ps)
    ps <- phyloseq::prune_taxa(rowSums(m > 0) >= 5, ps)
    m <- phyloseq::otu_table(ps)
    ps <- phyloseq::prune_samples(colSums(m > 0) >= 2, ps)
    ps
}

#' Convert from relative abundance (100) to counts
#' 
#' \code{relative_abundace_to_counts} converts from relative abundance to
#' counts.
#'
#' @param x A matrix. Taxa in rownames. Samples in colnames.
#' @param total_reads Numeric vector with total reads per sample. It must be in
#' the same order as in the colnames of the matrix.
#' @param total_sum Total scaling. Default 100 (percentage).
#'
#' @return A matrix of counts
#' @export
#'
relative_abundance_to_counts <- function(x, total_reads, total_sum = 100) {
    t(apply(x, 1, function(x) round(x * total_reads / total_sum)))
}


#' Calculate log Fold Change
#' 
#' \code{log_fold_change} calculates log fold change
#'
#' @param num vector with numerator
#' @param denom vector with denominator
#'
#' @return data transformed to log2 fold change
#' @export
#'
log_fold_change <- function(num, denom) {
    if (num >= denom) {
        return(log2(num / denom))
    } else if (num < denom) {
        return(-log2(denom / num))
    }
}