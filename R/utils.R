
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
