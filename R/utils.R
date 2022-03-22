
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

#' Apply CLR to matrix
#' 
#' \code{apply_clr} applies a centered-log ratio transformation to a matrix.
#' Features (e.g. taxa, OTUs) must be in the rows and samples in the columns.
#'
#' @param x A count matrix.
#' @param pseudocount A pseudocount to add. Default = 1.
#'
#' @return A matrix with CLR normalization
#'
norm_CLR <- function(x, pseudocount = 0) {
    ## Centered log ratio transformation of a vector
    ## Sources: 
    ## + https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6755255/
    ## + https://www.r-bloggers.com/2021/08/calculate-geometric-mean-in-r/
    
    mat <- x + pseudocount
    
    if (any(is.na(mat) | any(mat < 0))) {
        stop("Input vector must not contain NAs, or negative numbers.")
    } else if (any(mat == 0)) {
        warning(
            "0s found in matrix. 1 was added as pseudocount.", call. = FALSE
        )
        mat <- mat + 1
    }
        
    apply(mat, 2, function(x) log(x / exp(mean(log(x))))) # exp(mean(log(x))) is the geometric mean
}


#' Calculate log2 fold change
#' 
#' \code{log2_fold_change} calculates the log2 fold change of the features
#' of a matrix per sample. Denom is usually the treated condition.
#'
#' @param mat A matrix. Features in rows and samples in columns.
#' @param condition_vector A vector or factor with the names of the conditions.
#' Only two conditions are allowed and must be in the same order as the names
#' in the samples.
#' @param denom Condition used as reference. E.g. control condition.
#' @param log If log is TRUE, it's assumed that the matrix is already log
#' transoformed.
#'
#' @return
#' A named vector of log2 fold changes per feature.
#'
log2_fold_change <- function(mat, condition_vector, denom = NULL, log = FALSE) {
    
    condition_vector <- as.factor(condition_vector)
    conditions <- levels(condition_vector)
    
    if (length(conditions) != 2)
        stop('Two and only two levels are needed.', call. = FALSE)
    
    if (!is.null(denom)) {
        num_lvl <- conditions[conditions != denom] # treated
        denom_lvl <- conditions[conditions == denom] # control
        
    } else {
        num_lvl <- levels(condition_vector)[2] # treated
        denom_lvl <- levels(condition_vector)[1] # control
    }
    
    taxa <- rownames(mat)
    log2FoldChange <- vector("double", length(taxa))
    names(log2FoldChange) <- taxa
    
    for (i in seq_along(taxa)) {
        
        df <- cbind(condition_vector, as.data.frame(mat[i,]))
        colnames(df) <- c("condition", "value")
        
        num <- mean(df$value[df$condition == num_lvl])
        denom <- mean(df$value[df$condition == denom_lvl])
        
        if (log) {
            log2FoldChange[i] <- num - denom  # treated - control
        } else {
            if (num >= denom) {
                log2FoldChange[i] <- log2(num / denom) # treated / control
                
            } else if (num < denom) {
                log2FoldChange[i] <- -log2(denom / num) # treated / control
            }
        }
            
    }
    
    log2FoldChange
    
}


#' TSS normalization of matrix
#' 
#' \code{norm_TSS} Applies TSS normalization to a matrix of count data.
#' @param mat A numeric matrix.
#'
#' @return A TSS-normalized matrix
#'
norm_TSS <- function(mat) {
    apply(mat, 2, function(x) x / sum(x) * 1e6)
}


#' Plot enrichment object
#' 
#' \code{plot_enrichment} make plot from enrichment object.
#'
#' @param enrichment Enrichment object from benchdamic.
#' @param enrichment_col Column with enrichmnet annotations.
#' @param levels_to_plot Labels used for plotting. Default all.
#' @param conditions Named vector. condB first (i.e., reference, control).
#'
#' @return A ggplot object
#' @export
#' 
plot_enrichment <- function(
    enrichment, enrichment_col, levels_to_plot, conditions
) {
    
    enrichment_col_var <- rlang::sym(enrichment_col)
    
    if (length(conditions) != 2 || !is.character(conditions))
        stop("The conditions argument should be a character vector of length 2.", call. = FALSE)
    
    ## Create summary table
    summary_tbl <- purrr::map(enrichment, ~ {
        df <- .x$data
        df <- df %>%
            dplyr::filter(DA != "non-DA") %>% 
            dplyr::count(DA, !! enrichment_col_var) %>% 
            dplyr::filter(dplyr::if_all(.cols = dplyr::everything(), .fns = ~ !is.na(.x))) %>% 
            dplyr::mutate(
                n = ifelse(grepl("DOWN", DA), -n, n)
            )
        if (nrow(df) == 0) { # Thi is an empty dataframe
            df[1,] <- c("UP Abundant", NA, 0)
            df[2,] <- c("DOWN Abundant", NA, 0)
            df[[3]] <- as.numeric(df[[3]])
        }
        df
    }) %>%  
        dplyr::bind_rows(.id = "method")
    
    ## Change factors
    levels <- c("DOWN Abundant", "UP Abundant")
    names(levels) <- conditions
    summary_tbl[["DA"]] <- forcats::fct_recode(summary_tbl[["DA"]], !!!levels)
    
    ## Make plot
    summary_tbl %>% ggplot2::ggplot(
        mapping = ggplot2::aes(x = method, y = n)
    ) +
        ggplot2::geom_col(
            mapping = ggplot2::aes(fill = type),
            position = ggplot2::position_dodge(0.9, preserve = "single")
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::scale_fill_discrete(na.translate = F) +
        ggplot2::scale_y_continuous(labels = abs) +
        ggplot2::labs(
            y = "Number of features", x = "DA methods"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(hjust = 1, angle = 45)
        )
}



