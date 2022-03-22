
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
    enrichment, enrichment_col, levels_to_plot = NULL, conditions
) {
    
    enrichment_col_var <- rlang::sym(enrichment_col)
    
    if (length(conditions) != 2 || !is.character(conditions))
        stop("The conditions argument should be a character vector of length 2.", call. = FALSE)
    
    if (!all(names(conditions) == c("condB", "condA")))
        stop(paste0(
            "The conditions vector must be named. condB first (control, reference)",
            " and condA second (treatment, target)."
        ), call. = FALSE)
    
    ## Create summary table 1 - Number of features
    summary_tbl_1 <- purrr::map(enrichment, ~ {
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
   

    ## Create summary table 2 - with p values of the fischer test
    summary_tbl_2 <- purrr::map(enrichment, ~ {
       
       my_grid <- expand.grid(
           c("DOWN Abundant", "UP Abundant"), levels_to_plot
       )
       
       colnames(my_grid) <- c("DA", enrichment_col)
       
       list_vct <- .x$tests
       
       if (!length(list_vct)) {
           my_grid[["pval"]] <- NA 
           return(my_grid)
       } else {
           df <- tibble::tibble(names(list_vct), unlist(list_vct)) %>% 
               magrittr::set_colnames(c("names", "pval")) %>% 
               tidyr::separate(col = "names", into = c("DA", enrichment_col), sep = '-')
           df <- dplyr::left_join(my_grid, df)
           return(df)
       }
   }) %>% 
       dplyr::bind_rows(.id = "method")

    ## Merge the two tables
    summary_tbl <- dplyr::left_join(
       summary_tbl_1, summary_tbl_2 
    )
    
    ## Change factors
    levels <- c("DOWN Abundant", "UP Abundant")
    names(levels) <- conditions
    summary_tbl[["DA"]] <- forcats::fct_recode(summary_tbl[["DA"]], !!!levels)
    
    ## Select annotations to display
    if (!is.null(levels_to_plot)) {
        summary_tbl <- summary_tbl[summary_tbl[[enrichment_col]] %in% levels_to_plot, ]
    }
    
    ## Add text for the pval
    
    condB <- conditions[["condB"]]
    condA <- conditions[["condA"]]
    
    summary_tbl <- summary_tbl %>% 
        dplyr::mutate(
            symb = dplyr::case_when(
                pval <= 0.01 & pval > 0.05 ~ "+",
                pval <= 0.05 & pval > 0.01 ~ "*",
                pval <= 0.01 & pval > 0.001 ~ "**",
                pval <= 0.001 ~ "***"
            ),
            ypos = dplyr::case_when(
                DA == condB ~ n - 4,
                DA == condA ~ n + 4
            )
    )
    
    ## Add classification of methods
    
    summary_tbl <- dplyr::left_join(summary_tbl, method_classification())
    
    ## Make plot
    
    max_DA <- max(abs(summary_tbl$n))
    
    summary_tbl %>% ggplot2::ggplot(
        mapping = ggplot2::aes(x = method, y = n)
    ) +
        ggplot2::geom_col(
            mapping = ggplot2::aes(fill = !!enrichment_col_var),
            position = ggplot2::position_dodge(0.9, preserve = "single")
        ) +
        ggplot2::geom_text(
            mapping = ggplot2::aes(x = method, y = ypos, label = symb, color = !!enrichment_col_var),
            position = ggplot2::position_dodge(0.9, preserve = "single")
        ) +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::scale_fill_brewer(type = 'qual', palette = "Set2", na.translate = FALSE ) +
        ggplot2::scale_color_brewer(type = 'qual', palette = "Set2") +
        ggplot2::scale_y_continuous(
            labels = abs,
            limits = c(-max_DA, max_DA),
            sec.axis = ggplot2::sec_axis(
                trans = ~ . / max(abs(summary_tbl$n)),
                breaks = c(-0.5, 0.5),
                labels = c(condB, condA)
            )
        ) +
        ggplot2::labs(
            y = "Number of features", x = "DA methods"
        ) +
        # ggplot2::geom_text(
            # data = annotation_table(condB, condA), 
            # mapping = ggplot2::aes(
                # x = x, y = y, hjust = hjust, vjust = vjust, label = text
            # )
        # ) +
        ggplot2::facet_wrap(~method_class, nrow = 1, scales = "free_x") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
            legend.position = "bottom",
            legend.title = ggplot2::element_blank(),
            axis.text.y.right = ggplot2::element_text(angle = 270, hjust = 0.5),
            axis.ticks.y.right = ggplot2::element_blank()
        ) 
}


# annotation_table <- function(condB, condA) {
#     tibble::tribble(
#         ~x, ~y, ~text, ~hjust, ~vjust,
#         Inf, -Inf, condB, 1.05, -1,
#         Inf, Inf, condA, 1.05, 1.5
#     )
# }

annotation_table <- function(condB, condA) {
    tibble::tribble(
        ~x, ~y, ~text, ~hjust, ~vjust,
        Inf, -Inf, condB, 1.05, -1,
        Inf, Inf, condA, 1.05, 1.5
    )
}


method_classification <- function() {
    tibble::tribble(
        ~ method, ~ method_class,
        "ALDEx2.none.iqlr.wilcox", "Compositional",
        "DESeq2.poscounts", "RNA-Seq",
        "DESeq2.poscounts.weighted", "scRNA-Seq",
        "edgeR.TMM", "RNA-Seq",
        "edgeR.TMM.weighted", "scRNA-Seq",
        "limma.TMM", "RNA-Seq",
        "limma.TMM.weighted", "scRNA-Seq",
        "wilcox.CLR", "Compositional",
        "wilcox.none", "Classical",
        "wilcox.TSS", "Classical"
        
    )
}




