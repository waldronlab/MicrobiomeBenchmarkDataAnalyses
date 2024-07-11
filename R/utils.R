
#' Edit the taxa names by mia
#' 
#' \code{editMiaTaxaNames} corrects the taxonomic names of the 
#' `agglomerateByRank` function of the `mia` package. For some reason,the names
#' are coming bad for some datasets. This is a workaround.
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

#' Quiet
#' 
#' \code{quiet} hides messages. This was taken from stackOverflow (reference
#' must be added).
#'
#' @param x Expression, command, etc.
#'
#' @return Nothing. NULL.
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
        # m <- phyloseq::otu_table(ps)
        # return(phyloseq::prune_samples(colSums(m > 1) >= 2, ps)) ## apply higher filtering
    }
}


#' Get direction columns
#' 
#' \code{get_direction_cols} gets the names of the columns with the
#' effect sizes for each DA method. The output is suitable for
#' enrichment analyses with benchdamic.
#'
#' @param x Output of \code{\link{run_DA}}.
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
        
        mean_condB <- mean(mat[features[i], condition_vector == condB])
        mean_condA <- mean(mat[features[i], condition_vector == condA])
        
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

#' Plot enrichment object
#' 
#' \code{plot_enrichment} make plot from benchdamic enrichment object.
#'
#' @param enrichment Enrichment output from benchdamic.
#' @param enrichment_col Column with enrichment annotations.
#' @param levels_to_plot Levels used for plotting. Default all.
#' @param conditions Named vector. The names must be condB and condA. Example:
#' c(condB = 'control', condA = 'treatment')
#' 
#' @return A ggplot object
#' @export
#' 
plot_enrichment <- function(
    enrichment, enrichment_col, levels_to_plot = NULL, conditions
) {
    method <- n <- ypos <- symb <- NULL
    
    enrichment_col_var <- rlang::sym(enrichment_col)
    
    if (length(conditions) != 2 || !is.character(conditions))
        stop("The conditions argument should be a character vector of length 2.", call. = FALSE)
    
    if (is.null(names(conditions)) || !all(names(conditions) == c("condB", "condA")))
        stop(paste0(
            "The conditions vector must be named. condB first,",
            "Example: c(condB = 'control', condA = 'treatment')"
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
   

    ## Create summary table 2 - p values of the fischer test
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
    
    ## This step is necessary to include all of the methods
    meth_class <- get_meth_class()
    summary_tbl$method <- factor(
        summary_tbl$method, levels = as.character(meth_class$method)
    )
    
    ## Add classification of methods
    summary_tbl <- dplyr::left_join(summary_tbl, meth_class)
    
    ## Make plot
    
    max_DA <- max(abs(summary_tbl$n))
    max_DA <- max_DA + 10
    
    summary_tbl[[enrichment_col]] <- 
        factor(
            summary_tbl[[enrichment_col]],
            levels = levels_to_plot
        )
    
    return(summary_tbl)
    
    summary_tbl %>%
        ggplot2::ggplot(
        mapping = ggplot2::aes(x = method, y = n)
    ) +
        ggplot2::geom_col(
            mapping = ggplot2::aes(fill = !!enrichment_col_var),
            position = ggplot2::position_dodge(preserve = "single")
        ) +
        ggplot2::geom_text(
            mapping = ggplot2::aes(
                x = method, y = ypos, label = symb, 
                color = !!enrichment_col_var, #group = !!enrichment_col_var
            ),
            position = ggplot2::position_dodge(width = 0.9, preserve = "single") 
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
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::labs(
            y = "Number of features", x = "DA methods"
        ) +
        # ggplot2::geom_text(
            # data = annotation_table(condB, condA), 
            # mapping = ggplot2::aes(
                # x = x, y = y, hjust = hjust, vjust = vjust, label = text
            # )
        # ) +
        ggplot2::facet_grid(
            . ~method_class, scales = "free_x", space = "free_x"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(hjust = 1, angle = 45),
            legend.position = "bottom",
            legend.title = ggplot2::element_blank(),
            axis.text.y.right = ggplot2::element_text(angle = 270, hjust = 0.5),
            axis.ticks.y.right = ggplot2::element_blank()
        ) 
}

#' Methods classification
#' 
#' \code{method_classification} returns a tibble of the methods and their
#' classification.
#'
#' @return A tibble.
#' @export
#'
get_meth_class <- function() {
    
    effect_size_col <- method_class <- base_method <- method <- NULL
    
    methods_classification %>% 
        dplyr::select(-effect_size_col) %>% 
        dplyr::relocate(method_class, base_method, method) %>% 
        dplyr::arrange(method_class, base_method, method)
    
    ## create mappings for color
    # base_methods <- sort(unique(df$base_method))
    # set.seed(12345)
    # colors <- randomcoloR::distinctColorPalette(length(base_methods))
    # base_methods_colors <- 
    #     data.frame(base_method = base_methods, color = colors)
    
    ## Create mappings for shape
   #  norm <- sort(unique(df$norm))
   #  shapes <- seq_along(norm)
   #  norm_shapes <- 
   #      data.frame(norm = norm, shape = shapes)
   #  
   # dplyr::left_join(df, norm_shapes, by = 'norm') 
    
}

#' Plot positives
#' 
#' \code{plot_positives} is a version of plotPositives
#'
#' @param x A dataframe. Positives with additional columns.
#'
#' @return List of plots
#' @export
#'
plot_positives <- function(x) {
    
    if (any(!c('linetype', 'shape', 'color') %in% colnames(x)))
        stop(
            'Dataframe must include color, shape and linetype. Be sure to use the plot_positives function.', 
            call. = FALSE
        )
    
    max_top <- max(x$TP - x$FP, na.rm = TRUE)
    min_top <- min(x$TP - x$FP, na.rm = TRUE)
    
    list_of_tables <- split(x, x$method_class)
    
    list_of_plots <- purrr::map(list_of_tables, ~ {
        
        df <- .x |> 
            purrr::modify_if(
                .p = is.character, .f = function(.y) forcats::fct_inorder(.y)
            )
        
        # df <- .x %>% 
        #     dplyr::mutate(
        #         dplyr::across(.cols = tidyselect:::where(is.character), .fns = ~ forcats::fct_inorder(.x)),
        #         linetype = forcats::fct_inorder(as.character(linetype))
        #     )
        
        p1 <- df %>%
            ggplot2::ggplot(ggplot2::aes(x = top, y = TP - FP)) +
            ggplot2::geom_path(
                ggplot2::aes(
                    color = method, group = method, linetype = method
                )
            ) +
            ggplot2::geom_point(
                ggplot2::aes(color = method, shape = method), 
                size = 3
            ) +
            ggplot2::facet_wrap(.~ method_class) +
            ggplot2::scale_y_continuous(
                limits = c(min_top - 3, max_top + 3)
            ) +
            ggplot2::scale_color_manual(
                name = 'method',
                values = as.character(df$color2), 
                labels = as.character(df$method)
                
            ) +
            ggplot2::scale_shape_manual(
                name = 'method',
                values = as.integer(df$shape),
                labels = as.character(df$method)
                # values = as.integer(levels(df$shape)), labels = levels(df$norm)
            ) +
            ggplot2::scale_linetype_manual(
                name = 'method',
                labels = unique(as.character(df$method)),
                values = as.integer(df$linetype)
            ) +
            ggplot2::labs(
                x = 'Top', y = 'TP - FP'
            ) +
            # ggplot2::guides(color=guide_legend(override.aes=list(fill=NA))) +
            # ggplot2::guides(
            #     color = ggplot2::guide_legend(order = 1)
            # ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                legend.position = c(0.01, 0.97),
                legend.justification = c("left", "top"),
                # legend.title = ggplot2::element_blank(),
                # legend.box.background = element_rect(
                #     size = 2, color = "red", fill = 'blue'
                # ),
                legend.box.background = ggplot2::element_blank(),
                legend.key = ggplot2::element_blank(),
                legend.background = ggplot2::element_blank(),
                legend.title = ggplot2::element_blank()
            )
        p1
    })
    list_of_plots
}

#' Plot enrichment
#' 
#' \code{plot_entichment_2}
#' 
#' @param data A data.frame
#' @param dir Direction
#'
#' @return A ggplot object.
#' @export
#'
plot_enrichment_2 <- function(data, dir) {
    data |> 
        dplyr::mutate(
            direction = ifelse(n > 0, dir[['up']], dir[['down']]),
            direction = factor(direction, levels = dir),
            taxon_annotation = factor(taxon_annotation),
            n = ifelse(n < 0, n * -1, n)
        ) |> 
        ggplot2::ggplot(mapping = ggplot2::aes(method, n)) +
        ggplot2::geom_col(
            mapping = ggplot2::aes(method, n, fill = taxon_annotation),
            position = ggplot2::position_dodge2(0.9, preserve = 'single')
        ) +
        ggplot2::geom_text(
            mapping = ggplot2::aes(label = symb, color = taxon_annotation),
            position = ggplot2::position_dodge2(0.9, preserve = 'single'),
            show.legend = FALSE
        ) +
        ggplot2::facet_grid(
            direction ~ method_class, space = 'free', scales = 'free_x'
        ) +
        ggplot2::labs(
            x = 'Method',
            y = 'Number of taxa'
        ) +
        ggplot2::scale_fill_brewer(name = 'Taxon annotaion', type = 'qual', palette = 'Set1') +
        ggplot2::scale_color_brewer(name = 'Taxon annotation', type = 'qual', palette = 'Set1') +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            axis.title = ggplot2::element_text(face = 'bold'),
            legend.position = 'bottom', 
            # legend.title = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(face = 'bold', size = 10),
            strip.background = ggplot2::element_blank()
            
        )
}