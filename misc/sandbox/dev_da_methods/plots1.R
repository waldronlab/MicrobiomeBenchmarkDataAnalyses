

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
            y = "Number of features", x = "DA methods 2"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(hjust = 1, angle = 45)
        )
}


    