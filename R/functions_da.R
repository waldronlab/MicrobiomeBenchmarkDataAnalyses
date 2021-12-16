## This file contains functions for differential abundance (DA) analyses

deseq2Poscounts <- function(se, grp, ref = NULL) {
    se <- groupsAsFactors(se, grp = grp, ref = ref)
    design <- stats::as.formula("~ GROUP")
    suppressMessages({
        dds <- DESeq2::DESeqDataSet(se = se, design = design)
        dds <- DESeq2::DESeq(
            dds, sfType = "poscounts", test = "LRT", reduced = ~ 1,
            quiet = TRUE
        )
        dds_results <- DESeq2::results(dds, pAdjustMethod = "BH", alpha = 0.1)
    })
    dds_results_tbl <- dds_results %>% 
        as.data.frame() %>% 
        tibble::as_tibble(rownames = "rownames")
    select_cols <- c("rownames", "log2FoldChange", "pvalue", "padj")
    dds_results_tbl <- dds_results_tbl[, select_cols ]
    dds_results_tbl %>%
        magrittr::set_colnames(c("TAXA", "FC", "PVAL", "ADJPVAL"))
}

edgerTmm <- function(se, grp, ref = NULL) {
    
    se <- groupsAsFactors(se, grp = grp, ref = ref)
    dge <- dgeTmm(se = se) # Normalization step, adds norm.factors (TMM) to sample data
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    design <- stats::model.matrix(stats::as.formula("~ GROUP"), col_data)
    
    dge <- edgeR::estimateDisp(dge, design)
    glmFit <- edgeR::glmFit(dge, design)
    glmRes <- edgeR::glmLRT(glmFit, coef = 2)
    
    taxnames <- rownames(glmRes$table)
    pval <- glmRes$table$PValue
    padj <- stats::p.adjust(pval, "BH")
    fc <- glmRes$table$logFC
    
    tibble::tibble(
        TAXA = taxnames,
        FC = fc,
        PVAL = pval,
        ADJPVAL = padj
    )
}

limmaVoomTmm <- function(se, grp, ref = NULL) {
    
    se <- groupsAsFactors(se, grp = grp, ref = ref)
    dge <- dgeTmm(se = se) # Normalization step, adds norm.factors (TMM) to sample data
    sample_data <- as.data.frame(SummarizedExperiment::colData(se))
    design <- stats::model.matrix(stats::as.formula("~ GROUP"), sample_data)
    NFs <- dge$samples$norm.factors
    
    v <- limma::voom(dge, design, plot = FALSE, lib.size = colSums(SummarizedExperiment::assay(se)) * NFs)
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = 2, n = nrow(dge), sort.by = "none")
    
    taxnames <- rownames(tt)
    fc <- tt$logFC
    pval <- tt$P.Value
    padj <- stats::p.adjust(pval, method = "BH")
    
    tibble::tibble(
        TAXA = taxnames,
        FC = fc,
        PVAL = pval,
        ADJPVAL = padj
    )
}

aldex2 <- function(se, grp, ref = NULL) {
    
    se <- groupsAsFactors(se, grp = grp, ref = ref)
    reads <- SummarizedExperiment::assay(se) + 1L
    conditions <- as.character(SummarizedExperiment::colData(se)[["GROUP"]])
    suppressMessages({
        aldex_obj <- ALDEx2::aldex(reads = reads, conditions = conditions,  mc.samples = 128,
            test = "t", effect = TRUE, include.sample.summary = FALSE,
            verbose = FALSE, denom = "iqlr", iterate = FALSE)
    })
    
    taxanames <- rownames(aldex_obj)
    fc <- aldex_obj$effect
    pval <- aldex_obj$wi.ep
    padj <- aldex_obj$wi.eBH
    tibble::tibble(TAXA = taxanames, FC = fc, PVAL = pval, ADJPVAL = padj)
}

metagenomeSeq <- function(se, grp, ref = NULL) {
    
    se <- groupsAsFactors(se, grp = grp, ref = ref)
    counts <- SummarizedExperiment::assay(se)
    pheno_data <- Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::colData(se)))
    feature_data <- Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::rowData(se)))
    MRE <- metagenomeSeq::newMRexperiment(counts = counts, phenoData = pheno_data, featureData = feature_data)
    CSS <- metagenomeSeq::calcNormFactors(obj = MRE, p = metagenomeSeq::cumNormStatFast(MRE, rel = 0.1))
    normFactor <- CSS$normFactors
    normFactor <- log2(normFactor/stats::median(normFactor) + 1)
    metagenomeSeq::normFactors(MRE) <- normFactor
    conditions <- Biobase::pData(MRE)[["GROUP"]]
    design <- stats::model.matrix( ~ conditions + normFactor)
    suppressMessages({
        fit <- metagenomeSeq::fitZig(MRE, design, verbose = FALSE, useCSSoffset = FALSE,
            control = metagenomeSeq::zigControl(maxit = 1000))
    })
    coefficients <- tibble::as_tibble(fit@fit$coefficients, rownames = "taxonName") %>%
        magrittr::set_colnames(c("taxonName", "intercept", "FC", "normFactor"))
    pvalues <- metagenomeSeq::MRfulltable(fit, number = nrow(get("counts", Biobase::assayData(MRE)))) %>%
        tibble::as_tibble(rownames = "taxonName")
    coefficients
    dplyr::full_join(coefficients, pvalues) %>%
        dplyr::select(taxonName, FC, pvalues, adjPvalues) %>%
        magrittr::set_colnames(c("TAXA", "FC", "PVAL", "ADJPVAL"))
}


corncob <- function(se, grp, ref = NULL, bootstrap = FALSE, B = NULL, fdr_cutoff = 0.05) {
    
    se <- groupsAsFactors(se, grp = grp, ref = ref)
    
    if (methods::isClass(se, "TreeSummmarizedExperiment")) {
        
        ps <- mia::makePhyloseqFromTreeSummarizedExperiment(se)
        
    } else if (methods::isClass(se, "SummarizedExperiment")) {
        
        sam_tbl <- phyloseq::sample_data(as.data.frame(SummarizedExperiment::colData(se)))
        tax_tbl <- phyloseq::tax_table(as.data.frame(SummarizedExperiment::rowData(se)))
        otu_tbl <- phyloseq::otu_table(SummarizedExperiment::assay(se))
        
        ps <- phyloseq::phyloseq(sam_tbl, tax_tbl, otu_tbl)
        
    }
    
    design <- stats::as.formula(paste0("~ GROUP"))
    
    if (isFALSE(bootstrap)) {
        B <- NULL
    }
    
    corncob_res <- corncob::differentialTest(formula = design,
        phi.formula = design,
        formula_null = ~ 1,
        phi.formula_null = design,
        test = "Wald", boot = bootstrap,
        data = ps,
        fdr_cutoff = fdr_cutoff,
        B = B)
    
    taxonNames <- names(corncob_res[["p"]]) # Same names and order for the other vectors
    pvalues <- corncob_res[["p"]]
    adjp <- corncob_res[["p_fdr"]]
    FC <- vapply(corncob_res[["all_models"]], function(x) {
        
        if (any(is.na(x))) {
            NA
        } else {
            stats::coef(x)[2]
        }
    }, double(1))
    
    tibble::tibble(taxonNames, FC, pvalues, adjp) %>%
        magrittr::set_colnames(c("TAXA", "FC", "PVAL", "ADJPVAL"))
}


groupsAsFactors <- function(se, grp, ref = NULL) {
    col_data <- SummarizedExperiment::colData(se)
    if (!grp %in% colnames(col_data))
        stop("Can't find column '", grp, "' in the data.", call. = FALSE)
    if (any(is.na(col_data[[grp]]))) {
        na <- sum(is.na(col_data[[grp]]))
        warning(
            na, "rows with NAs in column '", grp, "' were dropped.", 
            call. = FALSE
        ) # dropping occurred in the line code below
    }
    col_data <- col_data[!is.na(col_data[[grp]]),]
    grps <- unique(sort(as.character(col_data[[grp]])))
    if (length(grps) != 2)
        stop(
            "Only two conditions must be compared. Column '", grp,
            "' doesn't have two conditions.", call. = FALSE
        )
    if (!is.null(ref)) {
        if (!ref %in% grps)
            stop(
                "Reference (", ref, ") not in the '", grp,
                "' column", call. = FALSE
            )
        first_grp <- grps[which(grps == ref)]
        col_data[["GROUP"]] <- 
            as.factor(ifelse(col_data[[grp]] == first_grp, 0, 1))
        message(ref, " was used as reference.")
    } else {
        col_data[["GROUP"]] <-
            as.factor(ifelse(col_data[[grp]] == grps[1], 0, 1))
        message(grps[1], " was used as reference.")
    }
    SummarizedExperiment::colData(se) <- col_data
    se 
}


dgeTmm <- function(se) {
    
    # NF stands for normalization factor
    counts = SummarizedExperiment::assay(se)
    NF <- edgeR::calcNormFactors(object = counts, method = "TMM")
    NF <- NF/exp(mean(log(NF)))
    dge <- edgeR::DGEList(counts = counts)
    dge$samples$norm.factors <- NF
    dge
}


## <- se[rowSums(SummarizedExperiment::assay(se)) > 0]
## se


# For plots ---------------------------------------------------------------

volcano_plot <- function(x) {
    x %>% 
        dplyr::mutate(
            sig = ifelse(
                .data$ADJPVAL > 0.1 | is.na(.data$ADJPVAL), "unchanged", "DA"
            )
        ) %>% 
        ggplot2::ggplot(ggplot2::aes(.data$FC, -log10(.data$PVAL))) +
        ggplot2::geom_point(
            ggplot2::aes(color = .data$sig), shape = 1, size = 2, stroke = 0.8
        ) +
        ggplot2::labs(y = "-log10(P-value)", x = "log2(Fold Change)") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.title = ggplot2::element_blank())
}


ancom_bc <- function(se, grp, ref = NULL) {
    se <- groupsAsFactors(se, grp, ref)
    ps <- mia::makePhyloseqFromTreeSummarizedExperiment(se)
    out = ANCOMBC::ancombc(phyloseq = ps, formula = grp, 
        p_adj_method = "bonferroni", zero_cut = 0.90, lib_cut = 1000, 
        group = grp, struc_zero = TRUE, neg_lb = TRUE, 
        tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, 
        global = TRUE)
    res = out$res
    tab_coef <- res$beta
    tab_p <- res$p_val
    tab_q <- res$q_val
    do.call(cbind, list(tab_coef, tab_p, tab_q)) %>% 
        magrittr::set_colnames(c("FC", "PVAL", "ADJPVAL")) %>%
        tibble::as_tibble(rownames = "TAXA")
}

quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
} 