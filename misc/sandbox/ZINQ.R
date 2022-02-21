library(ZINQ)
library(MicrobiomeBenchmarkData)
library(magrittr)

tse <- getDataset("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)[[1]]
ps <- mia::makePhyloseqFromTreeSummarizedExperiment(tse)

zinq <- function(x, grp, ref = NULL, pval_type = "Cauchy") {
    
    if (class(x) == "phyloseq") {
        m <- microbiome::abundances(x)
        taxa <- rownames(m)
        metadata <- microbiome::meta(x)
        
    } else if (class(x) %in% c("SummarizedExperiment", "TreeSummarizedExperiment")) {
        m <- SummarizedExperiment::assay(tse)
        taxa <- rownames(m)
        metadata <- SummarizedExperiment::colData(tse) |>
            as.data.frame()
    }
    
    list_of_abundances <- vector("list", length(taxa))
    names(list_of_abundances) <- taxa
    
    for (i in seq_along(taxa)) {
        abundance <- data.frame(m[taxa[i],])
        covariates <- metadata[rownames(abundance), grp]
        # return(covariates)
        df <- cbind(abundance = abundance, covariates)
        colnames(df) <- c("abundance", "covariate")
        
        if (!is.null(ref)) {
            df[["covariate"]] <- 
                stats::relevel(factor(df[["covariate"]]), ref = ref)
        } else {
            df[["covariate"]] <- factor(df[["covariate"]])
        }
        
        ## Calculate fold change
        num_lvl <- levels(df[["covariate"]])[1]
        denom_lvl <- levels(df[["covariate"]])[2]
        
        num <- mean(df$abundance[df$covariate == num_lvl])
        denom <- mean(df$abundance[df$covariate == denom_lvl])
    
        
        if (num >= denom) {
            log2FoldChange <- log2(num / denom)
        } else if (num < denom) {
            log2FoldChange <- -log2(denom / num)
        }
        
        output <- tryCatch(
            error = function(e) NULL, {
                
                ZINQ::ZINQ_tests(
                    formula.logistic = abundance ~ covariate, 
                    formula.quantile = abundance ~ covariate, 
                    C = "covariate", y_CorD = "D", data = df
                )
            })
        
            pvalues <- vector("double", 3)
            names(pvalues) <- c("log2FoldChange", "Cauchy", "MinP")
            
        if ( !is.null(output)) {
            pvalues[["log2FoldChange"]] <- log2FoldChange
            pvalues[["Cauchy"]] <- ZINQ::ZINQ_combination(output, method="Cauchy")
            pvalues[["MinP"]] <- ZINQ::ZINQ_combination(output, method="MinP")
            
        } else {
            pvalues[["log2FoldChange"]] <- log2FoldChange
            pvalues[["Cauchy"]] <- NA
            pvalues[["MinP"]] <- NA
        }
        
            list_of_abundances[[i]] <- pvalues
    }
    
    output <- do.call("rbind", list_of_abundances) |>
        as.data.frame()
    output[["adj_Cauchy"]] <- stats::p.adjust(output[["Cauchy"]], method = "fdr")
    output[["adj_MinP"]] <- stats::p.adjust(output[["MinP"]], method = "fdr")
    # return(output)
    statInfo <- output
    pValMat <- output[,c("Cauchy", "adj_Cauchy")]
    colnames(pValMat) <- c("rawP", "adjP")
    name <- "ZINQ"
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}

x = zinq(x = tse, grp = "hmp_body_subsite", ref = "subgingival_plaque")




values = assay(tse["OTU_97.34979",])



colData(tse)[["hmp_body_subsite"]]

values_df <- values %>%
    t() %>% 
    as.data.frame()



