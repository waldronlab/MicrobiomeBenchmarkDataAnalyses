
## This code is borrowed from the lefser package
## I need it here  because there are some non-exported functions,
## necessary for lefser2 to work.


fillPmatZmat <- function(group,
                         block,
                         expr_sub,
                         p.threshold)
{
    # creates a list of boolean vectors, each vector indicates
    # existance (TRUE) or absence (FALSE) of a class/sub-class combination
    combos <- apply(
        expand.grid(levels(group), levels(block)), 1L, paste0, collapse = "")
    combined <- paste0(as.character(group), as.character(block))
    logilist <- lapply(setNames(nm = sort(combos)), `==`, combined)
    
    ## uses Wilcoxon rank-sum test to test for significant differential abundances between
    ## subclasses of one class against subclasses of all othe classes; results are saved in
    ## "pval_mat" and "z_mat" matrices
    whichlist <- lapply(logilist, which)
    sblock <- seq_along(levels(block))
    texp_sub <- t(expr_sub)
    iters <- expand.grid(sblock, sblock + length(sblock))
    group_formats <- apply(iters, 1L, function(x) {
        ind <- unlist(whichlist[x])
        apply(texp_sub, 2L, function(g) {
            wx <- suppressWarnings(coin::wilcox_test(g ~ group, subset = ind))
            cbind.data.frame(
                p.value = coin::pvalue(wx), statistic = coin::statistic(wx)
            )
        })
    })
    
    res <- lapply(group_formats, function(x) do.call(rbind, x))
    pval_mat <- do.call(cbind, lapply(res, `[[`, "p.value"))
    z_mat <- do.call(cbind, lapply(res, `[[`, "statistic"))
    
    rownames(pval_mat) <- rownames(expr_sub)
    rownames(z_mat) <- rownames(expr_sub)
    
    ## converts "pval_mat" into boolean matrix "logical_pval_mat" where
    ## p-values <= wilcoxon.threshold
    logical_pval_mat <- pval_mat <= p.threshold * 2.0
    logical_pval_mat[is.na(logical_pval_mat)] <- FALSE
    
    ## determines which rows (features) have all p-values<=0.05
    ## and selects such rows from the matrix of z-statistics
    sub <- apply(logical_pval_mat, 1L, all)
    z_mat_sub <- z_mat[sub, , drop = FALSE]
    # confirms that z-statistics of a row all have the same sign
    sub <- abs(rowSums(z_mat_sub)) == rowSums(abs(z_mat_sub))
    expr_sub[names(sub[sub]), , drop = FALSE]
}

## ensures that more than half of the values in each for each feature are unique
## if that is not the case then a count value is altered by adding it to a small value
## generated via normal distribution with mean=0 and sd=5% of the count value
createUniqueValues <- function(df, group){
    orderedrows <- rownames(df)
    splitdf <- split(df, group)
    maxim <- vapply(table(group), function(x) max(x * 0.5, 4), numeric(1L))
    for (i in seq_along(splitdf)) {
        sdat <- splitdf[[i]]
        splitdf[[i]][] <- lapply(sdat, function(cols) {
            if (length(unique(cols)) > maxim[i])
                cols
            else
                abs(cols + rnorm(
                    length(cols), mean = 0, sd = max(cols * 0.05, 0.01))
                )
        })
    }
    df <- do.call(rbind, unname(splitdf))
    df[match(orderedrows, rownames(df)),, drop = FALSE]
}

contastWithinClassesOrFewPerClass <-
    function(expr_sub_t_df, rand_s, min_cl, ncl, groups) {
        cols <- expr_sub_t_df[rand_s, , drop = FALSE]
        cls <- expr_sub_t_df$class[rand_s]
        # if the number of classes is less than the actual number (typically two)
        # of classes in the dataframe then return TRUE
        if (length(unique(cls)) < ncl) {
            return (TRUE)
        }
        # detect if for each class there are not fewer than the minimum (min_cl) number of samples
        if (any(table(cls) < min_cl)) {
            return (TRUE)
        }
        # separate the randomly selected samples (cols) into a list of the two classes
        drops <- c("class")
        by_class <-
            lapply(seq_along(groups), function(x) {
                cols[cols[, "class"] == groups[x],!(names(cols) %in% drops)]
            })
        
        # makes sure that within each class all features have at least min_cl unique count values
        for (i in seq_along(groups)) {
            unique_counts_per_microb = apply(by_class[[i]], 2, function(x) {
                length(unique(x))
            })
            if ((any(unique_counts_per_microb <= min_cl) &
                 min_cl > 1) |
                (min_cl == 1 & any(unique_counts_per_microb <= 1))) {
                return (TRUE)
            }
        }
        return (FALSE)
        
    }

ldaFunction <- function (data, lfk, rfk, min_cl, ncl, groups) {
    # test 1000 samples for contrast within classes per feature
    # and that there is at least a minimum number of samples per class
    for (j in 1:1000) {
        rand_s <- sample(seq_len(lfk), rfk, replace = TRUE)
        if (!contastWithinClassesOrFewPerClass(data, rand_s, min_cl, ncl, groups)) {
            break
        }
    }
    # lda with rfk number of samples
    lda.fit <- MASS::lda(class ~ ., data = data, subset = rand_s)
    # coefficients that transform observations to discriminants
    w <- lda.fit$scaling[, 1]
    # scaling of lda coefficients
    w.unit <- w / sqrt(sum(w ^ 2))
    sub_d <- data[rand_s,]
    ss <- sub_d[,-match("class", colnames(sub_d))]
    xy.matrix <- as.matrix(ss)
    # the original matrix is transformed
    LD <- xy.matrix %*% w.unit
    # effect size is calculated as difference between averaged disciminants
    # of two classes
    effect_size <-
        abs(mean(LD[sub_d[, "class"] == 0]) - mean(LD[sub_d[, "class"] == 1]))
    # scaling lda coefficients by the efect size
    scal <- w.unit * effect_size
    # mean count values per fclass per feature
    rres <- lda.fit$means
    rowns <- rownames(rres)
    lenc <- length(colnames(rres))
    
    coeff <- vector("numeric", length(scal))
    for (v in seq_along(scal)) {
        if (!is.na(scal[v])) {
            coeff[v] <- abs(scal[v])
        } else{
            coeff[v] <- 0
        }
        
    }
    # count value differences between means of two classes for each feature
    lda.means.diff <- (lda.fit$means[2,] - lda.fit$means[1,])
    # difference between a feature's class means and effect size adjusted lda coefficient
    # are averaged for each feature
    (lda.means.diff + coeff) / 2
}

.numeric01 <- function(x) {
    x <- as.factor(x)
    uvals <- levels(x)
    ifelse(x == uvals[1L], 0L, 1L)
}

.trunc <- function(scores_df, trim.names){
    Names <- gsub("`", "", scores_df[["Names"]])
    if (trim.names) {
        listNames <- strsplit(Names, "\\||\\.")
        Names <- vapply(listNames, tail, character(1L), 1L)
    }
    scores_df[["Names"]] <- Names
    return(scores_df)
}
