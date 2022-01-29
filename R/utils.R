## Hide messages
quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
} 


## Log Fold Change
log_fold_change <- function(num, denom) {
    if (num >= denom) {
        return(log2(num / denom))
    } else if (num < denom) {
        return(-log2(denom / num))
    }
}