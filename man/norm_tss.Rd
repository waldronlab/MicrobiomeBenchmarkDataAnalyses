% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DA_methods.R
\name{norm_tss}
\alias{norm_tss}
\title{TSS normalization}
\usage{
norm_tss(mat, total_sum = 1e+06)
}
\arguments{
\item{mat}{A numeric matrix of counts.}

\item{total_sum}{The tolal sum of the scaling, e.g. 100 or 1e.}
}
\value{
A TSS-normalized matrix,
}
\description{
\code{norm_TSS} Applies TSS normalization to a matrix of count data.
}
