% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{plot_enrichment}
\alias{plot_enrichment}
\title{Plot enrichment object}
\usage{
plot_enrichment(enrichment, enrichment_col, levels_to_plot = NULL, conditions)
}
\arguments{
\item{enrichment}{Enrichment output from benchdamic.}

\item{enrichment_col}{Column with enrichment annotations.}

\item{levels_to_plot}{Levels used for plotting. Default all.}

\item{conditions}{Named vector. The names must be condB and condA. Example:
c(condB = 'control', condA = 'treatment')}
}
\value{
A ggplot object
}
\description{
\code{plot_enrichment} make plot from benchdamic enrichment object.
}
