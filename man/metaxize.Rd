% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaxa2_proc.R
\name{metaxize}
\alias{metaxize}
\title{\code{metaxize} converts Metaxa2 taxnomomy-reliability output files to a single tibble expressing proportional abundance of reads at a chosen taxonomic level and at a chosen reliability score cutoff}
\usage{
metaxize(x, tlev = "family", rscore = 60, regext = "CT\\\\d+")
}
\arguments{
\item{x}{A filepath to a directory containing (only) taxonomy-reliability.txt files}

\item{tlev}{A string specifying the taxonomic level at which results should be aggregated}

\item{rscore}{A number between 0 and 100 specifying the relaiability score cutoff; only reads having a reliability score >= \code{rscore} will pass filter}
}
\value{
A tibble of sample-wise taxonomic proportions
}
\description{
\code{metaxize} converts Metaxa2 taxnomomy-reliability output files to a single tibble expressing proportional abundance of reads at a chosen taxonomic level and at a chosen reliability score cutoff
}
