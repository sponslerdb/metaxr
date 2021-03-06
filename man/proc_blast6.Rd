% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vsearch_proc.R
\name{proc_blast6}
\alias{proc_blast6}
\title{\code{proc_blast6} processes blast6 input, assigning column names, adding a sample ID field, and extracting accession numbers into their own field}
\usage{
proc_blast6(x)
}
\arguments{
\item{x}{A VSEARCH blast6 output file}
}
\value{
A tibble
}
\description{
\code{proc_blast6} processes blast6 input, assigning column names, adding a sample ID field, and extracting accession numbers into their own field
}
