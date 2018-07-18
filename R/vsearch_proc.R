#' \code{load_blast6} imports blast6-formatted VSEARCH output files from a given directory and merges them into a single tibble
#'
#' @param x A file path to directory containing (only) VSEARCH output files
load_blast6 <- function(x) {
  map(list.files(path = x,
                 full.names = TRUE,
                 ignore.case = TRUE),
      proc_blast6) %>%
    bind_rows()
}


#' \code{proc_blast6} processes blast6 input, assigning column names, adding a sample ID field, and extracting accession numbers into their own field
#'
#' @param x A VSEARCH blast6 output file
#' @return A tibble
proc_blast6 <- function(x) {
  read_tsv(x,
           col_names = c("qseqid", "sseqid", "pident",
                         "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send",
                         "evalue", "bitscore"),
           guess_max = 100000) %>%
    mutate(sample = rep(str_extract(x, "CT\\d+"))) %>%
    separate("sseqid", c("a", "gi", "b", "sseqid", "z"), fill = "warn", sep = "\\|", remove = TRUE) %>%
    select(sample, qseqid, sseqid, pident, length, gapopen, qstart, qend, sstart, send, evalue, bitscore)
}


#' \code{top_hit} extracts the top ranked subject sequence for each query sequence per sample
#'
#' @param x  A tibble produced by \code{load_blast6}
#' @return A tibble
top_hit <- function(x) {
  group_by(x, sample, qseqid) %>% slice(1)
}


#' \code{prepNCBI} prepares the NCBI databases that \code{taxonomizr} needs; this takes some time, so it is better to do it once outside of the functions that use the databases later
#'
#' @return A list consisting of NCBI names and nodes
prepNCBI <- function() {
  list(taxonomizr::read.nodes('~/metaxr/inst/extdata/nodes.dmp'),
       taxonomizr::read.names('~/metaxr/inst/extdata/names.dmp'))
}


#'taxonomize
#'
#'
taxonomize <- function(x, db) {
  bind_cols(x, as.tibble(
    getTaxonomy(
      taxonomizr::accessionToTaxa(
        x$sseqid, "~/metaxr/inst/extdata/accessionTaxa.sql"),
      db[[1]], db[[2]])))
}

#'
#'
#'
#'
tallyGen <- function(x) {
  x %>% group_by(sample, genus) %>%
    tally() %>%
    mutate(prop = n/sum(n)) %>%
    left_join(x, by = c("sample", "genus")) %>%
    select(sample, family, genus, prop) %>%
    distinct()
}
