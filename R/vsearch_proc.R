#' \code{load_blast6} imports blast6-formatted VSEARCH output files from a given directory and merges them into a single tibble
#'
#' @param x A file path to directory containing (only) VSEARCH output files
#' @return A tibble
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
  read_tsv(x, col_names = c("qseqid", "sseqid", "pident", # read in file and add column names
                            "length", "mismatch", "gapopen",
                            "qstart", "qend", "sstart", "send",
                            "evalue", "bitscore"),
           guess_max = 100000) %>% # it is necessary to increase guess_max so that it will assign the correct data type and avoid parsing errors
    mutate(sample = rep(str_extract(x, "CT\\d+"))) %>% # add sample field
    separate("sseqid", c("a", "gi", "b", "sseqid", "z"), fill = "warn", sep = "\\|", remove = TRUE) %>% # pull out accession number from sseqid field
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
#' @return A list consisting of NCBI names and nodes for input into \code{taxonomize}
prepNCBI <- function() {
  list(taxonomizr::read.nodes('~/metaxr/inst/extdata/nodes.dmp'),
       taxonomizr::read.names('~/metaxr/inst/extdata/names.dmp'))
}


#' \code{taxonomize} is a wrapper for the taxonomizr::getTaxonomy and taxonomizr::accessionToTaxa; assigns taxonomic data to NCBI accession numbers using the data base creates by \code{prepNCBI}
#'
#' @param x A tibble produced by \code{load_blast6} and \code{top_hit}
#' @param db A data base created by \code{prepNCBI}
#' @return A tibble
taxonomize <- function(x, db) {
  bind_cols(x, as.tibble(
    taxonomizr::getTaxonomy(
      taxonomizr::accessionToTaxa(
        x$sseqid, "~/metaxr/inst/extdata/accessionTaxa.sql"),
      db[[1]], db[[2]])))
}

#' \code{tally_gen} sums reads by genus for each sample, calculates proportional abundance, and culls rare genera whose proportional abundance is less than min_prop
#'
#' @param x A tibble produced by \code{taxonomize}
#' @param min_prop A real number specifying the minumum proportional abundance below which a taxon will be dropped as a likely false positive
#' @return A tibble
tally_gen <- function(x, min_prop = 0.001) {
  x %>% group_by(sample, genus) %>%
    tally() %>% # sum reads by genus within each sample
    mutate(gen_prop = n/sum(n)) %>% # turn counts to proportions
    left_join(x, by = c("sample", "genus")) %>% # join summary data back into original data frame by sample and genus
    select(sample, family, genus, gen_prop) %>%
    filter(gen_prop >= min_prop) %>% # remove genera falling below min_prop in abundance
    distinct() %>% # remove duplicate rows created by the left_join
    na.omit # hits to accession numbers that do not have taxonomic data at family and genus level (these are rare) are omitted
}



vsearch_eval <- function(x) {
  x %>%
    group_by(sample) %>%
    summarize(ITS2.reads_to_gen = n()) %>%
    add_meta
}

# ITS2_gen_class_rate <- add_meta(ITS2_vsearch) %>%
#   group_by(sample) %>%
#   summarize(reads_classified = n()) %>%
#   full_join(filter(read_summary, marker == "ITS2_reads"), by = "sample") %>%
#   mutate(class_rate = reads_classified/reads) %>%
#   select(sample, site, date, reads, reads_classified, class_rate)
