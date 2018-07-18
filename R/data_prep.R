#' \code{metaxize} converts Metaxa2 taxnomomy-reliability output files to a single tibble expressing proportional abundance of reads at a chosen taxonomic level and at a chosen reliability score cutoff
#'
#' @param x A filepath to a directory containing (only) taxonomy-reliability.txt files
#' @param tlev A string specifying the taxonomic level at which results should be aggregated
#' @param rscore A number between 0 and 100 specifying the relaiability score cutoff; only reads having a reliability score >= \code{rscore} will pass filter
#' @return A tibble of sample-wise taxonomic proportions
metaxize <- function(x, tlev = "family", rscore = 60) {
  lev <- noquote(tlev)
  files <- list.files(path = x,
             full.names = TRUE,
             ignore.case = TRUE)
  tlist <- map(files, function(y) {
    read_tsv(y, col_names = c("ID", "tax_tree")) %>%
      separate(tax_tree, c("kingdom", "phylum", "class", "order",
                           "family", "genus", "species", "x", "y"),
               fill = "warn", sep = ";", remove = TRUE) %>%
      select(lev) %>%
      separate(lev, c("taxon", "rel_score"), fill = "warn", sep = "\\(", remove = TRUE) %>%
      mutate(rel_score = as.numeric(str_sub(rel_score, 1, -2))) %>%
      filter(rel_score >= rscore) %>%
      group_by(taxon) %>%
      summarise(reads = n()) %>%
      mutate(prop = reads/sum(reads),
             taxon = str_remove(taxon, ".__")) %>%
      mutate(sample = rep(str_extract(y, "CT\\d+"))) %>%
      select(sample, taxon, prop)
  }
  )
  bind_rows(tlist)
}

#' \code{consensus} returns taxa common to both dataset x and dataset y
#'
#' @param x,y Tibbles output from \code{metaxize}
consensus <- function(x, y) {
  inner_join(x, y, by = c("sample", "taxon"))
}

#'
#'
#'
load_blast6 <- function(x) {
  read_tsv(x,
           col_names = c("qseqid", "sseqid", "pident",
                         "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send",
                         "evalue", "bitscore"),
           guess_max = 100000) %>%
    mutate(sample = rep(str_extract(x, "CT\\d+"))) %>%
    separate("sseqid", c("a", "gi", "b", "sseqid", "z"), fill = "warn", sep = "\\|", remove = TRUE) %>%
    select(sample, qseqid, sseqid, pident, length, gapopen, qstart, qend, sstart, send, evalue, bitscore) %>%
    group_by(qseqid) %>%
    slice(1)
}

#'
#'
#'
load_all <- function(x) {
  map(list.files(path = x,
                 full.names = TRUE,
                 ignore.case = TRUE),
      load_blast6) %>%
    bind_rows()
}

#'
#'
#'
prepNCBI <- function() {
  list(taxonomizr::read.nodes('~/metaxr/data/nodes.dmp'),
       taxonomizr::read.names('~/metaxr/data/names.dmp'))
}


#'taxonomize
#'
#'
taxonomize <- function(x, ncbi) {
  bind_cols(x, as.tibble(
    getTaxonomy(
      taxonomizr::accessionToTaxa(
        x$sseqid, "~/metaxr/data/accessionTaxa.sql"),
      ncbi[[1]], ncbi[[2]])))
}
