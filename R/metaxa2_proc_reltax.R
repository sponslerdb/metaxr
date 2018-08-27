#' \code{metaxize_reltax} converts Metaxa2 taxnomomy-reliability output files to a single tibble expressing proportional abundance of reads at a chosen taxonomic level and at a chosen reliability score cutoff
#'
#' @param x A filepath to a directory containing (only) taxonomy-reliability.txt files
#' @param tlev A string specifying the taxonomic level at which results should be aggregated
#' @param rscore A number between 0 and 100 specifying the relaiability score cutoff; only reads having a reliability score >= \code{rscore} will pass filter
#' @return A tibble of sample-wise taxonomic proportions
metaxize_reltax <- function(x, tlev = "family", rscore = 80, regext = "CT\\d+") {
  rscore <- enquo(rscore)
  tlev <- enquo(tlev)
  files <- list.files(path = x,
             full.names = TRUE,
             ignore.case = TRUE)
  tlist <- map(files, function(y) {
    read_tsv(y, col_names = c("ID", "tax_tree")) %>%
      separate(tax_tree, c("kingdom", "phylum", "class", "order",
                           "family", "genus", "species", "x", "y"),
               fill = "warn", sep = ";", remove = TRUE) %>%
      select(!! tlev) %>%
      separate(!! tlev, c("taxon", "rel_score"), fill = "warn", sep = "\\(", remove = TRUE) %>%
      mutate(rel_score = as.numeric(str_sub(rel_score, 1, -2))) %>%
      filter(rel_score >= !! rscore) %>%
      group_by(taxon) %>%
      summarise(reads = n()) %>%
      mutate(prop = reads/sum(reads),
             taxon = str_remove(taxon, ".__")) %>%
      mutate(sample = rep(str_extract(y, regext))) %>%
      select(sample, taxon, prop)
  }
  )
  bind_rows(tlist)
}

#' \code{consensus_xyz} returns taxa common to at least two of datasets x, y, and z
#'
#' @param x,y,z Tibbles output from \code{metaxize}
#' @return
consensus_xyz_reltax <- function(x, y, z, tlev) {
  xy <- inner_join(x, y, by = c("sample", "taxon"))
  colnames(xy) <- c("sample", "taxon",
                    paste(deparse(substitute(x)), "prop", sep = "."),
                    paste(deparse(substitute(y)), "prop", sep = "."))
  xz <- inner_join(x, z, by = c("sample", "taxon"))
  colnames(xz) <- c("sample", "taxon",
                    paste(deparse(substitute(x)), "prop", sep = "."),
                    paste(deparse(substitute(z)), "prop", sep = "."))
  yz <- inner_join(y, z, by = c("sample", "taxon"))
  colnames(xz) <- c("sample", "taxon",
                    paste(deparse(substitute(y)), "prop", sep = "."),
                    paste(deparse(substitute(z)), "prop", sep = "."))

  sieve <- xy %>%
    full_join(xz, by = c("sample", "taxon")) %>%
    full_join(yz, by = c("sample", "taxon")) %>%
    select(sample, taxon)

  x_sieve <- semi_join(x, sieve, b = c("sample", "taxon"))

  y_sieve <- semi_join(y, sieve, b = c("sample", "taxon"))

  z_sieve <- semi_join(z, sieve, b = c("sample", "taxon"))

  out <- full_join(x_sieve, y_sieve, by = c("sample", "taxon")) %>%
    full_join(z_sieve, by = c("sample", "taxon"))

  colnames(out) <- c("sample", tlev,
                     "X.prop",
                     "Y.prop",
                     "Z.prop")


  out <- out %>%
    rowwise() %>%
    mutate(med_prop = median(c(X.prop,
                               Y.prop,
                               Z.prop), na.rm = TRUE),
           mean_prop = mean(c(X.prop,
                              Y.prop,
                              Z.prop), na.rm = TRUE),
           max_prop = max(c(X.prop,
                              Y.prop,
                              Z.prop), na.rm = TRUE))


  colnames(out) <- c("sample", tlev,
                     paste(deparse(substitute(x)), "prop", sep = "."),
                     paste(deparse(substitute(y)), "prop", sep = "."),
                     paste(deparse(substitute(z)), "prop", sep = "."),
                     "med_prop",
                     "mean_prop",
                     "max_prop")

  return(out)
}


#' \code{consensus_xy} returns taxa common to datasets x, y
#'
#' @param x,y Tibbles output from \code{metaxize}
#' @return
consensus_xy_reltax <- function(x, y) {
  xy <- inner_join(x, y, by = c("sample", "family", "genus"))
  colnames(xy) <- c("sample",
                    "family",
                    "genus",
                    paste(deparse(substitute(x)), "prop", sep = "."),
                    paste(deparse(substitute(y)), "prop", sep = "."))
  return(xy)
#
#   out <- xy %>%
#     rowwise() %>%
#     mutate(mean_prop = mean(c(X.prop,
#                               Y.prop), na.rm = TRUE),
#            max_prop = max(c(X.prop,
#                             Y.prop), na.rm = TRUE))
#
#   colnames(out) <- c("sample", tlev,
#                      paste(deparse(substitute(x)), "prop", sep = "."),
#                      paste(deparse(substitute(y)), "prop", sep = "."),
#                      "mean_prop",
#                      "max_prop")
#
  # return(out)
}

#' \code{consensus_fg} filters genus-level data set by family-level data set
#'
#' @param gen genus data set
#' @param fam family data set
#' @return left_join of family and genus such that all families are retained but only genera with family-level hits are retained
consensus_fg_reltax <- function(gen, fam, min) {
  left_join(fam, gen, by = c("sample", "family")) %>%
    filter(max_prop >= min) %>%
    select(-prop)
}


#' \code{add_meta} joins a set of metadata to data by sample field
#'
#' @param x a filepath to metadata file containing site field shared with data
#' @return
add_meta_reltax <- function(x) {
  key <- read_csv("~/metaxr/inst/extdata/CT_sample_key.csv") %>%
    select(sample, colony, site, date)
  full_join(x, key, by = "sample") %>%
    mutate(date = as.Date(date, "%m/%d/%Y")) %>%
    select(sample, colony, site, date, everything()) %>%
    arrange(site, date)
}


#' \code{metaxeval} reports classification rates
#'
#' @param x A filepath to a directory containing (only) taxonomy-reliability.txt files
#' @param tlev A string specifying the taxonomic level at which results should be aggregated
#' @param rscore A number between 0 and 100 specifying the relaiability score cutoff; only reads having a reliability score >= \code{rscore} will pass filter
#' @return a tibble with raw read count, classified read count (reads passing rel_score filter), and classification rate
metaxeval_reltax <- function(x, tlev = "family", rscore = 60, regext = "CT\\d+") {
  rscore <- enquo(rscore)
  tlev <- enquo(tlev)
  files <- list.files(path = x,
                      full.names = TRUE,
                      ignore.case = TRUE)
  raw <- map(files, function(y) {
    read_tsv(y, col_names = c("ID", "tax_tree")) %>%
      separate(tax_tree, c("kingdom", "phylum", "class", "order",
                           "family", "genus", "species", "x", "y"),
               fill = "warn", sep = ";", remove = TRUE) %>%
      select(!! tlev) %>%
      separate(!! tlev, c("taxon", "rel_score"), fill = "warn", sep = "\\(", remove = TRUE) %>%
      mutate(rel_score = as.numeric(str_sub(rel_score, 1, -2))) %>%
      mutate(sample = rep(str_extract(y, regext))) %>%
      select(sample, rel_score)
  }
  )
  raw_out <- raw %>%
    bind_rows() %>%
    group_by(sample) %>%
    nest() %>%
    mutate(raw_reads = as.integer(map(data, nrow))) %>%
    select(sample, raw_reads)

  pass_filter <- raw %>%
    bind_rows() %>%
    filter(rel_score >= !! rscore) %>%
    group_by(sample) %>%
    nest() %>%
    mutate(class_reads = as.integer(map(data, nrow))) %>%
    select(sample, class_reads)

  eval <- full_join(raw_out, pass_filter, by = "sample") %>%
    mutate(class_rate = class_reads / raw_reads)
}
