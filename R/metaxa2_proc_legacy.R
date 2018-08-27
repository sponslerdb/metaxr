#' \code{metaxize} converts Metaxa2 output files to a single tibble expressing proportional abundance of reads at the family level
#'
#' @param x A filepath to a directory containing (only) taxonomy.txt files
#' @return A tibble of sample-wise taxonomic proportions
metaxize <- function(x, regext = "CT\\d+") {
  # list files in source directory
  files <- list.files(path = x,
                      full.names = TRUE,
                      ignore.case = TRUE)

  # process by sample
  tlist <- map(files, function(y) {
    read_tsv(y, col_names = c("ID", "tax_tree", "perID", "length", "relscore")) %>%
      separate(tax_tree, c("kingdom", "phylum", "class", "order",
                           "family", "genus", "species", "x", "y"),
               fill = "warn", sep = ";", remove = TRUE) %>%
      select(ID, family) %>%
      group_by(family) %>%
      summarise(reads = n()) %>%
      mutate(prop = reads/sum(reads),
             family = str_remove(family, ".__"),
             family = replace_na(family, "UNDETERMINED"),
             sample = rep(str_extract(y, regext))) %>%
      select(sample, family, prop)
  }
  )

  # combine samples into single output tibble
  out <- bind_rows(tlist) %>%
    filter(family == "UNDETERMINED" | prop >= 0.001 )

  return(out)
}


#' \code{consensus_xyz_fam} returns taxa common to at least two of datasets x, y, and z
#'
#' @param x,y,z Tibbles output from \code{metaxize}
#' @return
consensus_xyz <- function(x, y, z) {
  xy <- inner_join(x, y, by = c("sample", "family"))
  colnames(xy) <- c("sample", "family",
                    paste(deparse(substitute(x)), "prop", sep = "."),
                    paste(deparse(substitute(y)), "prop", sep = "."))
  xz <- inner_join(x, z, by = c("sample", "family"))
  colnames(xz) <- c("sample", "family",
                    paste(deparse(substitute(x)), "prop", sep = "."),
                    paste(deparse(substitute(z)), "prop", sep = "."))
  yz <- inner_join(y, z, by = c("sample", "family"))
  colnames(xz) <- c("sample", "family",
                    paste(deparse(substitute(y)), "prop", sep = "."),
                    paste(deparse(substitute(z)), "prop", sep = "."))

  sieve <- xy %>%
    full_join(xz, by = c("sample", "family")) %>%
    full_join(yz, by = c("sample", "family")) %>%
    select(sample, family)

  x_sieve <- semi_join(x, sieve, b = c("sample", "family"))

  y_sieve <- semi_join(y, sieve, b = c("sample", "family"))

  z_sieve <- semi_join(z, sieve, b = c("sample", "family"))

  out <- full_join(x_sieve, y_sieve, by = c("sample", "family")) %>%
    full_join(z_sieve, by = c("sample", "family"))

  colnames(out) <- c("sample", "family",
                     "X.prop",
                     "Y.prop",
                     "Z.prop")


  out <- out %>%
    rowwise() %>%
    mutate(med_prop = median(c(X.prop,
                               Y.prop,
                               Z.prop), na.rm = TRUE)) %>%
    filter(family == "UNDETERMINED" | med_prop >= 0.001)


  colnames(out) <- c("sample", "family",
                     paste(deparse(substitute(x)), "prop", sep = "."),
                     paste(deparse(substitute(y)), "prop", sep = "."),
                     paste(deparse(substitute(z)), "prop", sep = "."),
                     "med_prop")

  return(out)
}


#' \code{add_meta} joins a set of metadata to data by sample field
#'
#' @param x a filepath to metadata file containing site field shared with data
#' @return
add_meta <- function(x) {
  key <- read_csv("~/metaxr/inst/extdata/CT_sample_key.csv") %>%
    select(sample, colony, site, date)
  full_join(x, key, by = "sample") %>%
    mutate(date = lubridate::as_date(date)) %>%
    select(sample, colony, site, date, everything()) %>%
    arrange(site, date)
}


#' \code{consensus_fg} filters genus-level data set by family-level data set
#'
#' @param gen genus data set
#' @param fam family data set
#' @return left_join of family and genus such that all families are retained but only genera with family-level hits are retained
consensus_fg <- function(gen, fam) {
  left_join(fam, gen, by = c("sample", "family")) %>%
    select(-prop)
}


#' \code{metaxeval_leg} reports classification rates
#'
#' @param x A filepath to a directory containing (only) taxonomy-reliability.txt files
#' @param tlev A string specifying the taxonomic level at which results should be aggregated
#' @param rscore A number between 0 and 100 specifying the relaiability score cutoff; only reads having a reliability score >= \code{rscore} will pass filter
#' @return a tibble with raw read count, classified read count (reads passing rel_score filter), and classification rate

metaxeval <- function(x, regext = "CT\\d+") {
  files <- list.files(path = x,
                      full.names = TRUE,
                      ignore.case = TRUE)
  tlist <- map(files, function(y) {
    read_tsv(y, col_names = c("ID", "tax_tree", "perID", "length", "relscore")) %>%
      separate(tax_tree, c("kingdom", "phylum", "class", "order",
                           "family", "genus", "species", "x", "y"),
               fill = "warn", sep = ";", remove = TRUE) %>%
      summarise(sequences = n()) %>%
      mutate(sample = rep(str_extract(y, regext))) %>%
      select(sample, sequences)
  }
  )
  out <- bind_rows(tlist)
  return(out)
}

