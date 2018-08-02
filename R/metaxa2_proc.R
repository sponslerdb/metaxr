#' \code{metaxize} converts Metaxa2 taxnomomy-reliability output files to a single tibble expressing proportional abundance of reads at a chosen taxonomic level and at a chosen reliability score cutoff
#'
#' @param x A filepath to a directory containing (only) taxonomy-reliability.txt files
#' @param tlev A string specifying the taxonomic level at which results should be aggregated
#' @param rscore A number between 0 and 100 specifying the relaiability score cutoff; only reads having a reliability score >= \code{rscore} will pass filter
#' @return A tibble of sample-wise taxonomic proportions
metaxize <- function(x, tlev = "family", rscore = 60, regext = "CT\\d+") {
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
      mutate(sample = rep(str_extract(y, regext))) %>%
      select(sample, taxon, prop)
  }
  )
  bind_rows(tlist)
}

#' \code{consensus} returns taxa common to both dataset x and dataset y
#'
#' @param x,y Tibbles output from \code{metaxize}
#' @return
consensus_ff <- function(x, y, z) {
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

  colnames(out) <- c("sample", "taxon",
                     paste(deparse(substitute(x)), "prop", sep = "."),
                     paste(deparse(substitute(y)), "prop", sep = "."),
                     paste(deparse(substitute(z)), "prop", sep = "."))

  out <- out %>%
    #replace_na(list(ITS2_fam.prop = 0, rbcL_fam.prop = 0, trnL_fam.prop = 0)) %>%
    rowwise() %>%
    mutate(med_prop = median(c(ITS2_fam.prop, rbcL_fam.prop, trnL_fam.prop), na.rm = TRUE),
           mean_prop = mean(c(ITS2_fam.prop, rbcL_fam.prop, trnL_fam.prop), na.rm = TRUE))

  return(out)
}

#' \code{}
#'
#'
#' @return
consensus_fg <- function(gen, fam) {
  inner_join(gen, fam, by = c("sample", c("family" = "taxon"))) %>%
    add_metadata() %>%
    arrange(site, date)
}


#' \code{}
#'
#'
#' @return
add_metadata <- function(x) {
  key <- read_csv("~/metaxr/inst/extdata/CT_sample_key.csv")
  full_join(x, key, by = "sample") %>%
    mutate(date = as.Date(date, "%m/%d/%Y")) %>%
    select(sample, colony, site, date, family, ITS2_fam.prop,
           rbcL_fam.prop, trnL_fam.prop, med_prop, mean_prop, ITS2_genus = genus)
}


#' \code{metaxeval} reports classification rates
#'
#' @param x A filepath to a directory containing (only) taxonomy-reliability.txt files
#' @param tlev A string specifying the taxonomic level at which results should be aggregated
#' @param rscore A number between 0 and 100 specifying the relaiability score cutoff; only reads having a reliability score >= \code{rscore} will pass filter
#' @return
metaxeval <- function(x, tlev = "family", rscore = 60, regext = "CT\\d+") {
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
