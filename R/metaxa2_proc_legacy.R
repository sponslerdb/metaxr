#' \code{metaxize} converts Metaxa2 output files to a single tibble expressing proportional abundance of reads at the family level
#'
#' @param x A filepath to a directory containing (only) taxonomy.txt files
#' @param min_prop A real number specifying the minumum proportional abundance below which a taxon will be dropped as a likely false positive
#' @return A tibble of sample-wise taxonomic proportions
metaxize <- function(x, regext = "CT\\d+", min_prop = 0.0001) {
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
             family = replace_na(family, "*undetermined*"),
             family = factor(family),
             sample = rep(str_extract(y, regext))) %>%
      select(sample, family, prop)
  }
  )
  # combine samples into single output tibble
  out <- bind_rows(tlist) %>% # collate sample-wise tibbles
    filter(family == "*undetermined*" | prop >= min_prop ) # trim to taxa comprising at least 0.1% of reads, but retain undetermined reads regardless of rarity

  return(out)
}


#' \code{consensus_xyz} returns taxa common to at least two of the marker datasets x, y, and z; reports median proportional read count
#'
#' @param x,y,z Tibbles output from \code{metaxize}
#' @param min_prop A real number specifying the minumum proportional abundance below which a taxon will be dropped as a likely false positive
#' @return a tibble
consensus_xyz <- function(x, y, z, min_prop) {
  # taxa common to markers x and y
  xy <- inner_join(x, y, by = c("sample", "family"))
  colnames(xy) <- c("sample", "family",
                    paste(deparse(substitute(x)), "prop", sep = "."),
                    paste(deparse(substitute(y)), "prop", sep = "."))
  # taxa common to markers x and z
  xz <- inner_join(x, z, by = c("sample", "family"))
  colnames(xz) <- c("sample", "family",
                    paste(deparse(substitute(x)), "prop", sep = "."),
                    paste(deparse(substitute(z)), "prop", sep = "."))
  # taxa common to markers y and z
  yz <- inner_join(y, z, by = c("sample", "family"))
  colnames(xz) <- c("sample", "family",
                    paste(deparse(substitute(y)), "prop", sep = "."),
                    paste(deparse(substitute(z)), "prop", sep = "."))
  # taxa common to any two of x, y, and z
  sieve <- xy %>%
    full_join(xz, by = c("sample", "family")) %>%
    full_join(yz, by = c("sample", "family")) %>%
    select(sample, family)
  # filter original marker datasets by the taxa in the sieve dataset created above
  x_sieve <- semi_join(x, sieve, b = c("sample", "family"))
  y_sieve <- semi_join(y, sieve, b = c("sample", "family"))
  z_sieve <- semi_join(z, sieve, b = c("sample", "family"))
  # join filtered marker datasets
  out <- full_join(x_sieve, y_sieve, by = c("sample", "family")) %>%
    full_join(z_sieve, by = c("sample", "family"))
  # rename columns
  colnames(out) <- c("sample", "family", "X.prop", "Y.prop", "Z.prop")
  # calculate median read count proportions
  out <- out %>%
    rowwise() %>%
    mutate(med_prop = median(c(X.prop,
                               Y.prop,
                               Z.prop), na.rm = TRUE)) %>%
    filter(family == "*undetermined*" | med_prop >= min_prop) %>% # filter to exclude taxa with a median proportional abundance of less than 0.1%
    group_by(sample) %>%
    mutate(scaled_prop = med_prop*(1/sum(med_prop))) # rescale median proportional abundance so that it totals to 1 for each sample
  # rename columns to bear marker names
  colnames(out) <- c("sample", "family",
                     paste(deparse(substitute(x)), "prop", sep = "."),
                     paste(deparse(substitute(y)), "prop", sep = "."),
                     paste(deparse(substitute(z)), "prop", sep = "."),
                     "med_prop",
                     "scaled_prop")
  return(out)
}


#' \code{add_meta} joins a set of metadata to data by sample field
#'
#' @param x a filepath to metadata file containing site field shared with data
#' @return a tibble containing site, hive, and date for each sample
add_meta <- function(x) {
  key <- read_csv("~/metaxr/inst/extdata/CT_sample_key.csv") %>%
    select(sample, colony, site, date)
  full_join(x, key, by = "sample") %>%
    mutate(date = lubridate::as_date(date),) %>%
    select(sample, colony, site, date, everything()) %>%
    arrange(site, date)
}


#' \code{consensus_fg} filters genus-level data set by family-level data set
#'
#' @param gen genus data set
#' @param fam family data set
#' @return left_join of family and genus such that all families are retained but only genera with family-level hits are retained
consensus_fg <- function(gen, fam) {
  left_join(fam, gen, by = c("sample", "family"))
}


#' \code{read_count} reports total number of reads for each library
#'
#' @param x A filepath to a directory containing (only) taxonomy-reliability.txt files
#' @return a tibble showing raw read count by library
read_count <- function(x, regext = "CT\\d+") {
  files <- list.files(path = x,
                      full.names = TRUE,
                      ignore.case = TRUE)
  tlist <- map(files, function(y) {
    read_tsv(y) %>%
      summarise(sequences = n()) %>%
      mutate(sample = rep(str_extract(y, regext))) %>%
      select(sample, sequences)
  }
  )
  out <- bind_rows(tlist) %>%
    add_meta()
  return(out)
}






metax_count <- function(x, regext = "CT\\d+", marker) {
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
      mutate(family = str_remove(family, ".__"),
             family = replace_na(family, "*undetermined*"),
             family = factor(family),
             sample = rep(str_extract(y, regext))) %>%
      select(sample, family)
  }
  )
  # combine samples into single output tibble
  det_reads <- bind_rows(tlist) %>% # collate sample-wise tibbles
    filter(family != "*undetermined*") %>%
    group_by(sample) %>%
    summarize(reads_to_fam = n())

  all_reads <- bind_rows(tlist) %>% # collate sample-wise tibbles
    group_by(sample) %>%
    summarize(reads = n())

  out <- full_join(all_reads, det_reads, by = "sample") %>%
    mutate(fam_class_rate = reads_to_fam/reads)

  colnames(out) <- c("sample",
                     paste(deparse(substitute(marker)), "reads", sep = "."),
                     paste(deparse(substitute(marker)), "reads_to_fam", sep = "."),
                     paste(deparse(substitute(marker)), "fam_class_rate", sep = "."))
  return(out)
}



tax_freq <- function(x, tax) {

  tax <- enquo(tax)

  samples <- x %>%
    select(site, sample) %>%
    group_by(site) %>%
    summarise(samples = n_distinct(sample))

  detections <- x %>%
    group_by(!!tax, site) %>%
    summarise(detections = n_distinct(sample)) %>%
    arrange(-detections)

  tax_freq <- full_join(samples, detections, by = "site") %>%
    mutate(freq = detections/samples) %>%
    arrange(-freq, site) %>%
    filter(!!tax != "*undetermined*") %>%
    select(site, !!tax, samples, detections, freq)
}

tax_abund <- function(x, tax, prop) {
  tax <- enquo(tax)
  prop <- enquo(prop)

  samples <- x %>%
    select(site, sample) %>%
    group_by(site) %>%
    summarise(samples = n_distinct(sample))

  abund <- x %>%
    select(sample, site, !!tax, !!prop) %>%
    distinct()

  abund_summary <- full_join(samples, abund, by = "site") %>%
    group_by(site, family) %>%
    summarise(mean_abund = sum(scaled_prop)/unique(samples),
              max_abund = max(scaled_prop)) %>%
    filter(!!tax != "*undetermined*") %>%
    arrange(-mean_abund)
}
