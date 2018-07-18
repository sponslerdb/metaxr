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
consensus <- function(x, y) {
  out <- inner_join(x, y, by = c("sample", "taxon"))
  colnames(out) <- c("sample", "taxon",
                     paste(deparse(substitute(x)), "prop", sep = "."),
                     paste(deparse(substitute(y)), "prop", sep = "."))
  return(out)
}

#'
#'
#'
#'
fam_filter_gen <- function(gen, fam) {
  inner_join(gen, fam, by = c("sample", c("family" = "taxon"))) %>%
    #select(sample, family, trnL_fam.prop, its2_fam.prop, genus, its2_genus.prop = prop) %>%
    add_metadata()
}


#'
#'
#'
add_metadata <- function(x) {
  key <- read_csv("~/metaxr/inst/extdata/CT_sample_key.csv")
  full_join(x, key, by = "sample") %>%
    select(sample, colony, site, date, family, trnL_fam.prop,
           its2_fam.prop, genus, prop)
}


