# Arun R. Chavan
# Started: 2026-08-26

###############################################################################

#' Get the taxonomic rank for a taxon
#'
#' @param .taxon Character string. e.g. "Rodentia"
#' @param .taxa Data frame containing taxonomic information
#' @param verbose logical.
#'
#' @returns Character string. e.g. "rank07" or "family".
#' @export
#'
#' @examples
get_taxonomic_rank <- function(.taxon, .taxa, verbose = FALSE) {
  tax_rank <- names(which(apply(.taxa, 2, function(x) any(grepl(.taxon, x)))))
  if (length(tax_rank) == 0) {
    stop("Taxon not found in data")
  }
  if (length(tax_rank) > 1) {
    stop("Taxon found in multiple taxonomic ranks.")
  }
  if (verbose) {
    message(paste0(.taxon, " was found in the taxonomic rank: ", tax_rank))
  }
  return(tax_rank)
}


###############################################################################

#' Get a list of species in a taxon
#'
#' @param .taxon Character string. e.g. "Rodentia"
#' @param .taxa Data frame containing the taxonomic information. Must contain
#'    species names in column `binomial`.
#' @param .rank Taxonomic rank. Computed if it's not provided.
#'
#' @returns Character vector.
#' @export
#'
#' @examples
get_species_in_taxon <- function(.taxon, .taxa, .rank = NULL) {
  if (is.null(.rank)) {
    .rank <- get_taxonomic_rank(.taxon, .taxa)
  }
  .taxa$binomial[which(.taxa[[.rank]] == .taxon)]
}

###############################################################################

#' Get node number for the MRCA of a taxon
#'
#' @param .phy Tree
#' @param .taxon Taxon. Character string.
#' @param .taxa Data frame containing taxonomic information. Must contain
#'    species names in column `binomial`.
#'
#' @returns Integer.
#' @export
#'
#' @examples
get_node <- function(.phy, .taxon, .taxa) {
  spp_list <- get_species_in_taxon(.taxon, .taxa)
  stopifnot(is.monophyletic(.phy, spp_list))
  getMRCA(.phy, spp_list)
}

###############################################################################
