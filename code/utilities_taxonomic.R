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
  tax_rank <- names(which(apply(.taxa, 2, function(x) any(x == .taxon))))
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
#' Get the taxonomic classification of a taxon
#'
#' @param .taxon Taxon. Character string.
#' @param .taxa Data frame containing taxonomic information.
#'
#' @returns Character vector.
#' @export
#'
#' @examples
get_classification <- function(.taxon, .taxa) {
  ranks <- c(
    "rank03",
    "rank04",
    "rank05",
    "rank06",
    "rank07",
    "rank08",
    "family"
  )
  tax_rank <- get_taxonomic_rank(.taxon, .taxa)

  # if the rank is low, we want to start from infraclass to keep things short
  if (tax_rank %in% c("rank01", "rank02", "rank03")) {
    shortpath <- tax_rank
  } else {
    shortpath <- ranks[seq_len(which(ranks == tax_rank))]
  }

  taxonomy <- .taxa |>
    select(all_of(shortpath)) |>
    filter(.data[[tax_rank]] == .taxon) |>
    distinct()

  stopifnot(nrow(taxonomy) == 1)

  x <- unlist(taxonomy)
  x <- x[!is.na(x)]
  return(x)
}

###############################################################################
#' Get a list of taxa that are monophyletic in all trees
#'
#' @param .taxa Data frame containing taxonomic information. Must contain
#'   species names in column `binomial`.
#' @param .treelist List of trees.
#'
#' @returns Character vector of monophyletic taxa.
#' @export
#'
#' @examples
get_monophyletic_taxa <- function(.taxa, .treelist) {
  ranks <- c(paste0("rank", str_pad(seq_len(8), 2, pad = "0")), "family")
  taxnames <- map(ranks, function(x) {
    table(.taxa[[x]]) |> (function(y) y[y > 1])()
  }) |>
    unlist() |>
    names() |>
    unique() |>
    set_names()

  taxranks <- taxnames |>
    map_chr(~ get_taxonomic_rank(.x, .taxa))

  taxspecies <- map2(taxnames, taxranks, ~ get_species_in_taxon(.x, .taxa, .y))

  monophyletic_in_all <- map_lgl(taxspecies, function(.spp) {
    map_lgl(.treelist, ~ is.monophyletic(.x, .spp)) |>
      all()
  })

  monophyletic_in_all[monophyletic_in_all] |> names()
}

###############################################################################
