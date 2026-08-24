# Arun Chavan
# Started: 2026-08-23

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(phytools)
library(ape)
library(conflicted)

conflicts_prefer(purrr::map, dplyr::filter)

clrs <- c(
  altricial = '#fc8d59',
  intermediate = '#ffffbf',
  precocial = '#91bfdb'
)

analysis <- "looo"
resdir <- here("results/scm", analysis)

###############################################################################
# data ========================================================================
###############################################################################

orders <- fs::dir_ls(here(resdir), type = "directory") |>
  fs::path_rel(here(resdir)) |>
  as.character()

df <- tibble(leave_out = orders)

df <- df |>
  mutate(
    tipdata = map(leave_out, ~ read_rds(here(resdir, .x, "tipdata.rds")))
  ) |>
  mutate(
    tree = map(
      leave_out,
      ~ read_rds(here(resdir, .x, "consensus", "tree_pruned.rds"))
    )
  ) |>
  mutate(taxa = map(leave_out, ~ read_csv(here(resdir, .x, "taxa.csv")))) |>
  mutate(
    ace = map(leave_out, ~ read_rds(here(resdir, .x, "consensus", "ace.rds")))
  )

###############################################################################
# plot PP at crown major node ==============================================
###############################################################################

get_taxonomic_rank <- function(.taxon, .taxa) {
  tax_rank <- names(which(apply(.taxa, 2, function(x) any(grepl(.taxon, x)))))
  if (length(tax_rank) == 0) {
    stop("Taxon not found in data")
  }
  message(paste0(.taxon, " was found in the taxonomic rank: ", tax_rank))
  return(tax_rank)
}

get_species_in_taxon <- function(.taxon, .taxa) {
  tax_rank <- get_taxonomic_rank(.taxon, .taxa)
  .taxa$binomial[which(.taxa[[tax_rank]] == .taxon)]
}

get_node <- function(.phy, .taxon, .taxa) {
  spp_list <- get_species_in_taxon(.taxon, .taxa)
  getMRCA(.phy, spp_list)
}

get_pp_at_node <- function(.ace, .node) .ace[which(rownames(.ace) == .node), ]


get_pp_df <- function(.taxon) {
  df |>
    mutate(node = map2_dbl(tree, taxa, ~ get_node(.x, .taxon, .y))) |>
    mutate(pp = map2(ace, node, get_pp_at_node)) |>
    select(leave_out, pp) |>
    unnest_longer(col = pp, indices_to = "precocity")
}

plot_pp_df <- function(.taxon) {
  get_pp_df(.taxon) |>
    ggplot(aes(pp, leave_out, fill = precocity)) +
    geom_col() +
    scale_fill_manual(values = clrs, breaks = names(clrs), name = NULL) +
    scale_y_discrete(
      name = "Order left out",
      expand = expansion(0.04, 0)
    ) +
    scale_x_continuous(
      name = glue("PP at node: {.taxon}"),
      expand = expansion(0.025, 0)
    ) +
    coord_cartesian() +
    theme_classic()
}


for (i in c("Eutheria", "Theria", "Boreoeutheria")) {
  get_pp_df(i) |>
    write_csv(here(resdir, glue("pp_", i, ".csv")))

  ggsave(
    here(resdir, glue("pp_", i, ".pdf")),
    width = 5.5,
    height = 5,
    plot = plot_pp_df(i)
  )
}

# end #########################################################################
