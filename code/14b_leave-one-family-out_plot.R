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

analysis <- "lofo"
resdir <- here("results/scm", analysis)

###############################################################################
# data ========================================================================
###############################################################################

df <- read_csv(here(resdir, "families.csv"))

df$classification <- vector("character", nrow(df))
for (i in seq_len(nrow(df))) {
  csfn_i <- paste0(
    na.omit(c(df$rank07[[i]], df$rank08[[i]], df$family[[i]])),
    collapse = " - "
  )
  df$classification[[i]] <- csfn_i
  rm(csfn_i)
}

df <- df |>
  rename(leave_out = family)

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
    left_join(select(df, leave_out, classification)) |>
    ggplot(aes(pp, classification, fill = precocity)) +
    geom_col() +
    scale_fill_manual(values = clrs, breaks = names(clrs), name = NULL) +
    scale_y_discrete(
      position = "right",
      name = "Family left out",
      expand = expansion(0.02, 0)
    ) +
    scale_x_continuous(
      name = glue("PP at node: {.taxon}"),
      expand = expansion(0.025, 0)
    ) +
    coord_cartesian() +
    theme_classic() +
    theme(legend.position = "top")
}


plot_pp_df("Eutheria")

for (i in c("Eutheria", "Theria", "Boreoeutheria")) {
  get_pp_df(i) |>
    write_csv(here(resdir, glue("pp_", i, ".csv")))

  ggsave(
    here(resdir, glue("pp_", i, ".pdf")),
    width = 9,
    height = 14,
    plot = plot_pp_df(i)
  )
}

# remove auto-saved full simmap summary. too large.
walk(df$leave_out, ~ unlink(here(resdir, .x, "consensus/simmap_summary.rds")))

# end #########################################################################
