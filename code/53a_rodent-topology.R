# Arun Chavan
# Started: 2026-08-25

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(ape)
library(ggtree)
library(conflicted)

conflicts_prefer(dplyr::filter(), purrr::map())

resdir <- here("results/53_rodent-topology")
fs::dir_create(resdir)

source(here("code/utilities_scm.R"))
source(here("code/utilities_scm_prep.R"))
source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################

taxa <- read_csv(here("data/taxa/species.csv")) |>
  filter(rank07 == "Rodentia") |>
  select(-binomial) |>
  rename(binomial = binomial_upham19) |>
  filter(!is.na(binomial)) |>
  mutate(precocity = "dummy") # to get build_treelist() to work

treelist <- build_treelist(taxa)

###############################################################################
# make plots ==================================================================
###############################################################################

plot_tree_with_labels <- function(.treename) {
  tr <- treelist[[.treename]]
  hl_nodes <- tibble(
    clade = c(
      "Anomaluromorpha",
      "Castorimorpha",
      "Hystricomorpha",
      "Myomorpha",
      "Sciuromorpha",
      "Caviidae",
      "Cricetidae",
      "Muridae",
      "Sciuridae",
      "Heteromyidae"
    )
  ) |>
    mutate(rank = map_chr(clade, ~ get_taxonomic_rank(.x, taxa))) |>
    mutate(
      monophyletic = map2_lgl(
        clade,
        rank,
        ~ is.monophyletic(tr, get_species_in_taxon(.x, taxa, .rank = .y))
      )
    ) |>
    filter(monophyletic) |>
    mutate(node = map_int(clade, ~ get_node(tr, .x, taxa)))

  max_tip_distance <- max(node.depth.edgelength(tr))

  ggtree(tr) +
    geom_highlight(
      data = hl_nodes |> filter(rank == "rank08"),
      mapping = aes(fill = clade, node = node),
      to.bottom = TRUE,
      type = "gradient",
      alpha = 0.2
    ) +
    geom_cladelab(
      data = hl_nodes |> filter(rank == "family"),
      mapping = aes(node = node, label = clade),
      offset = 0.2,
      barsize = 1,
      fontsize = 3
    ) +
    geom_cladelab(
      data = hl_nodes |> filter(rank == "rank08"),
      mapping = aes(node = node, label = clade, color = clade),
      offset = 9,
      barsize = 1,
      fontsize = 3
    ) +
    scale_x_continuous(limits = c(NA, max_tip_distance + 25)) +
    labs(title = .treename) +
    theme(
      legend.position = "none"
    )
}

quartz(
  type = "pdf",
  file = here(resdir, "trees_labelled.pdf"),
  width = 10,
  height = 6
)
for (i in names(treelist)) {
  print(plot_tree_with_labels(i))
}
dev.off()

###############################################################################
# extract topology ============================================================
###############################################################################

extract_triparition_topology <- function(.phy, .clades) {
  clades <- set_names(.clades)
  nodes <- map_int(clades, ~ get_node(.phy, .x, taxa))
  pairs <- combn(names(nodes), 2)
  pairwise_mrca <- apply(pairs, 2, function(x) {
    getMRCA(.phy, nodes[x])
  })
  anc_node <- names(table(pairwise_mrca)[which(table(pairwise_mrca) == 2)])
  pairs_with_outer <- pairs[, which(pairwise_mrca == anc_node)]
  outermost <- intersect(pairs_with_outer[, 1], pairs_with_outer[, 2])
  sisters <- sort(setdiff(clades, outermost))
  glue(
    r"[
  ({outermost}, ({paste0(sisters, collapse = ", ")}))
  ]"
  )
}

trees <- tibble(treename = names(treelist))

trees <- trees |>
  mutate(
    topology_mrc = map_chr(
      treename,
      ~ extract_triparition_topology(
        treelist[[.x]],
        c("Anomaluromorpha", "Castorimorpha", "Myomorpha")
      )
    )
  )


trees |> print(n = Inf)

# write
write_tsv(trees, here(resdir, "mrc-topology.tsv"))

# end #########################################################################
