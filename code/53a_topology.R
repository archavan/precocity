# Arun Chavan
# Started: 2026-08-28

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

resdir <- here("results/53a_topology")
fs::dir_create(resdir)

scmdir <- here("results/scm/main")

source(here("code/utilities_scm.R"))
source(here("code/utilities_scm_prep.R"))
source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################

# Extract topologies from the full tree to avoid lineages with a single species
# (which fail getMRCA) that exist in the smaller 432 species tree (e.g.
# Tubulidentata).

taxa <- read_csv(here("data/taxa/species.csv")) |>
  select(-binomial) |>
  rename(binomial = binomial_upham19) |>
  filter(!is.na(binomial)) |>
  mutate(precocity = "dummy") # to get build_treelist() to work

treelist <- build_treelist(taxa)

treeindices <- read_csv(here(scmdir, "treeindices.csv")) |>
  mutate(treetype = ifelse(treename == "consensus", "consensus", "sample"))

stopifnot(setequal(treeindices$treename, names(treelist)))

pp_eutheria <- read_csv(here(scmdir, "posterior-prob.csv")) |>
  filter(clade == "Eutheria") |>
  select(treename, altricial, intermediate, precocial)

###############################################################################
# tripartitions ===============================================================
###############################################################################

# Each element is one three-way split to resolve. A clade is one or more named
# taxa; where the tripartition needs a group with no name in the taxonomy
# (the mouse-related clade), it is given as the union of its named parts.

tripartitions <- list(
  eutheria = list(
    Xenarthra = "Xenarthra",
    Afrotheria = "Afrotheria",
    Boreoeutheria = "Boreoeutheria"
  ),
  afrotheria = list(
    Afroinsectivora = c("Afrosoricida", "Macroscelidea"),
    Tubulidentata = "Tubulidentata",
    Paenungulata = "Paenungulata"
  ),
  rodentia = list(
    Hystricomorpha = "Hystricomorpha",
    Sciuromorpha = "Sciuromorpha",
    mouse_related = c("Myomorpha", "Anomaluromorpha", "Castorimorpha")
  ),
  mouse_related = list(
    Anomaluromorpha = "Anomaluromorpha",
    Castorimorpha = "Castorimorpha",
    Myomorpha = "Myomorpha"
  )
)

clade_species <- function(.taxon_names, .taxa) {
  map(.taxon_names, ~ get_species_in_taxon(.x, .taxa)) |> list_c()
}

tripartition_species <- map(tripartitions, ~ map(.x, clade_species, taxa))

walk(tripartition_species, function(.clades) {
  stopifnot(lengths(.clades) > 0)
  stopifnot(!anyDuplicated(list_c(.clades)))
})

###############################################################################
# extract topology ============================================================
###############################################################################

# `getMRCA()` returns NULL for a clade represented by a single tip, so
# single-species clades are resolved to their tip index instead.
clade_node <- function(.phy, .spp) {
  if (length(.spp) == 1) {
    return(match(.spp, .phy$tip.label))
  }
  stopifnot(is.monophyletic(.phy, .spp))
  getMRCA(.phy, .spp)
}

is_clade <- function(.phy, .spp) {
  length(.spp) == 1 || is.monophyletic(.phy, .spp)
}

extract_tripartition_topology <- function(.phy, .clades) {
  if (!all(map_lgl(.clades, ~ is_clade(.phy, .x)))) {
    return("non-monophyletic")
  }

  pairs <- combn(names(.clades), 2)
  pairwise_mrca <- apply(pairs, 2, function(.x) {
    getMRCA(.phy, list_c(.clades[.x]))
  })

  if (n_distinct(pairwise_mrca) == 1) {
    return("polytomy")
  }

  anc_node <- names(which(table(pairwise_mrca) == 2))
  pairs_with_outer <- pairs[, which(pairwise_mrca == anc_node)]
  outermost <- intersect(pairs_with_outer[, 1], pairs_with_outer[, 2])
  sisters <- sort(setdiff(names(.clades), outermost))

  glue("({outermost}, ({paste0(sisters, collapse = ', ')}))") |> as.character()
}

topologies <- imap(
  tripartition_species,
  function(.clades, .name) {
    treeindices |>
      mutate(
        tripartition = .name,
        topology = map_chr(
          treename,
          ~ extract_tripartition_topology(treelist[[.x]], .clades)
        )
      ) |>
      select(treeindex, treename, treetype, tripartition, topology)
  }
)

###############################################################################
# posterior probability by topology ===========================================
###############################################################################

topology_pp <- topologies |>
  list_rbind() |>
  left_join(pp_eutheria, by = "treename") |>
  mutate(tripartition = fct(tripartition, names(tripartitions)))

p_topology_pp <- topology_pp |>
  filter(treetype == "sample") |>
  ggplot(aes(precocial, topology, fill = tripartition)) +
  geom_violin(scale = "width", linewidth = 0.25) +
  facet_grid(row = vars(tripartition), space = "free_y", scales = "free_y") +
  scale_y_discrete(name = "Topology", position = "right", labels = function(x) {
    str_replace_all(x, "_", "-") |>
      str_to_title()
  }) +
  scale_x_continuous(
    limits = c(0, 1),
    name = "PP(Preocial) at Eutherian node"
  ) +
  scale_fill_discrete(name = "Node", labels = function(x) {
    str_replace_all(x, "_", "-") |>
      str_to_title()
  }) +
  theme_classic() +
  theme(
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.widths = unit(2, "in"),
    panel.heights = unit(2, "in"),
    legend.position = "left"
  )

ggsave(
  here(resdir, "topology-pp-eutheria.pdf"),
  p_topology_pp,
  width = 7,
  height = 3.5
)

n_topologies <- topology_pp |>
  filter(treetype == "sample") |>
  count(tripartition, topology)

###############################################################################
# output ======================================================================
###############################################################################

write_tsv(topology_pp, here(resdir, "topology-vs-eutherian-pp.tsv"))
write_tsv(n_topologies, here(resdir, "topology-counts.tsv"))

# end #########################################################################
