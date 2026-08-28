# Arun Chavan
# Started: 2026-08-23

###############################################################################
# setup =======================================================================
###############################################################################

library(phytools)
library(tidyverse)
library(here)
library(conflicted)

conflicts_prefer(dplyr::filter, purrr::map)

clrs <- c(
  altricial = '#fc8d59',
  intermediate = '#ffffbf',
  precocial = '#91bfdb'
)

analysis_name <- "exclude-manual"
resdir <- here("results/scm", analysis_name)

source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################

# taxonomic information
spp <- read_csv(here("data/taxa/species.csv")) |>
  select(-binomial) |>
  rename(binomial = binomial_upham19)
tipdata <- read_rds(here(resdir, "tipdata.rds"))
taxa <- spp |> filter(binomial %in% names(tipdata))
stopifnot(setequal(taxa$binomial, names(tipdata)))

# consensus results
simmap_summary_consensus <- read_rds(here(
  resdir,
  "consensus/simmap_summary.rds"
))

# results from all trees including consensus
asr <- read_csv(here(resdir, "treeindices.csv"))
asr <- asr |>
  mutate(treetype = ifelse(treename == "consensus", "consensus", "sample")) |>
  mutate(
    respath = ifelse(
      treetype == "consensus",
      "consensus",
      paste0("sample/", treename)
    )
  ) |>
  mutate(tr = map(respath, ~ read_rds(here(resdir, .x, "tree_pruned.rds")))) |>
  mutate(fit = map(respath, ~ read_rds(here(resdir, .x, "model-fit.rds")))) |>
  mutate(ace = map(respath, ~ read_rds(here(resdir, .x, "ace.rds"))))

###############################################################################
# transition rates ============================================================
###############################################################################
get_qmatrices <- function(.fit) {
  fit_names <- attr(.fit, "row.names")
  map(seq_along(fit_names), ~ as.Qmatrix(attr(.fit, "models")[[.x]])) |>
    set_names(fit_names)
}

qtibbles <- function(.qmatlist) {
  qtib1 <- function(.qmat) {
    expand_grid(from = dimnames(.qmat)[[1]], to = dimnames(.qmat)[[2]]) |>
      filter(from != to) |>
      mutate(transition = paste0(from, " → ", to)) |>
      mutate(rate = map2_dbl(from, to, ~ .qmat[.x, .y]))
  }
  .qmatlist |>
    map(qtib1) |>
    bind_rows(.id = "model") |>
    mutate(model = str_remove(model, "fit_") |> str_to_upper())
}

asr <- asr |>
  mutate(rates = map(fit, ~ get_qmatrices(.x) |> qtibbles()))

qrates <- asr |>
  select(treeindex, treename, treetype, rates) |>
  unnest(col = rates)

p_qrates <- qrates |>
  filter(treetype != "consensus") |>
  ggplot(aes(rate, transition, fill = model)) +
  geom_violin(
    scale = "width",
    position = position_dodge(width = 0.65),
    linewidth = 0.25
  ) +
  scale_y_discrete(labels = str_to_title, name = "Transition") +
  scale_x_continuous(name = "Rate") +
  scale_fill_discrete(name = "Model") +
  theme_classic() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "grey"),
    panel.heights = unit(2.5, "in")
  )

quartz(
  type = "pdf",
  file = here(resdir, "fitted-rates.pdf"),
  width = 6,
  height = 3
)
print(p_qrates)
dev.off()

write_csv(qrates, here(resdir, "fitted-rates.csv"))

###############################################################################
# posterior probability distributions #########################################
###############################################################################

## get monophyletic groups for which to plot pp distributions =================
# should be monophyletic in all trees.
monophyletic <- get_monophyletic_taxa(taxa, asr$tr)

## extract PP values ==========================================================
get_pp_at_node <- function(.ace, .node) .ace[which(rownames(.ace) == .node), ]

get_pp_for_clades <- function(.ace, .phy, .taxa, .clades) {
  map(
    set_names(.clades),
    ~ get_pp_at_node(.ace, get_node(.phy, .x, .taxa))
  ) |>
    bind_rows(.id = "clade")
}

asr <- asr |>
  mutate(pp = map2(ace, tr, ~ get_pp_for_clades(.x, .y, taxa, monophyletic)))

pp_monophyletic <- asr |>
  select(treeindex, treename, treetype, pp) |>
  unnest(col = pp)

## write ======================================================================
write_csv(pp_monophyletic, here(resdir, "posterior-prob.csv"))

## plot =======================================================================

plot_pp_dist <- function(
  .taxon,
  .title = .taxon
) {
  df <- pp_monophyletic |>
    filter(treetype != "consensus") |>
    filter(clade == .taxon) |>
    select(treeindex, all_of(levels(tipdata))) |>
    pivot_longer(-treeindex, names_to = "precocity", values_to = "pp")

  ggplot(df, aes(treeindex, pp, fill = precocity)) +
    geom_col(width = 0.9, linewidth = 0) +
    scale_fill_manual(values = clrs) +
    scale_x_discrete(breaks = c(1, 25, 50, 75, 100), expand = c(0.04, 0.5)) +
    labs(
      y = "Posterior probability",
      x = "Sampled tree",
      title = .title,
      caption = paste0(get_classification(.taxon, taxa), collapse = " → ")
    ) +
    theme_bw(base_line_size = 0.25) +
    theme(
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      plot.caption = element_text(size = 5.5),
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(6, "pt"),
      legend.title = element_text(size = 6),
      legend.box.spacing = unit(0, "pt"),
      plot.title = element_text(size = 6),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      axis.ticks.length = unit(2, "pt")
    )
}


quartz(
  type = "pdf",
  file = here(resdir, "pp-distribution.pdf"),
  width = 3,
  height = 2
)
walk(monophyletic, possibly(function(x) print(plot_pp_dist(x)), quiet = FALSE))
dev.off()


###############################################################################
# ancestral states on the consensus phylogeny #################################
###############################################################################
tr_consensus <- asr$tr[[which(asr$treetype == "consensus")]]

add_cladelab <- function(.taxon, ln.offset, lab.offset, ...) {
  stopifnot(length(get_species_in_taxon(.taxon, taxa)) > 1)
  arc.cladelabels(
    text = .taxon,
    node = get_node(tr_consensus, .taxon, taxa),
    stretch = 1,
    cex = 0.5,
    mark.node = FALSE,
    ln.offset = ln.offset,
    lab.offset = lab.offset,
    ...
  )
}

quartz(
  type = "pdf",
  file = here(resdir, "consensus_asr.pdf"),
  width = 15,
  height = 9,
  pointsize = 14
)
par(oma = c(0, 1.5, 0, 1), xpd = NA)
plot(
  simmap_summary_consensus,
  cex = c(0.2, 0.1),
  colors = clrs,
  type = "arc",
  show.tip.label = FALSE,
  arc_height = 0.25,
  lwd = 1,
  fsize = 0.00001,
  xpd = NA
)
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
legend(
  x = min(pp$x.lim),
  y = max(pp$y.lim),
  xjust = 0,
  yjust = 1,
  levels(tipdata),
  pch = 21,
  pt.cex = 1,
  pt.bg = clrs[levels(tipdata)],
  bty = "n",
  cex = 0.6
)

add_cladelab("Metatheria", 1.02, 1.04, orientation = "horizontal")
add_cladelab("Prototheria", 1.02, 1.04, orientation = "horizontal")
# infraclass
add_cladelab("Laurasiatheria", 1.1, 1.12)
add_cladelab("Euarchontoglires", 1.1, 1.12)
add_cladelab("Afrotheria", 1.02, 1.04)
add_cladelab("Xenarthra", 1.02, 1.04, orientation = "horizontal")
# orders
add_cladelab("Primates", 1.06, 1.08)
add_cladelab("Artiodactyla", 1.06, 1.08)
add_cladelab("Rodentia", 1.06, 1.08)
add_cladelab("Carnivora", 1.06, 1.08)
add_cladelab("Chiroptera", 1.06, 1.08)
add_cladelab("Lagomorpha", 1.06, 1.08)
arc.cladelabels(
  tr_consensus,
  text = "Pinnipedia",
  node = getMRCA(
    tr_consensus,
    c(
      get_species_in_taxon("Odobenidae", taxa),
      get_species_in_taxon("Otariidae", taxa),
      get_species_in_taxon("Phocidae", taxa)
    )
  ),
  ln.offset = 1.02,
  lab.offset = 1.04,
  stretch = 1,
  cex = 0.5,
  mark.node = FALSE
)
arc.cladelabels(
  tr_consensus,
  text = "Eulipotyphla",
  node = getMRCA(
    tr_consensus,
    c(
      get_species_in_taxon("Soricidae", taxa),
      get_species_in_taxon("Erinaceidae", taxa)
      # Solenodontidae, Talpidae not in data
    )
  ),
  ln.offset = 1.06,
  lab.offset = 1.08,
  stretch = 1,
  cex = 0.5,
  mark.node = FALSE
)
arc.cladelabels(
  tr_consensus,
  text = "Herpestoidea",
  node = getMRCA(
    tr_consensus,
    c(
      get_species_in_taxon("Eupleridae", taxa),
      get_species_in_taxon("Herpestidae", taxa),
      get_species_in_taxon("Hyaenidae", taxa)
    )
  ),
  ln.offset = 1.02,
  lab.offset = 1.04,
  stretch = 1,
  cex = 0.5,
  mark.node = FALSE
)

# families
add_cladelab("Felidae", 1.02, 1.04)
add_cladelab("Canidae", 1.02, 1.04)
add_cladelab("Ursidae", 1.02, 1.04)
add_cladelab("Mustelidae", 1.02, 1.04)
add_cladelab("Vespertilionidae", 1.02, 1.04)
add_cladelab("Muridae", 1.02, 1.04)
add_cladelab("Cricetidae", 1.02, 1.04)
add_cladelab("Sciuridae", 1.02, 1.04)
arc.cladelabels(
  tr_consensus,
  text = "Caviomorpha",
  node = getMRCA(
    tr_consensus,
    c(
      get_species_in_taxon("Erethizontidae", taxa),
      get_species_in_taxon("Caviidae", taxa),
      get_species_in_taxon("Cuniculidae", taxa),
      get_species_in_taxon("Echimyidae", taxa),
      get_species_in_taxon("Octodontidae", taxa),
      get_species_in_taxon("Chinchillidae", taxa),
      get_species_in_taxon("Capromyidae", taxa)
      # Dasyproctidae, Dinomyidae, Ctenomyidae, Abrocomidae, Myocastoridae not
      # in data
    )
  ),
  ln.offset = 1.02,
  lab.offset = 1.04,
  stretch = 1,
  cex = 0.5,
  mark.node = FALSE
)
dev.off()

# end =========================================================================
