# Arun Chavan
# Started: 2024-11-25

# background ==================================================================

# Process the simmap output and plot the results.

###############################################################################
# setup =======================================================================
###############################################################################

library(phytools)
library(tidyverse)
library(here)
library(cowplot)
library(fs)
library(patchwork)

clrs <- c(
  altricial = '#fc8d59',
  intermediate = '#ffffbf',
  precocial = '#91bfdb'
)

analysis_name <- "main"
resdir <- here("results/scm", analysis_name)
figdir <- here("results/99_paper-figs")

source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################
# tipdata
prec_data <- read_csv(here("data/03_coded/precocity.csv"))
prec_data <- prec_data %>%
  mutate(
    precocity = fct(precocity, c("altricial", "intermediate", "precocial"))
  )
prec_tipdata <- set_names(prec_data$precocity, prec_data$binomial)

# consensus results
tr_consensus <- read_rds(here(resdir, "consensus/tree_pruned.rds"))
simmap_summary_consensus <- read_rds(here(
  resdir,
  "consensus/simmap_summary.rds"
))

# results from sampled trees
asr <- read_csv(here(resdir, "posterior-prob.csv")) |>
  filter(treetype == "sample")

# taxonomic information =======================================================
taxa <- prec_data |>
  select(all_of(c(paste("rank0", 1:8, sep = ""), "family", "binomial")))

main_nodes <- c("Mammalia", "Eutheria", "Theria", "Prototheria", "Metatheria")

###############################################################################
# posterior probability distributions #########################################
###############################################################################

pp_df_main_nodes <- asr |>
  filter(clade %in% main_nodes) |>
  mutate(clade = fct(clade, main_nodes)) |>
  pivot_longer(
    cols = levels(prec_data$precocity),
    names_to = "precocity",
    values_to = "pp"
  )

ppdist_main_nodes_faceted <- ggplot(
  pp_df_main_nodes,
  aes(treeindex, pp, fill = precocity)
) +
  geom_col(width = 0.9, linewidth = 0) +
  facet_grid(cols = vars(clade)) +
  scale_fill_manual(values = clrs, guide = "none") +
  labs(y = "Posterior probability", x = "Sampled tree") +
  theme_half_open(
    font_family = "Source Sans 3",
    line_size = 0.25,
    font_size = 6,
    rel_small = 5.5 / 6,
    rel_tiny = 5 / 6,
    rel_large = 7 / 6
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    plot.tag = element_text(face = 2, family = "Source Sans 3", size = 8),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6)
  )

###############################################################################
# ancestral states on the consensus phylogeny #################################
###############################################################################

add_cladelab <- function(.taxon, ln.offset, lab.offset, text = NULL, ...) {
  stopifnot(length(get_species_in_taxon(.taxon, taxa)) > 1)
  if (is.null(text)) {
    text <- .taxon
  }
  arc.cladelabels(
    text = text,
    node = get_node(tr_consensus, .taxon, taxa),
    stretch = 1,
    cex = 1,
    mark.node = FALSE,
    ln.offset = ln.offset,
    lab.offset = lab.offset,
    ...
  )
}

plot_asr_results <- function() {
  plot(
    simmap_summary_consensus,
    cex = c(0.2, 0.1),
    colors = clrs,
    type = "arc",
    show.tip.label = FALSE,
    arc_height = 0.05,
    lwd = 1,
    part = 0.5,
    fsize = 0.00001,
    xpd = NA
  )
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  legend(
    x = min(pp$x.lim),
    y = max(pp$y.lim),
    xjust = 0,
    yjust = 1,
    str_to_sentence(levels(prec_tipdata)),
    pch = 21,
    pt.cex = 2.5,
    pt.bg = clrs,
    bty = "n",
    cex = 1.2
  )

  nodelabels(
    node = get_node(tr_consensus, "Eutheria", taxa),
    text = "Eutheria",
    frame = "none",
    cex = 1.3,
    adj = c(-0.2, 0.5)
  )

  add_cladelab("Metatheria", 1.02, 1.04, orientation = "horizontal")
  add_cladelab("Prototheria", 1.02, 1.04, orientation = "horizontal")
  # infraclass
  add_cladelab("Laurasiatheria", 1.1, 1.12)
  add_cladelab("Euarchontoglires", 1.1, 1.12)
  add_cladelab("Afrotheria", 1.1, 1.12, text = "Afroth.")
  add_cladelab("Xenarthra", 1.1, 1.12, text = "Xen.")
  # orders
  add_cladelab("Primates", 1.06, 1.08)
  add_cladelab("Artiodactyla", 1.06, 1.08)
  add_cladelab("Rodentia", 1.06, 1.08)
  add_cladelab("Carnivora", 1.06, 1.08)
  add_cladelab("Chiroptera", 1.06, 1.08)
  add_cladelab("Lagomorpha", 1.02, 1.04)
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
    cex = 1,
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
        # Solenodontidae and Talpidae not in data
      )
    ),
    ln.offset = 1.06,
    lab.offset = 1.08,
    stretch = 1,
    cex = 1,
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
    cex = 1,
    mark.node = FALSE
  )

  # families
  add_cladelab("Felidae", 1.02, 1.04)
  add_cladelab("Canidae", 1.02, 1.04)
  add_cladelab("Ursidae", 1.02, 1.04, text = "Urs.")
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
        get_species_in_taxon("Dasyproctidae", taxa),
        get_species_in_taxon("Dinomyidae", taxa),
        get_species_in_taxon("Ctenomyidae", taxa),
        get_species_in_taxon("Echimyidae", taxa),
        get_species_in_taxon("Octodontidae", taxa),
        get_species_in_taxon("Chinchillidae", taxa),
        get_species_in_taxon("Capromyidae", taxa),
        get_species_in_taxon("Myocastoridae", taxa)
        # Abrocomidae not in data
      )
    ),
    ln.offset = 1.02,
    lab.offset = 1.04,
    stretch = 1,
    cex = 1,
    mark.node = FALSE
  )
}

quartz(
  type = "pdf",
  file = here(figdir, "fig_scm_a_consensus-asr.pdf"),
  width = 13.5,
  height = 7.5,
  family = "Source Sans 3",
  pointsize = 11
)
par(omi = c(0, 0.5, 0, 0.4), mai = c(0, 0, 0, 0), xpd = NA)
plot_asr_results()
dev.off()

png(
  here(figdir, "fig_scm_a_consensus-asr.png"),
  width = 13.5,
  height = 7.5,
  units = "in",
  res = 600,
  type = "cairo",
  family = "Source Sans 3",
  pointsize = 11
)
par(omi = c(0, 0.5, 0, 0.4), mai = c(0, 0, 0, 0), xpd = NA)
plot_asr_results()
dev.off()

###############################################################################
# put everything together #####################################################
###############################################################################

asr_png <- png::readPNG(
  here(figdir, "fig_scm_a_consensus-asr.png"),
  native = TRUE
)

composed <- wrap_plots(
  wrap_elements(full = asr_png),
  ppdist_main_nodes_faceted
) +
  plot_layout(ncol = 1, heights = c(4, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = 2, family = "Source Sans 3", size = 8),
    plot.title = element_text(hjust = 0.5)
  )

quartz(
  type = "pdf",
  file = here(figdir, "fig_scm_composed_v1.pdf"),
  width = 6.5,
  height = 5,
  family = "Source Sans 3"
)
print(composed)
dev.off()

quartz(
  type = "png",
  file = here(figdir, "fig_scm_composed_v1.png"),
  width = 6.5,
  height = 5,
  dpi = 600,
  family = "Source Sans 3"
)
print(composed)
dev.off()
# end =========================================================================
