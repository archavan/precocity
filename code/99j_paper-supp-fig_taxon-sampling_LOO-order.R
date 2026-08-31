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

analysis <- "loo_order"
resdir <- here("results/scm", analysis)
figdir <- here("results/99_paper-figs")

source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################

pp <- read_csv(here(resdir, "pp-distribution_focal-nodes.csv")) |>
  filter(treetype == "sample")

###############################################################################
# plot ========================================================================
###############################################################################

focal_nodes <- c(
  "Mammalia",
  "Theria",
  "Eutheria",
  "Prototheria",
  "Metatheria",
  "Xenarthra",
  "Afrotheria",
  "Boreoeutheria",
  "Laurasiatheria",
  "Euarchontoglires",
  "Rodentia",
  "Primates",
  "Carnivora",
  "Chiroptera",
  "Artiodactyla"
)


pp_eutheria <- pp |>
  filter(clade == "Eutheria") |>
  pivot_longer(
    cols = all_of(names(clrs)),
    values_to = "pp_at_node",
    names_to = "precocity"
  ) |>
  ggplot(aes(pp_at_node, taxonset, fill = precocity)) +
  geom_violin(scale = "width", linewidth = 0.2, color = "black") +
  scale_x_continuous(
    name = "Posterior Probability at Eutheria node",
    limits = c(0, 1)
  ) +
  scale_y_discrete(name = "Order Left Out") +
  scale_fill_manual(
    values = clrs,
    guide = "none"
  ) +
  facet_grid(
    cols = vars(precocity),
    labeller = labeller(precocity = str_to_sentence),
    axes = "all",
    axis.labels = "margins"
  ) +
  theme_classic(base_family = "Source Sans 3", base_line_size = 0.25) +
  theme(
    panel.grid.major.y = element_line(color = "grey"),
    strip.background = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

# write plots and files =======================================================
ggsave(
  here(figdir, "supp-fig_taxon-sampling_LOO-order.png"),
  pp_eutheria,
  width = 4.2,
  height = 3.25,
  units = "in",
  dpi = 600
)

# end #########################################################################
