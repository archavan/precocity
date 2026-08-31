# Arun Chavan
# Started: 2026-08-30

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(patchwork)
library(conflicted)

conflicts_prefer(dplyr::filter, purrr::map)

resdir <- here("results/53a_topology")
figdir <- here("results/99_paper-figs")

source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################

top <- read_tsv(here(resdir, "topology-vs-eutherian-pp.tsv"))

###############################################################################
# plot ========================================================================
###############################################################################

p_top <- top |>
  mutate(topology = str_replace(topology, "mouse_related", "Mouse-related")) |>
  mutate(
    tripartition = tripartition |> str_replace("_", "-") |> str_to_sentence()
  ) |>
  pivot_longer(
    cols = c(altricial, intermediate, precocial),
    names_to = "precocity",
    values_to = "pprob"
  ) |>
  ggplot(aes(pprob, topology)) +
  geom_violin(
    scale = "width",
    fill = "grey",
    color = "black",
    linewidth = 0.25
  ) +
  facet_grid(
    rows = vars(tripartition),
    cols = vars(precocity),
    labeller = as_labeller(str_to_sentence),
    scales = "free_y",
    space = "free_y",
    axes = "all",
    axis.labels = "margins",
    switch = "y"
  ) +
  # facet_wrap(
  #   vars(tripartition),
  #   ncol = 1,
  #   scales = "free_y",
  #   space = "free_y",
  #   axes = "all",
  #   axis.labels = "margins"
  # ) +
  scale_x_continuous(
    limits = c(0, 1),
    "Posterior Probability at Eutheria node"
  ) +
  scale_y_discrete(position = "right", name = "Topology") +
  theme_classic(
    ink = "black",
    base_family = "Source Sans 3",
    base_line_size = 0.25
  ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 7),
    strip.text.y.left = element_text(size = 7, angle = 0, hjust = 1),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7)
  )

p_top

ggsave(
  here(figdir, "supp-fig_topology.png"),
  p_top,
  width = 6,
  height = 2.2,
  units = "in",
  dpi = 600
)

# end =========================================================================
