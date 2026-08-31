# Arun Chavan
# Started: 2026-08-30

###############################################################################
# setup =======================================================================
###############################################################################

library(phytools)
library(tidyverse)
library(here)
library(patchwork)
library(conflicted)
library(glue)

conflicts_prefer(dplyr::filter, purrr::map)

clrs <- c(
  altricial = '#fc8d59',
  intermediate = '#ffffbf',
  precocial = '#91bfdb'
)

analysis_name <- "threshold-grid"
resdir <- here("results/scm", analysis_name)
figdir <- here("results/99_paper-figs")

source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################

pp <- read_csv(here(resdir, "pp-distribution_focal-nodes.csv"))

###############################################################################
# grid ========================================================================
###############################################################################

pp_long <- pp |>
  filter(clade == "Eutheria") |>
  pivot_longer(
    cols = all_of(c("altricial", "intermediate", "precocial")),
    names_to = "precocity",
    values_to = "pp_eutheria"
  ) |>
  separate_wider_delim(taxonset, delim = "_", names = c("eye", "litter")) |>
  mutate(eye = str_remove(eye, "eye")) |>
  mutate(litter = str_remove(litter, "litter"))

p_grid <- ggplot(pp_long, aes(pp_eutheria, fill = precocity)) +
  geom_histogram(
    position = "dodge",
    binwidth = 0.1,
    boundary = 0.5,
    color = "black",
    linewidth = 0.15,
  ) +
  facet_grid(
    rows = vars(eye),
    cols = vars(litter),
    scales = "free",
    labeller = labeller(
      litter = function(x) {
        glue("Litter size\n(1.5, {x})")
      },
      eye = function(x) {
        glue("Eyes open (days)\n(0, {x})")
      }
    )
  ) +
  scale_fill_manual(
    values = clrs,
    breaks = names(clrs),
    labels = str_to_sentence,
    name = NULL
  ) +
  scale_x_continuous("PP at ancestral node for Eutheria") +
  scale_y_continuous("Count") +
  theme_classic(
    ink = "black",
    base_family = "Source Sans 3",
    base_line_size = 0.25
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(8, "pt")
  )


p_grid
# write plots =================================================================

ggsave(
  here(figdir, "supp-fig_scm_threshold-grid.png"),
  p_grid,
  width = 4.5,
  height = 3.25,
  units = "in",
  dpi = 600
)

# end =========================================================================
