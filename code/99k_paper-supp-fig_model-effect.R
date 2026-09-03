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

conflicts_prefer(dplyr::filter, purrr::map)

resdir <- here("results/scm")
figdir <- here("results/99_paper-figs")

source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################

main <- read_csv(here(resdir, "main/posterior-prob.csv")) |>
  filter(treetype == "sample" & clade == "Eutheria") |>
  mutate(
    fit = map(
      treename,
      ~ read_rds(here(resdir, "main/sample", .x, "model-fit.rds"))
    )
  ) |>
  mutate(weight_ER = map_dbl(fit, ~ .x["fit_er", "weight"]))

case78 <- read_csv(here(resdir, "case78/posterior-prob.csv")) |>
  filter(treetype == "sample" & clade == "Eutheria") |>
  mutate(
    fit = map(
      treename,
      ~ read_rds(here(resdir, "case78/sample", .x, "model-fit.rds"))
    )
  ) |>
  mutate(weight_ER = map_dbl(fit, ~ .x["fit_er", "weight"]))

###############################################################################
# plot ========================================================================
###############################################################################

plot_corr <- function(.dataset, .title) {
  .dataset |>
    ggplot(aes(weight_ER, altricial)) +
    geom_point(size = 1, alpha = 0.5) +
    scale_x_continuous(name = "AIC weight of ER model", n.breaks = 3) +
    scale_y_continuous(name = "PP(Altricial) at Eutheria node") +
    labs(title = .title) +
    theme_classic(base_family = "Source Sans 3", base_line_size = 0.25) +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5),
      aspect.ratio = 1,
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7),
      plot.tag = element_text(size = 8, face = "bold")
    )
}

p_er <- wrap_plots(
  plot_corr(.dataset = case78, "Case (1978)"),
  plot_corr(.dataset = main, "Full dataset")
) +
  plot_annotation(tag_levels = "a")

ggsave(
  here(figdir, "supp-fig_model-effect_ER-correlation.png"),
  p_er,
  width = 4,
  height = 2,
  dpi = 600
)
# end =========================================================================
