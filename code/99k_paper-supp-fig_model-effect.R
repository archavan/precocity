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

bind_rows(
  main = select(main, treeindex, weight_ER, altricial),
  case78 = select(case78, treeindex, weight_ER, altricial),
  .id = "dataset"
) |>
  mutate(dataset = fct(dataset, c("main", "case78"))) |>
  ggplot(aes(weight_ER, altricial)) +
  geom_point(size = 1, alpha = 0.5) +
  facet_wrap(
    vars(dataset),
    scales = "free",
    labeller = as_labeller(c(case78 = "Case (1978)", main = "Full dataset"))
  ) +
  scale_x_continuous(name = "AIC weight of ER model") +
  scale_y_continuous(name = "PP(Altricial) at Eutheria node") +
  theme_classic(base_family = "Source Sans 3", base_line_size = 0.25) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    aspect.ratio = 1,
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7)
  )

ggsave(
  here(figdir, "supp-fig_model-effect_ER-correlation.png"),
  last_plot(),
  width = 4,
  height = 2,
  dpi = 600
)

# end =========================================================================
