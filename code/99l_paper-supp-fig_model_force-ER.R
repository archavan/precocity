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

analysis_name <- "main-ER"
resdir <- here("results/scm", analysis_name)
figdir <- here("results/99_paper-figs")

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
  mutate(fit = map(respath, ~ read_rds(here(resdir, .x, "model-fit.rds"))))

pp <- read_csv(here(resdir, "posterior-prob.csv"))
rates <- read_csv(here(resdir, "fitted-rates.csv"))

###############################################################################
# PP Precocial at select nodes ================================================
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

pp_focal <- pp |>
  filter(clade %in% focal_nodes) |>
  mutate(clade = fct(clade, focal_nodes)) |>
  mutate(fill = ifelse(clade == "Eutheria", "yellow", "grey")) |>
  pivot_longer(
    cols = c(altricial, intermediate, precocial),
    names_to = "precocity",
    values_to = "pprob"
  ) |>
  ggplot(aes(pprob, clade, fill = fill)) +
  geom_violin(
    scale = "width",
    color = "black",
    linewidth = 0.25
  ) +
  facet_grid(cols = vars(precocity), labeller = as_labeller(str_to_sentence)) +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0, 1),
    name = "Posterior Probability at node"
  ) +
  scale_y_discrete(name = "Node") +
  labs(tag = "a") +
  theme_classic(
    base_family = "Source Sans 3",
    ink = "black",
    base_line_size = 0.25
  ) +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    panel.grid.major = element_line(color = "grey", linewidth = 0.2),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 7),
    strip.background = element_blank()
  )

###############################################################################
# distribution of fitted rates ================================================
###############################################################################

p_dist_rates <- rates |>
  filter(treetype == "sample") |>
  ggplot(aes(rate, transition)) +
  geom_violin(
    scale = "width",
    position = position_dodge(width = 0.7),
    linewidth = 0.25,
    color = "black",
    fill = "grey",
  ) +
  facet_grid(rows = vars(model), scales = "free") +
  scale_fill_discrete(guide = "none") +
  scale_y_discrete(name = "Transition") +
  scale_x_continuous(name = "Fitted Rate") +
  labs(tag = "b") +
  theme_classic(base_family = "Source Sans 3", base_line_size = 0.25) +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    panel.grid.major = element_line(color = "grey", linewidth = 0.2),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(8, "pt"),
    strip.text = element_text(size = 7),
    strip.background = element_blank()
  )


###############################################################################
# combine plots ===============================================================
###############################################################################

composed <- wrap_plots(
  free(pp_focal),
  wrap_plots(p_dist_rates, plot_spacer(), ncol = 1, heights = c(1, 1)),
  nrow = 1,
  widths = c(3.5, 1)
)

composed
# write plots =================================================================

ggsave(
  here(figdir, "supp-fig_scm_force-ER_reporting.png"),
  composed,
  width = 6.5,
  height = 3,
  units = "in",
  dpi = 600
)

# end =========================================================================
