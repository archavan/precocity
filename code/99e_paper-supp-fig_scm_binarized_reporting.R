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

analysis_name <- "binarized"
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
    cols = c(altricial, precocial),
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
# delta AIC ===================================================================
###############################################################################

prettify_mw_anova <- function(x) {
  x |>
    rownames_to_column("model") |>
    as_tibble() |>
    janitor::clean_names() |>
    mutate(model = str_remove(model, "fit_")) |>
    mutate(model = str_to_upper(model))
}

mw_sample <- asr |>
  filter(treetype == "sample") |>
  mutate(fit = map(fit, prettify_mw_anova)) |>
  unnest(fit)

p_delta_aic_sample <- mw_sample |>
  mutate(delta_aic = aic - min(aic), .by = treeindex) |>
  ggplot(aes(treeindex, delta_aic, color = model)) +
  geom_segment(
    aes(x = treeindex, xend = treeindex, y = 0, yend = delta_aic),
    linewidth = 0.25
  ) +
  geom_point(size = 0.5) +
  facet_grid(cols = vars(model)) +
  guides(
    color = guide_legend(
      title = "Model",
      override.aes = list(size = 2),
      position = "bottom"
    )
  ) +
  labs(x = "Sampled tree", y = "∆ AIC", tag = "b") +
  theme_classic(base_family = "Source Sans 3", base_line_size = 0.25) +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "grey"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    strip.text = element_text(size = 8),
    strip.background = element_blank()
  )

###############################################################################
# model weights ===============================================================
###############################################################################

m_weights <- mw_sample |>
  ggplot(aes(treeindex, weight, fill = model)) +
  geom_col(color = "white", linewidth = 0.1) +
  scale_fill_discrete(name = "Model") +
  guides(fill = guide_legend(position = "bottom")) +
  labs(x = "Sampled tree", y = "Model AIC weight", tag = "c") +
  theme_classic(base_line_size = 0.25, base_family = "Source Sans 3") +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.tag = element_text(size = 8, face = "bold"),
    legend.key.size = unit(0.75, "line")
  )

###############################################################################
# distribution of fitted rates ================================================
###############################################################################

p_dist_rates <- rates |>
  filter(treetype == "sample") |>
  ggplot(aes(rate, transition, fill = model)) +
  geom_violin(
    scale = "width",
    position = position_dodge(width = 0.7),
    linewidth = 0.25,
    color = "black"
  ) +
  facet_grid(rows = vars(model), scales = "free") +
  scale_fill_discrete(guide = "none") +
  scale_y_discrete(name = "Transition") +
  scale_x_continuous(name = "Fitted Rate") +
  labs(tag = "d") +
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
  free(wrap_plots(pp_focal, plot_spacer(), widths = c(2, 1))),
  free(
    wrap_plots(p_delta_aic_sample, m_weights, widths = c(2.6, 1)) &
      theme(
        legend.box.spacing = unit(0, "pt"),
        legend.key.spacing = unit(2, "pt")
      )
  ),
  wrap_plots(
    p_dist_rates,
    plot_spacer(),
    widths = c(0.8, 1)
  ),
  nrow = 3,
  heights = c(1.1, 1, 0.5)
)

composed
# write plots =================================================================

ggsave(
  here(figdir, "supp-fig_scm_binarized_reporting.png"),
  composed,
  width = 6.2,
  height = 6.65,
  units = "in",
  dpi = 600
)

# end =========================================================================
