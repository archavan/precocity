# Arun Chavan
# Started: 2024-12-07

# background ==================================================================

# setup =======================================================================
library(phytools)
library(tidyverse)
library(here)
library(patchwork)

analysis_name <- "case78-plus-pantheria"
resdir <- here("results/scm", analysis_name)
figdir <- here("results/99_paper-figs")

# functions ===================================================================
prettify_mw_anova <- function(x) {
  x %>% 
    rownames_to_column("model") %>% 
    as_tibble() %>% 
    janitor::clean_names() %>% 
    mutate(model = str_remove(model, "fit_")) %>% 
    mutate(model = str_to_upper(model))
}

# data ========================================================================
# consensus results
mw_consensus <- read_rds(here(resdir, "consensus/model-weights.rds"))

# results from sampled trees
asr <- tibble(tr_id = 1:100, tr_name = dir(here(resdir, "sample")))
asr$model_weights <- lapply(asr$tr_name, \(x) read_rds(here(resdir, "sample", x, "model-weights.rds")))

mw_sample <- asr %>% 
  mutate(model_weights = purrr::map(model_weights, prettify_mw_anova)) %>% 
  unnest(model_weights)

# plot AIC values for sampled trees ===========================================
p_delta_aic_sample <- mw_sample %>% 
  mutate(delta_aic = aic - min(aic), .by = tr_id) %>% 
  ggplot(aes(tr_id, model, fill = delta_aic)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE) +
  labs(x = "Sampled tree") +
  guides(fill = guide_colorbar(direction = "horizontal", 
                               position = "top",
                               title = "∆ AIC")) +
  theme_minimal() +
  theme(legend.text = element_text(angle = 270, hjust = 0, vjust = 0.5))

p_model_weights_sample <- mw_sample %>% 
  ggplot(aes(tr_id, model, fill = weight)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE) +
  labs(x = "Sampled tree") +
  guides(fill = guide_colourbar(direction = "horizontal", 
                               position = "top", 
                               title = "Model Weight")) +
  theme_minimal() +
  theme(legend.text = element_text(angle = 270, hjust = 0, vjust = 0.5))
  

aic_mw_composed <- wrap_plots(p_delta_aic_sample, p_model_weights_sample, ncol = 1) +
  plot_annotation(tag_levels = "a") &
  theme(
    text = element_text(family = "Source Sans Pro", size = 6.5),
    plot.tag = element_text(size = 8, face = "bold"),
    legend.key.height = unit(5, "pt"),
    legend.box.spacing = unit(0, "pt"),
    legend.box.margin = margin(0, 0, 0, 0)
  )

ggsave(filename = here(figdir, "supp-fig_scm_aic-and-mw_case-plus-pantheria.png"),
       plot = aic_mw_composed, 
       width = 6.5, height = 2.5, 
       device = ragg::agg_png,
       units = "in", dpi = 600)

# table for consensus tree ====================================================
mw_consensus %>% 
  prettify_mw_anova() %>% 
  rename(`Model` = model,
         `Log Lik.` = log_l,
         `D.F.` = d_f,
         AIC = aic,
         Weight = weight) %>% 
  write_rds(here(figdir, "supp-tbl_mw_case-plus-pantheria_consensus.rds"))

# end =========================================================================
