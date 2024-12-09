# Arun Chavan
# Started: 2024-12-07

# setup =======================================================================
library(phytools)
library(tidyverse)
library(here)
library(patchwork)

analysis_name <- "case78-plus-pantheria"
resdir <- here("results/scm", analysis_name)
figdir <- here("results/99_paper-figs")

# data ========================================================================
prec_data <- read_csv(here("data/03_coded/case78-plus-pantheria.csv"))

# consensus results
mw_consensus <- read_rds(here(resdir, "consensus/model-weights.rds"))

# results from sampled trees
asr <- tibble(tr_id = 1:100, tr_name = dir(here(resdir, "sample")))
asr$tr_pruned <- lapply(asr$tr_name, \(x) read_rds(here(resdir, "sample", x, "tree_pruned.rds")))
asr$model_weights <- lapply(asr$tr_name, \(x) read_rds(here(resdir, "sample", x, "model-weights.rds")))
asr$ace <- lapply(asr$tr_name, \(x) read_rds(here(resdir, "sample", x, "ace.rds")))

# taxonomic information =======================================================
taxa <- prec_data %>% 
  select(rank01, rank02, rank03, rank04, rank05, rank06, rank07, rank08,
         family, binomial)

# helper functions ############################################################
prettify_mw_anova <- function(x) {
  x %>% 
    rownames_to_column("model") %>% 
    as_tibble() %>% 
    janitor::clean_names() %>% 
    mutate(model = str_remove(model, "fit_")) %>% 
    mutate(model = str_to_upper(model))
}


get_taxonomic_rank <- function(.taxon) {
  tax_rank <- names(which(apply(taxa, 2, \(x) any(grepl(.taxon, x)))))
  if (length(tax_rank) == 0) stop("Taxon not found in data")
  message(paste0(.taxon, " was found in the taxonomic rank: ", tax_rank))
  return(tax_rank)
}

get_species_in_taxon <- function(.taxon) {
  tax_rank <- get_taxonomic_rank(.taxon)
  taxa$binomial[which(taxa[[tax_rank]] == .taxon)]
}

get_classification <- function(.taxon) {
  ranks <- c("rank03", "rank05", "rank06", "rank07", "rank08", "family")
  tax_rank <- get_taxonomic_rank(.taxon)
  col_index <- which(names(taxa) == tax_rank)
  
  # if the rank is low, we want to start from infraclass to keep things short
  if (tax_rank %in% c("rank01", "rank02", "rank03")) {
    shortpath <- tax_rank
  } else {
    shortpath <- ranks[seq_len(which(ranks == tax_rank))]
  }
  
  taxonomy <- taxa %>% 
    select(all_of(shortpath)) %>% 
    filter(.data[[tax_rank]] == .taxon) %>% 
    distinct()
  
  stopifnot(nrow(taxonomy) == 1)
  
  x <- unlist(taxonomy)
  x <- x[!is.na(x)]
  return(x)
}

get_node <- function(.phy, .taxon) {
  spp_list <- get_species_in_taxon(.taxon)
  getMRCA(.phy, spp_list)
}

get_pp_at_node <- function(.ace, .node) .ace[which(rownames(.ace) == .node), ]

get_pp_df <- function(.taxon) {
  map2_df(asr$tr_pruned,
          asr$ace,
          function(x, y) {
            node <- get_node(x, .taxon)
            get_pp_at_node(y, node)
          }) %>% 
    mutate(tree_id = 1:n()) %>% 
    relocate(tree_id)
}

# delta AIC ###################################################################
mw_sample <- asr %>% 
  mutate(model_weights = purrr::map(model_weights, prettify_mw_anova)) %>% 
  unnest(model_weights)

p_delta_aic_sample <- mw_sample %>% 
  mutate(delta_aic = aic - min(aic), .by = tr_id) %>% 
  ggplot(aes(tr_id, model, fill = delta_aic)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  coord_cartesian(expand = FALSE) +
  labs(x = "Sampled tree", y = "Model") +
  guides(fill = guide_colorbar(direction = "horizontal", 
                               position = "top",
                               title = "∆ AIC")) +
  theme_classic(base_line_size = 0.25, base_family = "Source Sans Pro") +
  theme(legend.key.height = unit(5, "pt"),
        legend.box.spacing = unit(0, "pt"),
        legend.box.margin = margin(0, 0, 0, 0),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7, angle = 270, hjust = 0, vjust = 0.5),
        legend.title = element_text(size = 7))


# Effect of model choice ######################################################

m_weights <- mw_sample %>% 
  ggplot(aes(tr_id, weight, fill = model)) +
  geom_col() +
  scale_fill_discrete(name = "Model") +
  labs(x = "Sampled tree", y = "Model weight") +
  theme_classic(base_line_size = 0.25, base_family = "Source Sans Pro") +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.75, "line")
  )

model_fit_effect <- mw_sample %>% 
  filter(model == "ER") %>%
  mutate(pp_eutheria = map2(ace, tr_pruned, ~ get_pp_at_node(.x, get_node(.y, "Eutheria")))) %>% 
  unnest_wider(pp_eutheria) %>% 
  ggplot(aes(weight, altricial)) +
  geom_point(size = 1.5, color = "black", fill = "grey", shape = 21) +
  labs(x = "Model Weight for ER",
       y = "Posterior prob. of altricial\nancestral state for Eutheria") +
  theme_classic(base_family = "Source Sans Pro", base_line_size = 0.25) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8))

# combine plots ###############################################################

m_weights_composed <- p_delta_aic_sample / (m_weights | model_fit_effect) +
  plot_layout(heights = c(1, 7.5)) +
  plot_annotation(tag_levels = "a") &
  theme(
    text = element_text(family = "Source Sans Pro", size = 6.5),
    plot.tag = element_text(size = 8, face = "bold")
  )

ggsave(here(figdir, "supp-fig_scm_case78_model-weights.png"),
       m_weights_composed, 
       width = 6.5, height = 3.5, units = "in", 
       dpi = 600)











mw_sample %>% 
  filter(model == "ARD") %>%
  mutate(pp_eutheria = map2(ace, tr_pruned, ~ get_pp_at_node(.x, get_node(.y, "Eutheria")))) %>% 
  unnest_wider(pp_eutheria) %>% 
  ggplot(aes(weight, precocial)) +
  geom_point(size = 1.5, color = "black", fill = "grey", shape = 21) +
  labs(x = "Model Weight for ARD",
       y = "Posterior prob. of altricial\nancestral state for Eutheria") +
  theme_classic(base_family = "Source Sans Pro", base_line_size = 0.25) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8))


mw_sample %>% 
  mutate(delta_aic = aic - min(aic), .by = tr_id) %>% 
  ggplot(aes(tr_id, delta_aic, fill = model)) +
  geom_segment(aes(x = tr_id, xend = tr_id, y = 0, yend = delta_aic)) +
  geom_point(shape = 21) +
  facet_grid(rows = vars(model)) +
  guides(fill = guide_legend(title = "Model", override.aes = list(size = 3))) +
  theme_bw(base_family = "Source Sans Pro", base_line_size = 0.25) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        strip.background = element_blank())
  




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

ggsave(filename = here(figdir, "supp-fig_scm_case-plus-pantheria_model-weights.png"),
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
