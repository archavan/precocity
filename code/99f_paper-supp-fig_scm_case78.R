# Arun Chavan
# Started: 2024-12-07

# background ==================================================================

# Process the simmap output and plot the results.

# setup =======================================================================
library(phytools)
library(tidyverse)
library(here)
library(cowplot)
library(fs)
library(patchwork)

clrs <- c(altricial = '#fc8d59',
          intermediate = '#ffffbf',
          precocial = '#91bfdb')

analysis_name <- "case78"
resdir <- here("results/scm", analysis_name)
figdir <- here("results/99_paper-figs")

# data ========================================================================
# tipdata
prec_data <- read_csv(here("data/03_coded/case78/precocity_case78_v1.csv"))

prec_data <- prec_data %>% 
  mutate(precocity_recoded = case_when(
    precocity == "fetal" ~ "altricial",
    precocity == "A" ~ "altricial",
    precocity == "SA" ~ "altricial",
    precocity == "SP" ~ "intermediate",
    precocity == "P" ~ "precocial"
  )) %>% 
  mutate(precocity_recoded = fct(precocity_recoded, 
                                 c("altricial", "intermediate", "precocial")))

prec_tipdata <- set_names(prec_data$precocity_recoded, prec_data$binomial)

# consensus results
tr_consensus <- read_rds(here(resdir, "consensus/tree_pruned.rds"))
simmap_summary_consensus <- read_rds(here(resdir, "consensus/simmap_summary.rds"))

# results from sampled trees
asr <- tibble(tr_id = 1:100, tr_name = dir(here(resdir, "sample")))
asr$tr_pruned <- lapply(asr$tr_name, \(x) read_rds(here(resdir, "sample", x, "tree_pruned.rds")))
asr$model_weights <- lapply(asr$tr_name, \(x) read_rds(here(resdir, "sample", x, "model-weights.rds")))
asr$ace <- lapply(asr$tr_name, \(x) read_rds(here(resdir, "sample", x, "ace.rds")))

# taxonomic information =======================================================
taxa <- prec_data %>% 
  select(rank01, rank02, rank03, rank04, rank05, rank06, rank07, rank08,
         family, binomial)

###############################################################################
# helper functions ############################################################
###############################################################################

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


###############################################################################
# posterior probability distributions #########################################
###############################################################################

pp_df_main_nodes <- list(
  Mammalia = get_pp_df("Mammalia"),
  Eutheria = get_pp_df("Eutheria"),
  Theria = get_pp_df("Theria"),
  # Prototheria = get_pp_df("Prototheria"),
  Metatheria = get_pp_df("Metatheria")
) %>% 
  bind_rows(.id = "clade") %>%
  mutate(clade = fct(clade, c("Mammalia", "Eutheria", "Theria", "Prototheria", "Metatheria"))) %>% 
  pivot_longer(cols = c(altricial, intermediate, precocial), 
               names_to = "precocity", 
               values_to = "pp")


ppdist_main_nodes_faceted <- ggplot(pp_df_main_nodes, aes(tree_id, pp, fill = precocity)) +
  geom_col(width = 0.9, linewidth = 0) +
  facet_grid(cols = vars(clade)) +
  scale_fill_manual(values = clrs, guide = "none") +
  labs(y = "Posterior probability",
       x = "Sampled tree") +
  theme_half_open(font_family = "Source Sans Pro",
                  line_size = 0.25,
                  font_size = 6,
                  rel_small = 5.5/6,
                  rel_tiny = 5/6,
                  rel_large = 7/6) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7))

###############################################################################
# ancestral states on the consensus phylogeny #################################
###############################################################################

add_cladelab <- function(.taxon, 
                         ln.offset, 
                         lab.offset, 
                         text = NULL,
                         ...) {
  stopifnot(length(get_species_in_taxon(.taxon)) > 1)
  if (is.null(text)) text <- .taxon
  arc.cladelabels(text = text, 
                  node = get_node(tr_consensus, .taxon),
                  stretch = 1,
                  cex = 1,
                  mark.node = FALSE,
                  ln.offset = ln.offset,
                  lab.offset = lab.offset, 
                  ...)
}

plot_asr_results <- function() {
  plot(simmap_summary_consensus, 
       cex = c(0.2, 0.1),
       colors = clrs,
       type = "arc",
       show.tip.label = FALSE,
       arc_height = 0.05,
       lwd = 1,
       part = 0.5,
       fsize = 0.00001, 
       xpd = NA)
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  legend(x = min(pp$x.lim), y = max(pp$y.lim), xjust = 0, yjust = 1,
         levels(prec_tipdata), 
         pch = 21, 
         pt.cex = 1.5, 
         pt.bg = clrs, 
         bty = "n", 
         cex = 1)
  
  nodelabels(node = get_node(tr_consensus, "Eutheria"), 
             text = "Eutheria", 
             frame = "none", 
             cex = 1.2,
             adj = c(-0.2, 0.5))
  
  add_cladelab("Metatheria", 1.02, 1.04, orientation = "horizontal")
#  add_cladelab("Prototheria", 1.02, 1.04, orientation = "horizontal")
  # infraclass
  add_cladelab("Laurasiatheria", 1.1, 1.12)
  add_cladelab("Euarchontoglires", 1.1, 1.12)
  add_cladelab("Afrotheria", 1.1, 1.12, text = "Afroth.")
#  add_cladelab("Xenarthra", 1.1, 1.12, text = "Xen.")
  # orders
  add_cladelab("Primates", 1.06, 1.08)
  add_cladelab("Cetartiodactyla", 1.06, 1.08)
  add_cladelab("Rodentia", 1.06, 1.08)
  add_cladelab("Carnivora", 1.06, 1.08)
  add_cladelab("Chiroptera", 1.06, 1.08)
  add_cladelab("Lagomorpha", 1.02, 1.04)
  arc.cladelabels(tr_consensus, 
                  text = "Pinnipedia", 
                  node = getMRCA(tr_consensus, 
                                 c(get_species_in_taxon("Odobenidae"), 
                                   get_species_in_taxon("Otariidae"), 
                                   get_species_in_taxon("Phocidae"))),
                  ln.offset = 1.02, 
                  lab.offset = 1.04, 
                  stretch = 1,
                  cex = 1,
                  mark.node = FALSE)
  # arc.cladelabels(tr_consensus, 
  #                 text = "Eulipotyphla", 
  #                 node = getMRCA(tr_consensus, 
  #                                c(get_species_in_taxon("Soricidae"), # Solenodontidae and Talpidae not in data
  #                                  get_species_in_taxon("Erinaceidae"))),
  #                 ln.offset = 1.06, 
  #                 lab.offset = 1.08, 
  #                 stretch = 1,
  #                 cex = 1,
  #                 mark.node = FALSE)
  # arc.cladelabels(tr_consensus, 
  #                 text = "Herpestoidea", 
  #                 node = getMRCA(tr_consensus, 
  #                                c(get_species_in_taxon("Herpestidae"),
  #                                  get_species_in_taxon("Hyaenidae"))),
  #                 ln.offset = 1.02, 
  #                 lab.offset = 1.04, 
  #                 stretch = 1,
  #                 cex = 1,
  #                 mark.node = FALSE)
  
  # families
#  add_cladelab("Felidae", 1.02, 1.04)
  add_cladelab("Canidae", 1.02, 1.04)
  add_cladelab("Ursidae", 1.02, 1.04, text = "Urs.")
  add_cladelab("Mustelidae", 1.02, 1.04)
  add_cladelab("Vespertilionidae", 1.02, 1.04)
  add_cladelab("Muridae", 1.02, 1.04)
  add_cladelab("Cricetidae", 1.02, 1.04)
  add_cladelab("Sciuridae", 1.02, 1.04)
  # arc.cladelabels(tr_consensus, 
  #                 text = "Caviomorpha", 
  #                 node = getMRCA(tr_consensus, 
  #                                c(get_species_in_taxon("Caviidae"), 
  #                                  get_species_in_taxon("Cuniculidae"), 
  #                                  get_species_in_taxon("Ctenomyidae"),
  #                                  get_species_in_taxon("Dasyproctidae"),
  #                                  get_species_in_taxon("Chinchillidae"),
  #                                  get_species_in_taxon("Erethizontidae"),
  #                                  get_species_in_taxon("Echimyidae"),
  #                                  get_species_in_taxon("Octodontidae"))),
  #                 ln.offset = 1.02, 
  #                 lab.offset = 1.04, 
  #                 stretch = 1,
  #                 cex = 1,
  #                 mark.node = FALSE)
}

# cairo_pdf(here(figdir, "supp-fig_scm-case78_a_consensus-asr.pdf"), 
#           width = 10, 
#           height = 5.5, 
#           family = "Source Sans Pro", 
#           pointsize = 10.5)
# par(omi = c(0, 0.5, 0, 0.4),
#     mai = c(0, 0, 0, 0),
#     xpd = NA)
# plot_asr_results()
# dev.off()

png(here(figdir, "supp-fig_scm_case78_a_consensus-asr.png"), 
    width = 10, 
    height = 5.5, 
    units = "in",
    res = 600,
    type = "cairo",
    family = "Source Sans Pro", 
    pointsize = 9)
par(omi = c(0, 0.5, 0, 0.4),
    mai = c(0, 0, 0, 0),
    xpd = NA)
plot_asr_results()
dev.off()

###############################################################################
# ancestral state reconstruction: put together ################################
###############################################################################

asr_png <- png::readPNG(here(figdir, "supp-fig_scm_case78_a_consensus-asr.png"), 
                        native = TRUE)

composed <- wrap_plots(wrap_elements(full = asr_png), 
                       ppdist_main_nodes_faceted) +
  plot_layout(ncol = 1, heights = c(4, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 2, 
                                family = "Source Sans Pro", 
                                size = 8),
        plot.title = element_text(hjust = 0.5))

ggsave(here(figdir, "supp-fig_scm_case78_composed_v1.png"), 
       composed, 
       width = 6.5, 
       height = 5, 
       units = "in",
       dpi = 600)

###############################################################################
# model choice and weights ####################################################
###############################################################################

# delta AIC ###################################################################
mw_sample <- asr %>% 
  mutate(model_weights = purrr::map(model_weights, prettify_mw_anova)) %>% 
  unnest(model_weights)

p_delta_aic_sample <- mw_sample %>% 
  mutate(delta_aic = aic - min(aic), .by = tr_id) %>% 
  ggplot(aes(tr_id, delta_aic, color = model)) +
  geom_segment(aes(x = tr_id, xend = tr_id, y = 0, yend = delta_aic),
               linewidth = 0.25) +
  geom_point(size = 0.5) +
  facet_grid(cols = vars(model)) +
  guides(color = guide_legend(title = "Model", 
                              override.aes = list(size = 2), 
                              position = "bottom")) +
  labs(x = "Sampled tree", y = "∆ AIC") +
  theme_bw(base_family = "Source Sans Pro", base_line_size = 0.25) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        strip.background = element_blank())

# model weighsts ##############################################################
m_weights <- mw_sample %>% 
  ggplot(aes(tr_id, weight, fill = model)) +
  geom_col() +
  scale_fill_discrete(name = "Model") +
  guides(fill = guide_legend(position = "bottom")) +
  labs(x = "Sampled tree", y = "Model weight") +
  theme_bw(base_line_size = 0.25, base_family = "Source Sans Pro") +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.75, "line")
  )

# Effect of model choice ######################################################
model_fit_effect <- mw_sample %>% 
  filter(model == "ER") %>%
  mutate(pp_eutheria = map2(ace, tr_pruned, ~ get_pp_at_node(.x, get_node(.y, "Eutheria")))) %>% 
  unnest_wider(pp_eutheria) %>% 
  ggplot(aes(weight, altricial)) +
  geom_point(size = 1, color = "black", fill = "grey", shape = 21, stroke = 0.4) +
  labs(x = "Model Weight for ER",
       y = "Posterior prob. of altricial\nancestral state for Eutheria") +
  theme_bw(base_family = "Source Sans Pro", base_line_size = 0.25) +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8))

# combine plots ###############################################################

mw_aic_patch <- wrap_plots(p_delta_aic_sample, m_weights, ncol = 2) +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    text = element_text(family = "Source Sans Pro", size = 6.5),
    plot.tag = element_text(size = 8, face = "bold"),
    legend.key.spacing = unit(2, "pt"),
    legend.box.spacing = unit(0, "pt")
  )

# write plots #################################################################

ggsave(here(figdir, "supp-fig_scm_case78_mw-aic.png"),
       mw_aic_patch, 
       width = 6.5, height = 2.4, units = "in", 
       dpi = 600)

ggsave(here(figdir, "supp-fig_scm_case78_model-choice.png"),
       model_fit_effect, 
       width = 3, height = 2.25, units = "in", 
       dpi = 600)

# so we can combine it with the same plot from case+pt
write_rds(model_fit_effect,
          here(figdir, "supp-fig_scm_case78_model-choice.rds"))

# summary stats ===============================================================
mw_sample %>% 
  mutate(delta_aic = aic - min(aic), .by = tr_id) %>% 
  group_by(model) %>% 
  summarise(mean_weight = mean(weight),
            mean_delta_aic = mean(delta_aic)) %>% 
  ungroup()


get_pp_summary_stats <- function(.taxon) {
  get_pp_df(.taxon) %>% 
    summarise(altricial_mean = mean(altricial),
              altricial_sd = sd(altricial),
              intermediate_mean = mean(intermediate),
              intermediate_sd = sd(intermediate),
              precocial_mean = mean(precocial),
              precocial_sd = sd(precocial))
}

get_pp_summary_stats("Theria")
get_pp_summary_stats("Mammalia")
get_pp_summary_stats("Eutheria")
get_pp_summary_stats("Prototheria")
get_pp_summary_stats("Metatheria")

# end =========================================================================
