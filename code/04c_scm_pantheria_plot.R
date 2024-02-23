# Arun Chavan
# Started: 2024-02-09

# background ==================================================================

# Process the simmap output and plot the results.

# setup =======================================================================
library(phytools)
library(tidyverse)
library(here)
library(cowplot)

clrs <- c(altricial = '#fc8d59',
          intermediate = '#ffffbf',
          precocial = '#91bfdb')

# data ========================================================================
# precocity data
prec_data <- read_csv(here("data/03_coded/pt_precocity-data_v1.csv"))

# trees
tr_consensus_full <- read.nexus(here("data/trees/upham2019/consensus/DNA-only/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"))

upham_taxa <- read_csv(here("results/taxa/upham2019_all-taxa.csv"))

# taxonomic information =======================================================
taxa <- prec_data %>% 
  select(infraclass, superorder, order, family, binomial) %>% 
  distinct() %>% 
  mutate(class = "Mammalia") %>% 
  mutate(subclass = ifelse(infraclass %in% c("Eutheria", "Metatheria"),
                           "Theria",
                           NA_character_)) %>% 
  select(class, subclass, infraclass, superorder, order, family, binomial)

# prepare data ================================================================
## trait data -----------------------------------------------------------------
# prepare named vector of trait data
prec_tipdata <- set_names(prec_data$precocity, prec_data$binomial)
prec_tipdata <- fct(prec_tipdata, c("altricial", "intermediate", "precocial"))

## consensus tree for visualization -------------------------------------------
# clean up tip labels in the tree. The Upham tree has additional information
# after the binomial. We will replace the tiplabels with just the binomial.
tr_consensus_full$tip.label <- upham_taxa$binomial[match(tr_consensus_full$tip.label,
                                                         upham_taxa$tiplab)]

tr_consensus <- drop.tip(
  tr_consensus_full,
  which(!(tr_consensus_full$tip.label %in% prec_data$binomial))
)
tr_consensus <- ladderize(tr_consensus)

# simmap output ===============================================================
asr <- read_rds(here("results/simmap/out/simmap_100-trees.rds"))
simmap_consensus <- read_rds(here("results/simmap/out/simmap_consensus-tree.rds"))

###############################################################################
# posterior probability distributions #########################################
###############################################################################

# functions ===================================================================
get_taxonomic_rank <- function(.taxon) {
  tax_rank <- names(which(apply(taxa, 2, function(x) any(grepl(.taxon, x)))))
  if (length(tax_rank) == 0) stop("Taxon not found in data")
  message(paste0(.taxon, " was found in the taxonomic rank: ", tax_rank))
  return(tax_rank)
}

get_species_in_taxon <- function(.taxon) {
  tax_rank <- get_taxonomic_rank(.taxon)
  taxa$binomial[which(taxa[[tax_rank]] == .taxon)]
}

get_classification <- function(.taxon) {
  ranks <- c("infraclass", "superorder", "order", "family")
  tax_rank <- get_taxonomic_rank(.taxon)
  col_index <- which(names(taxa) == tax_rank)
  
  # if the rank is low, we want to start from infraclass to keep things short
  if (tax_rank %in% c("class", "subclass", "infraclass")) {
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
  tax_rank <- get_taxonomic_rank(.taxon)
  spp_list <- taxa %>% 
    filter(.data[[tax_rank]] == .taxon) %>% 
    pull(binomial)
  getMRCA(.phy, spp_list)
}

get_pp_at_node <- function(.ace, .node) .ace[which(rownames(.ace) == .node), ]

get_pp_df <- function(.taxon) {
  map2_df(asr$trees,
          asr$simmap_ace,
          function(x, y) {
            node <- get_node(x, .taxon)
            get_pp_at_node(y, node)
          }) %>% 
    rownames_to_column("tree_id") %>% 
    mutate(tree_id = fct(tree_id, as.character(1:100)))
}

plot_pp <- function(
    .taxon, 
    .title = paste0(get_classification(.taxon), collapse = " → ")
) {
  df <- get_pp_df(.taxon) %>% 
    pivot_longer(-tree_id, names_to = "precocity", values_to = "pp")
  
  ggplot(df, aes(tree_id, pp, fill = precocity)) +
    geom_col(width = 0.9, linewidth = 0) +
    scale_fill_manual(values = clrs) +
    scale_x_discrete(breaks = c(1, 25, 50, 75, 100),
                     expand = c(0.04, 0.5)) +
    labs(y = "Posterior probability",
         x = "Sampled tree", 
         title = .title) +
    theme_bw(base_family = "Source Sans Pro", base_line_size = 0.25) +
    theme(
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(6, "pt"),
      legend.title = element_text(size = 6),
      legend.box.spacing = unit(0, "pt"),
      plot.title = element_text(size = 6),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      axis.ticks.length = unit(2, "pt")
    ) 
}


# plot ========================================================================
plot_and_save_pp_for_all_taxa_in_rank <- function(.rank) {
  x <- unique(taxa[[.rank]][!is.na(taxa[[.rank]])])
  names(x) <- x
  fs::dir_create(here("results/simmap/plots/pp", .rank))
  map(x, possibly(plot_pp)) %>% 
    iwalk(
      ~ ggsave(filename = here("results/simmap/plots/pp", 
                               .rank,
                               paste0(.y, ".pdf")),
               plot = .x, width = 2.75, height = 1.75, device = cairo_pdf)
    )
}

plot_and_save_pp_for_all_taxa_in_rank("class")
plot_and_save_pp_for_all_taxa_in_rank("subclass")
plot_and_save_pp_for_all_taxa_in_rank("infraclass")
plot_and_save_pp_for_all_taxa_in_rank("superorder")
plot_and_save_pp_for_all_taxa_in_rank("order")
plot_and_save_pp_for_all_taxa_in_rank("family")

###############################################################################
# ancestral states on the consensus phylogeny #################################
###############################################################################

simmap_consensus_summary <- summary(simmap_consensus)

add_cladelab <- function(.taxon, 
                         ln.offset, 
                         lab.offset, 
                         ...) {
  stopifnot(length(get_species_in_taxon(.taxon)) > 1)
  arc.cladelabels(text = .taxon, 
                  node = get_node(tr_consensus, .taxon),
                  stretch = 1,
                  cex = 0.5,
                  mark.node = FALSE,
                  ln.offset = ln.offset,
                  lab.offset = lab.offset, 
                  ...)
}

cairo_pdf(here("results/simmap/plots/consensus_asr.pdf"), 
          width = 12, height = 7, 
          family = "Source Sans Pro", pointsize = 14)
par(oma = c(0, 1.5, 0, 1), xpd = NA)
plot(simmap_consensus_summary, 
     cex = c(0.2, 0.1),
     colors = clrs,
     type = "arc",
     show.tip.label = FALSE,
     arc_height = 0.25,
     lwd = 1,
     fsize = 0.00001, xpd = NA)
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
legend(x = min(pp$x.lim), y = max(pp$y.lim), xjust = 0, yjust = 1,
       levels(prec_tipdata), 
       pch = 21, pt.cex = 1, 
       pt.bg = clrs, 
       bty = "n", cex = 0.6)

add_cladelab("Metatheria", 1.02, 1.04, orientation = "horizontal")
add_cladelab("Prototheria", 1.02, 1.04, orientation = "horizontal")
# infraclass
add_cladelab("Laurasiatheria", 1.1, 1.12)
add_cladelab("Euarchontoglires", 1.1, 1.12)
add_cladelab("Afrotheria", 1.02, 1.04)
add_cladelab("Xenarthra", 1.02, 1.04, orientation = "horizontal")
# orders
add_cladelab("Rodentia", 1.06, 1.08)
add_cladelab("Carnivora", 1.06, 1.08)
add_cladelab("Chiroptera", 1.06, 1.08)
add_cladelab("Lagomorpha", 1.06, 1.08)
arc.cladelabels(tr_consensus, 
                text = "Pinnipedia", 
                node = getMRCA(tr_consensus, 
                               c(get_species_in_taxon("Odobenidae"), 
                                 get_species_in_taxon("Otariidae"), 
                                 get_species_in_taxon("Phocidae"))),
                ln.offset = 1.02, 
                lab.offset = 1.04, 
                stretch = 1,
                cex = 0.5,
                mark.node = FALSE)
arc.cladelabels(tr_consensus, 
                text = "Eulipotyphla", 
                node = getMRCA(tr_consensus, 
                               c(get_species_in_taxon("Soricidae"), # Solenodontidae and Talpidae not in data
                                 get_species_in_taxon("Erinaceidae"))),
                ln.offset = 1.06, 
                lab.offset = 1.08, 
                stretch = 1,
                cex = 0.5,
                mark.node = FALSE)
arc.cladelabels(tr_consensus, 
                text = "Herpestoidea", 
                node = getMRCA(tr_consensus, 
                               c(get_species_in_taxon("Herpestidae"),
                                 get_species_in_taxon("Hyaenidae"))),
                ln.offset = 1.02, 
                lab.offset = 1.04, 
                stretch = 1,
                cex = 0.5,
                mark.node = FALSE)

# families
add_cladelab("Felidae", 1.02, 1.04)
add_cladelab("Canidae", 1.02, 1.04)
add_cladelab("Ursidae", 1.02, 1.04)
add_cladelab("Mustelidae", 1.02, 1.04)
add_cladelab("Vespertilionidae", 1.02, 1.04)
add_cladelab("Muridae", 1.02, 1.04)
add_cladelab("Cricetidae", 1.02, 1.04)
add_cladelab("Sciuridae", 1.02, 1.04)
arc.cladelabels(tr_consensus, 
                text = "Caviomorpha", 
                node = getMRCA(tr_consensus, 
                               c(get_species_in_taxon("Caviidae"), 
                                 get_species_in_taxon("Cuniculidae"), 
                                 get_species_in_taxon("Ctenomyidae"),
                                 get_species_in_taxon("Dasyproctidae"),
                                 get_species_in_taxon("Chinchillidae"),
                                 get_species_in_taxon("Erethizontidae"),
                                 get_species_in_taxon("Echimyidae"),
                                 get_species_in_taxon("Octodontidae"))),
                ln.offset = 1.02, 
                lab.offset = 1.04, 
                stretch = 1,
                cex = 0.5,
                mark.node = FALSE)
dev.off()

###############################################################################
# rerun model selection on the consensus tree #################################
###############################################################################

# I already did this, in the last script along with all the sampled trees but I
# didn't save the intermediate steps. For the consensus tree, I would like to
# show the model AIC values and model weights. Rerunning that part just to get
# this information.

run_model_comparison <- function(.phy, .dat) {
  message("    checking data...")
  stopifnot(all(.phy$tip.label %in% names(.dat)))
  
  # sort trait data
  .dat_sorted <- .dat[.phy$tip.label]
  
  # fit models
  message("    fitting models...")
  fit_er  <- fitMk(.phy, .dat_sorted, model = "ER",  pi = "estimated")
  fit_sym <- fitMk(.phy, .dat_sorted, model = "SYM", pi = "estimated")
  fit_ard <- fitMk(.phy, .dat_sorted, model = "ARD", pi = "estimated")
  
  # model comparison
  message("    running model comparison...")
  fit_aov <- anova(fit_er, fit_sym, fit_ard)

  return(fit_aov)
}

models_consensus <- run_model_comparison(tr_consensus, prec_tipdata)
rownames_to_column(models_consensus, "model") %>% 
  write_tsv(here("results/simmap/out/simmap_consensus_model-weights.tsv"))

# end =========================================================================
