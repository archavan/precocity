# Arun Chavan
# Started: 2024-02-09

# background ==================================================================

# Process the simmap output and plot the results.

# setup =======================================================================
library(phytools)
library(tidyverse)
library(here)
library(cowplot)
library(fs)

clrs <- c(altricial = '#fc8d59',
          intermediate = '#ffffbf',
          precocial = '#91bfdb')

analysis_name <- "case78"
resdir <- here("results/scm", analysis_name)

# data ========================================================================
# tip data
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
# posterior probability distributions #########################################
###############################################################################

# functions ===================================================================
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

plot_pp <- function(
    .taxon, 
    .title = .taxon
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
         title = .title,
         caption = paste0(get_classification(.taxon), collapse = " → ")) +
    theme_bw(base_family = "Source Sans Pro", base_line_size = 0.25) +
    theme(
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      plot.caption = element_text(size = 5.5),
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
  fs::dir_create(here(resdir, "plots/pp", .rank))
  map(x, possibly(plot_pp)) %>% 
    iwalk(
      ~ ggsave(filename = here(resdir, 
                               "plots/pp", 
                               .rank,
                               paste0(.y, ".pdf")),
               plot = .x, width = 2.75, height = 1.75, device = cairo_pdf)
    )
}

plot_and_save_pp_for_all_taxa_in_rank("rank01")
plot_and_save_pp_for_all_taxa_in_rank("rank02")
plot_and_save_pp_for_all_taxa_in_rank("rank03")
plot_and_save_pp_for_all_taxa_in_rank("rank04")
plot_and_save_pp_for_all_taxa_in_rank("rank05")
plot_and_save_pp_for_all_taxa_in_rank("rank06")
plot_and_save_pp_for_all_taxa_in_rank("rank07")
plot_and_save_pp_for_all_taxa_in_rank("rank08")
plot_and_save_pp_for_all_taxa_in_rank("family")

###############################################################################
# ancestral states on the consensus phylogeny #################################
###############################################################################

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

cairo_pdf(here(resdir, "plots", "consensus_asr.pdf"), 
          width = 12, height = 7, 
          family = "Source Sans Pro", pointsize = 14)
par(oma = c(0, 1.5, 0, 1), xpd = NA)
plot(simmap_summary_consensus, 
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
# add_cladelab("Prototheria", 1.02, 1.04, orientation = "horizontal")
# infraclass
add_cladelab("Laurasiatheria", 1.1, 1.12)
add_cladelab("Euarchontoglires", 1.1, 1.12)
add_cladelab("Afrotheria", 1.02, 1.04)
# add_cladelab("Xenarthra", 1.02, 1.04, orientation = "horizontal")
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
# arc.cladelabels(tr_consensus, 
#                 text = "Eulipotyphla", 
#                 node = getMRCA(tr_consensus, 
#                                c(get_species_in_taxon("Soricidae"), # Solenodontidae and Talpidae not in data
#                                  get_species_in_taxon("Erinaceidae"))),
                # ln.offset = 1.06, 
                # lab.offset = 1.08, 
                # stretch = 1,
                # cex = 0.5,
                # mark.node = FALSE)
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
#add_cladelab("Felidae", 1.02, 1.04)
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

# end =========================================================================
