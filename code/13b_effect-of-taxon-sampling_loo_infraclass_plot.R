# Arun Chavan
# Started: 2024-11-26

# setup =======================================================================
library(phytools)
library(tidyverse)
library(here)
library(cowplot)

resdir <- here("results/loo_infraclass/")

clrs <- c(altricial = '#fc8d59',
          intermediate = '#ffffbf',
          precocial = '#91bfdb')

# tip data ====================================================================
prec_data <- read_csv(here("data/03_coded/case78-plus-pantheria.csv"))

prec_data <- prec_data %>% 
  mutate(precocity = fct(precocity, 
                         c("altricial", "intermediate", "precocial")))

taxa <- prec_data %>% 
  select(-precocity) %>% 
  distinct()

# trees =======================================================================
tr <- read.nexus(here("data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"))

# simmap output ===============================================================
df <- read_rds(here(resdir, "output_minimal.rds"))

# add taxonomic information
infraclass_txinfo <- prec_data %>% 
  select(rank01, rank02, rank03, rank04, rank05) %>% 
  distinct() 

df <- right_join(infraclass_txinfo, df, by = c("rank05" = "leave_out"))

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

get_node <- function(.phy, .taxon) {
  spp_list <- get_species_in_taxon(.taxon)
  if (!all(spp_list %in% .phy$tip.label)) {
    warning("Some species not found in the tree. Make sure you know what you are doing!")
    spp_list <- spp_list[spp_list %in% .phy$tip.label]
  }
  stopifnot(length(spp_list) > 0)
  getMRCA(.phy, spp_list)
}

get_pp_at_node <- function(.ace, .node) .ace[which(rownames(.ace) == .node), ]

# plot ========================================================================
ace_pp <- df %>% 
  mutate(eutheria_node = map_int(tr_pruned, ~ get_node(.x, "Eutheria"))) %>% 
  mutate(eutheria_ace = map2(simmap_ace, eutheria_node, get_pp_at_node))

ace_pp <- ace_pp %>% 
  select(rank01, rank02, rank03, rank04, rank05, eutheria_ace) %>%
  unnest_longer(eutheria_ace, values_to = "pp", indices_to = "precocity") 

ace_pp_plot <- ggplot(ace_pp, aes(pp, rank05, fill = precocity)) +
  geom_col(color = "black", linewidth = 0.25) +
  scale_fill_manual(values = clrs) +
  theme_half_open(font_family = "Source Sans Pro") +
  labs(x = "PP at node Eutheria", y = "Infraclass left out")
  
ggsave(here(resdir, "eutheria_ace_pp.pdf"), 
       ace_pp_plot,
       width = 6.5, height = 4, device = cairo_pdf)

# end =========================================================================
