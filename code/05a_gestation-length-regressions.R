# Arun Chavan
# Started: 2024-02-10

# background ==================================================================

# PGLS of gestation length against body mass and precocity

# setup =======================================================================
library(ape)
library(nlme)
library(tidyverse)
library(here)
library(patchwork)

resdir <- here("results/pgls_gestation-len")
fs::dir_create(resdir)

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

# subset to eutheria ==========================================================
# We will be doing regression for gestation length only for eutherians. 
prec_data_euth <- filter(prec_data, infraclass == "Eutheria")

###############################################################################
# PGLS using the consensus tree for visualization #############################
###############################################################################

# model with body mass and precocity ==========================================

# subset to eutheria
tr_consensus_euth <- drop.tip(
  tr_consensus,
  which(!(tr_consensus$tip.label %in% prec_data_euth$binomial))
)

# match order of species in tree and data
prec_data_euth <- prec_data_euth[match(tr_consensus_euth$tip.label, prec_data_euth$binomial), ]
waldo::compare(tr_consensus_euth$tip.label, prec_data_euth$binomial)

spp <- prec_data_euth$binomial

# initiate correlation structure
cor_lambda_consensus <- corPagel(value = 1, tr_consensus_euth, form = ~spp)

# model
pgls_consensus <- gls(log(gestation_len_d) ~ log(adult_body_mass_g) + precocity,
                      data = prec_data_euth, 
                      correlation = cor_lambda_consensus)

anova(pgls_consensus)
summary(pgls_consensus)

# residuals of body mass regression agaist precocity ==========================
pgls_consensus_bm <- gls(log(gestation_len_d) ~ log(adult_body_mass_g),
                         data = prec_data_euth, 
                         correlation = cor_lambda_consensus)

prec_vs_bm_residuals <- prec_data_euth %>% 
  mutate(bm_residuals = pgls_consensus_bm$residuals) %>% 
  ggplot(aes(precocity, bm_residuals, )) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey25") +
  ggforce::geom_sina(aes(fill = precocity), shape = 21, size = 1, stroke = 0.15) +
  scale_fill_manual(values = clrs, guide = "none") +
  labs(y = "residuals") +
  theme_bw(base_family = "Source Sans Pro") +
  theme(panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank())

# plot for visualization ======================================================
# to limit abline to range of data

calc_segment_coords <- function(.xmin, .xmax, .intercept, .slope) {
  y    <- unname(.intercept + (.slope * .xmin))
  yend <- unname(.intercept + (.slope * .xmax))
  return(list(x = .xmin, xend = .xmax, y = y, yend = yend))
}

segment_coords_pgls_bm <- calc_segment_coords(
  .xmin = min(log(prec_data_euth$adult_body_mass_g)),
  .xmax = max(log(prec_data_euth$adult_body_mass_g)),
  .intercept = coef(pgls_consensus_bm)[1],
  .slope = coef(pgls_consensus_bm)[2]
)

segment_coords_pgls_altrial <- calc_segment_coords(
  .xmin = min(log(prec_data_euth$adult_body_mass_g[which(prec_data_euth$precocity == "altricial")])),
  .xmax = max(log(prec_data_euth$adult_body_mass_g[which(prec_data_euth$precocity == "altricial")])),
  .intercept = coef(pgls_consensus)[1],
  .slope = coef(pgls_consensus)[2]
)

segment_coords_pgls_intermediate <- calc_segment_coords(
  .xmin = min(log(prec_data_euth$adult_body_mass_g[which(prec_data_euth$precocity == "intermediate")])),
  .xmax = max(log(prec_data_euth$adult_body_mass_g[which(prec_data_euth$precocity == "intermediate")])),
  .intercept = coef(pgls_consensus)[1] + coef(pgls_consensus)[3],
  .slope = coef(pgls_consensus)[2]
)

segment_coords_pgls_precocial <- calc_segment_coords(
  .xmin = min(log(prec_data_euth$adult_body_mass_g[which(prec_data_euth$precocity == "precocial")])),
  .xmax = max(log(prec_data_euth$adult_body_mass_g[which(prec_data_euth$precocity == "precocial")])),
  .intercept = coef(pgls_consensus)[1] + coef(pgls_consensus)[4],
  .slope = coef(pgls_consensus)[2]
)

gl_vs_bm <- ggplot(prec_data_euth, 
       aes(log(adult_body_mass_g), log(gestation_len_d), fill = precocity)) +
  geom_segment(x = segment_coords_pgls_bm$x,
               xend = segment_coords_pgls_bm$xend,
               y = segment_coords_pgls_bm$y,
               yend = segment_coords_pgls_bm$yend,
               color = "grey25", linetype = 2) +
  geom_segment(x = segment_coords_pgls_altrial$x,
               xend = segment_coords_pgls_altrial$xend,
               y = segment_coords_pgls_altrial$y,
               yend = segment_coords_pgls_altrial$yend,
               color = clrs["altricial"]) +
  geom_segment(x = segment_coords_pgls_intermediate$x,
               xend = segment_coords_pgls_intermediate$xend,
               y = segment_coords_pgls_intermediate$y,
               yend = segment_coords_pgls_intermediate$yend,
               color = clrs["intermediate"]) +
  geom_segment(x = segment_coords_pgls_precocial$x,
               xend = segment_coords_pgls_precocial$xend,
               y = segment_coords_pgls_precocial$y,
               yend = segment_coords_pgls_precocial$yend,
               color = clrs["precocial"]) +
  geom_point(shape = 21, stroke = 0.3, size = 1.75) +
  scale_fill_manual(values = clrs,
                    guide = guide_legend(override.aes = list(size = 3))) +
  # gl ~ bm
  annotate(geom = "segment",
           x = 5, xend = 6, y = 7, yend = 7,
           linetype = 2, color = "grey25") +
  annotate(geom = "text",
           x = 6.2, y = 7, 
           label = "log(gestation length) ~ log(adult body mass)",
           hjust = 0, family = "Source Sans Pro") +
  # gl ~ bm + prec
  annotate(geom = "segment",
           x = 5, xend = 6, y = 6.8, yend = 6.8,
           color = clrs["altricial"]) +
  annotate(geom = "segment",
           x = 5, xend = 6, y = 6.75, yend = 6.75,
           color = clrs["intermediate"]) +
  annotate(geom = "segment",
           x = 5, xend = 6, y = 6.7, yend = 6.7,
           color = clrs["precocial"]) +
  annotate(geom = "text",
           x = 6.2, y = 6.75, 
           label = "log(gestation length) ~ log(adult body mass) + precocity",
           hjust = 0, family = "Source Sans Pro") +
  labs(x = "log(adult body mass in grams)",
       y = "log(gestation length in days)") +
  theme_bw(base_family = "Source Sans Pro") +
  theme(panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.key.size = unit(12, "pt"),
        legend.background = element_blank())


gl_vs_bm_with_residuals_inset <- gl_vs_bm + 
  inset_element(prec_vs_bm_residuals, left = 0.65, right = 0.999, bottom = 0.001, top = 0.3)

ggsave(here(resdir, "gl-vs-bm_residuals-inset.pdf"), 
       gl_vs_bm_with_residuals_inset, device = cairo_pdf,
       width = 8, height = 6)

# end =========================================================================

