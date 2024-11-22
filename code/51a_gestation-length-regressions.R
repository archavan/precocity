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
prec_data <- read_csv(here("data/03_coded/case78-plus-pantheria.csv"))
pt <- read_csv(here("data/02_pruned/pantheria_upham2019.csv"))
prec_data <- left_join(prec_data, pt)
prec_data <- filter(prec_data, rank03 == "Eutheria") # regression only for eutheria

prec_data <- filter(prec_data, 
                    !is.na(adult_body_mass_g) & !is.na(gestation_len_d))

prec_data <- mutate(prec_data, precocity = fct(precocity, 
                                               c("altricial",
                                                 "intermediate",
                                                 "precocial")))

# trees
tr_consensus_full <- read.nexus(here("data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"))

# taxonomic information =======================================================
taxa <- prec_data %>% 
  select(all_of(c(paste0("rank0", 1:8), "family", "binomial"))) %>% 
  distinct()

# prepare data ================================================================
## consensus tree for visualization -------------------------------------------
tr_consensus <- drop.tip(
  tr_consensus_full,
  which(!(tr_consensus_full$tip.label %in% prec_data$binomial))
)
tr_consensus <- ladderize(tr_consensus)

###############################################################################
# PGLS using the consensus tree for visualization #############################
###############################################################################

# model with body mass and precocity ==========================================

# match order of species in tree and data
prec_data <- prec_data[match(tr_consensus$tip.label, prec_data$binomial), ]
waldo::compare(tr_consensus$tip.label, prec_data$binomial)

spp <- prec_data$binomial

# initiate correlation structure
cor_lambda_consensus <- corPagel(value = 1, tr_consensus, form = ~spp)

# model
pgls_consensus <- gls(log(gestation_len_d) ~ log(adult_body_mass_g) + precocity,
                      data = prec_data, 
                      correlation = cor_lambda_consensus,
                      na.action = na.exclude)

anova(pgls_consensus)
summary(pgls_consensus)

# residuals of body mass regression agaist precocity ==========================
pgls_consensus_bm <- gls(log(gestation_len_d) ~ log(adult_body_mass_g),
                         data = prec_data, 
                         correlation = cor_lambda_consensus,
                         na.action = na.exclude)

prec_vs_bm_residuals <- prec_data %>% 
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
  .xmin = min(log(prec_data$adult_body_mass_g)),
  .xmax = max(log(prec_data$adult_body_mass_g)),
  .intercept = coef(pgls_consensus_bm)[1],
  .slope = coef(pgls_consensus_bm)[2]
)

segment_coords_pgls_altrial <- calc_segment_coords(
  .xmin = min(log(prec_data$adult_body_mass_g[which(prec_data$precocity == "altricial")])),
  .xmax = max(log(prec_data$adult_body_mass_g[which(prec_data$precocity == "altricial")])),
  .intercept = coef(pgls_consensus)[1],
  .slope = coef(pgls_consensus)[2]
)

segment_coords_pgls_intermediate <- calc_segment_coords(
  .xmin = min(log(prec_data$adult_body_mass_g[which(prec_data$precocity == "intermediate")])),
  .xmax = max(log(prec_data$adult_body_mass_g[which(prec_data$precocity == "intermediate")])),
  .intercept = coef(pgls_consensus)[1] + coef(pgls_consensus)[3],
  .slope = coef(pgls_consensus)[2]
)

segment_coords_pgls_precocial <- calc_segment_coords(
  .xmin = min(log(prec_data$adult_body_mass_g[which(prec_data$precocity == "precocial")])),
  .xmax = max(log(prec_data$adult_body_mass_g[which(prec_data$precocity == "precocial")])),
  .intercept = coef(pgls_consensus)[1] + coef(pgls_consensus)[4],
  .slope = coef(pgls_consensus)[2]
)

gl_vs_bm <- ggplot(prec_data, 
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

