# Arun Chavan
# Started: 2024-11-21

# setup =======================================================================
library(ape)
library(nlme)
library(tidyverse)
library(here)
library(glue)
library(patchwork)


clrs <- c(altricial = '#fc8d59',
          intermediate = '#ffffbf',
          precocial = '#91bfdb')

resdir <- here("results/99_paper-figs")
fs::dir_create(resdir)

text_size_base <- 6
text_size_min <- 5
text_size_max <- 7
text_size_tag <- 8
segment_linewidth <- 0.25
geom_pt_size <- 1.25
use_font <- "Source Sans Pro"

# data ========================================================================
prec <- read_csv(here("data/03_coded/case78-plus-pantheria.csv"))  
pt <- read_csv(here("data/02_pruned/pantheria_upham2019.csv"))
tr_consensus_full <- read.nexus(here("data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"))

###############################################################################
# gestation length vs body mass ###############################################
###############################################################################

# trait data for regression ===================================================
prec_data_full <- left_join(prec, pt) %>% 
  filter(!is.na(adult_body_mass_g) & !is.na(gestation_len_d))

prec_data_full <- mutate(prec_data_full, precocity = fct(precocity, 
                                                         c("altricial",
                                                           "intermediate",
                                                           "precocial")))

prec_data <- filter(prec_data_full, rank03 == "Eutheria") # regression only for eutheria

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

# model with body mass and precocity ==========================================

# match order of species in tree and data
prec_data <- prec_data[match(tr_consensus$tip.label, prec_data$binomial), ]
stopifnot(identical(tr_consensus$tip.label, prec_data$binomial))

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
pgls_consensus_summary <- summary(pgls_consensus)
pgls_consensus_summary$tTable %>% as.data.frame()

# residuals of body mass regression agaist precocity ==========================
pgls_consensus_bm <- gls(log(gestation_len_d) ~ log(adult_body_mass_g),
                         data = prec_data, 
                         correlation = cor_lambda_consensus,
                         na.action = na.exclude)

prec_vs_bm_residuals <- prec_data %>% 
  mutate(bm_residuals = pgls_consensus_bm$residuals) %>% 
  ggplot(aes(precocity, bm_residuals)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey25", 
             linewidth = segment_linewidth) +
  ggforce::geom_sina(aes(fill = precocity), shape = 21, size = geom_pt_size/2, stroke = 0.15) +
  scale_fill_manual(values = clrs, guide = "none") +
  labs(y = "residuals") +
  theme_bw(base_family = use_font) +
  theme(plot.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(size = text_size_base),
        axis.title.y = element_text(size = text_size_max),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.25),
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
                   aes(log(adult_body_mass_g), log(gestation_len_d), 
                       fill = precocity)) +
  geom_segment(x = segment_coords_pgls_bm$x,
               xend = segment_coords_pgls_bm$xend,
               y = segment_coords_pgls_bm$y,
               yend = segment_coords_pgls_bm$yend,
               color = "grey25", linetype = 2, 
               linewidth = segment_linewidth) +
  geom_segment(x = segment_coords_pgls_altrial$x,
               xend = segment_coords_pgls_altrial$xend,
               y = segment_coords_pgls_altrial$y,
               yend = segment_coords_pgls_altrial$yend,
               color = clrs["altricial"],
               linewidth = segment_linewidth) +
  geom_segment(x = segment_coords_pgls_intermediate$x,
               xend = segment_coords_pgls_intermediate$xend,
               y = segment_coords_pgls_intermediate$y,
               yend = segment_coords_pgls_intermediate$yend,
               color = clrs["intermediate"],
               linewidth = segment_linewidth) +
  geom_segment(x = segment_coords_pgls_precocial$x,
               xend = segment_coords_pgls_precocial$xend,
               y = segment_coords_pgls_precocial$y,
               yend = segment_coords_pgls_precocial$yend,
               color = clrs["precocial"],
               linewidth = segment_linewidth) +
  geom_point(shape = 21, stroke = 0.15, size = geom_pt_size) +
  scale_fill_manual(values = clrs,
                    guide = guide_legend(override.aes = list(size = 2.5), 
                                         position = "inside")) +
  # gl ~ bm
  annotate(geom = "segment",
           x = 5, xend = 6, y = 7, yend = 7,
           linetype = 2, color = "grey25",
           linewidth = segment_linewidth) +
  annotate(geom = "text",
           x = 6.2, y = 7, 
           label = "log(gestation length) ~ log(adult body mass)",
           hjust = 0, 
           family = use_font,
           size = text_size_base/.pt) +
  # gl ~ bm + prec
  annotate(geom = "segment",
           x = 5, xend = 6, y = 6.8, yend = 6.8,
           color = clrs["altricial"],
           linewidth = segment_linewidth) +
  annotate(geom = "segment",
           x = 5, xend = 6, y = 6.75, yend = 6.75,
           color = clrs["intermediate"],
           linewidth = segment_linewidth) +
  annotate(geom = "segment",
           x = 5, xend = 6, y = 6.7, yend = 6.7,
           color = clrs["precocial"],
           linewidth = segment_linewidth) +
  annotate(geom = "text",
           x = 6.2, y = 6.75, 
           label = "log(gestation length) ~ log(adult body mass) + precocity",
           hjust = 0, 
           family = use_font,
           size = text_size_base/.pt) +
  labs(x = "log(adult body mass in grams)",
       y = "log(gestation length in days)") +
  theme_bw(base_family = use_font) +
  theme(panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = text_size_base),
        axis.title = element_text(size = text_size_max),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position.inside = c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = text_size_min),
        legend.key.size = unit(text_size_max, "pt"),
        legend.background = element_blank())


gl_vs_bm

###############################################################################
# combine subplots ############################################################
###############################################################################

gl_vs_bm_with_residuals_inset <- gl_vs_bm + 
  inset_element(prec_vs_bm_residuals, 
                left = 0.675, 
                right = 0.999, 
                bottom = 0.001, 
                top = 0.35) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = text_size_tag))

gl_vs_bm_with_residuals_inset

# save files ==================================================================

write_rds(pgls_consensus, here(resdir, "supp-tbl_pgls-coeffs-for-table.rds"))

ggsave(here(resdir, "supp-fig_pgls.png"), 
       plot = gl_vs_bm_with_residuals_inset, 
       width = 5, height = 4, dpi = 600)

ggsave(here(resdir, "supp-fig_pgls.pdf"), 
       plot = gl_vs_bm_with_residuals_inset, 
       device = cairo_pdf,
       width = 5, height = 4)

###############################################################################
# plot with only gl ~ bm regression and recreate the plot from Martin and MacLarnon ===
###############################################################################

# model with body mass and precocity ==========================================
bm_only_pgls <- ggplot(prec_data, 
       aes(log(adult_body_mass_g), log(gestation_len_d), 
           fill = precocity)) +
  geom_segment(x = segment_coords_pgls_bm$x,
               xend = segment_coords_pgls_bm$xend,
               y = segment_coords_pgls_bm$y,
               yend = segment_coords_pgls_bm$yend,
               color = "grey25", linetype = 2, 
               linewidth = segment_linewidth) +
  geom_point(shape = 21, stroke = 0.25, size = geom_pt_size) +
  scale_fill_manual(values = clrs,
                    guide = guide_legend(override.aes = list(size = 2.5), 
                                         position = "inside")) +
  # gl ~ bm
  annotate(geom = "segment",
           x = 6, xend = 7, y = 7, yend = 7,
           linetype = 2, color = "grey25",
           linewidth = segment_linewidth) +
  annotate(geom = "text",
           x = 7.2, y = 7, 
           label = "log(gestation length) ~ log(adult body mass)",
           hjust = 0, 
           family = use_font,
           size = text_size_base/.pt) +
  labs(x = "log(adult body mass in grams)",
       y = "log(gestation length in days)") +
  theme_bw(base_family = use_font) +
  theme(panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = text_size_base),
        axis.title = element_text(size = text_size_max),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position.inside = c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = text_size_base),
        legend.key.size = unit(text_size_max, "pt"),
        legend.background = element_blank())

# histogram of residuals ======================================================
residual_histogram <- prec_data %>% 
  mutate(bm_residuals = pgls_consensus_bm$residuals) %>% 
  ggplot(aes(bm_residuals, fill = precocity)) +
  geom_histogram(bins = 50, color = "black", size = 0.25) +
  scale_fill_manual(values = clrs,
                    guide = guide_legend(override.aes = list(size = 2.5), 
                                         position = "inside")) +
  labs(x = "Residuals", y = "Count") +
  theme_bw(base_family = use_font) +
  theme(panel.border = element_rect(linewidth = 0.5, fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = text_size_base),
        axis.title = element_text(size = text_size_max),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = text_size_base),
        legend.key.size = unit(text_size_max, "pt"),
        legend.background = element_blank())

# combine plots ===============================================================
bm_only_patch <- wrap_plots(bm_only_pgls, residual_histogram, nrow = 1) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = text_size_tag))

# save plots ==================================================================
ggsave(here(resdir, "supp-fig_pgls_v2.png"), 
       plot = bm_only_patch, 
       width = 6.5, height = 3, dpi = 600)

ggsave(here(resdir, "supp-fig_pgls_v2.pdf"), 
       plot = bm_only_patch, 
       device = cairo_pdf,
       width = 6.5, height = 3)

###############################################################################
