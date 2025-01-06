# Arun Chavan
# Started: 2024-11-25

# background ==================================================================

# Compare the extent of altriciality between eutherians and non-eutherians.

# setup =======================================================================
library(tidyverse)
library(here)
library(ggtext)
library(cowplot)

resdir <- here("results/99_paper-figs/")

clrs <- c(altricial = '#fc8d59',
          intermediate = '#ffffbf',
          precocial = '#91bfdb')

# data ========================================================================
prec <- read_csv("data/03_coded/case78-plus-pantheria.csv")
pt <- read_csv("data/02_pruned/pantheria_upham2019.csv")

prec <- left_join(prec, pt)

prec <- filter(prec, 
               !is.na(adult_body_mass_g) & !is.na(neonate_body_mass_g))

prec <- mutate(prec,
               bm_ratio = neonate_body_mass_g / adult_body_mass_g)

# comparison ==================================================================
alt_comp <- ggplot(prec, aes(log10(bm_ratio), fill = precocity)) +
  geom_histogram(position = position_dodge(preserve = "single"),
                 color = "black", 
                 linewidth = 0.1,
                 binwidth = 0.2) +
  facet_grid(rows = vars(rank03), scales = "free_y") +
  scale_fill_manual(values = clrs, name = NULL) +
  labs(x = "log<sub>10</sub>(neonate body mass [g] ÷ adult body mass [g])") +
  theme_bw(base_family = "Source Sans Pro", base_size = 6) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.title.x = element_markdown(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.box.background = element_blank(),
    legend.key.size = unit(0.12, "inch"), 
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 7)
  )


ggsave(here(resdir, "fig_altriciality-comparison.pdf"), 
       alt_comp, width = 3, height = 3.5, device = cairo_pdf)

ggsave(here(resdir, "fig_altriciality-comparison.png"), 
       alt_comp, width = 3, height = 3.5, dpi = 600, units = "in")

# by order ====================================================================
prec <- prec %>% 
  mutate(log_bm_ratio = log10(bm_ratio)) 

bm_ratio_by_order <- ggplot(prec, aes(log_bm_ratio, rank07, fill = precocity)) +
  geom_jitter(height = 0.1, shape = 21, stroke = 0.25, size = 1.2) +
  facet_grid(rows = vars(rank03), space = "free", scales = "free") +
  scale_fill_manual(values = clrs, name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 2, stroke = 0.3), 
                             position = "top")) +
  labs(x = "log<sub>10</sub>(neonate body mass ÷ adult body mass)",
       y = "Order") +
  theme_bw(base_family = "Source Sans Pro", 
           base_line_size = 0.25) +
  theme(strip.clip = "off",
        strip.text.y = element_text(angle = 0, hjust = 0, size = 7,
                                    color = "black"),
        strip.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25),
        panel.grid.minor = element_blank(),
        axis.title.x = element_markdown(size = 7, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7, 
                                   margin = margin(l = 0),
                                   color = "black"),
        legend.margin = margin(0, 0, 0, 0, "pt"),
        legend.box.spacing = unit(0, "pt"),
        legend.key.size = unit(10, "pt"),
        legend.spacing.y = unit(20, "pt")) 

ggsave(here(resdir, "fig_altriciality-comparison_by-order.pdf"), 
       bm_ratio_by_order, width = 3.5, height = 3.8, device = cairo_pdf)

ggsave(here(resdir, "fig_altriciality-comparison_by-order.png"), 
       bm_ratio_by_order, width = 3.5, height = 3.8, dpi = 600, units = "in")

# bm ratio % ==================================================================
bm_ratio_pct <- ggplot(prec, aes(bm_ratio*100, rank07, fill = precocity)) +
  geom_jitter(height = 0.1, shape = 21, stroke = 0.2, size = 1) +
  facet_grid(rows = vars(rank03), space = "free", scales = "free") +
  scale_fill_manual(values = clrs, name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 2, stroke = 0.3), 
                             position = "inside")) +
  labs(x = "Relative Neonatal Body Mass (%)",
       y = "Order") +
  theme_bw(base_family = "Source Sans Pro", base_line_size = 0.25) +
  theme(strip.clip = "off",
        strip.text.y = element_text(angle = 0, hjust = 0, size = 6),
        strip.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25),
        panel.grid.minor = element_blank(),
        axis.title.x = element_markdown(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.background = element_blank(),
        legend.position.inside = c(0.99, 0.99),
        legend.justification.inside = c(1, 1),
        legend.key.size = unit(0.075, "in")) 

ggsave(here(resdir, "fig_altriciality-comparison_by-order_pct.pdf"), 
       bm_ratio_pct, width = 3.1, height = 3.1, device = cairo_pdf)

ggsave(here(resdir, "fig_altriciality-comparison_by-order_pct.png"), 
       bm_ratio_pct, width = 3.1, height = 3.1, dpi = 600, units = "in")

# end =========================================================================
