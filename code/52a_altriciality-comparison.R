# Arun Chavan
# Started: 2024-03-05

# background ==================================================================

# Compare the extent of altriciality between eutherians and non-eutherians.

# setup =======================================================================
library(tidyverse)
library(here)
library(ggtext)
library(patchwork)

resdir <- here("results/altriciality-comparison")
fs::dir_create(resdir)

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
  geom_histogram(position = position_dodge(preserve = "single")) +
  facet_grid(rows = vars(rank03), scales = "free_y") +
  scale_fill_manual(values = clrs, name = NULL) +
  labs(x = "log<sub>10</sub>(neonate body mass ÷ adult body mass)") +
  theme_bw(base_family = "Source Sans Pro") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_markdown(),
        legend.box.background = element_rect(color = "black"),
        legend.key.size = unit(0.12, "inch"), 
        legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1))


ggsave(here(resdir, "body-mass-ratio-distribution.pdf"), 
       alt_comp, width = 5, height = 6, device = cairo_pdf)

# overlapping range of the metric of altriciality =============================
prec <- prec %>% 
  mutate(log_bm_ratio = log10(bm_ratio)) 

max_non_eutherian_metric <- max(prec$log_bm_ratio[prec$rank03 != "Eutheria"])

prec %>% 
  filter(log_bm_ratio <= max_non_eutherian_metric) %>% 
  filter(rank03 == "Eutheria") %>% 
  print(width = Inf)

prec %>% 
  arrange(log_bm_ratio) %>%
  write_csv(here(resdir, "bm-ratio-sorted.csv"))

# neonate vs adult body mass ==================================================
bm_scatter <- ggplot(prec, 
                     aes(adult_body_mass_g, 
                         neonate_body_mass_g, 
                         fill = precocity)) +
  annotate("rect", 
           xmin = min(prec$adult_body_mass_g[prec$rank03 != "Eutheria"]) - 10, 
           xmax = max(prec$adult_body_mass_g[prec$rank03 != "Eutheria"]) + 10000, 
           ymin = min(prec$neonate_body_mass_g[prec$rank03 != "Eutheria"]) - 0.003, 
           ymax = max(prec$neonate_body_mass_g[prec$rank03 != "Eutheria"]) + 0.2,
           alpha = .4,
           fill = "grey",
           color = NA) + 
  annotate(geom = "text",
           x = 20000,
           y = 0.05,
           label = "Non-eutherian", 
           family = "Source Sans Pro",
           hjust = 1, vjust = 1) +
  geom_point(shape = 21, stroke = 0.25) +
  scale_x_log10(name = "Adult Body Mass (g)") +
  scale_y_log10(name = "Neonate Body Mass (g)") +
  scale_fill_manual(values = clrs, name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 3, stroke = 0.5), 
                             position = "inside")) +
  coord_fixed() +
  theme_bw(base_family = "Source Sans Pro") +
  theme(legend.position.inside = c(0.01, 0.99), 
        legend.justification.inside = c(0, 1),
        legend.box.background = element_blank(),
        legend.key.height = unit(0.15, "in"))


ggsave(here(resdir, "neonate-vs-adult-body-mass.pdf"), 
       plot = bm_scatter,
       width = 5.5, height = 5, device = cairo_pdf)

# by order ====================================================================
bm_ratio_by_order <- ggplot(prec, aes(log_bm_ratio, rank07, fill = precocity)) +
  geom_jitter(height = 0.1, shape = 21, stroke = 0.2, size = 1) +
  facet_grid(rows = vars(rank03), space = "free", scales = "free") +
  scale_fill_manual(values = clrs, name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 2, stroke = 0.3), 
                             position = "inside")) +
  labs(x = "log<sub>10</sub>(neonate body mass ÷ adult body mass)",
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
        legend.position.inside = c(0.01, 0.99),
        legend.justification.inside = c(0, 1),
        legend.key.size = unit(0.075, "in")) 

ggsave(here(resdir, "log-bm-ratio-by-order.pdf"), 
       bm_ratio_by_order, width = 3.1, height = 3.1, device = cairo_pdf)

# difference in medians =======================================================
prec %>% 
  group_by(rank03) %>% 
  summarise(median_log_bm_ratio = median(log_bm_ratio)) %>% 
  ungroup() %>% 
  mutate(diff_to_eutheria = median_log_bm_ratio - median_log_bm_ratio[rank03 == "Eutheria"])

# end =========================================================================
