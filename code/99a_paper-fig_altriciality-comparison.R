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

# end =========================================================================
