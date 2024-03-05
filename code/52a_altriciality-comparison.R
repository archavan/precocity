# Arun Chavan
# Started: 2024-03-05

# background ==================================================================

# Compare the extent of altriciality between eutherians and non-eutherians.

# setup =======================================================================
library(tidyverse)
library(here)
library(ggtext)

resdir <- here("results/altriciality-comparison")
fs::dir_create(resdir)

clrs <- c(altricial = '#fc8d59',
          intermediate = '#ffffbf',
          precocial = '#91bfdb')

# data ========================================================================
prec <- read_csv("data/03_coded/pantheria-plus-underrep-taxa.csv")
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

 # end =========================================================================
