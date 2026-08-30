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

clrs <- c(
  altricial = '#fc8d59',
  intermediate = '#ffffbf',
  precocial = '#91bfdb'
)

# data ========================================================================
prec <- read_csv("data/03_coded/precocity.csv")
pt <- read_csv("data/02_pruned/pantheria_upham2019.csv")

prec <- left_join(prec, pt)

prec <- filter(prec, !is.na(adult_body_mass_g) & !is.na(neonate_body_mass_g))

prec <- mutate(prec, bm_ratio = neonate_body_mass_g / adult_body_mass_g)

prec <- filter(prec, rank03 != "Prototheria")

# by order ====================================================================
prec <- prec %>%
  mutate(log_bm_ratio = log10(bm_ratio))

bm_ratio_by_order <- ggplot(prec, aes(log_bm_ratio, rank07, fill = precocity)) +
  geom_jitter(height = 0.1, shape = 21, stroke = 0.25, size = 1.2) +
  facet_grid(rows = vars(rank03), space = "free", scales = "free") +
  scale_fill_manual(
    values = clrs,
    name = NULL,
    breaks = names(clrs),
    labels = str_to_sentence
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 2, stroke = 0.3),
      position = "top"
    )
  ) +
  labs(
    x = "log<sub>10</sub>(neonate body mass ÷ adult body mass)",
    y = "Order"
  ) +
  theme_bw(base_family = "Source Sans 3", base_line_size = 0.25) +
  theme(
    strip.clip = "off",
    strip.text.y = element_text(
      angle = 0,
      hjust = 0,
      size = 7,
      color = "black"
    ),
    strip.background = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_markdown(size = 7, color = "black"),
    axis.title.y = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.text = element_text(
      size = 7,
      margin = margin(l = 0),
      color = "black"
    ),
    legend.margin = margin(0, 0, 0, 0, "pt"),
    legend.box.spacing = unit(0, "pt"),
    legend.key.size = unit(10, "pt"),
    legend.spacing.y = unit(20, "pt")
  )


bm_ratio_by_order

quartz(
  type = "pdf",
  file = here(resdir, "fig_altriciality-comparison_by-order.pdf"),
  width = 3.5,
  height = 3.7,
  family = "Source Sans 3",
)
print(bm_ratio_by_order)
dev.off()

ggsave(
  here(resdir, "fig_altriciality-comparison_by-order.png"),
  bm_ratio_by_order,
  width = 3.5,
  height = 3.7,
  dpi = 600,
  units = "in"
)

# end =========================================================================
