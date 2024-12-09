library(tidyverse)
library(here)

figdir <- here("results/99_paper-figs/")

# Mean estimtes and 70% CI from Werneburg et al, Table 1
df <- tribble(
  ~Trait,          ~mean, ~lower, ~upper,
   "Birth",         125,   89,     161,
   "Fur developed", 149,   112,    186,
   "Eyelids open",  150,   95,     205
)

plac <- ggplot(df, aes(mean, Trait)) +
  geom_vline(xintercept = 125, linetype = 2, color = "red", linewidth = 0.25) +
  geom_point(size = 2) +
  geom_segment(aes(x = lower, xend = upper), linewidth = 0.5) +
  labs(x = "Days", y = NULL) +
  theme_bw(base_family = "Source Sans Pro", base_line_size = 0.25) +
  theme(
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    panel.grid.minor = element_blank()
    )

ggsave(here(figdir, "supp-fig_werneburg-estimates.png"), 
       plac, 
       width = 3, height = 1.5, units = "in", dpi = 600)

# end =========================================================================
