library(tidyverse)

df <- tribble(
  ~trait,    ~mean, ~lower, ~upper,
   "birth",   125,   89,     161,
   "fur",     149,   112,    186,
   "eyelids", 150,   95,     205
)

plac <- ggplot(df, aes(mean, trait)) +
  geom_vline(xintercept = 125, linetype = 2, color = "red") +
  geom_point(size = 3) +
  geom_segment(aes(x = lower, xend = upper)) +
  labs(x = "days", y = NULL) +
  theme_bw(base_family = "Source Sans Pro") +
  theme(panel.grid.minor = element_blank())

ggsave("results/werneburg-estimates.pdf", 
       plac, width = 5, height = 1.5, device = cairo_pdf)
