library(here)
library(tidyverse)
library(patchwork)

figdir <- here("results/99_paper-figs/")

# read pre-prepared plots =====================================================

read_plot_rds <- function(.analysis) {
  read_rds(here(figdir, 
                paste0("supp-fig_scm_", .analysis, "_model-choice.rds")))
}

mwe_case78 <- read_plot_rds("case78") +
  ggtitle("Case (1978)")
mwe_case_plus_pt <- read_plot_rds("case78-plus-pantheria") +
  ggtitle("Case (1978) + PanTHERIA")

# patchwork ===================================================================

mwe_patch <- wrap_plots(mwe_case78, mwe_case_plus_pt, nrow = 1) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.title = element_text(size = 8),
    plot.tag = element_text(size = 8, face = "bold")
  )

# write files =================================================================

ggsave(here(figdir, "supp-fig_model-fit-effect.png"),
       plot = mwe_patch,
       width = 5, height = 2.5, dpi = 600, units = "in")

# end =========================================================================
