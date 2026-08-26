# Arun Chavan
# Started: 2026-08-23

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(phytools)
library(ape)
library(conflicted)

conflicts_prefer(purrr::map, dplyr::filter)

clrs <- c(
  altricial = '#fc8d59',
  intermediate = '#ffffbf',
  precocial = '#91bfdb'
)

analysis <- "loo_order"
resdir <- here("results/scm", analysis)

source(here("code/utilities_taxonomic.R"))

###############################################################################
# data ========================================================================
###############################################################################

focal_clades <- c(
  "Mammalia",
  "Theria",
  "Eutheria",
  "Metatheria",
  "Boreoeutheria",
  "Xenarthra",
  "Afrotheria",
  "Laurasiatheria",
  "Euarchontoglires"
)

taxonsets <- read_csv(here(resdir, "taxonsets.csv"))
trees <- read_csv(here(resdir, "treeindices.csv")) |>
  mutate(treetype = ifelse(treename == "consensus", "consensus", "sample")) |>
  mutate(
    respath = ifelse(
      treetype == "consensus",
      "consensus",
      paste0("sample/", treename)
    )
  )

read_taxonset_res <- function(.taxonset) {
  list(
    treewise = trees |>
      mutate(
        ace = map(
          respath,
          ~ read_rds(here(resdir, .taxonset, .x, "ace.rds"))
        )
      ) |>
      mutate(
        fit = map(
          respath,
          ~ read_rds(here(resdir, .taxonset, .x, "model-fit.rds"))
        )
      ) |>
      mutate(
        tr = map(
          respath,
          ~ read_rds(here(resdir, .taxonset, .x, "tree_pruned.rds"))
        )
      ),
    tipdata = read_rds(here(resdir, .taxonset, "tipdata.rds")),
    taxa = read_csv(
      here(resdir, .taxonset, "taxa.csv"),
      show_col_types = FALSE
    )
  )
}


taxonsets <- taxonsets |>
  mutate(res = map(taxonset, read_taxonset_res, .progress = TRUE))


###############################################################################
# get PP at crown major node ==============================================
###############################################################################

get_pp_at_node <- function(.ace, .node) .ace[which(rownames(.ace) == .node), ]

get_pp_for_clades <- function(.ace, .phy, .taxa, .clades = focal_clades) {
  map(
    set_names(focal_clades),
    ~ get_pp_at_node(.ace, get_node(.phy, .x, .taxa))
  ) |>
    bind_rows(.id = "clade")
}

get_qmatrices <- function(.fit) {
  fit_names <- attr(.fit, "row.names")
  map(seq_along(fit_names), ~ as.Qmatrix(attr(.fit, "models")[[.x]])) |>
    set_names(fit_names)
}

for (i in seq_len(nrow(taxonsets))) {
  message(paste0("Working on taxonset: ", taxonsets$taxonset[i]))
  df_i <- taxonsets$res[[i]]$treewise
  df_i$pp_focal <- vector("list", nrow(df_i))
  df_i$q_matrix <- vector("list", nrow(df_i))
  for (j in seq_len(nrow(df_i))) {
    df_i$pp_focal[[j]] <- get_pp_for_clades(
      df_i$ace[[j]],
      df_i$tr[[j]],
      taxonsets$res[[i]]$taxa
    )
    df_i$q_matrix[[j]] <- get_qmatrices(df_i$fit[[j]])
  }
  taxonsets$res[[i]]$treewise <- df_i
  rm(df_i, i, j)
}

pp <- taxonsets |>
  mutate(
    pp = map(
      res,
      ~ .x$treewise |> select(treeindex, treename, treetype, pp_focal)
    )
  ) |>
  select(taxonset, pp) |>
  unnest(col = pp) |>
  unnest(col = pp_focal)

plot_pp_dist_for_clade <- function(.clade) {
  pp |>
    filter(treetype == "sample") |>
    filter(clade == .clade) |>
    pivot_longer(
      cols = all_of(names(clrs)),
      values_to = "pp_at_node",
      names_to = "precocity"
    ) |>
    ggplot(aes(pp_at_node, taxonset, fill = precocity)) +
    geom_violin(scale = "width") +
    scale_x_continuous(name = "Posterior Probability", limits = c(0, 1)) +
    scale_y_discrete(name = "Left Out") +
    scale_fill_manual(
      values = clrs,
      breaks = names(clrs),
      labels = str_to_sentence(names(clrs)),
      name = NULL
    ) +
    facet_grid(
      cols = vars(precocity),
      labeller = labeller(precocity = str_to_sentence)
    ) +
    labs(title = paste0("Ancestral state for ", .clade)) +
    theme_classic() +
    theme(
      panel.grid.major = element_line(color = "grey", linewidth = 0.25),
      panel.widths = unit(6, "in"),
      panel.heights = unit(4.5, "in"),
      strip.background = element_blank()
    )
}


# write plots and files =======================================================

pdf(
  here(resdir, "pp-distribution_focal-nodes.pdf"),
  onefile = TRUE,
  width = 9,
  height = 6
)
map(focal_clades, ~ plot_pp_dist_for_clade(.x))
dev.off()

write_csv(pp, here(resdir, "pp-distribution_focal-nodes.csv"))

write_rds(taxonsets, here(resdir, "processed.rds"), compress = "gz")

# end #########################################################################
