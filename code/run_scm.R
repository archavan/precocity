# Arun Chavan
# Started: 2026-08-21

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(phytools)
library(optparse)
library(conflicted)

conflicts_prefer(purrr::map)

source(here("code/utilities_scm.R"))

## parse arguments ============================================================

parser <- OptionParser() |>
  add_option(
    "--analysis",
    help = "Name of analysis. Used to create outdir",
    required = TRUE
  ) |>
  add_option(
    "--tipdata",
    help = "Path to tipdata file. Should be a named vector as .rds.",
    required = TRUE
  ) |>
  add_option(
    "--treelist",
    help = "Path to a named multiPhylo object .rds.",
    required = TRUE
  ) |>
  add_option(
    "--treeindex",
    help = "Index number of tree to use from treelist",
    type = "integer",
    required = TRUE
  ) |>
  add_option("--model", default = "all", help = "ER, SYM, ARD, or all") |>
  add_option(
    "--nsim",
    default = 1000,
    type = "integer",
    help = "Number of simmaps."
  )

opt <- parse_args(parser)
stopifnot(opt$model %in% c("all", "ER", "SYM", "ARD"))

## data and derived options ===================================================
seed <- opt$treeindex

treelist <- read_rds(opt$treelist)
tipdata <- read_rds(opt$tipdata)

stopifnot(opt$treeindex <= length(treelist), opt$treeindex >= 1)
stopifnot(!is.null(names(treelist)))
stopifnot(!is.null(names(tipdata)))

treename <- names(treelist)[opt$treeindex]
is_consensus <- treename == "consensus"

tr <- treelist[[treename]]
tipdata <- tipdata[tr$tip.label]
stopifnot(identical(tr$tip.label, names(tipdata)))

## setup directories ==========================================================
if (is_consensus) {
  resdir <- here("results/scm", opt$analysis, "consensus")
} else {
  resdir <- here("results/scm", opt$analysis, "sample", treename)
}
fs::dir_create(resdir)

###############################################################################
# run simmap ==================================================================
###############################################################################
set.seed(seed)
if (opt$model == "all") {
  fit <- fit_models(tr, tipdata)
} else {
  fit <- fit_single_model(tr, tipdata, opt$model)
}

set.seed(seed)
if (is_consensus) {
  message("Running simmap...")
  simmap <- run_simmap(fit, opt$nsim)
  simmap_summary <- summary(simmap)
  simmap_ace <- simmap_summary$ace
} else {
  message("Running simmap...")
  simmap_ace <- run_simmap_and_get_ace(fit, opt$nsim)
}

###############################################################################
# write output ================================================================
###############################################################################
write_rds(fit, here(resdir, "model-fit.rds"), compress = "gz")
write_rds(simmap_ace, here(resdir, "ace.rds"), compress = "gz")
write_rds(tr, here(resdir, "tree_pruned.rds"), compress = "gz")
if (is_consensus) {
  write_rds(simmap_summary, here(resdir, "simmap_summary.rds"), compress = "gz")
}

# end #########################################################################
