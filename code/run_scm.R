# Arun Chavan
# Started: 2026-08-21

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(optparse)
library(conflicted)

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
    help = "Path to a named list object .rds with trees.",
    required = TRUE
  ) |>
  add_option(
    "--treeindex",
    help = "Index number of tree to use from treelist",
    required = TRUE
  ) |>
  add_option("--model", default = "all", help = "ER, SYM, ARD, or all")

opt <- parse_args(parser)

## data =======================================================================
treelist <- read_rds(opt$treelist)
tipdata <- read_rds(opt$tipdata)


stopifnot(!is.null(names(treelist)))
stopifnot(!is.null(names(tipdata)))

tr <- treelist[[opt$treeindex]]
stopifnot(identical(tr$tip.label, names(tipdata)))

# derived options =============================================================
seed <- opt$treeindex
treename <- names(treelist)[treeindex]
treetype <- ifelse(treename == "consensus", "consensus", "sample")

savesimmap <- ifelse(treetype == "consensus", TRUE, FALSE)

## setup directories ==========================================================
if (treetype == "consensus") {
  resdir <- here("results/scm", opt$analysis, "consensus")
}
if (treetype == "sample") {
  resdir <- here("results/scm", opt$analysis, "sample", treename)
}

###############################################################################
# run simmap ==================================================================
###############################################################################

if (opt$model == "all") {
  fit <- fit_models(tr, tipdata)
} else {
  fit <- fit_single_model(tr, tipdata, opt$model)
}

if (treetype == "consensus") {
  simmap <- run_simmap(fit, 1000)
  simmap_summary <- summary(simmap)
  simmap_ace <- simmap_summary$ace
}

if (treetype == "model") {
  simmap_ace <- run_simmap_and_get_ace(fit, 1000)
}

###############################################################################
# write output ================================================================
###############################################################################
write_rds(fit, here(resdir, "model-fit.rds"))
write_rds(simmap_ace, here(resdir, "ace.rds"))
if (savesimmap) {
  write_rds(simmap_summary, here(resdir, "simmap_summary.rds"), compress = "gz")
}

# end #########################################################################
