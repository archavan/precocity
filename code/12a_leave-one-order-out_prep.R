# Arun Chavan
# Started: 2024-11-26

###############################################################################
# background ==================================================================
###############################################################################

# Run stochastic character mapping while leaving one order out at a time to
# test the effect of taxon sampling systematically. Each order is run over the
# consensus tree and all 100 sampled trees from Upham et al., so that every
# order gives a distribution of ancestral state estimates rather than a single
# one.

###############################################################################
# setup =======================================================================
###############################################################################
library(phytools)
library(tidyverse)
library(here)
library(glue)
library(ape)
library(conflicted)

conflicts_prefer(dplyr::filter, purrr::map)

source(here("code/utilities_scm.R"))
source(here("code/utilities_scm_prep.R"))

## parameters =================================================================
analysis <- "loo_order"
rank_column <- "rank07"
model_to_test <- "all"
nsim <- 1000

###############################################################################
# data ========================================================================
###############################################################################

prec_data <- read_precocity()

treelist <- build_treelist(prec_data)

###############################################################################
# prepare for analysis ========================================================
###############################################################################

taxonsets <- taxonsets_leave_one_out(prec_data, rank_column)

prep_scm_batch(
  taxonsets,
  .trait_data = prec_data,
  .treelist = treelist,
  .analysis = analysis,
  .model = model_to_test,
  .nsim = nsim
)

# end =========================================================================
