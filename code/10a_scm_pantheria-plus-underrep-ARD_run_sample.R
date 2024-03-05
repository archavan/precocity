# Arun Chavan
# Started: 2024-02-02

# background ==================================================================

# Run stochastic character mapping on a samples of 100 trees from Upham et al. 

# setup =======================================================================
library(tidyverse)
library(here)
library(phytools)

source(here("code/utilities_scm.R"))

# read and process tree_id ====================================================
args <- commandArgs(trailingOnly = TRUE)
tree_id <- as.integer(args[1])
tree_name <- paste0("tr", str_pad(tree_id, 3, "left", "0"))

message("Running analyses for tree number ", tree_id)

# define results directory ====================================================
model_to_test <- "ARD"
analysis_name <- paste0("pantheria-plus-underrep-", model_to_test)

resdir <- here("results/scm", analysis_name, "sample", tree_name)
fs::dir_create(resdir)

# tip data ====================================================================
prec_data <- read_csv(here("data/03_coded/pantheria-plus-underrep-taxa.csv"))

prec_data <- prec_data %>% 
  mutate(precocity = fct(precocity, 
                         c("altricial", "intermediate", "precocial")))

prec_tipdata <- set_names(prec_data$precocity, prec_data$binomial)

# trees =======================================================================
tr100 <- read.nexus(here("data/trees/upham2019/sample/2024-02-21/output.nex"))

# run simmap ==================================================================
tr <- tr100[[tree_id]]
tr_pruned <- prune_tree_for_fitmk(tr, prec_tipdata)
fit_mk <- fit_single_model(tr_pruned, prec_tipdata, model_to_test)
simmap_ace <- run_simmap_and_get_ace(fit_mk, 1000) 

# write output ================================================================
write_rds(tr_pruned, file = here(resdir, "tree_pruned.rds"))
write_rds(fit_mk, here(resdir, "model-fit.rds"))
write_rds(simmap_ace, here(resdir, "ace.rds"))

# end =========================================================================

