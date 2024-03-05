# Arun Chavan
# Started: 2024-02-22

# background ==================================================================

# Run stochastic character mapping on consensus from Upham et al. 

# setup =======================================================================
library(tidyverse)
library(here)
library(phytools)

source(here("code/utilities_scm.R"))

# define results directory ====================================================
model_to_test <- "ER"
analysis_name <- paste0("pantheria-plus-underrep-", model_to_test)

resdir <- here("results/scm", analysis_name, "consensus")
fs::dir_create(resdir)

# tip data ====================================================================
prec_data <- read_csv(here("data/03_coded/pantheria-plus-underrep-taxa.csv"))

prec_data <- prec_data %>% 
  mutate(precocity = fct(precocity, 
                         c("altricial", "intermediate", "precocial")))

prec_tipdata <- set_names(prec_data$precocity, prec_data$binomial)

# trees =======================================================================
tr <- read.nexus(here("data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"))

# run simmap ==================================================================
tr_pruned <- prune_tree_for_fitmk(tr, prec_tipdata)
fit_mk <- fit_single_model(tr_pruned, prec_tipdata, model_to_test)
simmap_consensus <- run_simmap(fit_mk, 1000)
simmap_consensus_summary <- summary(simmap_consensus)
simmap_ace <- simmap_consensus_summary$ace

# write output ================================================================
write_rds(tr_pruned, file = here(resdir, "tree_pruned.rds"))
write_rds(fit_mk, here(resdir, "model-fit.rds"))
write_rds(simmap_ace, here(resdir, "ace.rds"))
write_rds(simmap_consensus_summary, here(resdir, "simmap_summary.rds"), compress = "gz")

# end =========================================================================

