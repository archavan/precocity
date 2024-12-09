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
analysis_name <- "case78"

resdir <- here("results/scm", analysis_name, "consensus")
fs::dir_create(resdir)

# tip data ====================================================================
prec_data <- read_csv(here("data/03_coded/case78/precocity_case78_v1.csv"))

prec_data <- prec_data %>% 
  mutate(precocity_recoded = case_when(
    precocity == "fetal" ~ "altricial",
    precocity == "A" ~ "altricial",
    precocity == "SA" ~ "altricial",
    precocity == "SP" ~ "intermediate",
    precocity == "P" ~ "precocial"
  )) %>% 
  mutate(precocity_recoded = fct(precocity_recoded, 
                                 c("altricial", "intermediate", "precocial")))

prec_tipdata <- set_names(prec_data$precocity_recoded, prec_data$binomial)

# trees =======================================================================
tr <- read.nexus(here("data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"))

# run simmap ==================================================================
tr_pruned <- prune_tree_for_fitmk(tr, prec_tipdata)
fit_aov <- fit_models(tr_pruned, prec_tipdata)
simmap_consensus <- run_simmap(fit_aov, 1000)
simmap_consensus_summary <- summary(simmap_consensus)
simmap_ace <- simmap_consensus_summary$ace

# write output ================================================================
write_rds(tr_pruned, file = here(resdir, "tree_pruned.rds"))
write_rds(fit_aov, here(resdir, "model-weights.rds"))
write_rds(simmap_ace, here(resdir, "ace.rds"))
write_rds(simmap_consensus_summary, here(resdir, "simmap_summary.rds"), compress = "gz")

# end =========================================================================

