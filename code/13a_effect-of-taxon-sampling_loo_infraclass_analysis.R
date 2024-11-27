# Arun Chavan
# Started: 2024-11-26

# background ==================================================================

# Run stochastic character mapping on consensus from Upham et al. while leaving
# one order out at a time to test the effect of taxon sampling systematically.

# setup =======================================================================
library(phytools)
library(tidyverse)
library(here)

source(here("code/utilities_scm.R"))

# define results directory ====================================================

resdir <- here("results/loo_infraclass")
fs::dir_create(resdir)

# tip data ====================================================================
prec_data <- read_csv(here("data/03_coded/case78-plus-pantheria.csv"))

prec_data <- prec_data %>% 
  mutate(precocity = fct(precocity, 
                         c("altricial", "intermediate", "precocial")))

# trees =======================================================================
tr <- read.nexus(here("data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"))

# generate receptacle for results =============================================
df <- prec_data %>% 
  filter(rank03 == "Eutheria") %>% 
  select(rank05) %>% 
  rename(leave_out = rank05) %>% 
  distinct()

df <- df %>% 
  mutate(prec = map(leave_out, ~ filter(prec_data, rank05 != .x))) %>% 
  mutate(tipdata = map(prec, ~ set_names(.x$precocity, .x$binomial))) %>% 
  mutate(tr_pruned = map(tipdata, ~ prune_tree_for_fitmk(tr, .x))) 

# run simmap ==================================================================

df <- df %>% 
  mutate(fit_aov = map2(tr_pruned, tipdata, fit_models, .progress = TRUE)) %>% 
  mutate(simmap_consensus = map(fit_aov, ~ run_simmap(.x, 1000), .progress = TRUE)) %>% 
  mutate(simmap_consensus_summary = map(simmap_consensus, summary, .progress = TRUE)) %>% 
  mutate(simmap_ace = map(simmap_consensus_summary, ~ .x$ace, .progress = TRUE))

# write output ================================================================
# simmap summary object takes up too much space
df_without_summary <- select(df, -simmap_consensus_summary)
df_minimal <- select(df, leave_out, prec, tipdata, tr_pruned, fit_aov, simmap_ace)

# write_rds(df, file = here(resdir, "output.rds"), compress = "gz")
write_rds(df_minimal, 
          file = here(resdir, "output_minimal.rds"), 
          compress = "gz")

# end =========================================================================
