# Arun Chavan
# Started: 2024-02-02

# background ==================================================================

# Run stochastic character mapping on a samples of 100 trees from Upham et al. 

# setup =======================================================================
library(tidyverse)
library(here)
library(phytools)

resdir <- here("results/scm/pantheria/out")
fs::dir_create(resdir)

# data ========================================================================
# precocity data
prec_data <- read_csv(here("data/03_coded/pt_precocity-data_v1.csv"))

# trees
tr_consensus_full <- read.nexus(here("data/trees/upham2019/consensus/DNA-only/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"))

tr100 <- read.nexus(here("data/trees/upham2019/sample/2024-02-01/output.nex"))

taxa <- read_csv(here("results/taxa/upham2019_all-taxa.csv"))

# prepare data ================================================================

## trait data -----------------------------------------------------------------
# prepare named vector of trait data
prec_tipdata <- set_names(prec_data$precocity, prec_data$binomial)

## consensus tree for visualization -------------------------------------------

# clean up tip labels in the tree. The Upham tree has additional information
# after the binomial. We will replace the tiplabels with just the binomial.
tr_consensus_full$tip.label <- taxa$binomial[match(tr_consensus_full$tip.label,
                                                   taxa$tiplab)]

tr_consensus <- drop.tip(
  tr_consensus_full,
  which(!(tr_consensus_full$tip.label %in% prec_data$binomial))
)
tr_consensus <- ladderize(tr_consensus)

# set up receptacle ===========================================================
asr <- tibble(tree_id = 1:100)

asr$trees <- vector("list", nrow(asr))
for (i in seq_len(nrow(asr))) {
  asr$trees[[i]] <- tr100[[i]]
}

# function ====================================================================

run_simmap <- function(.phy, .dat) {
  message("    checking data...")
  stopifnot(all(.phy$tip.label %in% names(.dat)))
  
  # sort trait data
  .dat_sorted <- .dat[.phy$tip.label]
  
  # fit models
  message("    fitting models...")
  fit_er  <- fitMk(.phy, .dat_sorted, model = "ER",  pi = "estimated")
  fit_sym <- fitMk(.phy, .dat_sorted, model = "SYM", pi = "estimated")
  fit_ard <- fitMk(.phy, .dat_sorted, model = "ARD", pi = "estimated")
  
  # model comparison
  message("    running model comparison...")
  fit_aov <- anova(fit_er, fit_sym, fit_ard)
  
  # stochastic character mapping
  message("    running simmap....")
  scm <- simmap(fit_aov, nsim = 1000)
  message("    finished running simmap.")
  return(scm)
}

# run simmap on the consensus tree ============================================

# This is for visualization only. 
simmap_consensus <- run_simmap(tr_consensus, prec_tipdata)

write_rds(simmap_consensus, 
          here(resdir, "simmap_consensus-tree.rds"), 
          compress = "gz")

# Run simmap on the sample of trees ===========================================

# Each simmap output is ~ 300MB. For 100 trees, it will be 30 GB. Since we will
# only really use the ancestral state estimates from the sample of trees we
# will only save that part of the output (~70 KB) for all 100 trees.

get_ace_from_simmap <- function(...) {
  .simmap <- run_simmap(...)
  
  message("    summarizing simmap...")
  .simmap_summary <- summary(.simmap)
  message("    finished summarizing simmap")
  
  return(.simmap_summary$ace)
}

asr$simmap_ace <- vector("list", nrow(asr))

for (i in seq_len(nrow(asr))) {
  message("Running simmap for tree ", i)
  asr$simmap_ace[[i]] <- get_ace_from_simmap(.phy = asr$trees[[i]],
                                             .dat = prec_tipdata)
}

write_rds(asr, here(resdir, "simmap_100-trees.rds"), compress = "gz")

# end =========================================================================

