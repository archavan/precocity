# Arun Chavan
# Started: 2023-11-20

# background ==================================================================

# Here we get lists of taxa from the two mammal trees, dosReis2012 and upham2019

# setup =======================================================================
library(tidyverse)
library(here)
library(ape)

# data ========================================================================
dr12 <- read.nexus("data/trees/dosReis2012/mammalia_dosReis.tree")
up19 <- read.nexus("data/trees/upham2019/consensus/DNA-only/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre")

# get taxa ====================================================================
dr12_tips <- enframe(dr12$tip.label, name = NULL, value = "binomial")

up19_tips <- enframe(up19$tip.label, name = NULL, value = "tiplab")

up19_tips <- up19_tips %>% 
  mutate(binomial = str_replace(tiplab, ".*([A-Z][a-z]+_[a-z]+).*", "\\1"))

# write lists of taxa to file =================================================
write_csv(dr12_tips, "data/trees/dosReis2012/dosReis2012_all-taxa.csv")
write_csv(up19_tips, "data/trees/upham2019/upham2019_all-taxa.csv")

# end =========================================================================
