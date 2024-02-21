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

# clean up tip labels in the tree. The Upham tree has additional information
# after the binomial. We will replace the tiplabels with just the binomial.
up19_tips <- up19_tips %>% 
  mutate(binomial = str_replace(tiplab, ".*([A-Z][a-z]+_[a-z]+).*", "\\1"))

tr_consensus <- up19
tr_consensus$tip.label <- up19_tips$binomial[match(tr_consensus$tip.label,
                                                   up19_tips$tiplab)]

stopifnot(all(map2_lgl(tr_consensus$tip.label, up19$tip.label, ~ str_detect(.y, .x))))

# write lists of taxa to file =================================================
write_csv(dr12_tips, "data/trees/dosReis2012/dosReis2012_all-taxa.csv")
write_csv(up19_tips, "data/trees/upham2019/upham2019_all-taxa.csv")

write.nexus(tr_consensus, 
            file = "data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree")

# end =========================================================================

