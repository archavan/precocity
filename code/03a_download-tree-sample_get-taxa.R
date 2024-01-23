# Arun Chavan
# Started: 2024-01-22

# background ==================================================================

# Here we generate a list of taxa from PanTHERIA for which we have precocity
# data. We will use that list to download a sample of 100 trees from upham2019
# on which we can perform all further analyses.

# setup =======================================================================
library(tidyverse)

# data ========================================================================
up19_tips <- read_csv("results/taxa/upham2019_all-taxa.csv") 
pt_coded <- read_csv("data/03_coded/pt_precocity-data_v1.csv")

# get taxa ====================================================================
pt_coded_tips <- up19_tips %>% 
  filter(binomial %in% pt_coded$binomial) %>% 
  select(binomial)

# write lists of taxa to file =================================================
write_tsv(pt_coded_tips, "results/taxa/upham2019_pt-precocity-data-taxa.tsv")

# paste this list into http://vertlife.org/phylosubsets/ to download a sample of
# trees. I used DNA-only, node-dated set of trees (with a total of 4098 species)

# end =========================================================================
