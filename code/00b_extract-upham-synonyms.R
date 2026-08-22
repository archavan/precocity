# Arun Chavan
# Started: 2026-08-22

# background ==================================================================

# Extract the curated Upham 2019 synonyms into their own lookup table.
#
# Case 1978 uses older binomials than the Upham 2019 tree. The mapping between
# them was curated by hand in the synonym_upham19 / synonym_ref columns of
# data/traits/case-1978/appendix-2.csv, which is a hard place to find it.
# This script lifts those columns into data/taxa/, alongside the other
# taxonomic lookups. appendix-2.csv remains as is -- this file is derived
# from it, so re-run this script if the curation there changes.
#
# Coverage is limited to species Case 1978 covers. It is not a general
# MSW3-to-Upham synonym map, and trait species outside Case 1978 with outdated
# names are not reconciled by it.

# setup =======================================================================
library(tidyverse)
library(here)

# data ========================================================================
appendix2 <- read_csv(here("data/traits/case-1978/appendix-2.csv"))

upham_spp <- read_csv(here("data/trees/upham2019/upham2019_all-taxa.csv"))

# extract synonyms ============================================================
synonyms <- appendix2 |>
  filter(!is.na(synonym_upham19)) |>
  select(
    binomial_synonym = binomial_case78,
    binomial_upham19 = synonym_upham19,
    ref = synonym_ref
  ) |>
  arrange(binomial_upham19)

# check =======================================================================

# every target must be a tip in the tree we analyse on
stopifnot(all(synonyms$binomial_upham19 %in% upham_spp$binomial))

# the mapping must be one-to-one in both directions
stopifnot(!any(duplicated(synonyms$binomial_synonym)))
stopifnot(!any(duplicated(synonyms$binomial_upham19)))

# write files =================================================================
write_csv(synonyms, here("data/taxa/synonyms_upham2019.csv"))

# end =========================================================================
