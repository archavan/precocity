# Arun Chavan
# Started: 2026-08-22

###############################################################################
# background ==================================================================
###############################################################################

# Trying to gather data for more xenarthrans and afrotherians since their tip
# states have a strong influence on the crown eutherian ancestral state.

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)

resdir <- here("data/03_coded/xenafr")
fs::dir_create(resdir)

###############################################################################
# data ========================================================================
###############################################################################

spp <- read_csv(here("data/taxa/species.csv"))

case78 <- read_csv(here("data/03_coded/case78/precocity_case78_v1.csv"))
pantheria <- read_csv(here(
  "data/03_coded/pantheria/precocity_pantheria_v1.csv"
))

already_coded <- union(case78$binomial, pantheria$binomial)

###############################################################################
# species with no data ========================================================
###############################################################################

## write species to be coded to file ==========================================
to_code <- spp |>
  filter(!binomial_upham19 %in% already_coded) |>
  filter(rank05 %in% c("Xenarthra", "Afrotheria")) |>
  filter(!is.na(binomial_upham19))

write_csv(to_code, here(resdir, "to-code.csv"))

## read manually coded data ===================================================
coded <- read_csv(here(resdir, "coded.csv")) |>
  filter(!is.na(precocity))

stopifnot(all(coded$precocity %in% c("altricial", "intermediate", "precocial")))
stopifnot(all(with(coded, binomial == binomial_upham19)))


###############################################################################
# combine data ================================================================
###############################################################################

write_csv(coded, here(resdir, "precocity_xenafr_v1.csv"))

# end #########################################################################
