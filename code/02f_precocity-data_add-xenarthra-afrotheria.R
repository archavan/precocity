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

resdir <- here("data/03_coded/xen-afr")
fs::dir_create(resdir)

###############################################################################
# data ========================================================================
###############################################################################

prec <- read_csv(here("data/03_coded/case78-plus-pantheria.csv"))
spp <- read_csv(here("data/taxa/species.csv"))

###############################################################################
# species with no data ========================================================
###############################################################################

## write species to be coded to file ==========================================
to_code <- spp |>
  filter(!binomial_upham19 %in% prec$binomial) |>
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

full <- bind_rows(
  prec,
  select(coded, -source, -notes, -binomial_upham19)
) |>
  arrange(
    rank01,
    rank02,
    rank03,
    rank04,
    rank05,
    rank06,
    rank07,
    rank08,
    family,
    binomial
  )

stopifnot(all(full$binomial %in% spp$binomial_upham19))

write_csv(full, here("data/03_coded/case78-plus-pantheria-plus-xenafr.csv"))

# end #########################################################################
