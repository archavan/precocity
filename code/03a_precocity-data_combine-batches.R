# Arun Chavan
# Started: 2024-11-21

###############################################################################
# background ==================================================================
###############################################################################

# Combine all batches of precocity coding into one master dataset.

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)

###############################################################################
# data ========================================================================
###############################################################################
spp <- read_csv(here("data/taxa/species.csv"))

case78 <- read_csv(here("data/03_coded/case78/precocity_case78_v1.csv"))
pantheria <- read_csv(here(
  "data/03_coded/pantheria/precocity_pantheria_v1.csv"
))
xenafr <- read_csv(here("data/03_coded/xenafr/precocity_xenafr_v1.csv"))

case78 <- case78 |>
  mutate(
    precocity = case_when(
      precocity == "fetal" ~ "altricial",
      precocity == "A" ~ "altricial",
      precocity == "SA" ~ "altricial",
      precocity == "SP" ~ "intermediate",
      precocity == "P" ~ "precocial"
    )
  ) |>
  mutate(coding_method = "inherited")

pantheria <- pantheria |>
  select(
    rank01,
    rank02,
    rank03,
    rank04,
    rank05,
    rank06,
    rank07,
    rank08,
    family,
    binomial,
    precocity,
    coding_method
  )

xenafr <- xenafr |>
  select(-source, -notes, -binomial_upham19) |>
  mutate(coding_method = "manual")


###############################################################################
# combined ####################################################################
###############################################################################
preference <- c("xenafr", "case78", "pantheria")

full <- lst(xenafr, case78, pantheria) |>
  bind_rows(.id = "batch")

conflicting_coding <- full |>
  distinct(binomial, precocity) |>
  filter(n() > 1, .by = binomial) |>
  arrange(binomial)

print(conflicting_coding)

# combine by prefering the batch in the order of preference if a species appears
# in multiple batches.
combined <- full |>
  mutate(batch = fct(batch, preference)) |>
  group_by(binomial) |>
  arrange(batch, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  relocate(batch, .after = precocity) |>
  mutate(batch = as.character(batch)) |>
  rename(batch_used = batch)


stopifnot(
  !anyNA(combined$precocity),
  all(combined$precocity %in% c("altricial", "intermediate", "precocial")),
  !any(duplicated(combined$binomial)),
  all(combined$binomial %in% spp$binomial_upham19)
)

###############################################################################
# write data ##################################################################
###############################################################################

arrange_taxa <- function(x) {
  x |>
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
}

arrange_taxa(combined) |>
  write_csv(here("data/03_coded/precocity.csv"))

# end =========================================================================
