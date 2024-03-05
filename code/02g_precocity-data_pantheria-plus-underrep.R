# Arun Chavan
# Started: 2024-02-25

# background ==================================================================

# precocity data for case78 supplemented with data from underrepresented taxa

# setup =======================================================================
library(tidyverse)
library(here)

# data ========================================================================
case78 <- read_csv("data/03_coded/case78/precocity_case78_v1.csv")
pt <- read_csv("data/03_coded/pantheria/precocity_pantheria_v1.csv")

case78 <- case78 %>% 
  mutate(precocity = case_when(
    precocity == "fetal" ~ "altricial",
    precocity == "A" ~ "altricial",
    precocity == "SA" ~ "altricial",
    precocity == "SP" ~ "intermediate",
    precocity == "P" ~ "precocial"
  ))

pt_tax_cvg <- read_csv(here("data/03_coded/pantheria/taxonomic-coverage.csv"))
spp_underrep <- read_csv(here("data/03_coded/pantheria/spp-from-underrepresented-orders.csv"))

# get data for underrep spp
underrep_spp_data_from_case78 <- case78 %>% 
  filter(rank07 %in% pt_tax_cvg$order[pt_tax_cvg$enrichment < 0.5]) %>% 
  filter(!binomial %in% pt$binomial)

pt <- pt %>% 
  select(rank01, rank02, rank03, rank04, rank05, rank06, rank07, rank08,
         family, binomial, precocity)

pt_plus_underrep <- full_join(pt, underrep_spp_data_from_case78)

# write =======================================================================
write_csv(pt_plus_underrep, here("data/03_coded/pantheria-plus-underrep-taxa.csv"))

# end =========================================================================
