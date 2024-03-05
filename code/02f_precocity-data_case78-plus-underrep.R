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

case78_tax_cvg <- read_csv(here("data/03_coded/case78/taxonomic-coverage.csv"))
spp_underrep <- read_csv(here("data/03_coded/case78/spp-from-underrepresented-orders.csv"))

# get data for underrep spp
underrep_spp_data_from_pt <- pt %>% 
  filter(rank07 %in% case78_tax_cvg$order[case78_tax_cvg$enrichment < 0.5]) %>% 
  filter(!binomial %in% case78$binomial)

underrep_spp_data_from_pt <- underrep_spp_data_from_pt %>% 
  select(rank01, rank02, rank03, rank04, rank05, rank06, rank07, rank08, 
         family, binomial, precocity)


case78_plus_underrep <- full_join(case78, underrep_spp_data_from_pt)

# write =======================================================================
write_csv(case78_plus_underrep, here("data/03_coded/case78-plus-underrep-taxa.csv"))

# end =========================================================================
