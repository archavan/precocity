# Arun Chavan
# Started: 2024-02-25

# background ==================================================================
# precocity data for species in the intersection of case78 and panteria datasets

# setup =======================================================================
library(tidyverse)
library(here)

# data ========================================================================
pt <- read_csv(here("data/03_coded/pantheria/precocity_pantheria_v1.csv"))
case78 <- read_csv(here("data/03_coded/case78/precocity_case78_v1.csv"))

pt <- pt %>% 
  select(-adult_body_mass_g, 
         -gestation_len_d,
         -litter_size,
         -age_at_eye_opening_d,
         -weaning_age_d,
         -neonate_body_mass_g,
         -litter_size_discrete,
         -age_at_eye_opening_discrete,
         -coding_method)

case78 <- case78 %>% 
  mutate(precocity = case_when(
    precocity == "fetal" ~ "altricial",
    precocity == "A" ~ "altricial",
    precocity %in% c("SA", "SP") ~ "intermediate",
    precocity == "P" ~ "precocial"
  ))

shared <- inner_join(pt, 
                     case78, 
                     by = c("rank01", "rank02", "rank03", "rank04", "rank05",
                            "rank06", "rank07", "rank08", 
                            "family", "binomial"),
                     suffix = c("_pt", "_case78"))

# write data ==================================================================
write_csv(shared, here("data/03_coded/pt_case78_intersect-spp.csv"))

# end =========================================================================
