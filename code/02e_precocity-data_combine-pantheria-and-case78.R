# Arun Chavan
# Started: 2024-11-21

# background ==================================================================

# Combined pantheria and case78 coded data in the following ways:

# 1. intersect
# 2. pantheria + underrepresented taxa from Case78
# 3. case78 + underrepresented taxa from pantheria
# 4. case78 + pantheria (not stricly a union because we will prefer case78 coding if both exist)

# setup =======================================================================
library(tidyverse)
library(here)

# data ========================================================================
case78 <- read_csv("data/03_coded/case78/precocity_case78_v1.csv")
pt <- read_csv("data/03_coded/pantheria/precocity_pantheria_v1.csv")

case78_tax_cvg <- read_csv(here("data/03_coded/case78/taxonomic-coverage.csv"))
pt_tax_cvg <- read_csv(here("data/03_coded/pantheria/taxonomic-coverage.csv"))

case78 <- case78 %>% 
  mutate(precocity = case_when(
    precocity == "fetal" ~ "altricial",
    precocity == "A" ~ "altricial",
    precocity == "SA" ~ "altricial",
    precocity == "SP" ~ "intermediate",
    precocity == "P" ~ "precocial"
  ))

pt <- pt %>% 
  select(rank01, rank02, rank03, rank04, rank05, rank06, rank07, rank08, 
         family, binomial, precocity)

###############################################################################
# intersect ###################################################################
###############################################################################

shared <- inner_join(pt, 
                     case78, 
                     by = c("rank01", "rank02", "rank03", "rank04", "rank05",
                            "rank06", "rank07", "rank08", 
                            "family", "binomial"),
                     suffix = c("_pt", "_case78"))


###############################################################################
# pantheria + underrep taxa from case78 #######################################
###############################################################################

pt_underrep_spp_data_from_case78 <- case78 %>% 
  filter(rank07 %in% pt_tax_cvg$order[pt_tax_cvg$enrichment < 0.5]) %>% 
  filter(!binomial %in% pt$binomial)

pt_plus_underrep <- full_join(pt, pt_underrep_spp_data_from_case78)

###############################################################################
# case78 + underrep taxa from pantheria #######################################
###############################################################################

case78_underrep_spp_data_from_pt <- pt %>% 
  filter(rank07 %in% case78_tax_cvg$order[case78_tax_cvg$enrichment < 0.5]) %>% 
  filter(!binomial %in% case78$binomial)

case78_plus_underrep <- full_join(case78, case78_underrep_spp_data_from_pt)

###############################################################################
# case78 + pantheria ##########################################################
###############################################################################

case78_plus_pt <- full_join(pt, 
                            case78, 
                            by = c("rank01", "rank02", "rank03", "rank04", "rank05",
                                   "rank06", "rank07", "rank08", 
                                   "family", "binomial"),
                            suffix = c("_pt", "_case78")) %>% 
  mutate(precocity = coalesce(precocity_case78, precocity_pt)) %>% # prefer case78
  select(-precocity_pt, 
         -precocity_case78) %>% 
  distinct()

###############################################################################
# write data ##################################################################
###############################################################################

arrange_taxa <- function(x) {
  x %>% 
    arrange(rank01, rank02, rank03, rank04, rank05, rank06, rank07, rank08,
            family, binomial)
}

write_csv(arrange_taxa(shared), here("data/03_coded/pt_case78_intersect-spp.csv"))
write_csv(arrange_taxa(pt_plus_underrep), here("data/03_coded/pantheria-plus-underrep-taxa.csv"))
write_csv(arrange_taxa(case78_plus_underrep), here("data/03_coded/case78-plus-underrep-taxa.csv"))
write_csv(arrange_taxa(case78_plus_pt), here("data/03_coded/case78-plus-pantheria.csv"))

# end =========================================================================
