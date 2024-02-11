# Arun Chavan
# Started: 2024-02-11

# background ==================================================================

# Organize taxonomic information:
# Some of the higher order taxa in pantheria are old, and now known not to be
# monophyletic groups. We need to get higher order clades from upham 2019 and
# merge that information with pantheria dataset so we have clade-based, updated,
# taxonomic information alongside the pantheria dataset.

# setup =======================================================================
library(tidyverse)
library(here)

# data ========================================================================
clades <- read_csv(here("data/taxa/higher-order-clades.csv"))
msw <- read_csv(here("data/taxa/msw3-all.csv"), 
                col_select = all_of(c("Order", 
                                      "Suborder", 
                                      "Infraorder", 
                                      "Superfamily", 
                                      "Family",
                                      "TaxonLevel",
                                      "Extinct?"))) %>% 
  janitor::clean_names() %>% 
  filter(taxon_level == "FAMILY") %>% 
  filter(!extinct) %>% 
  distinct() %>% 
  select(-taxon_level, -extinct)

# convert all columns to sentence case for ease of regexing
clades <- mutate(clades, across(everything(), str_to_sentence))
msw <- mutate(msw, across(everything(), str_to_sentence))

# add families to higer order clades from MSW data ============================
get_taxonomic_rank <- function(.taxon) {
  txrank <- names(which(apply(msw, 2, function(x) any(grepl(.taxon, x)))))
  stopifnot(length(txrank) == 1)
  return(txrank)
}

get_families <- function(.taxon) {
  txrank <- get_taxonomic_rank(.taxon)
  msw$family[which(msw[[txrank]] == .taxon)]
}

clades$family <- vector("list", nrow(clades))

# First fill in families for which we have a more specific taxonomic rank than
# order that we want to use.
clades$family[[which(clades$rank08 == "Vermilingua")]]     <- get_families("Vermilingua")
clades$family[[which(clades$rank08 == "Folivora")]]        <- get_families("Folivora")
clades$family[[which(clades$rank08 == "Caniformia")]]      <- get_families("Caniformia")
clades$family[[which(clades$rank08 == "Feliformia")]]      <- get_families("Feliformia")
clades$family[[which(clades$rank08 == "Castorimorpha")]]   <- get_families("Castorimorpha")
clades$family[[which(clades$rank08 == "Myomorpha")]]       <- get_families("Myomorpha")
clades$family[[which(clades$rank08 == "Anomaluromorpha")]] <- get_families("Anomaluromorpha")
clades$family[[which(clades$rank08 == "Hystricomorpha")]]  <- get_families("Hystricomorpha")
clades$family[[which(clades$rank08 == "Sciuromorpha")]]    <- get_families("Sciuromorpha")

# There are two mode clades from Chiroptera for which we want to use the more
# specific taxon than order, but these clades are not recognized in MSW. Add
# their families manually.

clades$family[[which(clades$rank08 == "Yangochiroptera")]] <- c("Emballonuridae",
                                                                "Furipteridae",
                                                                "Miniopteridae",
                                                                "Molossidae",
                                                                "Mormoopidae",
                                                                "Mystacinidae",
                                                                "Myzopodidae",
                                                                "Natalidae",
                                                                "Noctilionidae",
                                                                "Nycteridae",
                                                                "Phyllostomidae",
                                                                "Thyropteridae",
                                                                "Cistugidae",
                                                                "Vespertilionidae")
clades$family[[which(clades$rank08 == "Yinpterochiroptera")]] <- c("Craseonycteridae",
                                                                   "Hipposideridae",
                                                                   "Megadermatidae",
                                                                   "Pteropodidae",
                                                                   "Rhinolophidae",
                                                                   "Rhinopomatidae")

# Check if all rank07 (orders) are present in MSW
unique(clades$rank07[which(!is.na(clades$rank07))]) %>% 
  set_names(.) %>% 
  map_chr(possibly(get_taxonomic_rank, NA_character_)) %>% 
  is.na() %>% 
  which()
#> Eulipotyphla and Cetartiodactyla

# We need to manually add families for Eulipotyphla and Catardtiodactyla
clades$family[[which(clades$rank07 == "Eulipotyphla")]] <- c("Erinaceidae",
                                                             "Soricidae",
                                                             "Talpidae",
                                                             "Solenodontidae")

clades$family[[which(clades$rank07 == "Cetartiodactyla")]] <- c(get_families("Artiodactyla"),
                                                                get_families("Cetacea"))


# For everything else, get all families at the order level
for (i in seq_len(nrow(clades))) {
  if (is.null(clades$family[[i]])) {
    clades$family[[i]] <- get_families(clades$rank07[[i]])
  } else {
    clades$family[[i]] <- clades$family[[i]]
  }
}


# unnest ======================================================================
families <- unnest_longer(clades, family)

# check that we have all families =============================================

# We have almost all families covered now. Except for Neophontidae but now sure
# about its status
anti_join(msw, families, by = "family")
#> # A tibble: 1 × 5
#>   order        suborder infraorder superfamily family       
#> 1 Soricomorpha <NA>     <NA>       <NA>        Nesophontidae

anti_join(families, msw, by = "family")
#> # A tibble: 2 × 9
#>   rank01   rank02 rank03   rank04        rank05         rank06 rank07     rank08          family       
#> 1 Mammalia Theria Eutheria Boreoeutheria Laurasiatheria <NA>   Chiroptera Yangochiroptera Miniopteridae
#> 2 Mammalia Theria Eutheria Boreoeutheria Laurasiatheria <NA>   Chiroptera Yangochiroptera Cistugidae

# write files =================================================================
write_csv(families, here("data/taxa/families.csv"))

# end =========================================================================

