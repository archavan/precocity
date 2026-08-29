# Arun Chavan
# Started: 2026-08-28

###############################################################################
# background ==================================================================
###############################################################################

# Redo the algorithmic PanTHERIA precocity coding over a 3x3 grid of the two
# upper cut-points (age at eye opening, litter size) and run SCM on each cell,
# over the consensus tree and all 100 sampled trees. Only the rows coded
# algorithmically at baseline are recomputed per cell; manual and inherited
# (case78/xenafr) rows remain unchanged.

###############################################################################
# setup =======================================================================
###############################################################################
library(phytools)
library(tidyverse)
library(here)
library(glue)
library(ape)
library(conflicted)

conflicts_prefer(dplyr::filter, purrr::map)

source(here("code/utilities_scm.R"))
source(here("code/utilities_scm_prep.R"))

## parameters =================================================================
analysis <- "threshold-grid"
model_to_test <- "all"
nsim <- 1000

grid <- expand_grid(eye_upper = c(3, 5, 7), litter_upper = c(2, 3, 4)) |>
  mutate(cell = glue("eye{eye_upper}_litter{litter_upper}"))

###############################################################################
# data ========================================================================
###############################################################################

precocity <- read_csv(here("data/03_coded/precocity.csv"))
pantheria <- read_csv(here(
  "data/03_coded/pantheria/precocity_pantheria_v1.csv"
))

###############################################################################
# recompute algorithmic rows for one grid cell ================================
###############################################################################

#' Recode the algorithmically-coded rows of the precocity table at a given
#'   pair of upper cut-points, leaving manual/inherited rows untouched
#'
#' @param .precocity full precocity table
#' @param .pantheria pantheria table with litter_size, age_at_eye_opening_d
#' @param .eye_upper upper cut-point for age at eye opening (days)
#' @param .litter_upper upper cut-point for litter size
#' @param .eye_lower lower cut-point for age at eye opening, held fixed
#' @param .litter_lower lower cut-point for litter size, held fixed
#'
#' @return .precocity with algorithmic rows recoded, unresolved rows dropped
#' @export
#'
#' @examples
recode_algorithmic_rows <- function(
  .precocity,
  .pantheria,
  .eye_upper,
  .litter_upper,
  .eye_lower = 0,
  .litter_lower = 1.5
) {
  stopifnot(.eye_upper > .eye_lower, .litter_upper > .litter_lower)

  algorithmic <- .precocity |>
    filter(coding_method == "algorithmic") |>
    select(binomial) |>
    left_join(
      select(.pantheria, binomial, litter_size, age_at_eye_opening_d),
      by = "binomial"
    )

  recoded <- algorithmic |>
    mutate(
      litter_size_discrete = case_when(
        litter_size < .litter_lower ~ "low",
        litter_size >= .litter_lower & litter_size < .litter_upper ~ "medium",
        litter_size >= .litter_upper ~ "high"
      ),
      age_at_eye_opening_discrete = case_when(
        age_at_eye_opening_d == .eye_lower ~ "low",
        age_at_eye_opening_d > .eye_lower &
          age_at_eye_opening_d < .eye_upper ~ "medium",
        age_at_eye_opening_d >= .eye_upper ~ "high"
      )
    ) |>
    mutate(
      precocity = case_when(
        age_at_eye_opening_discrete == "high" ~ "altricial",
        age_at_eye_opening_discrete == "low" &
          litter_size_discrete == "low" ~ "precocial",
        age_at_eye_opening_discrete == "medium" &
          litter_size_discrete == "medium" ~ "intermediate",
        .default = NA_character_
      )
    )

  recoded <- recoded |>
    select(binomial, precocity) |>
    filter(!is.na(precocity))

  non_algorithmic <- .precocity |> filter(coding_method != "algorithmic")

  recoded_full <- .precocity |>
    select(-precocity) |>
    inner_join(recoded, by = "binomial")

  bind_rows(non_algorithmic, recoded_full)
}

###############################################################################
# baseline fidelity check ======================================================
###############################################################################

baseline <- recode_algorithmic_rows(
  precocity,
  pantheria,
  .eye_upper = 5,
  .litter_upper = 3
)

fidelity <- precocity |>
  select(binomial, precocity_baseline = precocity) |>
  full_join(
    select(baseline, binomial, precocity_reproduced = precocity),
    by = "binomial"
  )

n_mismatch <- sum(
  fidelity$precocity_baseline != fidelity$precocity_reproduced,
  na.rm = TRUE
) +
  sum(xor(
    is.na(fidelity$precocity_baseline),
    is.na(fidelity$precocity_reproduced)
  ))

stopifnot(n_mismatch == 0)

###############################################################################
# build the 3 x 3 grid of trait tables ========================================
###############################################################################

cell_data <- map2(grid$eye_upper, grid$litter_upper, function(x, y) {
  recode_algorithmic_rows(
    precocity,
    pantheria,
    .eye_upper = x,
    .litter_upper = y
  ) |>
    mutate(
      precocity = fct(precocity, c("altricial", "intermediate", "precocial"))
    )
})
names(cell_data) <- grid$cell

taxonsets <- tibble(taxonset = grid$cell, trait_data = cell_data) |>
  left_join(grid |> rename(taxonset = cell), by = "taxonset")

###############################################################################
# shared treelist ==============================================================
###############################################################################

all_binomials <- reduce(map(cell_data, "binomial"), union)
union_trait_data <- precocity |>
  filter(binomial %in% all_binomials) |>
  select(binomial) |>
  left_join(
    bind_rows(cell_data) |> distinct(binomial, precocity),
    by = "binomial"
  )

treelist <- build_treelist(union_trait_data)

###############################################################################
# write one directory per grid cell ============================================
###############################################################################

#' Write one SCM run directory per grid cell, sharing one treelist
#'
#' @param .taxonsets tibble with `taxonset`, `trait_data` (full per-cell
#'   table with binomial/precocity), and grid metadata columns
#' @param .treelist named multiPhylo, shared by all cells
#' @param .analysis analysis name, used as the results directory
#' @param .slurm named list of sbatch settings
#' @param .model model argument passed to run_scm.R
#' @param .nsim number of simmaps
#'
#' @return the taxonset summary table, invisibly
#' @export
#'
#' @examples
prep_scm_grid <- function(
  .taxonsets,
  .treelist,
  .analysis,
  .slurm = list(
    partition = "pi_medzhitov,day",
    mem_per_cpu = "2G",
    time_limit = "06:00:00",
    mail_user = "arun.chavan@yale.edu"
  ),
  .model = "all",
  .nsim = 1000
) {
  resdir <- here("results/scm", .analysis)
  fs::dir_create(resdir)

  treelist_path <- here(resdir, "treelist.rds")
  write_rds(.treelist, treelist_path, compress = "gz")
  write_csv(
    tibble(treeindex = seq_along(.treelist), treename = names(.treelist)),
    here(resdir, "treeindices.csv")
  )

  count_state <- function(.dat, .state) sum(.dat$precocity == .state)

  taxonset_table <- .taxonsets |>
    mutate(
      n_species = map_int(trait_data, nrow),
      n_altricial = map_int(trait_data, count_state, .state = "altricial"),
      n_intermediate = map_int(
        trait_data,
        count_state,
        .state = "intermediate"
      ),
      n_precocial = map_int(trait_data, count_state, .state = "precocial")
    )

  map2(.taxonsets$taxonset, .taxonsets$trait_data, function(x, y) {
    setdir <- here(resdir, x)
    fs::dir_create(setdir)
    fs::dir_create(here(setdir, "logs"))

    tipdata <- set_names(y$precocity, y$binomial)

    write_rds(tipdata, here(setdir, "tipdata.rds"), compress = "gz")
    write_csv(
      y |> select(-any_of(c("precocity", "batch_used", "coding_method"))),
      here(setdir, "taxa.csv")
    )
    write_lines(
      make_sbatch_script(
        .taxonset = x,
        .analysis = .analysis,
        .n_trees = length(.treelist),
        .tipdata_path = here(setdir, "tipdata.rds"),
        .treelist_path = treelist_path,
        .logdir = here(setdir, "logs"),
        .slurm = .slurm,
        .model = .model,
        .nsim = .nsim
      ),
      here(setdir, "run.sh")
    )
  })

  write_csv(
    taxonset_table |> select(-trait_data),
    here(resdir, "taxonsets.csv")
  )
  write_lines(
    glue(r"[sbatch {here(resdir, .taxonsets$taxonset, "run.sh")}]"),
    here(resdir, "run-all.sh")
  )

  message(glue(
    "{nrow(.taxonsets)} grid cells x {length(.treelist)} trees = \\
     {nrow(.taxonsets) * length(.treelist)} jobs"
  ))

  invisible(taxonset_table)
}

###############################################################################
# state changes vs baseline, for the taxonsets table ===========================
###############################################################################

baseline_states <- baseline |> select(binomial, baseline_precocity = precocity)

taxonsets <- taxonsets |>
  mutate(
    n_pt_dropped = map_int(trait_data, function(x) {
      length(setdiff(baseline$binomial, x$binomial))
    }),
    n_changed_state_vs_baseline = map_int(trait_data, function(x) {
      x |>
        select(binomial, precocity) |>
        inner_join(baseline_states, by = "binomial") |>
        summarize(n = sum(as.character(precocity) != baseline_precocity)) |>
        pull(n)
    })
  )

###############################################################################
# prepare for analysis =========================================================
###############################################################################

prep_scm_grid(
  taxonsets,
  .treelist = treelist,
  .analysis = analysis,
  .model = model_to_test,
  .nsim = nsim
)

###############################################################################
# per-species diagnostics: dropped and state-changed species per cell =========
###############################################################################

dropped_species <- map2(
  taxonsets$taxonset,
  taxonsets$trait_data,
  function(x, y) {
    tibble(taxonset = x, binomial = setdiff(baseline$binomial, y$binomial))
  }
) |>
  bind_rows() |>
  left_join(
    select(pantheria, binomial, litter_size, age_at_eye_opening_d),
    by = "binomial"
  ) |>
  left_join(baseline_states, by = "binomial")

changed_state_species <- map2(
  taxonsets$taxonset,
  taxonsets$trait_data,
  function(x, y) {
    y |>
      select(binomial, precocity) |>
      inner_join(baseline_states, by = "binomial") |>
      filter(as.character(precocity) != baseline_precocity) |>
      mutate(taxonset = x)
  }
) |>
  bind_rows() |>
  rename(cell_precocity = precocity) |>
  left_join(
    select(pantheria, binomial, litter_size, age_at_eye_opening_d),
    by = "binomial"
  ) |>
  relocate(taxonset)

write_csv(dropped_species, here("results/scm", analysis, "dropped-species.csv"))
write_csv(
  changed_state_species,
  here("results/scm", analysis, "changed-state-species.csv")
)

sessionInfo()

# end =========================================================================
