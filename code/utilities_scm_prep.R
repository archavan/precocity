# functions for preparing batched stochastic character mapping runs.

# A taxonset is a named subset of species. `prep_scm_batch()` runs every
# taxonset over every tree in a shared treelist, giving
# results/scm/<analysis>/<taxonset>/{consensus,sample/<treename>}.
# The `taxonsets_*()` functions generate taxonsets; `prep_scm_batch()` does not
# care how they were made.

# =============================================================================

#' Read the coded precocity data with ordered states
#'
#' @return tibble
#' @export
#'
#' @examples
read_precocity <- function() {
  read_csv(here("data/03_coded/precocity.csv")) |>
    mutate(
      precocity = fct(precocity, c("altricial", "intermediate", "precocial"))
    )
}

# =============================================================================

#' Build the named treelist used by all SCM analyses
#'
#' Trees are pruned once to the species in `.trait_data`. Every taxonset is a
#'   subset of those, and run_scm.R prunes further to the tipdata it is given,
#'   so one treelist serves every taxonset. Requires utilities_scm.R.
#'
#' @param .trait_data full trait table, used to prune the trees
#' @param .include_sample whether to include the 100 sampled trees
#'
#' @return multiPhylo, named, with "consensus" last
#' @export
#'
#' @examples
build_treelist <- function(.trait_data, .include_sample = TRUE) {
  tr_consensus <- read.nexus(here(
    "data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"
  ))

  if (.include_sample) {
    tr_sample <- read.nexus(here(
      "data/trees/upham2019/sample/2024-02-21/output.nex"
    ))

    treelist <- append(tr_sample, tr_consensus)
    stopifnot(length(treelist) == (length(tr_sample) + 1))
    names(treelist) <- append(names(tr_sample), "consensus")
  }

  if (!.include_sample) {
    treelist <- list(tr_consensus)
    class(treelist) <- "multiPhylo" # run_scm.R needs a list
    stopifnot(length(treelist) == 1)
    names(treelist) <- "consensus"
  }

  stopifnot(!anyDuplicated(names(treelist)))

  # This pruning is to avoid storing trees with all ~5000 tips. Pruning to the
  # actual taxonset happens later in run_scm.R
  tipdata <- set_names(.trait_data$precocity, .trait_data$binomial)
  treelist <- map(treelist, ~ prune_tree_for_fitmk(.x, tipdata))
  class(treelist) <- "multiPhylo"

  return(treelist)
}

# =============================================================================

#' Taxonsets that each drop one taxon at a given rank
#'
#' @param .trait_data full trait table
#' @param .rank_column column holding the rank to iterate over
#'
#' @return tibble with `taxonset` and `keep`
#' @export
#'
#' @examples
taxonsets_leave_one_out <- function(.trait_data, .rank_column) {
  taxa <- unique(na.omit(.trait_data[[.rank_column]]))
  rank_values <- .trait_data[[.rank_column]]

  out <- tibble(taxonset = taxa) |>
    mutate(
      rank_column = .rank_column,
      n_dropped = map_int(taxonset, ~ sum(rank_values == .x, na.rm = TRUE)),
      keep = map(taxonset, function(taxon) {
        # `%in%` over `!=` because it keeps rows whose rank is NA
        .trait_data$binomial[!(rank_values %in% taxon)]
      })
    )

  stopifnot(all(out$n_dropped + map_int(out$keep, length) == nrow(.trait_data)))

  return(out)
}

# =============================================================================

#' Taxonsets that retain a graded fraction of one clade
#'
#' With `.nested = TRUE` each smaller sample is drawn from the next larger one
#'   within a replicate, so the series is monotone and the grading is not
#'   confounded with independent sampling noise. With `.within` set, sampling is
#'   stratified on that column, holding the clade's composition fixed as it
#'   shrinks.
#'
#' @param .trait_data full trait table
#' @param .clade name of the clade to downsample
#' @param .rank_column column the clade name lives in
#' @param .fractions retained fractions of the clade
#' @param .n_replicates number of independent replicate series
#' @param .nested whether each fraction is drawn from the previous one
#' @param .within column to stratify sampling on, or NULL
#'
#' @return tibble with `taxonset` and `keep`
#' @export
#'
#' @examples
taxonsets_downsample <- function(
  .trait_data,
  .clade,
  .rank_column,
  .fractions = c(0.8, 0.6, 0.4, 0.2, 0.1),
  .n_replicates = 10,
  .nested = TRUE,
  .within = NULL
) {
  stopifnot(all(.fractions > 0), all(.fractions < 1))

  in_clade <- which(.trait_data[[.rank_column]] == .clade)
  stopifnot(length(in_clade) > 0)

  clade_species <- .trait_data$binomial[in_clade]
  other_species <- .trait_data$binomial[-in_clade]

  strata <- if (is.null(.within)) {
    set_names(rep("all", length(clade_species)), clade_species)
  } else {
    set_names(as.character(.trait_data[[.within]][in_clade]), clade_species)
  }
  n_full <- table(strata)

  draw <- function(.pool, .fraction) {
    tibble(binomial = .pool, stratum = unname(strata[.pool])) |>
      nest_by(stratum) |>
      mutate(
        kept = list({
          target <- max(1, round(.fraction * n_full[[stratum]]))
          stopifnot(target <= nrow(data))
          sample(data$binomial, size = target)
        })
      ) |>
      pull(kept) |>
      list_c()
  }

  fractions <- sort(.fractions, decreasing = TRUE)

  sample_series <- function() {
    if (.nested) {
      accumulate(fractions, draw, .init = clade_species)[-1]
    } else {
      map(fractions, ~ draw(clade_species, .x))
    }
  }

  tibble(replicate = seq_len(.n_replicates)) |>
    mutate(sampled = map(replicate, ~ sample_series())) |>
    unnest_longer(sampled, indices_to = "step") |>
    mutate(
      clade = .clade,
      fraction = fractions[step],
      taxonset = glue(
        "{.clade}_{sprintf('%03d', round(fraction * 100))}pct",
        "_rep{sprintf('%02d', replicate)}"
      ),
      keep = map(sampled, ~ c(other_species, .x))
    ) |>
    select(taxonset, clade, fraction, replicate, keep)
}

# =============================================================================

#' Write the sbatch script for one taxonset
#'
#' @param .taxonset taxonset name
#' @param .analysis analysis name
#' @param .n_trees length of the treelist, becomes the array size
#' @param .tipdata_path path to that taxonset's tipdata
#' @param .treelist_path path to the shared treelist
#' @param .logdir directory for slurm logs
#' @param .slurm named list of sbatch settings
#' @param .model model argument passed to run_scm.R
#' @param .nsim number of simmaps
#'
#' @return character
#' @export
#'
#' @examples
make_sbatch_script <- function(
  .taxonset,
  .analysis,
  .n_trees,
  .tipdata_path,
  .treelist_path,
  .logdir,
  .slurm,
  .model,
  .nsim
) {
  glue(
    r"[
#!/usr/bin/env bash
#SBATCH --array=1-{.n_trees}
#SBATCH --partition={.slurm$partition}
#SBATCH --job-name={.analysis}_{.taxonset}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={.slurm$mem_per_cpu}
#SBATCH --time={.slurm$time_limit}
#SBATCH --output={.logdir}/%A_%3a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={.slurm$mail_user}

module load R/4.4.1-foss-2022b

# to fix slow renv startup
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

cd {here()} # must start r session in project root dir to use renv

Rscript --no-save code/run_scm.R \
--analysis {.analysis}/{.taxonset} \
--tipdata {fs::path_rel(.tipdata_path, here())} \
--treelist {fs::path_rel(.treelist_path, here())} \
--treeindex $SLURM_ARRAY_TASK_ID \
--model {.model} \
--nsim {.nsim}
  ]"
  )
}

# =============================================================================

#' Prepare a batch of SCM runs, one per taxonset, over a shared set of trees
#'
#' Writes one directory per taxonset containing its tipdata, taxa table and
#'   sbatch script, one shared treelist at the analysis root, a table of the
#'   taxonsets, and a run-all.sh that submits every array job.
#'
#' Assumes run_scm.R prunes the tree it is given down to the tipdata it is
#'   given, so that all taxonsets can share one treelist.
#'
#' @param .taxonsets tibble with `taxonset` (unique, filesystem-safe) and `keep`
#'   (list-column of binomials). Any further columns are carried into
#'   taxonsets.csv as metadata.
#' @param .trait_data full trait table; needs `binomial` and `precocity`
#' @param .treelist named multiPhylo, shared by all taxonsets
#' @param .analysis analysis name, used as the results directory
#' @param .slurm named list of sbatch settings
#' @param .model model argument passed to run_scm.R
#' @param .nsim number of simmaps
#' @param .min_species smallest taxonset allowed through
#'
#' @return the taxonset table, invisibly
#' @export
#'
#' @examples
prep_scm_batch <- function(
  .taxonsets,
  .trait_data,
  .treelist,
  .analysis,
  .slurm = list(
    partition = "pi_medzhitov,day",
    mem_per_cpu = "2G",
    time_limit = "06:00:00",
    mail_user = "arun.chavan@yale.edu"
  ),
  .model = "all",
  .nsim = 1000,
  .min_species = 10
) {
  stopifnot(all(c("taxonset", "keep") %in% names(.taxonsets)))
  stopifnot(nrow(.taxonsets) > 0)
  stopifnot(!anyDuplicated(.taxonsets$taxonset))
  stopifnot(all(str_detect(.taxonsets$taxonset, "^[A-Za-z0-9._-]+$")))
  stopifnot(!is.null(names(.treelist)), !anyDuplicated(names(.treelist)))

  all_kept <- list_c(.taxonsets$keep)
  stopifnot(all(all_kept %in% .trait_data$binomial))

  tips_in_all_trees <- reduce(map(.treelist, ~ .x$tip.label), intersect)
  stopifnot(all(all_kept %in% tips_in_all_trees))

  resdir <- here("results/scm", .analysis)
  fs::dir_create(resdir)

  treelist_path <- here(resdir, "treelist.rds")
  write_rds(.treelist, treelist_path, compress = "gz")
  write_csv(
    tibble(treeindex = seq_along(.treelist), treename = names(.treelist)),
    here(resdir, "treeindices.csv")
  )

  count_state <- function(.keep, .state) {
    sum(.trait_data$precocity[.trait_data$binomial %in% .keep] == .state)
  }

  taxonset_table <- .taxonsets |>
    mutate(
      n_species = map_int(keep, length),
      n_altricial = map_int(keep, ~ count_state(.x, "altricial")),
      n_intermediate = map_int(keep, ~ count_state(.x, "intermediate")),
      n_precocial = map_int(keep, ~ count_state(.x, "precocial"))
    )

  stopifnot(all(taxonset_table$n_species >= .min_species))

  walk2(.taxonsets$taxonset, .taxonsets$keep, function(taxonset, keep_species) {
    setdir <- here(resdir, taxonset)
    fs::dir_create(setdir)
    fs::dir_create(here(setdir, "logs"))

    set_data <- .trait_data |> filter(binomial %in% keep_species)
    tipdata <- set_names(set_data$precocity, set_data$binomial)

    write_rds(tipdata, here(setdir, "tipdata.rds"), compress = "gz")
    write_csv(
      set_data |>
        select(-any_of(c("precocity", "batch_used", "coding_method"))),
      here(setdir, "taxa.csv")
    )
    write_lines(
      make_sbatch_script(
        .taxonset = taxonset,
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

  write_csv(taxonset_table |> select(-keep), here(resdir, "taxonsets.csv"))
  write_lines(
    glue(r"[sbatch {here(resdir, .taxonsets$taxonset, "run.sh")}]"),
    here(resdir, "run-all.sh")
  )

  message(glue(
    "{nrow(.taxonsets)} taxonsets x {length(.treelist)} trees = \\
     {nrow(.taxonsets) * length(.treelist)} jobs"
  ))

  invisible(taxonset_table)
}

###############################################################################
