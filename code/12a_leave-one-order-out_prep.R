# Arun Chavan
# Started: 2024-11-26

###############################################################################
# background ==================================================================
###############################################################################

# Run stochastic character mapping on consensus from Upham et al. while leaving
# one order out at a time to test the effect of taxon sampling systematically.

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

## parameters =================================================================
analysis <- "looo"
model_to_test <- "all"
nsim <- 1000
partition <- "pi_medzhitov,day"
mem_per_cpu <- "2G"
time_limit <- "06:00:00"

resdir <- here("results/scm", analysis)
fs::dir_create(resdir)

###############################################################################
# data ========================================================================
###############################################################################

prec_data <- read_csv(here("data/03_coded/precocity.csv"))

prec_data <- prec_data %>%
  mutate(
    precocity = fct(precocity, c("altricial", "intermediate", "precocial"))
  )

## trees ======================================================================
tr <- read.nexus(here(
  "data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"
))

###############################################################################
# prepare for analysis ========================================================
###############################################################################

# We have to adapt the run_scm.R architecture, which was written to run the
# analysis on same tipdata and taxa with a sample of 100 trees, to run just the
# consenus tree but with different combination of taxa and tipdata. Hence the
# gymnastics below.

## generate receptacle for results ============================================
orders <- prec_data$rank07 |> unique()

df <- tibble(leave_out = orders)

df <- df |>
  mutate(prec = map(leave_out, ~ filter(prec_data, rank07 != .x))) |>
  mutate(tipdata = map(prec, ~ set_names(.x$precocity, .x$binomial))) |>
  mutate(tr_pruned = map(tipdata, ~ prune_tree_for_fitmk(tr, .x))) |>
  mutate(
    tr_pruned = map(tr_pruned, function(x) {
      x <- list(x)
      class(x) <- "multiPhylo" # run_scm.R needs a list
      stopifnot(length(x) == 1)
      names(x) <- "consensus"
      return(x)
    })
  )

for (i in seq_len(nrow(df))) {
  fs::dir_create(here(resdir, df$leave_out[[i]]))
  fs::dir_create(here(resdir, df$leave_out[[i]], "logs"))
  write_rds(df$tipdata[[i]], here(resdir, df$leave_out[[i]], "tipdata.rds"))
  write_rds(df$tr_pruned[[i]], here(resdir, df$leave_out[[i]], "treelist.rds"))
  write_csv(
    df$prec[[i]] |> select(-precocity, -batch_used, -coding_method),
    here(resdir, df$leave_out[[i]], "taxa.csv")
  )
}

## generate bash scripts ======================================================

make_bash <- function(.leave_out) {
  glue(
    r"[
#!/usr/bin/env bash
#SBATCH --array=1-{length(df$tr_pruned[[which(df$leave_out == .leave_out)]])}
#SBATCH --partition={partition}
#SBATCH --job-name={.leave_out}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={time_limit}
#SBATCH --output={here(resdir, .leave_out, "logs")}/%A_%3a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R/4.4.1-foss-2022b

# to fix slow renv startup
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

cd {here()} # must start r session in project root dir to use renv

Rscript --no-save code/run_scm.R \
--analysis {analysis}/{.leave_out} \
--tipdata {fs::path_rel(here(resdir, .leave_out, "tipdata.rds"), here())} \
--treelist {fs::path_rel(here(resdir, .leave_out, "treelist.rds"), here())} \
--treeindex $SLURM_ARRAY_TASK_ID \
--model {model_to_test} \
--nsim {nsim}
  ]"
  )
}

df <- df |>
  mutate(bash = map_chr(leave_out, make_bash))

# df$bash[[1]] |> writeLines()

walk2(df$leave_out, df$bash, ~ write_lines(.y, here(resdir, .x, "run.sh")))

# script to launch all scripts
glue(r"[sbatch {here(resdir, df$leave_out, "run.sh")}]") |>
  write_lines(here(resdir, "run-all.sh"))

# end =========================================================================
