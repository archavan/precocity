# Arun Chavan
# Started: 2026-08-21

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(ape)
library(conflicted)

conflicts_prefer(purrr::map)

source(here("code/utilities_scm.R"))

## parameters =================================================================
analysis <- "case78-plus-pantheria"
model_to_test <- "all"
nsim <- 1000
partition <- "pi_medzhitov,day"
mem_per_cpu <- "2G"
time_limit <- "06:00:00"

resdir <- here("results/scm", analysis)
fs::dir_create(resdir)

logdir <- here(resdir, "logs")
fs::dir_create(logdir)

###############################################################################
# data ========================================================================
###############################################################################

## tipdata ====================================================================

prec <- read_csv(here("data/03_coded/precocity.csv"))
prec <- prec |> filter(batch_used %in% c("case78", "pantheria"))
prec <- prec |>
  mutate(
    precocity = fct(precocity, c("altricial", "intermediate", "precocial"))
  )

tipdata <- set_names(prec$precocity, prec$binomial)

## trees ======================================================================
tr_sample <- read.nexus(here(
  "data/trees/upham2019/sample/2024-02-21/output.nex"
))
tr_consensus <- read.nexus(here(
  "data/trees/upham2019/consensus/DNA-only/consensus_full_with-binomial-tiplabels.tree"
))


treelist <- append(tr_sample, tr_consensus)
stopifnot(length(treelist) == (length(tr_sample) + 1))
names(treelist) <- append(names(tr_sample), "consensus")
stopifnot(!anyDuplicated(names(treelist)))

treelist <- map(treelist, ~ prune_tree_for_fitmk(.x, tipdata))
class(treelist) <- "multiPhylo"

treeindices <- tibble(
  treeindex = seq_along(treelist),
  treename = names(treelist)
)


###############################################################################
# write data ==================================================================
###############################################################################

tipdata_path <- here(resdir, "tipdata.rds")
treelist_path <- here(resdir, "treelist.rds")

write_rds(tipdata, tipdata_path, compress = "gz")
write_rds(treelist, treelist_path, compress = "gz")
write_csv(treeindices, here(resdir, "treeindices.csv"))

###############################################################################
# write bash script ===========================================================
###############################################################################

bash <- glue(
  r"[
#!/usr/bin/env bash
#SBATCH --array=1-{length(treelist)}
#SBATCH --partition={partition}
#SBATCH --job-name=scm_{analysis}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={time_limit}
#SBATCH --output={logdir}/%A_%3a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R/4.4.1-foss-2022b

# to fix slow renv startup
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

cd {here()} # must start r session in project root dir to use renv

Rscript --no-save code/run_scm.R \
--analysis {analysis} \
--tipdata {fs::path_rel(tipdata_path, here())} \
--treelist {fs::path_rel(treelist_path, here())} \
--treeindex $SLURM_ARRAY_TASK_ID \
--model {model_to_test} \
--nsim {nsim}
  ]"
)


write_lines(bash, here(resdir, "run.sh"))

# end #########################################################################
