#!/usr/bin/env bash
#SBATCH --array=1-100
#SBATCH --partition=day
#SBATCH --job-name=scm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R/4.2.0-foss-2020b

cd /home/arc78/precocity # must start r session in project root dir to use renv

Rscript --no-save code/09a_scm_pantheria-plus-underrep-ER_run_sample.R $SLURM_ARRAY_TASK_ID

