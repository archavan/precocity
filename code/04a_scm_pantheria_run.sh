#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=simmap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=60G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R/4.2.0-foss-2020b

cd /home/arc78/precocity # must start r session in project root dir to use renv

Rscript --no-save code/04a_run-stochastic-character-mapping.R

