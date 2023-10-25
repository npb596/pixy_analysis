#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=15gb
#SBATCH -J pegas_run
#SBATCH -o pegas_run-%j.out
#SBATCH -t 50:00:00
#SBATCH -p jro0014_amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=npb0015@auburn.edu

module load R/4.2.2

Rscript 01_run_pegas.R
