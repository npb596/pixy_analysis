#!/bin/bash

#SBATCH -J missing_genos
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=8G
#SBATCH -t 10:00:00
#SBATCH -p normal
#SBATCH -e missing_genos.err
#SBATCH -o missing_genos.out

#apptainer exec --bind /beegfs/home/nbailey/PixyUpdate/pixy_analysis/data_generation/01_simulating-test-data /beegfs/data/soft/containers/R-4_4_1.sif sh run_simulation_partial_missing.sh

apptainer exec --bind /beegfs/data/nbailey/PixyUpdate/pixy_analysis/data_generation/01_simulating-test-data /beegfs/data/soft/containers/R-4_4_1.sif Rscript 02_simulate_partial_missing.R
