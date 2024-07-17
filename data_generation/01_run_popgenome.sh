#!/bin/bash

#SBATCH -J run_popgenome
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=8G
#SBATCH -t 1-00:00:00
#SBATCH -p normal
#SBATCH -e run_popgenome.err
#SBATCH -o run_popgenome.out

apptainer exec --bind /beegfs/data/nbailey/PixyUpdate/pixy_analysis/data_generation/ /beegfs/data/soft/containers/R-4_4_1.sif Rscript popgenome/02_compute_pi_vcf_folder.R $1 $2
