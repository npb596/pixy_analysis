#!/bin/bash

#SBATCH -J scikit_allel
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH -t 6-00:00:00
#SBATCH -p normal
#SBATCH -e scikit_allel.err
#SBATCH -o scikit_allel.out

rm data/pi_out/*
rm data/dxy_out/*
rm data/watterson_theta_out/*
rm data/tajima_d_out/*

ls -d ../01_simulating-test-data/data/simulated_missing_genos/* > tmp/missing_genos_dir_list.txt
ls -d ../01_simulating-test-data/data/simulated_missing_sites/* > tmp/missing_sites_dir_list.txt
#
while read vcffolder 
do
# 
./02_scikit_dxy_pi_vcfolder.sh $vcffolder
# 
done < tmp/missing_genos_dir_list.txt
#
#
while read vcffolder 
do
#
./02_scikit_dxy_pi_vcfolder.sh $vcffolder
#
done < tmp/missing_sites_dir_list.txt
#
#
#sbatch 02_scikit_dxy_pi_vcfolder.sh ../01_simulating-test-data/data/simulated_invar


#sbatch 02_scikit_dxy_pi_vcfolder.sh ../simulating-test-data/data/accuracy_invar
