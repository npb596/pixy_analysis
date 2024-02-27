#!/bin/bash

#SBATCH -J 1e6_data_sim
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH -t 2-00:00:00
#SBATCH -p jro0014_amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=npb0015@auburn.edu

module load python/anaconda/3.8.6

# simulate test data for evaluating the effect of missing data

mkdir -p data/simulated_var_only

# generate 100 datasets
python scripts/msprime_simulate_data.py --prefix data/simulated_var_only/pi_sim_ --datasets 100

# for each dataset, inject invariant sites 
# (automatically loops through all files in a folder)
sh scripts/inject_invariant_sites.sh data/simulated_var_only data/simulated_invar

# generate missing datasets from invariant data

# get list of simulated data
ls data/simulated_invar/*100000*vcf.gz > invar_list.txt

# apply missingness script to each file for each level of missingness
while read invar_vcf

do

	for i in $(seq 0 1000 100000)
	do
		echo $i
	
		#((i = $i * 100))

		sh scripts/create_missing_sites.sh $invar_vcf $i data/simulated_missing_sites/$i
	
	done
	
done < invar_list.txt
