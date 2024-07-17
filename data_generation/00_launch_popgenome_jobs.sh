#!/bin/bash

#sbatch 01_run_popgenome.sh 01_simulating-test-data/data/simulated_missing_genos 8
#sbatch 01_run_popgenome.sh 01_simulating-test-data/data/simulated_missing_sites 8
sbatch 01_run_popgenome.sh 01_simulating-test-data/data/simulated_invar 8
#sbatch 01_run_popgenome.sh ../01_simulating-test-data/data/accuracy_invar 12
