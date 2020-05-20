#!/bin/bash

#SBATCH --job-name e3NEdyson
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition=super

export OMP_NUM_THREADS=2
./build/programs/hubb_chain_2b_omp.x /data/tblommel/NEdysondata/omp/2/
