#!/bin/bash

#SBATCH --job-name tti_NEdyson
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --partition=super
#SBATCH --exclusive

export OMP_NUM_THREADS=64
./build/programs/molNEdyson ./data/h2sto.params.ini
