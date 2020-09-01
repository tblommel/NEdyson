#!/bin/bash

#SBATCH --job-name tti_NEdyson
#SBATCH --nodes=1
#SBATCH --time=00:15:00
#SBATCH --partition=batch

./build/programs/molNEdyson ./data/h2sto.params.ini
