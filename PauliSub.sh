#!/bin/bash

#SBATCH --job-name NEdyson
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=batch,super,debug

./build/programs/molNEdyson ./inputs/NEdyson-input0.ini
