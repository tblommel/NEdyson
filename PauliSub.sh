#!/bin/bash

#SBATCH --job-name NEdyson
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=batch,super

./build/programs/molNEdyson ./inputs/NEdyson-input0.ini
