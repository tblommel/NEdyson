#!/bin/bash

#SBATCH --job-name NEdyson
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --partition=ludicrous
#SBATCH --exclusive

./build/programs/molNEdyson ./inputs/NEdyson-input-H2-2.ini
