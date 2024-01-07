#!/bin/bash

#SBATCH --job-name NEdyson
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --partition=ludicrous

./build/programs/molNEdyson ./inputs/NEdyson-input-H2.ini
