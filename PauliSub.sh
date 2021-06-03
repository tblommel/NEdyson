#!/bin/bash

#SBATCH --job-name NEdyson_PP_hf
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=debug,super,batch

./build/programs/molNEdyson ./data/NEdyson-input.ini
