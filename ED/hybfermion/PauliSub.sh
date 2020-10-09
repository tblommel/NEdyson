#!/bin/bash

#SBATCH --job-name hybfermion
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=super

./build/fermion 10
