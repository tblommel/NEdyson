#!/bin/bash

#SBATCH --job-name NEdyson_PP_hf
#SBATCH --nodes=1
#SBATCH --time=36:00:00
#SBATCH --partition=super
#SBATCH --exclusive
#SBATCH --array=0-0

CASE_NUM=`printf %01d $SLURM_ARRAY_TASK_ID`

./build/programs/molNEdyson ./data/NEdyson-input-hf-tti.ini
