#!/bin/bash

#SBATCH --job-name NEdyson_PP_hf
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --partition=super,batch
#SBATCH --array=0-7

CASE_NUM=`printf %01d $SLURM_ARRAY_TASK_ID`

./build/programs/molNEdyson ./data/NEdyson-input${CASE_NUM}.ini
