#!/bin/bash

#SBATCH --job-name hodlr
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=super
#SBATCH --exclusive

./build/programs/SC /pauli-storage/tblommel/SC_mat_Legendre_points.h5 /pauli-storage/tblommel/hodlr_SC/omega0.2/w0.2gc16gs6.5h0.02.h5 6000 80 18 4e-2
