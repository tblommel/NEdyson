import os

def makefile(nt, dt, k, hfbool, number):
  file1 = open("inputs/NEdyson-input" + str(number) + ".ini", "w")
  file1.write("mode = GF2\n")
  file1.write("tti = 1\n")
  file1.write("hfbool = " + hfbool + "\n")
  file1.write("\n")
  file1.write("unrestricted = 0\n")
  file1.write("decomposed = 0\n")
  file1.write("repr = cheb\n")
  file1.write("repr-file = /home/tblommel/Libraries/gfmol/data/repn/cheby/512.hdf5\n")
  file1.write("\n")
  file1.write("hf-input = /pauli-storage/tblommel/NEdyson_data/He-VB2PP/He-VB2PPinp_new.h5\n")
  file1.write("\n")
  file1.write("boolPumpProbe = 0\n")
  file1.write("PumpProbe-file = /pauli-storage/tblommel/DNE/DNE.h5\n")
  file1.write("\n")
  file1.write("boolOutput = 1\n")
  if(hfbool == "1"):
    file1.write("output = /pauli-storage/tblommel/NEdyson_data/He-VB2PP/dt_conv_diff_k/tti_k" + k + "_dt" + dt + "_HF_full.h5\n")
  else:
    file1.write("output = /pauli-storage/tblommel/NEdyson_data/He-VB2PP/dt_conv_diff_k/tti_k" + k + "_dt" + dt + "_2B_full.h5\n")
  file1.write("\n")
  file1.write("boolOutputPP = 0\n")
  file1.write("outputPP = /pauli-storage/tblommel/DNE/DNE.h5\n")
  file1.write("\n")
  file1.write("nPumpProbe = 3.556433e-8\n")
  file1.write("lPumpProbe = 3.0235617994e6\n")
  file1.write("\n")
  file1.write("nt = " + nt + "\n")
  file1.write("ntau = 110\n")
  file1.write("dt = 0." + dt + "\n")
  file1.write("k = " + k + "\n")
  file1.write("\n")
  file1.write("beta = 100.0\n")
  file1.write("\n")
  file1.write("BootMaxIter = 50\n")
  file1.write("BootMaxErr = 1e-13\n")
  file1.write("CorrSteps = 15\n")
  file1.write("maxiter = 50\n")
  file1.write("\n")
  file1.write("etol = 1e-12\n")
  file1.write("decomp-prec = 1e-16\n")

  file1.close()

def makePS(number):
  file1 = open("PauliSub" + str(number) + ".sh", "w")
  file1.write("#!/bin/bash\n")
  file1.write("\n")
  file1.write("#SBATCH --job-name NEdyson_PP_hf\n")
  file1.write("#SBATCH --nodes=1\n")
  file1.write("#SBATCH --time=08:00:00\n")
  file1.write("#SBATCH --partition=super,batch\n")
  file1.write("#SBATCH --exclusive\n")
  file1.write("#SBATCH --array=0-0\n")
  file1.write("\n")
  file1.write("CASE_NUM=`printf %01d $SLURM_ARRAY_TASK_ID`\n")
  file1.write("#!/bin/bash\n")
  file1.write("./build/programs/molNEdyson ./inputs/NEdyson-input" + str(number) + ".ini\n")


nt_options = ["50","50","50","50","50","50"]
dt_options = ["12","08","04","02","01","005"]
k_options = ["5","4","3","2","1"]
hfbool_options = ["0"]
number = 0

for i,nt in enumerate(nt_options):
  for k in k_options:
    for hfbool in hfbool_options:
      makefile(nt, dt_options[i], k, hfbool, number)
      makePS(number)
      os.system("sbatch PauliSub" + str(number) + ".sh")
      number += 1
