#!/usr/bin/env python3

import argparse
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import warnings
import spectral
import greens

if __name__=="__main__":
  parser = argparse.ArgumentParser(description="Arguments") 
  parser.add_argument("--Adir",type=str)
  parser.add_argument("--AdirMat",type=str)
  args = parser.parse_args()

  #initialize the greens functions
  A = spectral.Spectral.blankA()

  #read in the ret functions
  A.read_from_file(args.Adir)

  AM = spectral.Spectral(1,201,4,10,0.1,1)

  infile = open(args.AdirMat,'r')
  cw=0
  for line in infile:
    AM.A[0,cw,0,0] = float(line)
    cw+=1

  fig,ax1 = plt.subplots(figsize=(13,13),ncols=1,nrows=1)
  for i in range(1):
    for j in range(1):
      ax1.plot(np.arange(-10,10.1,0.1),AM.A[0,:,i,j],label="Mathematica")
      ax1.plot(np.arange(-10,10.1,A.dw),A.A[1600,:,i,j],label="My Code")
      ax1.set_xlabel("omega")
      ax1.set_title("("+str(i+1)+","+str(j+1)+") component")
      ax1.set_xlim(left=-5, right =5)
  plt.legend()
  plt.show()
