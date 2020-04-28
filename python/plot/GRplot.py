#!/usr/bin/env python3

import argparse
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import warnings
import greens

if __name__=="__main__":
  parser = argparse.ArgumentParser(description="Arguments") 
  parser.add_argument("--dir",type=str)
  args = parser.parse_args()

  #initialize the greens functions
  G = greens.Greens.blankG()

  #read in the ret functions
  G.read_from_file_ret(args.dir)

  xvals = np.arange(0,1601*G.dt,G.dt)
  fig,axs=plt.subplots(1,1,figsize=(13,13))
  axs.plot(xvals,G.ret[:,0,0,0].real,label = "real 00")
  axs.plot(xvals,G.ret[:,0,0,0].imag,label = "imag 00")
  axs.plot(xvals,G.ret[:,0,1,0].real,label = "real 10")
  axs.plot(xvals,G.ret[:,0,1,0].imag,label = "imag 10")
  axs.plot(xvals,G.ret[:,0,0,1].real,label = "real 01")
  axs.plot(xvals,G.ret[:,0,0,1].imag,label = "imag 01")
  axs.plot(xvals,G.ret[:,0,1,1].real,label = "real 11")
  axs.plot(xvals,G.ret[:,0,1,1].imag,label = "imag 11")
  axs.set_ylim(bottom=-1.5,top=1.5)
  axs.set_title("Non interacting Ret Greens Func GR(t,0)")
  axs.legend()
  plt.show()
