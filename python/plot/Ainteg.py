#!/usr/bin/env python3

import argparse
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import warnings
import spectral

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="Arguments") 
	parser.add_argument("--dir",type=str)
	args = parser.parse_args()

	#initialize the greens functions
	A = spectral.Spectral.blankA()

	#read in the ret functions
	A.read_from_file(args.dir)
	integ=np.sum(A.A,1)
	print(integ[1600,:,:]*A.dw)

	A.plot("Non interacting spectral function")
