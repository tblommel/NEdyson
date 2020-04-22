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
	G.read_from_file_mat(args.dir)


	xvals = np.arange(0,G.dtau*(G.ntau+0.5),G.dtau)
	fig,axs=plt.subplots(1,1,figsize=(13,13))
	for i in range(G.size1):
		for j in range(G.size1):
			axs.plot(xvals, G.mat[:,i,j].real, label = str(i)+str(j)+"real")
	axs.set_title("Non interacting Mat Greens Func GM(tau)")
	axs.legend()
	plt.show()		
