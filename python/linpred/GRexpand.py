#!/usr/bin/env python3
import os.path
import argparse
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import warnings
import greens

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="Arguments") 
	parser.add_argument("--dir",type=str, required = True)
	parser.add_argument("-p",type=int, required = True)
	parser.add_argument("--nfit",type=int, required = True)
	parser.add_argument("--ntp",type=int, required = True)
	args = parser.parse_args()

	#initialize the greens functions
	G = greens.Greens.blankG()

	#read in the ret functions
	G.read_from_file_ret(args.dir)

	#expand
	G.expand_ret(args.p,args.nfit,args.ntp)

	#output
	actfile = args.dir + "_GRExp.dat"
	outfile = open(actfile,'w')
	outfile.write(str(args.nfit)+" "+str(args.ntp)+"\n")
	for tp in range(G.nt-args.nfit-args.ntp):
		for t in range(G.nt-tp):
			for i in range(G.size1):
				for j in range(G.size1):
					outfile.write(str(G.retexp[t,tp,i,j].real)+" "+str(G.retexp[t,tp,i,j].imag)+" ")
			outfile.write("\n")
	
	outfile.close()
