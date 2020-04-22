import os.path
import math
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.linalg as linalg
import warnings

def is_stable(a):
    """
    check if the autoregressive process
    described by coefficients a is stable
    (i.e. will decay to zero)
    """
    alpha = np.roots(np.concatenate(([1], a)))
    return all(np.abs(alpha) < 1.0)

def autocorrelations(x, p):
    """
    Calculate autocorrelation vector/matrix
    (slow but matrices are small so ok)
    """
    n = x.size

    r = np.zeros(p, dtype=complex)
    for i in range(p):
        for k in range(p, n):
            r[i] += x[k] * np.conj(x[k-i-1])

    R = np.zeros((p,p), dtype=complex)
    for i in range(p):
        for j in range(p):
            for k in range(p, n):
                R[i,j] += x[k-j-1] * np.conj(x[k-i-1])

    r /= n
    R /= n
    return r, R

def lpc(x, p, eps=None):
    """
    find linear prediction coefficients a
    """
    r, R = autocorrelations(x, p)
    lambda_min = np.min(linalg.eig(R)[0])

    if eps is None:
        eps = np.power(10, np.log10(np.abs(lambda_min)) + 1.0)

    a = linalg.solve(R + eps*np.identity(p), -r)

    if not is_stable(a):
        warnings.warn("solution unstable")
    return a

def lpredict(x, *, n, p, nfit):
    """
    extend sequence x by n elements using linear prediction with p lags over last nfit values
    """
    assert len(x.shape) == 1
    m = x.shape[0]

    a = lpc(x[-nfit:], p)

    y = np.zeros((m + n,), complex)
    y[:m] = x

    mse = 0.0
    for i in range(p, m):
        mse += (1.0/(m - p)) * np.abs(np.dot(-a[::-1], y[i-p:i]) - y[i])**2
    rmse = np.sqrt(mse)

    for i in range(m, m+n):
        y[i] = np.dot(-a[::-1], y[i-p:i])

    if(abs(y[-1].real)>10e3 or math.isnan(y[-1].real)):
        print("PUT OUT ZEROS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        return np.zeros((y[m:].shape),dtype=complex), rmse
    else:
        return y[m:], rmse



class Greens:
	@classmethod
	def blankG(cls):
		return cls(0,0,0,0,0,0,0,0,0,0,0)

	def read_from_file_mat(self,filepath):
		actfile = filepath+"_GM.dat"
		infile = open(actfile,'r')
		line = infile.readline()
		a = line.split()
		self.nt = int(a[0])
		self.ntau = int(a[1])
		self.size1 = int(a[2])
		self.dt = float(a[3])
		self.dtau = float(a[4])
		self.es = int(a[2])*int(a[2])
		self.sig = int(a[5])
		mat = np.zeros((self.ntau+1,self.size1,self.size1),dtype=complex)
		for i,line in enumerate(infile):
			a=line.split()
			for j in range(self.size1):
				for k in range(self.size1):
					mat[i,j,k] = complex(float(a[(j*self.size1+k)*2]),float(a[(j*self.size1+k)*2+1]))
		self.mat=mat

	def read_from_file_ret(self,filepath):
		#open file
		actfile = filepath+"_GR.dat"
		infile = open(actfile,'r')

		#read in first line
		line = infile.readline()
		a = line.split()

		#if it contains expanded information then read in next line		
		expbool=False
		if(a[0] == "expanded"):
			expbool=True
			line = infile.readline()
			a = line.split()

		#read in gfunc params
		self.nt = int(a[0])
		self.ntau = int(a[1])
		self.size1 = int(a[2])
		self.dt = float(a[3])
		self.dtau = float(a[4])
		self.es = int(a[2])*int(a[2])
		self.sig = int(a[5])
		ret = np.zeros((self.nt+1,self.nt+1,self.size1,self.size1),dtype=complex)

		#if we have expanded information, read in the parameters
		nfit=ntp=0
		if(expbool):
			line = infile.readline()
			a = line.split()
			nfit = a[0]
			ntp = a[1]
		retexp = np.zeros((self.nt,self.nt-nfit-ntp),dtype=complex)

		#read in the gfunc values
		for t in range(self.nt+1):
			for tp in range(t+1):
				line = infile.readline()	
				a=line.split()
				for j in range(self.size1):
					for k in range(self.size1):
						ret[t,tp,j,k] = complex(float(a[(j*self.size1+k)*2]),float(a[(j*self.size1+k)*2+1]))
		self.ret=ret
		
		#read in expanded values
		if(expbool):
			for tp in range(self.nt-nfit-ntp):
				for t in range(self.nt-tp):
					line = infile.readline()
					a = line.split()
					for i in range(self.size1):
						for j in range(self.size1):
							retexp[t,tp,i,j] = complex(float(a[(j*self.size1+k)*2]),float(a[(j*self.size1+k)*2+1]))
			self.retexp = retexp

	def read_from_file_les(self,filepath):
		actfile = filepath+"_GL.dat"
		infile = open(actfile,'r')
		line = infile.readline()
		a = line.split()
		self.nt = int(a[0])
		self.ntau = int(a[1])
		self.size1 = int(a[2])
		self.dt = float(a[3])
		self.dtau = float(a[4])
		self.es = int(a[2])*int(a[2])
		self.sig = int(a[5])
		les = np.zeros((self.nt+1,self.nt+1,self.size1,self.size1),dtype=complex)
		count1=0
		count2=0
		for i,line in enumerate(infile):
			a=line.split()
			for j in range(self.size1):
				for k in range(self.size1):
					les[count2,count1,j,k] = complex(-float(a[(k*self.size1+j)*2]),float(a[(k*self.size1+j)*2+1]))
					les[count1,count2,j,k] = complex(float(a[(j*self.size1+k)*2]),float(a[(j*self.size1+k)*2+1]))
			count1+=1
			if(count1>count2):
				count1=0
				count2+=1
		self.les=les

	def read_from_file_tv(self,filepath):
		actfile = filepath+"_GTV.dat"
		infile = open(actfile,'r')
		line = infile.readline()
		a = line.split()
		self.nt = int(a[0])
		self.ntau = int(a[1])
		self.size1 = int(a[2])
		self.dt = float(a[3])
		self.dtau = float(a[4])
		self.es = int(a[2])*int(a[2])
		self.sig = int(a[5])
		tv = np.zeros((self.nt+1,self.ntau+1,self.size1,self.size1),dtype=complex)
		count1=0
		count2=0
		for i,line in enumerate(infile):
			a=line.split()
			for j in range(self.size1):
				for k in range(self.size1):
					tv[count1,count2,j,k] = complex(float(a[(j*self.size1+k)*2]),float(a[(j*self.size1+k)*2+1]))
			count2+=1
			if(count2>self.ntau):
				count2=0
				count1+=1
		self.tv=tv

	def read_from_file(self,filepath):
		self.read_from_file_mat(filepath)
		self.read_from_file_ret(filepath)
		self.read_from_file_les(filepath)
		self.read_from_file_tv(filepath)


	def expand_ret(self,p,nfit,ntp):
		newdat = np.zeros((self.nt,self.nt-nfit-ntp,self.size1,self.size1),dtype=complex)
		for tp in range(self.nt-nfit-ntp):
			print(tp)
			for i in range(self.size1):
				for j in range(self.size1):
					newdat[:self.nt-tp,tp,i,j], err = lpredict(self.ret[:,tp,i,j],n=self.nt-tp,p=p,nfit=nfit)
					if(tp==273):
						for tt in range(newdat[:self.nt-tp,tp,i,j].shape[0]):
							print(i,j,tt,newdat[tt,tp,i,j])
		self.retexp = newdat



	def __init__(self,mat,les,ret,retexp,tv,sig,nt,ntau,size1,dt,dtau):
		self.mat = mat
		self.les = les
		self.ret = ret
		self.retexp = retexp
		self.tv = tv
		self.sig = sig
		self.nt = nt
		self.ntau = ntau
		self.dt = dt
		self.dtau = dtau
		self.size1 = size1
		self.es = size1*size1
