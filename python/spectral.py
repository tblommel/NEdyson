import os.path
import numpy as np
import matplotlib.pyplot as plt
import greens

class Spectral:
  @classmethod
  def blankA(cls):
    return cls(0,0,0,0,0,0)

  def read_from_file(self,filepath):
    actfile = filepath+"_A.dat"
    infile = open(actfile,'r')
    line = infile.readline()
    a = line.split()
    self.nt = int(a[0])
    self.nw = int(a[1])
    self.size1 = int(a[2])
    self.es = int(a[2])*int(a[2])
    self.wmax = float(a[3])
    self.dw = float(a[4])
    self.dt = float(a[5])
    
    A = np.zeros((self.nt+1,self.nw,self.size1,self.size1),dtype=float)
    count1=0
    count2=0
    for i,line in enumerate(infile):
      a=line.split()
      for j in range(self.size1):
        for k in range(self.size1):
          A[count1,count2,j,k] = float(a[j*self.size1+k])
      count2+=1
      if(count2>=self.nw):
        count2=0
        count1+=1
    self.A=A


  def plot(self,title):
    fig, ax1 = plt.subplots(figsize=(13,13),ncols=self.size1,nrows=self.size1)
    for i in range(self.size1):
      for j in range(self.size1):
        img=ax1[i,j].imshow(self.A[:,:,i,j], origin = 'lower', extent=[-10,10,0,24])
        fig.colorbar(img,ax=ax1[i,j])
        ax1[i,j].set_xlabel("omega")
        ax1[i,j].set_title("("+str(i)+","+str(j)+") component")
        ax1[i,j].set_ylabel("t_avg ((t+t')/2)")
    fig.suptitle(title,fontsize=16)
    plt.show()


  def AfromG(self, G, nw, wmax):
    self.nt = 2*G.nt
    self.nw = nw
    self.size1 = G.size1
    self.es = self.size1*self.size1
    self.wmax = wmax
    self.dw = 2*wmax/(nw-1)
    self.dt = G.dt
    
    A=np.zeros((self.nt+1,nw,self.size1,self.size1))
    
    cplxi=complex(0,1)
    for ita in range(self.nt+1):
      print(ita)
      maxx = min(ita,self.nt-ita)
      for itr in np.arange(ita%2,maxx+1,2):
        tmp = G.ret[(ita+itr)//2,(ita-itr)//2,:,:]
        for iw in range(nw):
          arg = (iw-(nw-1)/2)*self.dw*itr*self.dt*cplxi
          weight = np.exp(arg)
          A[ita,iw,:,:] += (weight*tmp).imag

    self.A=A
    
    if(G.retexp != 0):
      finishthispart=True
      


  def print_to_file(self,filepath):
    actfile = filepath+"_A.dat"
    outfile = open(actfile,'w')
    outfile.write(str(self.nt)+" "+str(self.nw)+" "+str(self.size1)+" "+str(self.wmax)+" "+str(self.dw)+" "+str(self.dt)+"\n")
    for t in range(self.nt+1):
      for j in range(self.nw):
        for i in range(self.size1):
          for k in range(self.size1):
            outfile.write(str(self.A[t,j,i,k])+" ")
        outfile.write("\n")

    

  def __init__(self,nt,nw,size1,wmax,dw,dt):
    self.nt = nt
    self.nw = nw
    self.size1 = size1
    self.es = size1*size1
    self.wmax = wmax
    self.dw = dw
    self.dt = dt
    self.A = np.zeros((self.nt+1,self.nw,self.size1,self.size1),dtype=float)
