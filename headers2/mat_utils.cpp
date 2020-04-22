#ifndef MAT_UTILS
#define MAT_UTILS

#include <complex>
#include "greens.hpp"
#include "elementops.hpp"
#include "mat_utils.hpp"


#define PI 3.1415926535897932384626433832795028841971693

namespace NEdyson{

void get_dftcorr_cubic(double th, double *corfac, cplx *endcor){
  double ai[4],ar[4],t;
  double t2,t4,t6;
  double cth, ctth, spth2,sth,sth4i,stth,th2,th4,tmth2,tth4i;
  int j;

  if (fabs(th) < 5.0e-2) {
    t=th;
    t2=t*t;
    t4=t2*t2;
    t6=t4*t2;
    *corfac=1.0-(11.0/720.0)*t4+(23.0/15120.0)*t6;
    ar[0]=(-2.0/3.0)+t2/45.0+(103.0/15120.0)*t4-(169.0/226800.0)*t6;
    ar[1]=(7.0/24.0)-(7.0/180.0)*t2+(5.0/3456.0)*t4-(7.0/259200.0)*t6;
    ar[2]=(-1.0/6.0)+t2/45.0-(5.0/6048.0)*t4+t6/64800.0;
    ar[3]=(1.0/24.0)-t2/180.0+(5.0/24192.0)*t4-t6/259200.0;
    ai[0]=t*(2.0/45.0+(2.0/105.0)*t2-(8.0/2835.0)*t4+(86.0/467775.0)*t6);
    ai[1]=t*(7.0/72.0-t2/168.0+(11.0/72576.0)*t4-(13.0/5987520.0)*t6);
    ai[2]=t*(-7.0/90.0+t2/210.0-(11.0/90720.0)*t4+(13.0/7484400.0)*t6);
    ai[3]=t*(7.0/360.0-t2/840.0+(11.0/362880.0)*t4-(13.0/29937600.0)*t6);
  } else {
    cth=cos(th);
    sth=sin(th);
    ctth=cth*cth-sth*sth;
    stth=2.0e0*sth*cth;
    th2=th*th;
    th4=th2*th2;
    tmth2=3.0e0-th2;
    spth2=6.0e0+th2;
    sth4i=1.0/(6.0e0*th4);
    tth4i=2.0e0*sth4i;
    *corfac=tth4i*spth2*(3.0e0-4.0e0*cth+ctth);
    ar[0]=sth4i*(-42.0e0+5.0e0*th2+spth2*(8.0e0*cth-ctth));
    ai[0]=sth4i*(th*(-12.0e0+6.0e0*th2)+spth2*stth);
    ar[1]=sth4i*(14.0e0*tmth2-7.0e0*spth2*cth);
    ai[1]=sth4i*(30.0e0*th-5.0e0*spth2*sth);
    ar[2]=tth4i*(-4.0e0*tmth2+2.0e0*spth2*cth);
    ai[2]=tth4i*(-12.0e0*th+2.0e0*spth2*sth);
    ar[3]=sth4i*(2.0e0*tmth2-spth2*cth);
    ai[3]=sth4i*(6.0e0*th-spth2*sth);
  }
  for(j=0;j<=3;j++) endcor[j]=cplx(ar[j],ai[j]);
}


void matsubara_ft(cplx *res, int m, green_func &Sig, cplx *sigdft, int sig, double beta){
  double corfac, arg, dtau, one;
  int ntau, m1, r, l, sg, size1=Sig.size1();
  cplx *z1, *z2, *dft, bcorr[4];

  one=(sig==-1?1.0:0.0);
  ntau=Sig.ntau();
  sg=Sig.element_size();
  dtau=beta/ntau;
  z1 = new cplx[sg];
  z2 = new cplx[sg];
  dft= new cplx[sg];
  arg=(2*m+one)*PI/ntau;
  get_dftcorr_cubic(arg,&corfac,bcorr);
  m1=m;
  while(m1<0) m1+=ntau;
  while(m1>=ntau) m1-=ntau;
  for(l=0;l<sg;l++) dft[l]=sigdft[m1*sg+l]*corfac;
  element_set_zero(size1,z1);
  for(int r=0;r<=3;r++){
    element_set(size1,z2,Sig.matptr(r));
    element_smul(size1,z2,bcorr[r]);
    element_incr(size1,z1,z2);
    element_set(size1,z2,Sig.matptr(ntau-r));
    element_smul(size1,z2,(1.*sig)*cplx(bcorr[r].real(),-bcorr[r].imag()));
    element_incr(size1,z1,z2);
  }
  for(l=0;l<sg;l++) res[l]=dtau*(z1[l]+dft[l]);
  delete[] z1;
  delete[] z2;
  delete[] dft;
}


double get_tau(int r, double beta, int ntau){
  return r*beta/ntau;
}


double get_omega(int m, double beta, int sig){
  int i=(sig==-1?1:0);
  return PI/beta*(2*m+i);
}

void matsubara_dft(cplx *mdft, green_func &G, int sig){
  int ntau, r, m, sg, size1=G.size1();
  double arg, one;
  cplx *z, *z1, expfac;
  sg=G.element_size();
  ntau=G.ntau();
  z= new cplx[sg];
  z1=new cplx[sg];
  one=(sig==-1?1.0:0.0);
  for(m=0;m<=ntau;m++){
    element_set_zero(size1,z);
    for(r=0;r<=ntau;r++){
      arg=((2*m+one)*r*PI)/ntau;
      expfac=cplx(cos(arg),sin(arg));
      element_set(size1,z1,G.matptr(r));
      element_smul(size1,z1,expfac);
      element_incr(size1,z,z1);
    }
    element_set(size1,mdft+m*sg,z);
  }
  delete[] z;
  delete[] z1;
}


void set_first_order_tail(cplx *gmat, cplx *one, double beta, int sg, int ntau, int sig, int size1){
  cplx *z1 = new cplx[sg];
  for(int r=0;r<=ntau;r++){
    element_set(size1,gmat+r*sg,one);
    element_smul(size1,gmat+r*sg,-0.5);
    if(sig==1){
      double tau=r*beta/ntau;
      element_set(size1,z1,one);
      element_smul(size1,z1,tau/beta);
      element_incr(size1,gmat+r*sg,z1);
    }
  }
  delete[] z1;
}

void force_mat_herm(green_func &G){
	int ntau=G.ntau(), m;
	cdmatrix Gmat,Gmat_herm;
	for(m=0;m<=ntau;m++){
		G.get_mat(m,Gmat);
		Gmat_herm=0.5*(Gmat+Gmat.adjoint());
		G.set_mat(m,Gmat_herm);
	}
}

}//namespace
#endif
