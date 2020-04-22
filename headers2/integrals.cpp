#ifndef INTEGRALS_IMPL
#define INTEGRALS_IMPL

#include "integrals.hpp"

namespace NEdyson{


// This one does integral for every tau, and puts result into B by incrementing it!!!
// Does not include the j==n term
void CTV1(integration::Integrator &I, green_func &A, green_func &Acc, green_func &B, int n, double dt){
  int k=I.k(), ntau=A.ntau(), size1=A.size1(), es=A.element_size(), i, j, m;
  double weight;
  cplx *Atmp = new cplx[es];
  cplx *resptr;
  cplx *Gptr;

  for(j=0;j<n;j++){
    weight = I.gregory_weights(n,j);
    if(j>n){ // Dont have AR
      element_conj(size1,Atmp,Acc.retptr(j,n));
      element_smul(size1,Atmp,-1);
    }
    else{
      element_set(size1,Atmp,A.retptr(n,j));
    }
    element_smul(size1,Atmp,dt);

    resptr=B.tvptr(n,0);
    Gptr=B.tvptr(j,0);
    if(weight != 1){
      for(m=0;m<=ntau;m++){
        element_incr(size1,resptr,weight,Atmp,Gptr);
        Gptr+=es;
        resptr+=es;
      }
    }
    else{
      for(m=0;m<=ntau;m++){
        element_incr(size1,resptr,Atmp,Gptr);
        Gptr+=es;
        resptr+=es;
      }
    }
  }

  delete[] Atmp;
}


void CTV2(integration::Integrator &I, green_func &A, green_func &B, int n, int m, double beta, cplx *res){
  int k=I.k(), ntau=A.ntau(), size1=A.size1(), es=A.element_size(), sig=B.sig(),j,l;
  double dtau=beta/ntau;
  cplx *BM;
  cplx *ARM;
  for(int j=0;j<es;j++) res[j]=0;

  if(m==0){} //do nothing
  else if(m<k){//use starting weights
    BM=B.matptr(ntau);
    ARM=A.tvptr(n,0);
    for(j=0;j<=k;j++){
      for(l=0;l<=k;l++){
        element_incr(size1,res,I.rcorr(m,j,l),ARM+l*es,BM-j*es);
      }
    }
  }
  else if(m<k+k+1){//use gregory weights
    BM=B.matptr(ntau);
    ARM=A.tvptr(n,m);
    for(l=0;l<=m;l++){//Omegas
      element_incr(size1,res,I.gregory_weights(m,l),ARM-l*es,BM-l*es);
    }
  }
  else{//use greg and ones
    BM=B.matptr(ntau);
    ARM=A.tvptr(n,m);
    for(l=0;l<=k;l++){//omegas
      element_incr(size1,res,I.omega(l),ARM-l*es,BM-l*es);
    }
    for(l=k+1;l<m-k;l++){//ones
      element_incr(size1,res,ARM-l*es,BM-l*es);
    }
    for(l=m-k;l<=m;l++){//omegas
      element_incr(size1,res,I.omega(m-l),ARM-l*es,BM-l*es);
    }
  }
  for(l=0;l<es;l++) res[l]*=(sig*dtau);
}


void CTV3(integration::Integrator &I, green_func &A, green_func &B, int n, int m, double beta, cplx *res){
  int k=I.k(), ntau=A.ntau(), size1=A.size1(), es=A.element_size(), sig=A.sig(),j,l;
  double dtau=beta/ntau;
  cplx *BM;
  cplx *ARM;
  for(j=0;j<es;j++) res[j]=0;

  if(m==ntau){}//do nothing
  else if(m>ntau-k){//use starting weights
    BM=B.matptr(0);
    ARM=A.tvptr(n,ntau);
    for(j=0;j<=k;j++){
      for(l=0;l<=k;l++){
        element_incr(size1,res,I.rcorr(ntau-m,j,l),ARM-l*es,BM+j*es);
      }
    }
  }
  else if(ntau-m>=k+k+1){//use greg and ones
    BM=B.matptr(0);
    ARM=A.tvptr(n,m);
    for(l=0;l<=k;l++){
      element_incr(size1,res,I.omega(l),ARM+l*es,BM+l*es);
    }
    for(l=k+1;l<ntau-m-k;l++){
      element_incr(size1,res,ARM+l*es,BM+l*es);
    }
    for(l=ntau-m-k;l<=ntau-m;l++){
      element_incr(size1,res,I.omega(ntau-m-l),ARM+l*es,BM+l*es);
    }
  }
  else{//use greg weights
    BM=B.matptr(0);
    ARM=A.tvptr(n,m);
		for(l=0;l<=ntau-m;l++){
			element_incr(size1,res,I.gregory_weights(ntau-m,l),ARM+l*es,BM+l*es);
		}
	}
  for(int l=0;l<es;l++) res[l]*=dtau;
}


// Does the les_lesadv integration for every n, puts into Q via increment!!!
// dt \sum_{j=0}^m w_{mj} SL_{nj} GA_{jm}     for n=0...m
void Cles2_tstp(integration::Integrator &I, green_func &A, green_func &Acc, green_func &B, int m, double dt, cplx *Q){
  int k=I.k(), size1=A.size1(), es=A.element_size(), sig=A.sig(),j,l,n,top=(m>k) ? (m):(k);
  double w;
  cplx *BA=new cplx [(top+1)*es];
  cplx *AL=new cplx [es];

  // First fill BA, since they are the same for each n
  for(j=0;j<=top;j++){//fill BA from B.retptr.  BA_jm=(BR_mj)^T*
    w=I.gregory_weights(m,j);
    if(m>=j){//we have BR
      element_conj(size1,BA+j*es,B.retptr(m,j));
      element_smul(size1,BA+j*es,dt*w);
    }
    else{//dont have BR
      element_set(size1,BA+j*es,B.retptr(j,m));
      element_smul(size1,BA+j*es,-dt*w);
    }
  }

  cplx *res = Q;
  //sum from j=0 to n-1, where we don't have AL
  for(n=0;n<=top;n++){
    for(j=0;j<n;j++){
      element_conj(size1,AL,Acc.lesptr(j,n));
      element_incr(size1,res,-1.,AL,BA+j*es);
    }
    res+=es;
  }

  res = Q;
  //sum from j=n to max(k,m), here we have AL
  for(n=0;n<=top;n++){
    for(j=n;j<=top;j++){
      element_incr(size1,res,A.lesptr(n,j),BA+j*es);
    }
    res+=es;
  }

  delete[] BA;
  delete[] AL;
}


// Does the les_tvvt integral for each n=0...m
void Cles3_tstp(integration::Integrator &I,green_func &A, green_func &B, int m, double beta, cplx *Q){
  int k=I.k(), size1=A.size1(), es=A.element_size(), sig=B.sig(),j,l,n, ntau=A.ntau();
  int top = (m>=k)?(m):(k);
  double dtau=beta/ntau;
  cplx weight;
  cplx cplxi = cplx(0.,1.);
  cplx *BVT = new cplx[(ntau+1)*es];
  cplx *ATV;
  cplx *bvttmp;

  // Fill BVT first
  weight = (double)sig*cplxi*dtau;
  for(j=0;j<=ntau;j++) element_conj(size1, BVT+j*es, B.tvptr(m,ntau-j));
  for(l=0;l<(ntau+1)*es;l++) BVT[l]*=weight;

  //do the integral
  cplx *res=Q;
  for(n=0;n<=top;n++){
    bvttmp = BVT;
    ATV = A.tvptr(n,0);
    if(ntau < 2*(k+1)-1){
      weight = I.gregory_weights(ntau,m);
      element_incr(size1,res,weight,ATV,bvttmp);
      ATV+=es;
      bvttmp+=es;
    }
    else{
      for(l=0;l<=k;l++){
        weight = I.omega(l);
        element_incr(size1, res,weight,ATV,bvttmp);
        ATV+=es;
        bvttmp+=es;
      }
      for(l=k+1;l<ntau-k;l++){
        element_incr(size1, res,ATV,bvttmp);
        ATV+=es;
        bvttmp+=es;
      }
      for(l=ntau-k;l<=ntau;l++){
        weight = I.omega(ntau-l);
        element_incr(size1, res,weight,ATV,bvttmp);
        ATV+=es;
        bvttmp+=es;
      }
    }
    res+=es;
  }

  delete[] BVT;
}


}//namespace
#endif
