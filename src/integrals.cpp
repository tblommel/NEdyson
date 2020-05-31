#ifndef INTEGRALS_IMPL
#define INTEGRALS_IMPL

#include "integrals.h"

namespace NEdyson{


// This one does integral for every tau, and puts result into C by incrementing it!!!
// For every \tau=m...
//   C[n,m] += \int_0^n dt A^R(n,t) B^{TV}(t,m)
//          += \sum_{j=0}^{n-1} dt w_{n,j} A^R(n,j) B^{TV}(j,m)
void CTV1(const INTEG &I, cplx *ctv, const GREEN &A, const GREEN &Acc, const GREEN &B, int n, double dt){
  assert(n>=I.k());
  assert(n<=A.nt());
  assert(A.nt() == B.nt());
  assert(A.nt() == Acc.nt());
  assert(A.size1() == Acc.size1());
  assert(A.size1() == B.size1());
  assert(A.ntau() == B.ntau());
  assert(A.ntau() == Acc.ntau());
  assert(A.sig()==Acc.sig());

  int k=I.k(), ntau=A.ntau(), size1=A.size1(), es=A.element_size(), i, j, m;
  double weight;
  cplx *Atmp = new cplx[es];
  cplx *resptr;
  cplx *Gptr;
  int ntop = (n>k)?n:k;

  for(j=0;j<=ntop;j++){
    weight = I.gregory_weights(n,j);
    if(j>n){ // Dont have AR
      element_conj(size1,Atmp,Acc.retptr(j,n));
      element_smul(size1,Atmp,-1);
    }
    else{
      element_set(size1,Atmp,A.retptr(n,j));
    }
    element_smul(size1,Atmp,dt);

    resptr=ctv;
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

// Does the integral C^{TV}_2[A,B](n,m) = \int_0^m d\tau A^{TV(n,\tau) B^M(tau-m)
//           m<k                        = \sum_{j,l=0}^k d\tau R_{m,j,l} A^{TV}(n,l) s B^M(ntau-j)
//           m>=k                       = \sum_{l=0}^m  d\tau w_{m,l} A^{TV}(n,m-l) s B^M(ntau-l)
// Puts the result into res, which should be size1*size1
void CTV2(const INTEG &I, const GREEN &A, const GREEN &B, int n, int m, double beta, cplx *res){
  assert(A.ntau()==B.ntau());
  assert(A.ntau()>=I.k());
  assert(A.ntau()>=m);
  assert(m>=0);
  assert(n>=0);
  assert(A.nt()>=n);
  assert(A.nt() == B.nt());
  assert(A.size1() == B.size1());

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


// Does the integral C^{TV}_3[A,B](n,m) = \int_m^\beta d\tau A^{TV(n,\tau) B^M(tau-m)
//           m>ntau-k                   = \sum_{j,l=0}^k d\tau R_{ntau-m,j,l} A^{TV}(n,ntau-l) s B^M(j)
//           m<=k                       = \sum_{l=0}^{ntau-m}  d\tau w_{ntau-m,l} A^{TV}(n,m+l) s B^M(l)
// Puts the result into res, which should be size1*size1
void CTV3(const INTEG &I, const GREEN &A, const GREEN &B, int n, int m, double beta, cplx *res){
  assert(A.ntau()==B.ntau());
  assert(A.ntau()>=I.k());
  assert(A.ntau()>=m);
  assert(m>=0);
  assert(n>=0);
  assert(A.nt()>=n);
  assert(A.nt() == B.nt());
  assert(A.size1() == B.size1());
  
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



// This function calls CTV1,2,3 to do the convolutions necessary for the timestepping of the TV component
// Results get put into C.tv(tstp,m) m=0...ntau.  Not by increment
void Ctv_tstp(int tstp, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I, double beta, double dt) {

  int ntau = C.ntau(), size1 = C.size1(), es = size1*size1;
  cplx *ctv = new cplx[(ntau+1)*es];
  cplx *tmp = new cplx[es];
  
  for(int m=0;m<=ntau;m++){
    CTV2(I,A,B,tstp,m,beta,ctv+m*es);
    CTV3(I,A,B,tstp,m,beta,tmp);
    element_incr(size1, ctv+m*es, tmp);
  }

  CTV1(I,ctv,A,Acc,B,tstp,dt);

  for(int m=0;m<=ntau;m++){
    element_set(size1, C.tvptr(tstp,m), ctv+m*es);
  }

  delete[] ctv;
  delete[] tmp;
}






// Does the les_lesadv integration for every n, puts into Q via increment!!!
// for n=j1...j2
//   computes the integral C_2^<[A,B](n,m) = \int_0^m dt A^<(n,t) B^A(t,m)
//                                         = \sum_{j=0}^{max(k,m)} w_{m,j} A^<(n,j) B^A(j,m)
// places into Q which should be size size1*size1*(j2-j1+1)
void Cles2_tstp(int j1, int j2, const INTEG &I, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double dt, cplx *Q){
  assert(j1<=j2);
  assert(j1>=0);
  assert(j2<=A.nt());
  assert(m<=A.nt());
  assert(A.nt()==Acc.nt());
  assert(A.nt()==B.nt());
  assert(A.nt()==Bcc.nt());
  assert(A.size1()==Acc.size1());
  assert(A.size1()==B.size1());
  assert(A.size1()==Bcc.size1());
  assert(A.sig()==Acc.sig());
  assert(B.sig()==Bcc.sig());
  
  int k=I.k(), size1=A.size1(), es=A.element_size(), sig=A.sig(),j,l,n,top=(m>k) ? (m):(k);
  double w;
  cplx *BA=new cplx [(top+1)*es];
  cplx *AL=new cplx [es];

  // First fill BA, since they are the same for each n
  for(j=0;j<=top;j++){//fill BA from B.retptr.  BA_jm=(BR_mj)^T*
    w=I.gregory_weights(m,j);
    if(m>=j){//we have BR
      element_conj(size1,BA+j*es,Bcc.retptr(m,j));
      element_smul(size1,BA+j*es,dt*w);
    }
    else{//dont have BR
      element_set(size1,BA+j*es,B.retptr(j,m));
      element_smul(size1,BA+j*es,-dt*w);
    }
  }

  cplx *res = Q;
  //sum from j=0 to n-1, where we don't have AL
  for(n=j1;n<=j2;n++){
    for(j=0;j<n;j++){
      element_conj(size1,AL,Acc.lesptr(j,n));
      element_incr(size1,res,-1.,AL,BA+j*es);
    }
    res+=es;
  }

  res = Q;
  //sum from j=n to max(k,m), here we have AL
  for(n=j1;n<=j2;n++){
    for(j=n;j<=top;j++){
      element_incr(size1,res,A.lesptr(n,j),BA+j*es);
    }
    res+=es;
  }

  delete[] BA;
  delete[] AL;
}


void Cles2_tstp(const INTEG &I, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double dt, cplx *Q){
  assert(m<=A.nt());
  assert(A.nt()==Acc.nt());
  assert(A.nt()==B.nt());
  assert(A.nt()==Bcc.nt());
  assert(A.size1()==Acc.size1());
  assert(A.size1()==B.size1());
  assert(A.size1()==Bcc.size1());
  assert(A.sig()==Acc.sig());
  assert(B.sig()==Bcc.sig());
  
  int k=I.k();
  int num = m>=k?m:k;
  return Cles2_tstp(0,num, I, A, Acc, B, Bcc, m, dt, Q);
}


// Does the les_tvvt integral for each n=j1...j2 and puts into Q via increment!!!!
// computes C_3^<[A,B](n,m) = -i \int_0^\beta d\tau A^{TV}(n,\tau) B^{VT}(\tau,m)
//                          = -i \dau \sum_{j=0}^{ntau} w_{ntau,j} A^{TV}_{n,j} B^{VT}_{j,m}
// places into Q via increment.  Q should be size size1*size1*(j2-j1+1)
void Cles3_tstp(int j1, int j2, const INTEG &I, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double beta, cplx *Q){
  assert(j1<=j2);
  assert(j1>=0);
  assert(j2<=A.nt());
  assert(m<=A.nt());
  assert(A.nt()==Acc.nt());
  assert(A.nt()==B.nt());
  assert(A.nt()==Bcc.nt());
  assert(A.size1()==Acc.size1());
  assert(A.size1()==B.size1());
  assert(A.size1()==Bcc.size1());
  assert(A.sig()==Acc.sig());
  assert(B.sig()==Bcc.sig());
  assert(A.ntau()>=I.k());
  assert(A.ntau()==B.ntau());
  assert(A.ntau()==Bcc.ntau());
  assert(A.ntau()==Acc.ntau());

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
  for(j=0;j<=ntau;j++) element_conj(size1, BVT+j*es, Bcc.tvptr(m,ntau-j));
  for(l=0;l<(ntau+1)*es;l++) BVT[l]*=weight;

  //do the integral
  cplx *res=Q;
  for(n=j1;n<=j2;n++){
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



void Cles3_tstp(const INTEG &I, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double beta, cplx *Q){
  assert(m<=A.nt());
  assert(A.nt()==Acc.nt());
  assert(A.nt()==B.nt());
  assert(A.nt()==Bcc.nt());
  assert(A.size1()==Acc.size1());
  assert(A.size1()==B.size1());
  assert(A.size1()==Bcc.size1());
  assert(A.sig()==Acc.sig());
  assert(B.sig()==Bcc.sig());
  assert(A.ntau()>=I.k());
  assert(A.ntau()==B.ntau());
  assert(A.ntau()==Bcc.ntau());
  assert(A.ntau()==Acc.ntau());
  
  int k=I.k();
  int top = m<k?k:m;
  return Cles3_tstp(0, top, I, A, Acc, B, Bcc, m, beta, Q);
}


// This does the Integral \int_t'^tstp d\bar{t} A^R(tstp,\bar{t}) B^R(\bar{t},t')
// For every value t' from 0 to tstp
// Places result into C(tstp,t') via increment!!
void incr_convolution_ret(int tstp, const std::vector<bool> &mask_ret, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I , double dt){
  assert(A.nt()==Acc.nt());
  assert(A.nt()==B.nt());
  assert(A.nt()==Bcc.nt());
  assert(A.nt()==C.nt());
  assert(tstp <= A.nt());
  assert(A.sig()==Acc.sig());
  assert(B.sig()==Bcc.sig());
  assert(A.size1()==Acc.size1());
  assert(A.size1()==B.size1());
  assert(A.size1()==Bcc.size1());
  assert(A.size1()==C.size1());
  assert(static_cast<int>(mask_ret.size())==tstp+1);

  int size1 = A.size1(), es = size1*size1, k = I.k(), ntop = (tstp >= k ? tstp : k);
  int n,m,j;
  cplx *ctmp = new cplx[es];
  cplx *aret;
  cplx *btmp = new cplx[es];
  cplx *aretf = 0;
  
  if(ntop == tstp){
    aret = A.retptr(tstp,0);
  }
  else{ // Integration not causal, need to conjugate some values
    aretf = new cplx[(ntop+1)*es];
    aret = aretf;
    for(j=0;j<=tstp;j++){
      element_set(size1, aretf+j*es, A.retptr(tstp, j));
    }
    for(j=tstp+1;j<=ntop;j++){
      element_conj(size1, aretf+j*es, Acc.retptr(j, tstp));
      element_smul(size1, aretf+j*es, -1);
    }
  }

  // Begin the integration. There are four cases
  for(m=0;m<=tstp;m++){
    if(mask_ret[m]){
      element_set_zero(size1, ctmp);
      if(tstp<k){ // use the I weights
        for(j=0;j<=k;j++){
          if(j>=m) element_incr(size1, ctmp, I.poly_integ(m,tstp,j), aret+j*es, B.retptr(j,m));
          else{
            element_conj(size1, btmp, Bcc.retptr(m,j));
            element_incr(size1, ctmp, -I.poly_integ(m,tstp,j), aret+j*es, btmp);
          }
        }
      }
      else if(tstp-m<k){ // use information from behind m
        for(j=0;j<=k;j++){
          if(tstp-j>=m) element_incr(size1, ctmp, I.gregory_weights(tstp-m,j), aret+(tstp-j)*es, B.retptr(tstp-j,m));
          else{
            element_conj(size1, btmp, Bcc.retptr(m,tstp-j));
            element_incr(size1, ctmp, -I.gregory_weights(tstp-m,j), aret+(tstp-j)*es, btmp);
          }
        }
      }
      else if(tstp-m <= 2*k + 2){ // Everything fine, just use Sigma weights
        for(j=0;j<=tstp-m;j++) element_incr(size1, ctmp, I.gregory_weights(tstp-m,j), aret+(j+m)*es, B.retptr(j+m,m));
      }
      else{ // Everything fine, use omega and 1
        for(j=m;j<=m+k;j++) element_incr(size1, ctmp, I.omega(j-m), aret+j*es, B.retptr(j,m));        for(j=m+k+1;j<tstp-k;j++) element_incr(size1, ctmp, aret+j*es, B.retptr(j,m));
        for(j=tstp-k;j<=tstp;j++) element_incr(size1, ctmp, I.omega(tstp-j), aret+j*es, B.retptr(j,m));
      }
      element_incr(size1, C.retptr(tstp, m), dt, ctmp);
    }
  }
  delete[] ctmp;
  delete[] btmp;
  if(aretf != 0) delete[] aretf;
  return;
}


// This function does the integrals required by the TV propagation equation and then places them into C^{TV}(tstp,m) for m=0...ntau via increment!!!
void incr_convolution_tv(int tstp, const std::vector<bool> &mask_tv, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I, double beta, double dt){
  assert(A.nt()==Acc.nt());
  assert(A.nt()==B.nt());
  assert(A.nt()==Bcc.nt());
  assert(A.nt()==C.nt());
  assert(A.ntau()==Acc.ntau());
  assert(A.ntau()==B.ntau());
  assert(A.ntau()==Bcc.ntau());
  assert(A.ntau()==C.ntau());
  assert(A.ntau()>=I.k());
  assert(tstp <= A.nt());
  assert(A.sig()==Acc.sig());
  assert(B.sig()==Bcc.sig());
  assert(A.size1()==Acc.size1());
  assert(A.size1()==B.size1());
  assert(A.size1()==Bcc.size1());
  assert(A.size1()==C.size1());
  assert(static_cast<int>(mask_tv.size())==C.ntau()+1);

  int ntau = A.ntau(), size1 = A.size1(), es = size1*size1, k = I.k(), ntop = tstp>=k?tstp:k;
  
  int j,m,n;
  cplx *ctmp1 = new cplx[es];
  cplx *ctmp2 = new cplx[es];
  cplx *bmat = B.matptr(0);
  cplx *aret;
  cplx *aretf=0;
  if(ntop == tstp) aret = A.retptr(tstp, 0);
  else{
    aretf = new cplx[(ntop+1)*es];
    aret = aretf;
    for(j=0;j<=tstp;j++) element_set(size1, aretf+j*es, A.retptr(tstp, j));
    for(j=tstp+1;j<=ntop;j++){
      element_conj(size1, aretf+j*es, Acc.retptr(j,tstp));
      element_smul(size1, aretf+j*es, -1);
    }
  }

  for(m=0;m<=ntau;m++){
    if(mask_tv[m]){
      CTV2(I, A, B, tstp, m, beta, ctmp1);
      CTV3(I, A, B, tstp, m, beta, ctmp2);
      element_incr(size1, ctmp1, ctmp2);
      element_set_zero(size1,ctmp2);
      if(tstp<=2*k+2){
        for(n=0;n<=ntop;n++){
          element_incr(size1, ctmp2, I.gregory_weights(tstp, n), aret+n*es, B.tvptr(n,m));
        }
      }
      else{
        for(n=0;n<=k;n++) element_incr(size1, ctmp2, I.omega(n), aret+n*es, B.tvptr(n,m));
        for(n=k+1;n<tstp-k;n++) element_incr(size1, ctmp2, aret+n*es, B.tvptr(n,m));
        for(n=tstp-k;n<=tstp;n++) element_incr(size1, ctmp2, I.omega(tstp-n), aret+n*es, B.tvptr(n,m));
      }
      element_smul(size1, ctmp2, dt);
      element_incr(size1, C.tvptr(tstp,m),ctmp1);
      element_incr(size1, C.tvptr(tstp,m),ctmp2);
    }
  }

  delete[] ctmp1;
  delete[] ctmp2;
  if(aretf!=0) delete[] aretf;
  return;
}


void incr_convolution_les(int tstp, const std::vector<bool> &mask_les, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I, double beta, double dt){
  assert(A.nt()==Acc.nt());
  assert(A.nt()==B.nt());
  assert(A.nt()==Bcc.nt());
  assert(A.nt()==C.nt());
  assert(A.ntau()==Acc.ntau());
  assert(A.ntau()==B.ntau());
  assert(A.ntau()==Bcc.ntau());
  assert(A.ntau()==C.ntau());
  assert(A.ntau()>=I.k());
  assert(tstp <= A.nt());
  assert(A.sig()==Acc.sig());
  assert(B.sig()==Bcc.sig());
  assert(A.size1()==Acc.size1());
  assert(A.size1()==B.size1());
  assert(A.size1()==Bcc.size1());
  assert(A.size1()==C.size1());
  assert(static_cast<int>(mask_les.size())==tstp+1);

  int ntau = A.ntau(), size1 = A.size1(), es = size1*size1, k = I.k(), ntop = tstp>=k?tstp:k;
  double dtau = beta/ntau;
  cplx cplxi = cplx(0,1);
  int n,j,m,l;

  // Fill list of Bs to use in integration
  cplx *badv = new cplx[(ntop+1)*es];
  cplx *bles = new cplx[(ntop+1)*es];
  cplx *bvt = new cplx[(ntau+1)*es];
  cplx *ctmp1 = new cplx[es];
  cplx *ctmp2 = new cplx[es];
  cplx *ctmp3 = new cplx[es];

  for(j=0;j<=tstp;j++){
    element_conj(size1,badv+j*es,Bcc.retptr(tstp,j));
    element_set(size1,bles+j*es,B.lesptr(j, tstp));
  }
  for(j=tstp+1;j<=ntop;j++){
    element_set(size1, badv+j*es,B.retptr(j,tstp));
    element_smul(size1, badv+j*es,-1);
    element_conj(size1, bles+j*es,Bcc.lesptr(tstp, j));
    element_smul(size1, bles+j*es,-1);
  }
  if(B.sig() == -1){
    for(m=0;m<=ntau;m++) element_conj(size1, bvt+m*es, Bcc.tvptr(tstp, ntau-m));
  }
  else{
    for(m=0;m<=ntau;m++){
      element_conj(size1, bvt+m*es, Bcc.tvptr(tstp, ntau-m));
      element_smul(size1, bvt+m*es, -1);
    }
  }

  // Do the integrals for every t = 0...tstp
  for(n=0;n<=tstp;n++){
    if(mask_les[n]){
      // integral from 0 to t AR BL
      int nup = n>k?n:k;
      element_set_zero(size1, ctmp1);
      if(n<=2*k+2){ // start and sig weights
        for(j=0;j<=nup;j++){
          if(j<=n) element_incr(size1, ctmp1, I.gregory_weights(n,j), A.retptr(n,j), bles+j*es);
          else{
            element_conj(size1, ctmp2, Acc.retptr(j,n));
            element_incr(size1, ctmp1, -I.gregory_weights(n,j), ctmp2, bles+j*es);
          }
        }
      }
      else{ // omegas and 1s
        for(j=0;j<=k;j++) element_incr(size1, ctmp1, I.omega(j), A.retptr(n,j), bles+j*es);
        for(j=k+1;j<n-k;j++) element_incr(size1, ctmp1, A.retptr(n,j), bles+j*es);
        for(j=n-k;j<=n;j++) element_incr(size1, ctmp1, I.omega(n-j), A.retptr(n,j),bles+j*es);
      }


      // integral from 0 to tstp AL BA
      element_set_zero(size1, ctmp3);
      if(tstp<=2*k+2){ // start and sig weights
        for(j=0;j<=n;j++){
          element_conj(size1, ctmp2, Acc.lesptr(j,n));
          element_incr(size1, ctmp3, -I.gregory_weights(tstp, j), ctmp2, badv+j*es);
        }
        for(j=n+1;j<=ntop;j++) element_incr(size1, ctmp3, I.gregory_weights(tstp, j), A.lesptr(n,j), badv+j*es);
      }
      else{ // omegas and 1s
        for(j=0;j<=k;j++){
          if(j>=n) element_incr(size1, ctmp3, I.omega(j), A.lesptr(n,j), badv+j*es);
          else{
            element_conj(size1, ctmp2, Acc.lesptr(j,n));
            element_incr(size1, ctmp3, -I.omega(j), ctmp2, badv+j*es);
          }
        }

        if(n<=k){ // already encountered cut
          for(j=k+1;j<tstp-k;j++) element_incr(size1, ctmp3, A.lesptr(n,j), badv+j*es);
        }
        else if(n<tstp-k){ // cut is in this range
          for(j=k+1;j<n;j++){
            element_conj(size1, ctmp2, Acc.lesptr(j,n));
            element_incr(size1, ctmp3, -1, ctmp2, badv+j*es);
          }
          for(j=n;j<tstp-k;j++) element_incr(size1, ctmp3, A.lesptr(n,j), badv+j*es);
        }
        else{ // cut is past this range
          for(j=k+1;j<tstp-k;j++){
            element_conj(size1, ctmp2, Acc.lesptr(j,n));
            element_incr(size1, ctmp3, -1, ctmp2, badv+j*es);
          }
        }

        for(j=tstp-k;j<=tstp;j++){
          if(j>=n) element_incr(size1, ctmp3, I.omega(tstp-j), A.lesptr(n,j), badv+j*es);
          else{
            element_conj(size1, ctmp2, Acc.lesptr(j,n));
            element_incr(size1, ctmp3, -I.omega(tstp-j), ctmp2, badv+j*es);
          }
        }
      }

      
      // integral from 0 to beta ATV BVT
      element_set_zero(size1, ctmp2);
      for(m=0;m<=k;m++) element_incr(size1, ctmp2, I.omega(m), A.tvptr(n, m), bvt+m*es);
      for(m=k+1;m<ntau-k;m++) element_incr(size1, ctmp2, A.tvptr(n,m), bvt+m*es);
      for(m=ntau-k;m<=ntau;m++) element_incr(size1, ctmp2, I.omega(ntau-m), A.tvptr(n,m), bvt+m*es);

      // Put into C
      element_incr(size1, C.lesptr(n,tstp), dt, ctmp1);
      element_incr(size1, C.lesptr(n,tstp), dt, ctmp3);
      element_incr(size1, C.lesptr(n,tstp), dtau*cplx(0,-1), ctmp2);
    }
  }
  delete[] badv;
  delete[] bles;
  delete[] bvt;
  delete[] ctmp1;
  delete[] ctmp2;
  delete[] ctmp3;
  return;
}


}//namespace
#endif
