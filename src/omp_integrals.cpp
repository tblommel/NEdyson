#ifndef OMP_INTEGRALS_IMPL
#define OMP_INTEGRALS_IMPL

#include "omp_integrals.h"

namespace NEdyson{

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
      ZMatrixMap(aretf+j*es, size1, size1).noalias() = ZMatrixMap(A.retptr(tstp,j), size1, size1);
    }
    for(j=tstp+1;j<=ntop;j++){
      ZMatrixMap(aretf+j*es, size1, size1).noalias() = -ZMatrixMap(Acc.retptr(j,tstp), size1, size1).adjoint();
    }
  }

  // Begin the integration. There are four cases
  for(m=0;m<=tstp;m++){
    if(mask_ret[m]){
      element_set_zero(size1, ctmp);
      if(tstp<k){ // use the I weights
        for(j=0;j<=k;j++){
          if(j>=m) ZMatrixMap(ctmp, size1, size1).noalias() += I.poly_integ(m,tstp,j) * ZMatrixMap(aret+j*es, size1, size1) * ZMatrixMap(B.retptr(j,m), size1, size1);
          else{
            ZMatrixMap(ctmp, size1, size1).noalias() -= I.poly_integ(m,tstp,j) * ZMatrixMap(aret+j*es, size1, size1) * ZMatrixMap(Bcc.retptr(m,j), size1, size1).adjoint();
          }
        }
      }
      else if(tstp-m<k){ // use information from behind m
        for(j=0;j<=k;j++){
          if(tstp-j>=m)  ZMatrixMap(ctmp, size1, size1).noalias() += I.gregory_weights(tstp-m,j) * ZMatrixMap(aret+(tstp-j)*es, size1, size1) * ZMatrixMap(B.retptr(tstp-j,m), size1, size1);
          else{
            ZMatrixMap(ctmp, size1, size1).noalias() -= I.gregory_weights(tstp-m,j) * ZMatrixMap(aret+(tstp-j)*es, size1, size1) * ZMatrixMap(Bcc.retptr(tstp-j,m), size1, size1).adjoint();
          }
        }
      }
      else if(tstp-m <= 2*k + 2){ // Everything fine, just use Sigma weights
        for(j=0;j<=tstp-m;j++)  ZMatrixMap(ctmp, size1, size1).noalias() += I.gregory_weights(tstp-m,j) * ZMatrixMap(aret+(j+m)*es, size1, size1) * ZMatrixMap(B.retptr(j+m,m), size1, size1);
      }
      else{ // Everything fine, use omega and 1
        for(j=m;j<=m+k;j++)  ZMatrixMap(ctmp, size1, size1).noalias() += I.omega(j-m) * ZMatrixMap(aret+j*es, size1, size1) * ZMatrixMap(B.retptr(j,m), size1, size1);
        for(j=m+k+1;j<tstp-k;j++)  ZMatrixMap(ctmp, size1, size1).noalias() += ZMatrixMap(aret+j*es, size1, size1) * ZMatrixMap(B.retptr(j,m), size1, size1);
        for(j=tstp-k;j<=tstp;j++) ZMatrixMap(ctmp, size1, size1).noalias() += I.omega(tstp-j) * ZMatrixMap(aret+j*es, size1, size1) * ZMatrixMap(B.retptr(j,m), size1, size1);
      }
      ZMatrixMap(C.retptr(tstp, m), size1, size1).noalias() += dt*ZMatrixMap(ctmp, size1, size1);
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
    for(j=0;j<=tstp;j++) ZMatrixMap(aretf+j*es, size1, size1).noalias() = ZMatrixMap(A.retptr(tstp,j), size1, size1);
    for(j=tstp+1;j<=ntop;j++){
      ZMatrixMap(aretf+j*es, size1, size1).noalias() = -ZMatrixMap(Acc.retptr(j,tstp), size1, size1).adjoint();
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
