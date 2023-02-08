#ifndef TTI_DYSON_INTEGRALS_IMPL
#define TTI_DYSON_INTEGRALS_IMPL

namespace NEdyson{

// returns (Sigma*G)^<(t,t)
double dyson::energy_conv(int tstp, const TTI_GREEN &Sig, const TTI_GREEN &G, double beta, double dt) const {
  assert(tstp <= Sig.nt());
  assert(G.nt() == Sig.nt());  
  assert(I.k() <= G.nt());
  assert(I.k() <= G.ntau());
  assert(G.sig() == Sig.sig());
  assert(G.ntau() == Sig.ntau());
  assert(G.size1() == Sig.size1());
  assert(G.size1() == nao_);
  assert(G.nt() == nt_);
  assert(G.ntau() == ntau_);

  int top = tstp >= k_ ? tstp : k_, sig = G.sig();
  cplx res1 = 0;
  cplx res2 = 0;
  
  // Sig^R * G^<
  for(int i=0; i<=tstp; i++) {
    res1 += I.gregory_weights(tstp, i) * ZMatrixConstMap(Sig.retptr(tstp-i), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.lesptr(i-tstp), nao_, nao_).transpose()).sum();
  }
  for(int i=tstp; i<=top; i++) {
    res1 += std::conj(I.gregory_weights(tstp, i) * ZMatrixConstMap(Sig.retptr(i-tstp), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.lesptr(tstp-i), nao_, nao_).transpose()).sum());
  }

  // Sig^< * G^A
  for(int i=0; i<=tstp; i++) {
    res1 -= std::conj(I.gregory_weights(tstp, i) * ZMatrixConstMap(Sig.lesptr(i-tstp), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.retptr(tstp-i), nao_, nao_).transpose()).sum());
  }
  for(int i=tstp; i<=top; i++) {
    res1 -= I.gregory_weights(tstp, i) * ZMatrixConstMap(Sig.lesptr(tstp-i), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.retptr(i-tstp), nao_, nao_).transpose()).sum();
  }

  res1 *= dt;

  // Sig^rm * G^lm
  res2 = Conv.energy(ZTensorView<3>(Sig.tvptr(tstp,0), ntau_+1, nao_, nao_),
                     ZTensorView<3>(G.tvptr(tstp,0), ntau_+1, nao_, nao_),
                     beta, (double) G.sig());

  res2 *= cplx(0., -1.);

  return (cplx(0.,-0.5)*(res1+res2)).real();
}


// This one does integral for every tau, and puts result into C by incrementing it!!!
// For every \tau=m...
//   C[n,m] += \int_0^n dt A^R(n,t) B^{TV}(t,m)
//          += \sum_{j=0}^{n-1} dt w_{n,j} A^R(n,j) B^{TV}(j,m)
void dyson::CTV1(cplx *ctv, const TTI_GREEN &A, const TTI_GREEN &Acc, const TTI_GREEN &B, int n, double dt) const {
  assert(A.nt() == Acc.nt());
  assert(A.nt() == B.nt());
  assert(A.nt() == nt_);
  assert(A.ntau() == Acc.ntau());
  assert(A.ntau() == B.ntau());
  assert(A.ntau() == ntau_);
  assert(A.size1() == Acc.size1());
  assert(A.size1() == B.size1());
  assert(A.size1() == nao_);
  assert(A.sig() == Acc.sig());
  assert(n >= 0);
  assert(n <= nt_);
  assert(nt_ >= k_);

  int i, j, m;
  int ntop = n > k_? n : k_;

  ZMatrixMap tmpMap = ZMatrixMap(tmp.data(), nao_, nao_);
  for(j=0; j<=ntop; j++) {

    if(j>n) { // here we do not have AR
      tmpMap = -dt * I.gregory_weights(n,j) * ZMatrixMap(Acc.retptr(j-n), nao_, nao_).adjoint();
      for(m=0; m<=ntau_; m++) {
        ZMatrixMap(ctv + m*es_, nao_, nao_).noalias() += tmpMap * ZMatrixMap(B.tvptr(j,m), nao_, nao_);
      }
    }

    else { // here we do
      ZMatrixMap AMap = ZMatrixMap(A.retptr(n-j), nao_, nao_);
      for(m=0; m<=ntau_; m++) {
        ZMatrixMap(ctv + m*es_, nao_, nao_).noalias() += dt*I.gregory_weights(n,j) * AMap * ZMatrixMap(B.tvptr(j,m), nao_, nao_);
      }
    }

  }
}
/*
//DEPRECIATED
// Does the integral C^{TV}_2[A,B](n,m) = \int_0^m d\tau A^{TV(n,\tau) B^M(tau-m)
//           m<k                        = \sum_{j,l=0}^k d\tau R_{m,j,l} A^{TV}(n,l) s B^M(ntau-j)
//           m>=k                       = \sum_{l=0}^m  d\tau w_{m,l} A^{TV}(n,m-l) s B^M(ntau-l)
// Puts the result into res, which should be size1*size1
void dyson::CTV2(const TTI_GREEN &A, const TTI_GREEN &B, int n, int m, double beta, cplx *res) const {
  assert(A.size1() == B.size1());
  assert(A.size1() == nao_);
  assert(A.nt() == B.nt());
  assert(A.nt() == nt_);
  assert(A.ntau() == B.ntau());
  assert(A.ntau() == ntau_);
  assert(n >= 0);
  assert(n <= nt_);
  assert(m >= 0);
  assert(m <= ntau_);
  assert(k_ <= ntau_);

  int j, l;
  memset(res, 0, es_*sizeof(cplx));

  cplx *BMptr = B.matptr(ntau_);
  cplx *ATVptr = A.tvptr(n,m);
  
  ZMatrixMap resMap = ZMatrixMap(res, nao_, nao_);
  
  if(m == 0) {}          // Do nothing
  else if(m < k_){       // Use starting weights
    for(j=0; j<=k_; j++) {
      ZMatrixMap MatMap = ZMatrixMap(B.matptr(ntau_-j), nao_, nao_);
      for(l=0; l<=k_; l++){
        resMap.noalias() += I.rcorr(m,j,l) * ZMatrixMap(A.tvptr(n,l), nao_, nao_) * MatMap;
      }
    }
  }
  else if(m < k_+k_+1){  // Use gregory weights
    for(l=0; l<=m; l++) {
      resMap.noalias() += I.gregory_weights(m,l) * ZMatrixMap(ATVptr, nao_, nao_) * ZMatrixMap(BMptr, nao_, nao_);
      ATVptr -= es_;
      BMptr -= es_;
    }
  }
  else{                  // Use greg and ones
    for(l=0; l<=k_; l++){// Omegas
      resMap.noalias() += I.omega(l) * ZMatrixMap(ATVptr, nao_, nao_) * ZMatrixMap(BMptr, nao_, nao_);
      ATVptr -= es_;
      BMptr -= es_;
    }
    for(l=k_+1; l<m-k_; l++){// Ones
      resMap.noalias() += ZMatrixMap(ATVptr, nao_, nao_) * ZMatrixMap(BMptr, nao_, nao_);
      ATVptr -= es_;
      BMptr -= es_;
    }
    for(l=m-k_; l<=m; l++){// Omegas
      resMap.noalias() += I.omega(m-l) * ZMatrixMap(ATVptr, nao_, nao_) * ZMatrixMap(BMptr, nao_, nao_);
      ATVptr -= es_;
      BMptr -= es_;
    }
  }

  resMap *= B.sig()*(beta/ntau_);
}

//DEPRECIATED
// Does the integral C^{TV}_3[A,B](n,m) = \int_m^\beta d\tau A^{TV(n,\tau) B^M(tau-m)
//           m>ntau-k                   = \sum_{j,l=0}^k d\tau R_{ntau-m,j,l} A^{TV}(n,ntau-l) s B^M(j)
//           m<=k                       = \sum_{l=0}^{ntau-m}  d\tau w_{ntau-m,l} A^{TV}(n,m+l) s B^M(l)
// Puts the result into res, which should be size1*size1
void dyson::CTV3(const TTI_GREEN &A, const TTI_GREEN &B, int n, int m, double beta, cplx *res) const {
  assert(A.size1() == B.size1());
  assert(A.size1() == nao_);
  assert(A.nt() == B.nt());
  assert(A.nt() == nt_);
  assert(A.ntau() == B.ntau());
  assert(A.ntau() == ntau_);
  assert(n >= 0);
  assert(n <= nt_);
  assert(m >= 0);
  assert(m <= ntau_);
  assert(k_ <= ntau_);
  
  int j, l;
  memset(res, 0, es_*sizeof(cplx));
  
  ZMatrixMap resMap = ZMatrixMap(res, nao_, nao_);

  cplx *BMptr;
  cplx *ATVptr;

  if(m == ntau_){}                 // Do nothing
  else if(m > ntau_-k_) {          // Use starting weights
    BMptr = B.matptr(0);
    ATVptr = A.tvptr(n,ntau_);
    for(j=0; j<=k_; j++) {
      ZMatrixMap BMMap = ZMatrixMap(BMptr+j*es_, nao_, nao_);
      for(l=0; l<=k_; l++) {
        resMap.noalias() += I.rcorr(ntau_-m,j,l) * ZMatrixMap(ATVptr-l*es_, nao_, nao_) * BMMap;
      }
    }
  }
  else if(ntau_-m >= k_+k_+1){     // Use greg and ones
    BMptr = B.matptr(0);
    ATVptr = A.tvptr(n,m);
    for(l=0; l<=k_; l++) {
      resMap.noalias() += I.omega(l) * ZMatrixMap(ATVptr, nao_, nao_) * ZMatrixMap(BMptr, nao_, nao_);
      ATVptr += es_;
      BMptr += es_;
    }
    for(l=k_+1; l<ntau_-m-k_; l++) {
      resMap.noalias() += ZMatrixMap(ATVptr, nao_, nao_) * ZMatrixMap(BMptr, nao_, nao_);
      ATVptr += es_;
      BMptr += es_;
    }
    for(l=ntau_-m-k_; l<=ntau_-m; l++) {
      resMap.noalias() += I.omega(ntau_-m-l) * ZMatrixMap(ATVptr, nao_, nao_) * ZMatrixMap(BMptr, nao_, nao_);
      ATVptr += es_;
      BMptr += es_;
    }
  }
  else {                           // Use greg weights
    BMptr = B.matptr(0);
    ATVptr = A.tvptr(n,m);
    for(l=0; l<=ntau_-m; l++) {
      resMap.noalias() += I.gregory_weights(ntau_-m,l) * ZMatrixMap(ATVptr, nao_, nao_) * ZMatrixMap(BMptr, nao_, nao_);
      ATVptr += es_;
      BMptr += es_;
    }
  }
  
  resMap *= (beta/ntau_);
}
*/


// This function calls CTV1,2,3 to do the convolutions necessary for the timestepping of the TV component
// Results get put into C.tv(tstp,m) m=0...ntau.  Not by increment
void dyson::Ctv_tstp(int tstp, TTI_GREEN &C, const TTI_GREEN &A, const TTI_GREEN &Acc, const TTI_GREEN &B, const TTI_GREEN &Bcc, double beta, double dt) const {
  assert(C.nt() == A.nt());
  assert(C.nt() == Acc.nt());
  assert(C.nt() == B.nt());
  assert(C.nt() == Bcc.nt());
  assert(C.nt() == nt_);
  assert(C.ntau() == A.ntau());
  assert(C.ntau() == Acc.ntau());
  assert(C.ntau() == B.ntau());
  assert(C.ntau() == Bcc.ntau());
  assert(C.ntau() == ntau_);
  assert(C.size1() == A.size1());
  assert(C.size1() == Acc.size1());
  assert(C.size1() == B.size1());
  assert(C.size1() == Bcc.size1());
  assert(C.size1() == nao_);
  assert(A.sig() == Acc.sig());
  assert(B.sig() == Bcc.sig());
  assert(tstp >= 0);
  assert(tstp <= nt_);
  assert(A.nt() >= k_);
  assert(A.ntau() >= k_);

  ZMatrixMap tmpMap = ZMatrixMap(tmp.data(), nao_, nao_);

//  for(int m=0; m<=ntau_; m++) {
//    CTV2(A, B, tstp, m, beta, NTauTmp.data() + m*es_);
//    CTV3(A, B, tstp, m, beta, tmp.data());
//    ZMatrixMap(NTauTmp.data() + m*es_, nao_, nao_).noalias() += tmpMap;
//  }
  auto ZTVNTT = ZTensorView<3>(NTauTmp.data(), ntau_+1, nao_, nao_);
  Conv.mixing(ZTVNTT,
              ZTensorView<3>(A.tvptr(tstp, 0), ntau_+1, nao_, nao_),
              ZTensorView<3>(B.matptr(0), ntau_+1, nao_, nao_),
              beta, (double)B.sig());

  // This does the A^R B^TV convolution and puts into first argument via increment
  CTV1(NTauTmp.data(), A, Acc, B, tstp, dt);

  // Copy over the results into C.tv(tstp,m)
  memcpy(C.tvptr(tstp,0), NTauTmp.data(), (ntau_+1)*es_*sizeof(cplx));
}



}//namespace
#endif
