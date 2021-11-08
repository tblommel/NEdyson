#ifndef DYSON_INTEGRALS_IMPL
#define DYSON_INTEGRALS_IMPL

namespace NEdyson{

// returns E_{int}(t) = \pm \frac{i}{2} * (
//                       \int_0^t du S^R(t,u) G^<(u,t)
//                    +  \int_0^t du S^<(t,u) G^A(u,t)
//                    -i \int_0^\beta d\tau S^\rciel(t,\tau) G^\lceil(\tau,t) )
double dyson::energy_conv(int tstp, const GREEN &Sig, const GREEN &G, double beta, double dt) const {
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
    res1 += I.gregory_weights(tstp, i) * ZMatrixConstMap(Sig.retptr(tstp,i), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.lesptr(i,tstp), nao_, nao_).transpose()).sum();
  }
  for(int i=tstp+1; i<=top; i++) {
    res1 += std::conj(I.gregory_weights(tstp, i) * ZMatrixConstMap(Sig.retptr(i,tstp), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.lesptr(tstp,i), nao_, nao_).transpose()).sum());
  }

  // Sig^< * G^A
  for(int i=0; i<=tstp; i++) {
    res1 -= std::conj(I.gregory_weights(tstp, i) * ZMatrixConstMap(Sig.lesptr(i,tstp), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.retptr(tstp, i), nao_, nao_).transpose()).sum());
  }
  for(int i=tstp+1; i<=top; i++) {
    res1 -= I.gregory_weights(tstp, i) * ZMatrixConstMap(Sig.lesptr(tstp, i), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.retptr(i, tstp), nao_, nao_).transpose()).sum();
  }

  res1 *= dt;

  // Sig^rm * G^lm
//  OLD CODE FOR UNIFORM MESH
//  for(int i=0; i <= ntau_; i++) {
//    res2 += I.gregory_weights(ntau_, i) * ZMatrixConstMap(Sig.tvptr(tstp, i), nao_, nao_).cwiseProduct(ZMatrixConstMap(G.tvptr(tstp, ntau_-i), nao_, nao_).conjugate()).sum();
//  }
  res2 = Conv.energy(ZTensorView<3>(Sig.tvptr(tstp,0), ntau_+1, nao_, nao_),
                     ZTensorView<3>(G.tvptr(tstp,0), ntau_+1, nao_, nao_),
                     beta, (double) G.sig());

  res2 *= cplx(0., -1.);

  return (cplx(0., sig * 0.5)*(res1+res2)).real();
}

// This one does integral for every tau, and puts result into C by incrementing it!!!
// For every \tau \equiv m...
//   C[n,m] += \int_0^n dt A^R(n,t) B^{TV}(t,m)
//          += \sum_{j=0}^{n} dt w_{n,j} A^R(n,j) B^{TV}(j,m)
// When time evolving the dyson eqn. We set B^{TV}(tstp, ...) = 0 in tstp_tv function 
// because this is what we are solving for.  This missing term in the integral is accounted
// for on the LHS of the matrix equation used to solve for B^{TV}(tstp, ...)
void dyson::CTV1(cplx *ctv, const GREEN &A, const GREEN &Acc, const GREEN &B, int n, double dt) const {
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
  int ntop = (n>k_) ? n : k_;

  ZMatrixMap tmpMap = ZMatrixMap(tmp.data(), nao_, nao_);
  for(j=0; j<=ntop; j++) {

    if(j>n) { // here we do not have AR
      tmpMap = -dt * I.gregory_weights(n,j) * ZMatrixMap(Acc.retptr(j,n), nao_, nao_).adjoint();
      for(m=0; m<=ntau_; m++) {
        ZMatrixMap(ctv + m*es_, nao_, nao_).noalias() += tmpMap * ZMatrixMap(B.tvptr(j,m), nao_, nao_);
      }
    }

    else { // here we do
      ZMatrixMap AMap = ZMatrixMap(A.retptr(n,j), nao_, nao_);
      for(m=0; m<=ntau_; m++) {
        ZMatrixMap(ctv + m*es_, nao_, nao_).noalias() += dt*I.gregory_weights(n,j) * AMap * ZMatrixMap(B.tvptr(j,m), nao_, nao_);
      }
    }

  }
}

// DEPRECIATED. WE NOW USE LEGENDRE GRID
// Does the integral C^{TV}_2[A,B](n,m) = \int_0^m d\tau A^{TV(n,\tau) B^M(tau-m)
//           m<k                        = \sum_{j,l=0}^k d\tau R_{m,j,l} A^{TV}(n,l) s B^M(ntau-j)
//           m>=k                       = \sum_{l=0}^m  d\tau w_{m,l} A^{TV}(n,m-l) s B^M(ntau-l)
// Puts the result into res, which should be size1*size1
void dyson::CTV2(const GREEN &A, const GREEN &B, int n, int m, double beta, cplx *res) const {
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

// DEPRECIATED.  WE NOW USE LEGENDRE GRID
// Does the integral C^{TV}_3[A,B](n,m) = \int_m^\beta d\tau A^{TV(n,\tau) B^M(tau-m)
//           m>ntau-k                   = \sum_{j,l=0}^k d\tau R_{ntau-m,j,l} A^{TV}(n,ntau-l) s B^M(j)
//           m<=k                       = \sum_{l=0}^{ntau-m}  d\tau w_{ntau-m,l} A^{TV}(n,m+l) s B^M(l)
// Puts the result into res, which should be size1*size1
void dyson::CTV3(const GREEN &A, const GREEN &B, int n, int m, double beta, cplx *res) const {
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



// This function calls CTV1,2,3 to do the convolutions necessary for the timestepping of the TV component
// Results get put into C.tv(tstp,m) m=0...ntau.  Not by increment
void dyson::Ctv_tstp(int tstp, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, double beta, double dt) const {
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

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;



  // First do the A^TV B^M convolution.  This gets put into temporary storage since CTV1 may need access to C^TV(tstp, m)
  start = std::chrono::system_clock::now();
//  This is the old implementation with the equidistant grid
//  Now use Xinyang Legendre convolutions
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
  end = std::chrono::system_clock::now();
  //TIMING
  elapsed_seconds = end-start;
  std::ofstream out1;
  std::string timing_data_dir = std::string(TIMING_DATA_DIR);
  out1.open(timing_data_dir + "Nao" + std::to_string(C.size1()) + "Nt" + std::to_string(C.nt()) + "Ntau" + std::to_string(C.ntau()) + "tv_int_tvm.dat", std::ofstream::app);
  out1 << elapsed_seconds.count() << "\n" ;
  out1.close();
  // TIMING


  // This does the A^R B^TV convolution and puts into first argument via increment
  start = std::chrono::system_clock::now();
  CTV1(NTauTmp.data(), A, Acc, B, tstp, dt);
  end = std::chrono::system_clock::now();
  //TIMING
  elapsed_seconds = end-start;
  std::ofstream out2;
  out2.open(timing_data_dir + "Nao" + std::to_string(C.size1()) + "Nt" + std::to_string(C.nt()) + "Ntau" + std::to_string(C.ntau()) + "tv_int_rtv.dat", std::ofstream::app);
  out2 << elapsed_seconds.count() << "\n" ;
  out2.close();
  // TIMING

  // Copy over the results into C.tv(tstp,m)
  memcpy(C.tvptr(tstp,0), NTauTmp.data(), (ntau_+1)*es_*sizeof(cplx));
}






// Does the les_lesadv integration for every n, puts into Q via increment!!!
// for n=j1...j2
//   computes the integral C_2^<[A,B](n,m) = \int_0^m dt A^<(n,t) B^A(t,m)
//                                         = \sum_{j=0}^{max(k,m)} w_{m,j} A^<(n,j) B^A(j,m)
// places into res which should be size size1*size1*(j2-j1+1)
void dyson::Cles2_tstp(int j1, int j2, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double dt, cplx *res) const {
  assert(j1 <= j2);
  assert(j1 >= 0);
  assert(j2 <= A.nt());
  assert(m <= A.nt());
  assert(A.nt() == Acc.nt());
  assert(A.nt() == B.nt());
  assert(A.nt() == Bcc.nt());
  assert(A.nt() == nt_);
  assert(A.nt() >= k_);
  assert(A.size1() == Acc.size1());
  assert(A.size1() == B.size1());
  assert(A.size1() == Bcc.size1());
  assert(A.size1() == nao_);
  assert(A.sig() == Acc.sig());
  assert(B.sig() == Bcc.sig());
  
  int j, l, n, top = m>k_ ? m : k_;

  // First fill X with BA, since they are the same for each n
  for(j=0; j<=top; j++) { // fill X from B.retptr.  BA_jm=(BccR_mj)*
    if(m >= j) { // We have BR
      ZMatrixMap(X.data() + j*es_, nao_, nao_).noalias() = dt*I.gregory_weights(m,j) * ZMatrixMap(Bcc.retptr(m,j), nao_, nao_).adjoint();
    }
    else {// Dont have BR
      ZMatrixMap(X.data() + j*es_, nao_, nao_).noalias() = -dt*I.gregory_weights(m,j) * ZMatrixMap(B.retptr(j,m), nao_, nao_);
    }
  }

  //sum from j=0 to n-1, where we don't have AL
  for(n=j1; n<=j2; n++) {
    ZMatrixMap resMap = ZMatrixMap(res + (n-j1)*es_, nao_, nao_);
    for(j=0; j<n; j++) {
      resMap.noalias() -= ZMatrixMap(Acc.lesptr(j,n), nao_, nao_).adjoint() * ZMatrixMap(X.data() + j*es_, nao_, nao_);
    }
  }

  //sum from j=n to max(k,m), here we have AL
  for(n=j1; n<=j2; n++) {
    ZMatrixMap resMap = ZMatrixMap(res + (n-j1)*es_, nao_, nao_);
    for(j=n; j<=top; j++) {
      resMap.noalias() += ZMatrixMap(A.lesptr(n,j), nao_, nao_) * ZMatrixMap(X.data() + j*es_, nao_, nao_);
    }
  }
}


void dyson::Cles2_tstp(const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double dt, cplx *res) const {
  int num = m>=k_ ? m : k_;
  return Cles2_tstp(0, num, A, Acc, B, Bcc, m, dt, res);
}


// Does the les_tvvt integral for each n=j1...j2 and puts into Q via increment!!!!
// computes C_3^<[A,B](n,m) = -i \int_0^\beta d\tau A^{TV}(n,\tau) B^{VT}(\tau,m)
//                          = -i \dau \sum_{j=0}^{ntau} w_{ntau,j} A^{TV}_{n,j} B^{VT}_{j,m}
// places into Q via increment.  Q should be size size1*size1*(j2-j1+1)
void dyson::Cles3_tstp(int j1, int j2, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double beta, cplx *res) const {
  assert(j1 <= j2);
  assert(j1 >= 0);
  assert(j2 <= A.nt());
  assert(m <= A.nt());
  assert(A.nt() == Acc.nt());
  assert(A.nt() == B.nt());
  assert(A.nt() == Bcc.nt());
  assert(A.nt() == nt_);
  assert(A.nt() >= k_);
  assert(A.size1() == Acc.size1());
  assert(A.size1() == B.size1());
  assert(A.size1() == Bcc.size1());
  assert(A.size1() == nao_);
  assert(A.sig() == Acc.sig());
  assert(B.sig() == Bcc.sig());
  assert(A.ntau() >= k_);
  assert(A.ntau() == B.ntau());
  assert(A.ntau() == Bcc.ntau());
  assert(A.ntau() == Acc.ntau());

  int sig=B.sig(), j, l, n;
  int top = m>=k_ ? m : k_;
  cplx cplxi = cplx(0.,1.);

  // Fill NTauTmp first
  for(j=0; j<=ntau_; j++) ZMatrixMap(NTauTmp.data() + j*es_, nao_, nao_).noalias() = (double)sig * cplxi * (ZMatrixMap(Bcc.tvptr(m,ntau_-j), nao_, nao_).adjoint());

  // Do the integral
  for(n=j1; n<=j2; n++) {
    auto ZTVRES = ZTensorView<2>(res + (n-j1)*es_, nao_, nao_);
    Conv.lesser(ZTVRES,
                ZTensorView<3>(A.tvptr(n,0), ntau_+1, nao_, nao_),
                ZTensorView<3>(NTauTmp.data(), ntau_+1, nao_, nao_), beta);
  }

/*
 * Old implementation with equidistant grid
 * Now we use Xinyang Legendre convolution
  //do the integral
  if(ntau_ < 2*(k_+1)-1) {
    for(n=j1; n<=j2; n++) {
      ZMatrixMap resMap = ZMatrixMap(res + (n-j1)*es_, nao_, nao_);
      for(m=0; m<=ntau_; m++) {
        resMap.noalias() += I.gregory_weights(ntau_, m) * ZMatrixMap(A.tvptr(n, m), nao_, nao_) * ZMatrixMap(NTauTmp.data() + m*es_, nao_, nao_);
      }
    }
  }
  else {
    for(n=j1; n<=j2; n++) {
      ZMatrixMap resMap = ZMatrixMap(res + (n-j1)*es_, nao_, nao_);
      for(m=0; m<=k_; m++) {
        resMap.noalias() += I.omega(m) * ZMatrixMap(A.tvptr(n,m), nao_, nao_) * ZMatrixMap(NTauTmp.data() + m*es_, nao_, nao_);
      }
      for(m=k_+1; m<ntau_-k_; m++) {
        resMap.noalias() += ZMatrixMap(A.tvptr(n,m), nao_, nao_) * ZMatrixMap(NTauTmp.data() + m*es_, nao_, nao_);
      }
      for(m=ntau_-k_; m<=ntau_; m++) {
        resMap.noalias() += I.omega(ntau_-m) * ZMatrixMap(A.tvptr(n,m), nao_, nao_) * ZMatrixMap(NTauTmp.data() + m*es_, nao_, nao_);
      }
    }
  }
*/
}

void dyson::Cles3_tstp(const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double beta, cplx *res) const {
  int top = m <= k_ ? k_ : m;
  return Cles3_tstp(0, top, A, Acc, B, Bcc, m, beta, res);
}

}//namespace
#endif
