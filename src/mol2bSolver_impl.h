//////////////////////////////////////////////////////////////////////////
//                           Full   Functions                           //
//////////////////////////////////////////////////////////////////////////

void molGF2Solver::solve_les(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;

  for(int t=0; t<=tstp; ++t) {
    ZMatrixMap(A1_aaa.data(), nao_, nao_) = (ZMatrixMap(G.retptr(tstp,t), nao_, nao_) - ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).adjoint());
    ZMatrixMap(Sigma.lesptr(t,tstp), nao_, nao_) = U_*U_*ZMatrixMap(G.lesptr(t,tstp), nao_, nao_)
                                     .cwiseProduct(ZMatrixMap(G.lesptr(t,tstp), nao_, nao_))
                                     .cwiseProduct(ZMatrixMap(A1_aaa.data(),    nao_, nao_).transpose());
  }
}


void molGF2Solver::solve_tv(int tstp, GREEN &Sigma, GREEN &G) const {
//  int nao2 = nao_ * nao_;
//  int nao3 = nao2 * nao_;
//  int sig = G.sig();
//  int ntau = G.ntau();
//
//  auto A1_aA = ZMatrixMap(A1_aaa.data(), nao_, nao2);
//  auto A1_Aa = ZMatrixMap(A1_aaa.data(), nao2, nao_);
//  auto A1_Q =  ZColVectorMap(A1_aaa.data(), nao3);
//
//  auto B1_Aa = ZMatrixMap(B1_aaa.data(), nao2, nao_);
//  auto B1_aA = ZMatrixMap(B1_aaa.data(), nao_, nao2);
//
//  auto Uexch_j_qkm = DMatrixConstMap(Uijkl_exch_.data(), nao_, nao2*nao_);
//
//  for(int t=0; t<=ntau; ++t){
//    auto gtv = ZMatrixMap(G.tvptr(tstp,t), nao_, nao_);
//    auto gvt = ZMatrixMap(G.tvptr(tstp,ntau-t), nao_, nao_);
//
//    for(int i=0; i<nao_; ++i){
//      // A^\rceil(t,t')knp = (G^\rceilT(t,t'))_kl Ui_lnp
//      A1_aA = gtv.transpose() * DMatrixConstMap(Uijkl_.data() + i * nao3, nao_, nao2);
//
//      // B^\rceil(t,t')_qkn = [A^\rceil(t,t')_knp * G^\rceil(t,t')_pq]^T
//      B1_aA = (A1_Aa * gtv).transpose();
//
//      // A^\rceil(t,t')_qkm = B^\rceil(t,t')_qkn * (G^lceil(t',t)^T)_nm
//      //               =                         * -sig (G^rceil(t,ntau-t')^* 
//      A1_Aa = -sig * B1_Aa * gvt.conjugate();
//
//      // Sigma^<(t,t')i_j = Uexch_jqkm A^<(t,t')_qkm
//      ZRowVectorMap(Sigma.tvptr(tstp,t) + i * nao_, nao_) += Uexch_j_qkm * A1_Q;
//    }
//  }
}


void molGF2Solver::solve_ret(int tstp, GREEN &Sigma, GREEN &G) const {
/*  ZMatrixMap GGTt = ZMatrixMap(A1_aaa.data(), nao_, nao_);
  ZMatrixMap GGtT = ZMatrixMap(B1_aaa.data(), nao_, nao_);
  ZMatrixMap GLTt = ZMatrixMap(B2_aaa.data(), nao_, nao_);
  ZMatrixMap GLtT = ZMatrixMap(rho_T .data(), nao_, nao_);

  for(int t=0; t<=tstp; ++t) {
    GGTt = ZMatrixMap(G.retptr(tstp,t), nao_, nao_) - ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).adjoint();
    GGtT = -GGTt.adjoint();
    GLtT = ZMatrixMap(G.lesptr(t,tstp), nao_, nao_);
    GLTt = -GLtT.adjoint();

    ZMatrixMap(Sigma.retptr(tstp,t), nao_, nao_) = U_*U_* GGTt
                                            .cwiseProduct(GGTt)
                                            .cwiseProduct(GLtT.transpose());
    ZMatrixMap(Sigma.retptr(tstp,t), nao_, nao_) -= U_*U_* GLTt
                                             .cwiseProduct(GLTt)
                                             .cwiseProduct(GGtT.transpose());
  }
*/
}


void molGF2Solver::solve(int tstp, GREEN &Sigma, GREEN &G) const {
  assert(G.sig() == Sigma.sig());

  assert(G.ntau() == Sigma.ntau());

  assert(tstp <= G.nt());
  assert(tstp <= Sigma.nt());

  assert(G.size1() == nao_);
  assert(Sigma.size1() == nao_);

  // Set self energy to be zero
  Sigma.set_tstp_zero(tstp);
  
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;

  // Perform contractions
  start = std::chrono::system_clock::now();
  solve_les(tstp, Sigma, G);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;

  start = std::chrono::system_clock::now();
  solve_ret(tstp, Sigma, G);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
}
  

  


void molGF2Solver::solve_loop(int tstp, GREEN &Sigma, GREEN &G) const {
//  assert(G.sig() == G.sig());
//  assert(Sigma.sig() == Sigma.sig());
//  assert(G.sig() == Sigma.sig());
//
//  assert(G.ntau() == G.ntau());
//  assert(Sigma.ntau() == Sigma.ntau());
//  assert(G.ntau() == Sigma.ntau());
//
//  assert(tstp <= G.nt());
//  assert(tstp <= G.nt());
//  assert(tstp <= Sigma.nt());
//  assert(tstp <= Sigma.nt());
//
//  assert(G.size1() == nao_);
//  assert(G.size1() == nao_);
//  assert(Sigma.size1() == nao_);
//  assert(Sigma.size1() == nao_);
//
//  Sigma.set_tstp_zero(tstp);
//
//  int sig = G.sig();
//  int ntau= G.ntau();
//
//  // Lesser
//  for(int t=0; t<=tstp; ++t){
//    for(int i=0; i<nao_; ++i){
//      for(int j=0; j<nao_; ++j){          
//        for(int l=0; l<nao_; ++l){
//          for(int n=0; n<nao_; ++n){
//            for(int p=0; p<nao_; ++p){
//              for(int k=0; k<nao_; ++k){
//                for(int q=0; q<nao_; ++q){
//                  for(int m=0; m<nao_; ++m){
//                    Sigma.lesptr(t,tstp)[i*nao_+j] += Uijkl_(i,l,n,p) * (2*Uijkl_(j,k,q,m)-Uijkl_(j,q,k,m)) * G.lesptr(t,tstp)[l*nao_+k] * (G.retptr(tstp,t)[m*nao_+n]-std::conj(G.lesptr(t,tstp)[n*nao_+m])) * G.lesptr(t,tstp)[p*nao_+q];
//                  }
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//  }
//  // TV
//  for(int t=0; t<=ntau; ++t){
//    for(int i=0; i<nao_; ++i){
//      for(int j=0; j<nao_; ++j){
//        for(int l=0; l<nao_; ++l){
//          for(int n=0; n<nao_; ++n){
//            for(int p=0; p<nao_; ++p){
//              for(int k=0; k<nao_; ++k){
//                for(int q=0; q<nao_; ++q){
//                  for(int m=0; m<nao_; ++m){
//                    Sigma.tvptr(tstp,t)[i*nao_+j] += -sig * Uijkl_(i,l,n,p) * (2*Uijkl_(j,k,q,m)-Uijkl_(j,q,k,m)) * G.tvptr(tstp,t)[l*nao_+k] * std::conj(G.tvptr(tstp,ntau-t)[n*nao_+m]) * G.tvptr(tstp,t)[p*nao_+q];
//                  }
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//  }
//  // Retarded
//  for(int t=0; t<=tstp; ++t){
//    cplx *sigRptr = Sigma.retptr(tstp,t);
//    cplx *GLptr = G.lesptr(t,tstp);
//    cplx *GRptr = G.retptr(tstp,t);
//    for(int i=0; i<nao_; ++i){
//      for(int j=0; j<nao_; ++j){
//        cplx *sigRptrij = sigRptr + i*nao_ + j;
//        for(int l=0; l<nao_; ++l){
//          for(int n=0; n<nao_; ++n){
//            for(int p=0; p<nao_; ++p){
//              for(int k=0; k<nao_; ++k){
//                cplx GRlk = GRptr[l*nao_+k];
//                cplx GLkl = std::conj(GLptr[k*nao_+l]);
//                for(int q=0; q<nao_; ++q){
//                  cplx GRpq = GRptr[p*nao_+q];
//                  cplx GLqp = std::conj(GLptr[q*nao_+p]);
//                  for(int m=0; m<nao_; ++m){
//                    cplx GLmn = GLptr[m*nao_+n];
//                    cplx GRnm = std::conj(GRptr[n*nao_+m]);
//                    *sigRptrij += Uijkl_(i,l,n,p)
//                      * (2* Uijkl_(j,k,q,m) - Uijkl_(j,q,k,m))
//                      *( GRlk * GLmn * (GRpq-GLqp)
//                          - GLkl * (-GRnm * GLqp + GLmn * GRpq));
//                  }
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//  }
}



