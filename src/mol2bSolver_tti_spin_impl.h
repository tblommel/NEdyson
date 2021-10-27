//////////////////////////////////////////////////////////////////////////
//                      Spin Functions                                  //
//////////////////////////////////////////////////////////////////////////


void tti_molGF2SolverSpin::solve_les(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  ZMatrixMap(Sigma[0].get().lesptr(-tstp), nao_, nao_) = -ZMatrixMap(Sigma[0].get().tvptr(tstp,0), nao_, nao_).adjoint();
  ZMatrixMap(Sigma[1].get().lesptr(-tstp), nao_, nao_) = -ZMatrixMap(Sigma[1].get().tvptr(tstp,0), nao_, nao_).adjoint();
}


void tti_molGF2SolverSpin::solve_tv(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao_ * nao2;
  int sig = G[0].get().sig();
  int ntau = G[0].get().ntau();

  auto A1_Aa = ZMatrixMap(A1_aaa.data(), nao2, nao_);
  auto A1_Q = ZColVectorMap(A1_aaa.data(), nao3);

  auto C1_aA = ZMatrixMap(C1_aaa.data(), nao_, nao2);
  auto C1_Aa = ZMatrixMap(C1_aaa.data(), nao2, nao_);
  auto B1_aA = ZMatrixMap(B1_aaa.data(), nao_, nao2);
  auto B1_Aa = ZMatrixMap(B1_aaa.data(), nao2, nao_);

  auto Uexch_i_jkl = DMatrixConstMap(Uijkl_exch_.data(), nao_, nao3);
  auto U_i_jkl = DMatrixConstMap(Uijkl_.data(), nao_, nao3);

  for(int t=0; t<=ntau; ++t){
    for(int s=0; s<ns_; ++s){
      auto gTV_s_lk = ZMatrixMap(G[s].get().tvptr(tstp,t), nao_, nao_);

      for(int i=0; i<nao_; ++i){
        cplx *sigptr = Sigma[s].get().tvptr(tstp,t) + i*nao_;
        auto Ui_l_np = DMatrixConstMap(Uijkl_.data() + i*nao3, nao_, nao2);
        C1_aA.noalias() = gTV_s_lk.transpose() * Ui_l_np;

        for(int sp=0; sp<ns_; ++sp){
          B1_aA.noalias() = (C1_Aa * ZMatrixMap(G[sp].get().tvptr(tstp,t),nao_,nao_)).transpose();
          A1_Aa.noalias() = -sig * B1_Aa *
                              ZMatrixMap(G[sp].get().tvptr(tstp,ntau-t),nao_,nao_).conjugate();
          ZColVectorMap(sigptr, nao_).noalias() += Uexch_i_jkl * A1_Q;

          if(s==sp) {
            ZColVectorMap(sigptr, nao_).noalias() -= U_i_jkl * A1_Q;
          }
        }
      }
    }
  }
}


void tti_molGF2SolverSpin::solve_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  ZMatrixMap(Sigma[0].get().retptr(tstp), nao_, nao_) = G[0].get().sig() * ZMatrixMap(Sigma[0].get().tvptr(tstp, G[0].get().ntau()), nao_, nao_) - ZMatrixMap(Sigma[0].get().tvptr(tstp, 0), nao_, nao_);
  ZMatrixMap(Sigma[1].get().retptr(tstp), nao_, nao_) = G[0].get().sig() * ZMatrixMap(Sigma[1].get().tvptr(tstp, G[0].get().ntau()), nao_, nao_) - ZMatrixMap(Sigma[1].get().tvptr(tstp, 0), nao_, nao_);
}

/*
void tti_molGF2SolverSpin::solve_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao_ * nao2;
  int sig = G[0].get().sig();
  int ntau = G[0].get().ntau();

  auto A1_Aa = ZMatrixMap(A1_aaa.data(), nao2, nao_);
  auto A1_Q = ZColVectorMap(A1_aaa.data(), nao3);

  auto C1_aA = ZMatrixMap(C1_aaa.data(), nao_, nao2);
  auto C1_Aa = ZMatrixMap(C1_aaa.data(), nao2, nao_);
  auto C2_aA = ZMatrixMap(C2_aaa.data(), nao_, nao2);
  auto C2_Aa = ZMatrixMap(C2_aaa.data(), nao2, nao_);
  
  auto B1_aA = ZMatrixMap(B1_aaa.data(), nao_, nao2);
  auto B1_Aa = ZMatrixMap(B1_aaa.data(), nao2, nao_);
  auto B2_aA = ZMatrixMap(B2_aaa.data(), nao_, nao2);
  auto B2_Aa = ZMatrixMap(B2_aaa.data(), nao2, nao_);

  auto Uexch_i_jkl = DMatrixConstMap(Uijkl_exch_.data(), nao_, nao3);
  auto U_i_jkl = DMatrixConstMap(Uijkl_.data(), nao_, nao3);

  for(int s=0; s<ns_; ++s){
    auto gR_s_lk = ZMatrixMap(G[s].get().retptr(tstp), nao_, nao_);
    auto gL_s_lk = ZMatrixMap(G[s].get().lesptr(-tstp), nao_, nao_);

    for(int i=0; i<nao_; ++i){
      cplx *sigptr = Sigma[s].get().retptr(tstp) + i*nao_;
      auto Ui_l_np = DMatrixConstMap(Uijkl_.data() + i*nao3, nao_, nao2);
      C1_aA.noalias() = gR_s_lk.transpose() * Ui_l_np;
      C2_aA.noalias() = -gL_s_lk.conjugate() * Ui_l_np;

      for(int sp=0; sp<ns_; ++sp){
        auto gL_sp_lk = ZMatrixMap(G[sp].get().lesptr(-tstp), nao_, nao_);
        auto gR_sp_lk = ZMatrixMap(G[sp].get().retptr(tstp), nao_, nao_);
        
        B2_aA.noalias() = - (C2_Aa * gL_sp_lk.adjoint()).transpose();
        B1_aA.noalias() = (C1_Aa * gR_sp_lk - C1_Aa * gL_sp_lk.adjoint() + C2_Aa * gR_sp_lk).transpose();

        A1_Aa.noalias() = B1_Aa * ZMatrixMap(G[sp].get().lesptr(-tstp),nao_,nao_).transpose()
                        + B2_Aa * ZMatrixMap(G[sp].get().retptr(tstp),nao_,nao_).conjugate();
        ZColVectorMap(sigptr, nao_).noalias() += Uexch_i_jkl * A1_Q;

        if(s==sp) {
          ZColVectorMap(sigptr, nao_).noalias() -= U_i_jkl * A1_Q;
        }
      }
    }
  }
}
*/

void tti_molGF2SolverSpin::solve(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  size_t ref_size = size_t(2);
  assert(G.size() == ref_size);
  assert(Sigma.size() == ref_size);

  assert(G[0].get().sig() == G[1].get().sig());
  assert(Sigma[0].get().sig() == Sigma[1].get().sig());
  assert(G[0].get().sig() == Sigma[1].get().sig());

  assert(G[0].get().ntau() == G[1].get().ntau());
  assert(Sigma[0].get().ntau() == Sigma[1].get().ntau());
  assert(G[0].get().ntau() == Sigma[1].get().ntau());

  assert(tstp <= G[0].get().nt());
  assert(tstp <= G[1].get().nt());
  assert(tstp <= Sigma[0].get().nt());
  assert(tstp <= Sigma[1].get().nt());

  assert(G[0].get().size1() == nao_);
  assert(G[1].get().size1() == nao_);
  assert(Sigma[0].get().size1() == nao_);
  assert(Sigma[1].get().size1() == nao_);

  // Set self-energy zero
  Sigma[0].get().set_tstp_zero(tstp);
  Sigma[1].get().set_tstp_zero(tstp);

  // tv needs to be first since ret,les are just BCs
  solve_tv(tstp, Sigma, G);
  solve_ret(tstp, Sigma, G);
  solve_les(tstp, Sigma, G);
}


void tti_molGF2SolverSpin::solve_loop(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  size_t ref_size = size_t(2);
  assert(G.size() == ref_size);
  assert(Sigma.size() == ref_size);

  assert(G[0].get().sig() == G[1].get().sig());
  assert(Sigma[0].get().sig() == Sigma[1].get().sig());
  assert(G[0].get().sig() == Sigma[1].get().sig());

  assert(G[0].get().ntau() == G[1].get().ntau());
  assert(Sigma[0].get().ntau() == Sigma[1].get().ntau());
  assert(G[0].get().ntau() == Sigma[1].get().ntau());

  assert(tstp <= G[0].get().nt());
  assert(tstp <= G[1].get().nt());
  assert(tstp <= Sigma[0].get().nt());
  assert(tstp <= Sigma[1].get().nt());

  assert(G[0].get().size1() == nao_);
  assert(G[1].get().size1() == nao_);
  assert(Sigma[0].get().size1() == nao_);
  assert(Sigma[1].get().size1() == nao_);

  Sigma[0].get().set_tstp_zero(tstp);
  Sigma[1].get().set_tstp_zero(tstp);

  int sig = G[0].get().sig();
  int ntau= G[0].get().ntau();

  // Lesser
  for(int s=0; s<ns_; ++s){
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        for(int sp=0; sp<ns_; ++sp){
          for(int l=0; l<nao_; ++l){
            for(int n=0; n<nao_; ++n){
              for(int p=0; p<nao_; ++p){
                for(int k=0; k<nao_; ++k){
                  for(int q=0; q<nao_; ++q){
                    for(int m=0; m<nao_; ++m){
                      Sigma[s].get().lesptr(-tstp)[i*nao_+j] += Uijkl_(i,l,n,p) * Uijkl_(j,k,q,m) * G[s].get().lesptr(-tstp)[l*nao_+k] * (G[sp].get().retptr(tstp)[m*nao_+n]-std::conj(G[sp].get().lesptr(-tstp)[n*nao_+m])) * G[sp].get().lesptr(-tstp)[p*nao_+q];
                      if(s==sp){
                        Sigma[s].get().lesptr(-tstp)[i*nao_+j] -= Uijkl_(i,l,n,p) * Uijkl_(j,q,k,m) * G[s].get().lesptr(-tstp)[l*nao_+k] * (G[sp].get().retptr(tstp)[m*nao_+n]-std::conj(G[sp].get().lesptr(-tstp)[n*nao_+m])) * G[sp].get().lesptr(-tstp)[p*nao_+q];
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // TV
  for(int t=0; t<=ntau; ++t){
    for(int s=0; s<ns_; ++s){
      for(int i=0; i<nao_; ++i){
        for(int j=0; j<nao_; ++j){
          for(int sp=0; sp<ns_; ++sp){
            for(int l=0; l<nao_; ++l){
              for(int n=0; n<nao_; ++n){
                for(int p=0; p<nao_; ++p){
                  for(int k=0; k<nao_; ++k){
                    for(int q=0; q<nao_; ++q){
                      for(int m=0; m<nao_; ++m){
                        Sigma[s].get().tvptr(tstp,t)[i*nao_+j] += -sig*Uijkl_(i,l,n,p) * Uijkl_(j,k,q,m) * G[s].get().tvptr(tstp,t)[l*nao_+k] * std::conj(G[sp].get().tvptr(tstp,ntau-t)[n*nao_+m]) * G[sp].get().tvptr(tstp,t)[p*nao_+q];
                        if(s==sp){
                          Sigma[s].get().tvptr(tstp,t)[i*nao_+j] -= -sig*Uijkl_(i,l,n,p) * Uijkl_(j,q,k,m) * G[s].get().tvptr(tstp,t)[l*nao_+k] * std::conj(G[sp].get().tvptr(tstp,ntau-t)[n*nao_+m]) * G[sp].get().tvptr(tstp,t)[p*nao_+q];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // Retarded
  for(int s=0; s<ns_; ++s){
  cplx *sigRptr = Sigma[s].get().retptr(tstp);
  cplx *GsLptr = G[s].get().lesptr(-tstp);
  cplx *GsRptr = G[s].get().retptr(tstp);
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        cplx *sigRptrij = sigRptr+i*nao_+j;
        for(int sp=0; sp<ns_; ++sp){
          cplx *GspLptr = G[sp].get().lesptr(-tstp);
          cplx *GspRptr = G[sp].get().retptr(tstp);
          for(int l=0; l<nao_; ++l){
            for(int n=0; n<nao_; ++n){
              for(int p=0; p<nao_; ++p){
                for(int k=0; k<nao_; ++k){
                  cplx GsRlk = GsRptr[l*nao_+k];
                  cplx GsLkl = std::conj(GsLptr[k*nao_+l]);
                  cplx GspRlk = GspRptr[l*nao_+k];
                  cplx GspLkl = std::conj(GspLptr[k*nao_+l]);
                  for(int q=0; q<nao_; ++q){
                    cplx GspRpq = GspRptr[p*nao_+q];
                    cplx GspLqp = std::conj(GspLptr[q*nao_+p]);
                    for(int m=0; m<nao_; ++m){
                      cplx GspLmn = GspLptr[m*nao_+n];
                      cplx GspRnm = std::conj(GspRptr[n*nao_+m]);

                      *sigRptrij += Uijkl_(i,l,n,p) * Uijkl_(j,k,q,m) *
                            (GsRlk * GspLmn * (GspRpq - GspLqp)
                          - (GsLkl * (GspLmn * GspRpq - GspRnm * GspLqp)));
                      if(s==sp){
                        *sigRptrij -= Uijkl_(i,l,n,p) * Uijkl_(j,q,k,m) *
                            (GsRlk * GspLmn * (GspRpq - GspLqp)
                          - (GsLkl * (GspLmn * GspRpq - GspRnm * GspLqp)));
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}



