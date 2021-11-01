//////////////////////////////////////////////////////////////////////////
//                           Full   Functions                           //
//////////////////////////////////////////////////////////////////////////

void molGF2Solver::solve_les(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;

  auto A1_aA = ZMatrixMap(A1_aaa.data(), nao_, nao2);
  auto A1_Aa = ZMatrixMap(A1_aaa.data(), nao2, nao_);
  auto A1_Q =  ZColVectorMap(A1_aaa.data(), nao3);

  auto B1_Aa = ZMatrixMap(B1_aaa.data(), nao2, nao_);
  auto B1_aA = ZMatrixMap(B1_aaa.data(), nao_, nao2);

  auto Uexch_j_qkm = DMatrixConstMap(Uijkl_exch_.data(), nao_, nao2*nao_);

  for(int t=0; t<=tstp; ++t){
    auto gles = ZMatrixMap(G.lesptr(t,tstp), nao_, nao_);
    auto gret = ZMatrixMap(G.retptr(tstp,t), nao_, nao_);

    for(int i=0; i<nao_; ++i){
      // A^<(t,t')knp = (G^<T(t,t'))_kl Ui_lnp
      A1_aA = gles.transpose() * DMatrixConstMap(Uijkl_.data() + i * nao3, nao_, nao2);

      // B^<(t,t')_qkn = [A^<(t,t')_knp * G^<(t,t')_pq]^T
      B1_aA = (A1_Aa * gles).transpose();

      // A^<(t,t')_qkm = B^<(t,t')_qkn * (G^>(t',t)^T)_nm
      //               =               * (G^R(t',t)^T + G^<(t',t)^T)_nm
      //               =               * (G^R(t',t)^T)_nm - G^<(t,t')^*_nm
      A1_Aa = B1_Aa * (gret.transpose() - gles.conjugate());

      // Sigma^<(t,t')i_j = Uexch_jqkm A^<(t,t')_qkm
      ZRowVectorMap(Sigma.lesptr(t,tstp) + i * nao_, nao_) += Uexch_j_qkm * A1_Q;
    }
  }
}


void molGF2Solver::solve_tv(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int sig = G.sig();
  int ntau = G.ntau();

  auto A1_aA = ZMatrixMap(A1_aaa.data(), nao_, nao2);
  auto A1_Aa = ZMatrixMap(A1_aaa.data(), nao2, nao_);
  auto A1_Q =  ZColVectorMap(A1_aaa.data(), nao3);

  auto B1_Aa = ZMatrixMap(B1_aaa.data(), nao2, nao_);
  auto B1_aA = ZMatrixMap(B1_aaa.data(), nao_, nao2);

  auto Uexch_j_qkm = DMatrixConstMap(Uijkl_exch_.data(), nao_, nao2*nao_);

  for(int t=0; t<=ntau; ++t){
    auto gtv = ZMatrixMap(G.tvptr(tstp,t), nao_, nao_);
    auto gvt = ZMatrixMap(G.tvptr(tstp,ntau-t), nao_, nao_);

    for(int i=0; i<nao_; ++i){
      // A^\rceil(t,t')knp = (G^\rceilT(t,t'))_kl Ui_lnp
      A1_aA = gtv.transpose() * DMatrixConstMap(Uijkl_.data() + i * nao3, nao_, nao2);

      // B^\rceil(t,t')_qkn = [A^\rceil(t,t')_knp * G^\rceil(t,t')_pq]^T
      B1_aA = (A1_Aa * gtv).transpose();

      // A^\rceil(t,t')_qkm = B^\rceil(t,t')_qkn * (G^lceil(t',t)^T)_nm
      //               =                         * -sig (G^rceil(t,ntau-t')^* 
      A1_Aa = -sig * B1_Aa * gvt.conjugate();

      // Sigma^<(t,t')i_j = Uexch_jqkm A^<(t,t')_qkm
      ZRowVectorMap(Sigma.tvptr(tstp,t) + i * nao_, nao_) += Uexch_j_qkm * A1_Q;
    }
  }
}


void molGF2Solver::solve_ret(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;

  auto A1_aA = ZMatrixMap(A1_aaa.data(), nao_, nao2);
  auto A1_Aa = ZMatrixMap(A1_aaa.data(), nao2, nao_);
  auto A1_Q =  ZColVectorMap(A1_aaa.data(), nao3);

  auto B1_Aa = ZMatrixMap(B1_aaa.data(), nao2, nao_);
  auto B1_aA = ZMatrixMap(B1_aaa.data(), nao_, nao2);

  auto B2_aA = ZMatrixMap(B2_aaa.data(), nao_, nao2);
  auto B2_Aa = ZMatrixMap(B2_aaa.data(), nao2, nao_);

  auto Uexch_j_qkm = DMatrixConstMap(Uijkl_exch_.data(), nao_, nao2*nao_);

  for(int t=0; t<=tstp; ++t){
    auto gles = ZMatrixMap(G.lesptr(t,tstp), nao_, nao_);
    auto gret = ZMatrixMap(G.retptr(tstp,t), nao_, nao_);

    for(int i=0; i<nao_; ++i){
      // A^<(t,t')knp = (G^<T(t,t'))_kl Ui_lnp
      A1_aA = -gles.conjugate() * DMatrixConstMap(Uijkl_.data() + i * nao3, nao_, nao2);

      // B^<(t,t')_qkn = [A^<(t,t')_knp * G^<(t,t')_pq]^T
      B1_aA = -(A1_Aa * gles.adjoint()).transpose();
      // B^R(t,t')_qkn = [A^<(t,t')_knp * G^R(t,t')_pq]^T
      B2_aA = (A1_Aa * gret).transpose();
      
      // A^R(t,t')knp = (G^RT(t,t'))_kl Ui_lnp
      A1_aA = gret.transpose() * DMatrixConstMap(Uijkl_.data() + i * nao3, nao_, nao2);

      // B^R(t,t')_qkn += [A^R(t,t')_knp * G^R(t,t')_pq + A^R(t,t')_knp * G^<(t,t')_pq]^T
      B2_aA += (A1_Aa * (gret-gles.adjoint())).transpose();

      // A^R(t,t')_qkm = B^R(t,t')_qkn * (G^<(t',t)^T)_nm + B^< G^A
      A1_Aa = B2_Aa * gles.transpose() + B1_Aa * gret.conjugate();

      // Sigma^R(t,t')i_j = Uexch_jqkm A^<(t,t')_qkm
      ZRowVectorMap(Sigma.retptr(tstp,t) + i * nao_, nao_) += Uexch_j_qkm * A1_Q;
    }
  }
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
  //TIMING
  if(tstp > 5) {
    std::ofstream out1;
    std::string timing_data_dir = std::string(TIMING_DATA_DIR);
    out1.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "les_cont.dat", std::ofstream::app);
    out1 << elapsed_seconds.count() << "\n" ;
    out1.close();
  }
  // TIMING

  start = std::chrono::system_clock::now();
  solve_tv(tstp, Sigma, G);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  //TIMING
  if(tstp > 5) {
    std::ofstream out2;
    std::string timing_data_dir = std::string(TIMING_DATA_DIR);
    out2.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "tv_cont.dat", std::ofstream::app);
    out2 << elapsed_seconds.count() << "\n" ;
    out2.close();
  }
  // TIMING


  start = std::chrono::system_clock::now();
  solve_ret(tstp, Sigma, G);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  if(tstp > 5) {
    std::ofstream out3;
    std::string timing_data_dir = std::string(TIMING_DATA_DIR);
    out3.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "ret_cont.dat", std::ofstream::app);
    out3 << elapsed_seconds.count() << "\n" ;
    out3.close();
  }
  // TIMING
}
  

  


void molGF2Solver::solve_loop(int tstp, GREEN &Sigma, GREEN &G) const {
  assert(G.sig() == G.sig());
  assert(Sigma.sig() == Sigma.sig());
  assert(G.sig() == Sigma.sig());

  assert(G.ntau() == G.ntau());
  assert(Sigma.ntau() == Sigma.ntau());
  assert(G.ntau() == Sigma.ntau());

  assert(tstp <= G.nt());
  assert(tstp <= G.nt());
  assert(tstp <= Sigma.nt());
  assert(tstp <= Sigma.nt());

  assert(G.size1() == nao_);
  assert(G.size1() == nao_);
  assert(Sigma.size1() == nao_);
  assert(Sigma.size1() == nao_);

  Sigma.set_tstp_zero(tstp);

  int sig = G.sig();
  int ntau= G.ntau();

  // Lesser
  for(int t=0; t<=tstp; ++t){
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){          
        for(int l=0; l<nao_; ++l){
          for(int n=0; n<nao_; ++n){
            for(int p=0; p<nao_; ++p){
              for(int k=0; k<nao_; ++k){
                for(int q=0; q<nao_; ++q){
                  for(int m=0; m<nao_; ++m){
                    Sigma.lesptr(t,tstp)[i*nao_+j] += Uijkl_(i,l,n,p) * (2*Uijkl_(j,k,q,m)-Uijkl_(j,q,k,m)) * G.lesptr(t,tstp)[l*nao_+k] * (G.retptr(tstp,t)[m*nao_+n]-std::conj(G.lesptr(t,tstp)[n*nao_+m])) * G.lesptr(t,tstp)[p*nao_+q];
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
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        for(int l=0; l<nao_; ++l){
          for(int n=0; n<nao_; ++n){
            for(int p=0; p<nao_; ++p){
              for(int k=0; k<nao_; ++k){
                for(int q=0; q<nao_; ++q){
                  for(int m=0; m<nao_; ++m){
                    Sigma.tvptr(tstp,t)[i*nao_+j] += -sig * Uijkl_(i,l,n,p) * (2*Uijkl_(j,k,q,m)-Uijkl_(j,q,k,m)) * G.tvptr(tstp,t)[l*nao_+k] * std::conj(G.tvptr(tstp,ntau-t)[n*nao_+m]) * G.tvptr(tstp,t)[p*nao_+q];
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
  for(int t=0; t<=tstp; ++t){
    cplx *sigRptr = Sigma.retptr(tstp,t);
    cplx *GLptr = G.lesptr(t,tstp);
    cplx *GRptr = G.retptr(tstp,t);
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        cplx *sigRptrij = sigRptr + i*nao_ + j;
        for(int l=0; l<nao_; ++l){
          for(int n=0; n<nao_; ++n){
            for(int p=0; p<nao_; ++p){
              for(int k=0; k<nao_; ++k){
                cplx GRlk = GRptr[l*nao_+k];
                cplx GLkl = std::conj(GLptr[k*nao_+l]);
                for(int q=0; q<nao_; ++q){
                  cplx GRpq = GRptr[p*nao_+q];
                  cplx GLqp = std::conj(GLptr[q*nao_+p]);
                  for(int m=0; m<nao_; ++m){
                    cplx GLmn = GLptr[m*nao_+n];
                    cplx GRnm = std::conj(GRptr[n*nao_+m]);
                    *sigRptrij += Uijkl_(i,l,n,p)
                      * (2* Uijkl_(j,k,q,m) - Uijkl_(j,q,k,m))
                      *( GRlk * GLmn * (GRpq-GLqp)
                          - GLkl * (-GRnm * GLqp + GLmn * GRpq));
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



