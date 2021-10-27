//////////////////////////////////////////////////////////////////////////
//                      TTI Spin Decomp Functions                       //
//////////////////////////////////////////////////////////////////////////


void tti_molGF2SolverSpinDecomp::solve_bubble_les(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  ZMatrixMap(Sigma[0].get().lesptr(-tstp), nao_, nao_) = -ZMatrixMap(Sigma[0].get().tvptr(tstp,0), nao_, nao_).adjoint();
  ZMatrixMap(Sigma[1].get().lesptr(-tstp), nao_, nao_) = -ZMatrixMap(Sigma[1].get().tvptr(tstp,0), nao_, nao_).adjoint();
}

/*
void tti_molGF2SolverSpinDecomp::solve_bubble_les(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;

  // X^<_s'naq = V_nap * G^<_s'pq
  for(int sp = 0; sp<ns_; ++sp){
    ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, naa, nao_).noalias() = DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G[sp].get().lesptr(-tstp),nao_,nao_);
  }
  // Y^<_s'nqa = X^<_s'naq
  for(int sp=0; sp<ns_; ++sp){
    for(int n=0; n<nao_; n++){
      ZMatrixMap(Y1sija.data() + sp*nao2*nalpha_ + n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X1sija.data() + sp*nao2*nalpha_ + n*naa, nalpha_, nao_).transpose();
    }
  }

  // X^>_s'nqb(tstp,t) = (G^>_s')^T(tstp,t) V^T_m(qb)
  //                   = (G^R_s')^T(tstp,t) + (G^<_s')^T(tstp,t) ...
  //                   = (G^R_s')^T(tstp,t) - (G^<_s')^*(t,tstp) ...
  for(int sp = 0; sp<ns_; ++sp){
    ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_, naa).noalias() =
      (ZMatrixMap(G[sp].get().retptr(tstp),nao_,nao_).transpose()
     - ZMatrixMap(G[sp].get().lesptr(-tstp),nao_,nao_).conjugate())
     * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
  }


  // P^<_ab(t, tstp) = (Y^T)_a(s'nq) * X_(s'nq)b
  ZMatrixMap(P1ab.data(), nalpha_, nalpha_).noalias() =
      ZMatrixMap(Y1sija.data(), 2*nao2, nalpha_).transpose()
    * ZMatrixMap(X1sija.data(), 2*nao2, nalpha_);


  // Z^<_kja = V_kjb * (P^<_ab)^T
  ZMatrixMap(X1sija.data(), nao2, nalpha_).noalias() =
      DMatrixConstMap(Vija_.data(), nao2, nalpha_)
    * ZMatrixMap(P1ab.data(), nalpha_, nalpha_).transpose();


  // Z^<_kaj = Z^<_kja
  for(int k=0; k<nao_; ++k){
    ZMatrixMap(X2sija.data() + k*naa, nalpha_, nao_).noalias() =
        ZMatrixMap(X1sija.data() + k*naa, nao_, nalpha_).transpose();
  }

  // Sigma_sij = Y_sika * Z_kaj
  for(int s=0; s<ns_; ++s){
    ZMatrixMap(Sigma[s].get().lesptr(-tstp), nao_, nao_).noalias() +=
        ZMatrixMap(Y1sija.data() + s*nao2*nalpha_, nao_, naa)
      * ZMatrixMap(X2sija.data(), naa, nao_);
  }
}
*/

void tti_molGF2SolverSpinDecomp::solve_bubble_tv(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;
  int ntau = G[0].get().ntau();
  int sig = G[0].get().sig();

  for(int t=0; t<=ntau; ++t){
    // X^\rceil_s'naq = V_nap * G^\rceil_s'pq
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, naa, nao_).noalias() = DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G[sp].get().tvptr(tstp,t),nao_,nao_);
    }
    // Y^\rceil_s'nqa = X^\rceil_s'naq
    for(int sp=0; sp<ns_; ++sp){
      for(int n=0; n<nao_; ++n){
        ZMatrixMap(Y1sija.data() + sp*nao2*nalpha_ + n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X1sija.data() + sp*nao2*nalpha_ + n*naa, nalpha_, nao_).transpose();
      }
    }

    // X^\lceil_s'nqb(tstp,t) = (G^\lceil_s')^T(t,tstp) V^T_m(qb)
    //                    = -sig(G^\rceil_s')^*(tstp,ntau-t) ...
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_, naa).noalias() =
       -(double)sig * ZMatrixMap(G[sp].get().tvptr(tstp,ntau-t),nao_,nao_).conjugate()
       * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
    }

    // P^\rceil_ab(t, tstp) = (Y^T)_a(s'nq) * X_(s'nq)b
    ZMatrixMap(P1ab.data(), nalpha_, nalpha_).noalias() =
        ZMatrixMap(Y1sija.data(), 2*nao2, nalpha_).transpose()
      * ZMatrixMap(X1sija.data(), 2*nao2, nalpha_);

    // Z^\rceil_kja = V_kjb * (P^<_ab)^T
    ZMatrixMap(X1sija.data(), nao2, nalpha_).noalias() =
        DMatrixConstMap(Vija_.data(), nao2, nalpha_)
      * ZMatrixMap(P1ab.data(), nalpha_, nalpha_).transpose();

    // Z^\rceil_kaj = Z^\rceil_kja
    for(int k=0; k<nao_; ++k){
      ZMatrixMap(X2sija.data() + k*naa, nalpha_, nao_).noalias() =
          ZMatrixMap(X1sija.data() + k*naa, nao_, nalpha_).transpose();
    }

    // Sigma^\rceil_sij = Y^\rceil_sika * Z^\rceil_kaj
    for(int s=0; s<ns_; ++s){
      ZMatrixMap(Sigma[s].get().tvptr(tstp, t), nao_, nao_).noalias() +=
          ZMatrixMap(Y1sija.data() + s*nao2*nalpha_, nao_, naa)
        * ZMatrixMap(X2sija.data(), naa, nao_);
    }
  }
}

/*
void tti_molGF2SolverSpinDecomp::solve_bubble_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;

  //  X^<_s'naq(tstp,t) = V_nap * G^<_s'pq(tstp,t)
  //                    = V_nap * -(G^<_s'pq(t,tstp))^\dagger
  //
  //  X^R_s'naq(tstp,t) = V_nap * G^R_s'pq(tstp,t)
  for(int sp = 0; sp<ns_; ++sp){
    ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, naa, nao_).noalias() = -1 * DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G[sp].get().lesptr(-tstp),nao_,nao_).adjoint();
    ZMatrixMap(X2sija.data() + sp*nao2*nalpha_, naa, nao_).noalias() =  DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G[sp].get().retptr(tstp),nao_,nao_);
  }
  // Y^<_s'nqa = X^<_s'naq
  for(int sp=0; sp<ns_; ++sp){
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y1sija.data() + sp*nao2*nalpha_ + n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X1sija.data() + sp*nao2*nalpha_ + n*naa, nalpha_, nao_).transpose();
      ZMatrixMap(Y2sija.data() + sp*nao2*nalpha_ + n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X2sija.data() + sp*nao2*nalpha_ + n*naa, nalpha_, nao_).transpose();
    }
  }

  
  // X^A_s'nqb(t,tstp) = (G^A_s'mn)^T(t,tstp) * V^T_m(qb)
  //                   = (G^R_s'mn)^*(tstp,t) * V^T_m(qb)
  //
  // X^<_s'nqb(t,tstp) = (G^<_s'mn)^T(t,tstp) * V^T_m(qb)
  //                   
  // X^>_s'nqb(t,tstp) = (G^>_s'mn)^T(t,tstp) * V^T_m(qb)
  //                   = (G^R_s'mn(t,tstp) + G^<_s'mn(t,tstp))^T ...
  //                   = (-G^R_s'mn(tstp,t)^* + G^<_s'mn(t,tstp)^T) ...
  //
  for(int sp = 0; sp<ns_; ++sp){
    ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_, naa).noalias() = 
       ZMatrixMap(G[sp].get().retptr(tstp),nao_,nao_).conjugate() 
     * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

    ZMatrixMap(X2sija.data() + sp*nao2*nalpha_, nao_, naa).noalias() = 
       ZMatrixMap(G[sp].get().lesptr(-tstp),nao_,nao_).transpose() 
     * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

    ZMatrixMap(X3sija.data() + sp*nao2*nalpha_, nao_, naa).noalias() = 
      (ZMatrixMap(G[sp].get().lesptr(-tstp),nao_,nao_).transpose() 
     - ZMatrixMap(G[sp].get().retptr(tstp),nao_,nao_).conjugate())
     * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
  }


  // P^<_ab(tstp, t) = (Y^<T)_a(s'nq) * X^>_(s'nq)b
  // P^R_ab(tstp, t) = (Y^RT)_a(s'nq) * X^<_(s'nq)b + (Y^<T)_a(s'nq) * X^A_(s'nq)b
  ZMatrixMap(P1ab.data(), nalpha_, nalpha_).noalias() = 
      ZMatrixMap(Y1sija.data(), 2*nao2, nalpha_).transpose()
    * ZMatrixMap(X3sija.data(), 2*nao2, nalpha_);
  ZMatrixMap(P2ab.data(), nalpha_, nalpha_).noalias() = 
      ZMatrixMap(Y2sija.data(), 2*nao2, nalpha_).transpose()
    * ZMatrixMap(X2sija.data(), 2*nao2, nalpha_) 
    + ZMatrixMap(Y1sija.data(), 2*nao2, nalpha_).transpose()
    * ZMatrixMap(X1sija.data(), 2*nao2, nalpha_);


  // Z^<_kja = V_kjb * (P^<_ab)^T
  // Z^R_kja = V_kjb * (P^R_ab)^T
  ZMatrixMap(X1sija.data(), nao2, nalpha_).noalias() = 
      DMatrixConstMap(Vija_.data(), nao2, nalpha_)
    * ZMatrixMap(P1ab.data(), nalpha_, nalpha_).transpose();
  ZMatrixMap(X2sija.data(), nao2, nalpha_).noalias() = 
      DMatrixConstMap(Vija_.data(), nao2, nalpha_)
    * ZMatrixMap(P2ab.data(), nalpha_, nalpha_).transpose();

  // Z^<_kaj = Z^<_kja
  // Z^R_kaj = Z^R_kja
  for(int k=0; k<nao_; ++k){
    ZMatrixMap(X1sija.data() + nao2*nalpha_ + k*naa, nalpha_, nao_).noalias() =
        ZMatrixMap(X1sija.data() + k*naa, nao_, nalpha_).transpose();
    ZMatrixMap(X2sija.data() + nao2*nalpha_ + k*naa, nalpha_, nao_).noalias() =
        ZMatrixMap(X2sija.data() + k*naa, nao_, nalpha_).transpose();
  }

  // Sigma^R_sij = Y^R_sika * Z^R_kaj + Y^< * Z^R + Y^R * Z^<
  for(int s=0; s<ns_; ++s){
    ZMatrixMap(Sigma[s].get().retptr(tstp), nao_, nao_).noalias() += 
        ZMatrixMap(Y2sija.data() + s*nao2*nalpha_, nao_, naa)
      * ZMatrixMap(X2sija.data() + nao2*nalpha_, naa, nao_)
      + ZMatrixMap(Y1sija.data() + s*nao2*nalpha_, nao_, naa)
      * ZMatrixMap(X2sija.data() + nao2*nalpha_, naa, nao_)
      + ZMatrixMap(Y2sija.data() + s*nao2*nalpha_, nao_, naa)
      * ZMatrixMap(X1sija.data() + nao2*nalpha_, naa, nao_);
  }
}
*/

void tti_molGF2SolverSpinDecomp::solve_bubble_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  ZMatrixMap(Sigma[0].get().retptr(tstp), nao_, nao_) = G[0].get().sig() * ZMatrixMap(Sigma[0].get().tvptr(tstp, G[0].get().ntau()), nao_, nao_) - ZMatrixMap(Sigma[0].get().tvptr(tstp, 0), nao_, nao_);
  ZMatrixMap(Sigma[1].get().retptr(tstp), nao_, nao_) = G[0].get().sig() * ZMatrixMap(Sigma[1].get().tvptr(tstp, G[0].get().ntau()), nao_, nao_) - ZMatrixMap(Sigma[1].get().tvptr(tstp, 0), nao_, nao_);
}

void tti_molGF2SolverSpinDecomp::solve_bubble(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  // TV needs to be first since ret,les are BCs
  solve_bubble_tv(tstp, Sigma, G);
  solve_bubble_ret(tstp, Sigma, G);
  solve_bubble_les(tstp, Sigma, G);
}


void tti_molGF2SolverSpinDecomp::solve_exch_les(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  ZMatrixMap(Sigma[0].get().lesptr(-tstp), nao_, nao_) = -ZMatrixMap(Sigma[0].get().tvptr(tstp,0), nao_, nao_).adjoint();
  ZMatrixMap(Sigma[1].get().lesptr(-tstp), nao_, nao_) = -ZMatrixMap(Sigma[1].get().tvptr(tstp,0), nao_, nao_).adjoint();
}


void tti_molGF2SolverSpinDecomp::solve_exch_tv(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int naa = nao_ * nalpha_;
  int ntau = G[0].get().ntau();
  int sig = G[0].get().sig();

  for(int t=0; t<=ntau; ++t){
    for(int s=0; s<ns_; ++s){
      // X^\lceil_nkb = (G^\lceilT_s)_nm (V^T)_m(kb)
      //              = -sig (G^\rceil*_s_nm(ntau-t)  V^T_m(kb)
      ZMatrixMap(X1sija.data(), nao_, naa) = 
        -(double)sig * ZMatrixMap(G[s].get().tvptr(tstp,ntau-t), nao_, nao_).conjugate()
        * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

      // Transpose last two arguments of X
      for(int n=0; n<nao_; ++n){
        ZMatrixMap(X2sija.data() + n*naa, nalpha_, nao_) = 
            ZMatrixMap(X1sija.data() + n*naa, nao_, nalpha_).transpose();
      }

      // X^\rceil_nbl = X^lceil_nbk * G^rceilT_kl
      ZMatrixMap(X1sija.data(), naa, nao_) = 
          ZMatrixMap(X2sija.data(), naa, nao_)
        * ZMatrixMap(G[s].get().tvptr(tstp,t), nao_, nao_).transpose();
      
      // Z^\rceil_nqjl = V_qjb * X^\rceil_nbl
      for(int n=0; n<nao_; ++n){
        ZMatrixMap(Zijkl.data() + n*nao3, nao2, nao_) = 
            DMatrixConstMap(Vija_.data(), nao2, nalpha_)
          * ZMatrixMap(X1sija.data() + n*naa, nalpha_, nao_);
      }

      // Y^\rceil_naq = V_nap G^\rceil_pq
      ZMatrixMap(Y1sija.data(), naa, nao_) = 
          DMatrixConstMap(Viaj_.data(), naa, nao_)
        * ZMatrixMap(G[s].get().tvptr(tstp,t), nao_, nao_);

      // Transpose last two arguments of Y
      for(int n=0; n<nao_; ++n){
        ZMatrixMap(Y2sija.data() + n*naa, nao_, nalpha_) = 
            ZMatrixMap(Y1sija.data() + n*naa, nalpha_, nao_).transpose();
      }

      // X_jla = Z^T_(jl)(nq) * Y_nqa
      ZMatrixMap(X1sija.data(), nao2, nalpha_) =
          ZMatrixMap(Zijkl.data(), nao2, nao2).transpose()
        * ZMatrixMap(Y2sija.data(), nao2, nalpha_);
      
      // Sigma_ij += V_ila X^T(la)j
      ZMatrixMap(Sigma[s].get().tvptr(tstp,t), nao_, nao_) -=
          DMatrixConstMap(Vija_.data(), nao_, naa)
        * ZMatrixMap(X1sija.data(), nao_, naa).transpose();
    }
  }
}

void tti_molGF2SolverSpinDecomp::solve_exch_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  ZMatrixMap(Sigma[0].get().retptr(tstp), nao_, nao_) = G[0].get().sig() * ZMatrixMap(Sigma[0].get().tvptr(tstp, G[0].get().ntau()), nao_, nao_) - ZMatrixMap(Sigma[0].get().tvptr(tstp, 0), nao_, nao_);
  ZMatrixMap(Sigma[1].get().retptr(tstp), nao_, nao_) = G[0].get().sig() * ZMatrixMap(Sigma[1].get().tvptr(tstp, G[0].get().ntau()), nao_, nao_) - ZMatrixMap(Sigma[1].get().tvptr(tstp, 0), nao_, nao_);
}


/*
void tti_molGF2SolverSpinDecomp::solve_exch_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int naa = nao_ * nalpha_;

  for(int s=0; s<ns_; ++s){
    // X^>_nkb = (G^>T_s)_nm (V^T)_m(kb)         
    //         = (-G^R*_snm + G^LT_smn) V^T_m(kb) -> X10
    // X^<_nkb = (G^<T_s)_nm (V^T)_m(kb)          -> X20
    // X^A_nkb = (G^AT_s)_nm (V^T)_m(kb)
    //         = (G^T*_s)_nm (V^T)_m(kb)          -> X30
    ZMatrixMap(X1sija.data(), nao_, naa) = 
        (-ZMatrixMap(G[s].get().retptr(tstp), nao_, nao_).conjugate()
         +ZMatrixMap(G[s].get().lesptr(-tstp), nao_, nao_).transpose())
      *  DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
    ZMatrixMap(X2sija.data(), nao_, naa) = 
        ZMatrixMap(G[s].get().lesptr(-tstp), nao_, nao_).transpose()
      * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
    ZMatrixMap(X3sija.data(), nao_, naa) = 
        ZMatrixMap(G[s].get().retptr(tstp), nao_, nao_).conjugate()
      * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

    // Transpose last two arguments of X
    // > -> X11
    // < -> X21
    // A -> X31
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(X1sija.data() + naa*nao_ + n*naa, nalpha_, nao_) = 
          ZMatrixMap(X1sija.data() + n*naa, nao_, nalpha_).transpose();
      ZMatrixMap(X2sija.data() + naa*nao_ + n*naa, nalpha_, nao_) = 
          ZMatrixMap(X2sija.data() + n*naa, nao_, nalpha_).transpose();
      ZMatrixMap(X3sija.data() + naa*nao_ + n*naa, nalpha_, nao_) = 
          ZMatrixMap(X3sija.data() + n*naa, nao_, nalpha_).transpose();
    }

    // ~X~^<_nbl = X^>_nbk * G^<T_kl
    //           = - X^>_nbk * G^<*_kl                    -> X10
    // ~X~^R_nbl = X^<_nbk * G^RT_kl + X^A_nbk * G^<T_kl
    //           = X^<_nbk * G^RT_kl - X^A_nbk * G^<*_kl  -> X20
    ZMatrixMap(X1sija.data(), naa, nao_) = 
      - ZMatrixMap(X1sija.data() + naa*nao_, naa, nao_)
      * ZMatrixMap(G[s].get().lesptr(-tstp), nao_, nao_).conjugate();
    ZMatrixMap(X2sija.data(), naa, nao_) =
        ZMatrixMap(X2sija.data() + naa*nao_, naa, nao_)
      * ZMatrixMap(G[s].get().retptr(tstp), nao_, nao_).transpose()
      - ZMatrixMap(X3sija.data() + naa*nao_, naa, nao_)
      * ZMatrixMap(G[s].get().lesptr(-tstp), nao_, nao_).conjugate();


    // Y^<_naq = V_nap G^<_pq
    //         = - V_nap G^<\dagger_pq -> Y10
    // Y^R_naq = V_nap G^R_pq          -> Y20
    ZMatrixMap(Y1sija.data(), naa, nao_) =
      - DMatrixConstMap(Viaj_.data(), naa, nao_)
      * ZMatrixMap(G[s].get().lesptr(-tstp), nao_, nao_).adjoint();
    ZMatrixMap(Y2sija.data(), naa, nao_) = 
        DMatrixConstMap(Viaj_.data(), naa, nao_)
      * ZMatrixMap(G[s].get().retptr(tstp), nao_, nao_);

    // Transpose last two arguments of Y
    // < -> Y11
    // R -> Y21
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y1sija.data() + naa*nao_ + n*naa, nao_, nalpha_) =
          ZMatrixMap(Y1sija.data() + n*naa, nalpha_, nao_).transpose();
      ZMatrixMap(Y2sija.data() + naa*nao_ + n*naa, nao_, nalpha_) =
          ZMatrixMap(Y2sija.data() + n*naa, nalpha_, nao_).transpose();
    }
    
    // Z^<_nqjl = V_qjb * ~X~^<_nbl
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Zijkl.data() + n*nao3, nao2, nao_) =
          DMatrixConstMap(Vija_.data(), nao2, nalpha_)
        * ZMatrixMap(X1sija.data() + n*naa, nalpha_, nao_);
    }

    // ~Y~^R_jla = Z^T^<_(jl)(nq) * Y^R_nqa   -> X10
    ZMatrixMap(X1sija.data(), nao2, nalpha_) =
        ZMatrixMap(Zijkl.data(), nao2, nao2).transpose()
      * ZMatrixMap(Y2sija.data() + naa*nao_, nao2, nalpha_);

    // Z^R_nqjl = V_qjb * ~X~^R_nbl
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Zijkl.data() + n*nao3, nao2, nao_) =
          DMatrixConstMap(Vija_.data(), nao2, nalpha_)
        * ZMatrixMap(X2sija.data() + n*naa, nalpha_, nao_);
    }

    // ~Y~^R_jla += Z^T^R_(jl)(nq) * Y^R_nqa
    //            + Z^T^R_(jl)(nq) * Y^<_nqa  +-> X10
    ZMatrixMap(X1sija.data(), nao2, nalpha_) +=
        ZMatrixMap(Zijkl.data(), nao2, nao2).transpose()
      * ( ZMatrixMap(Y2sija.data() + naa*nao_, nao2, nalpha_)
         +ZMatrixMap(Y1sija.data() + naa*nao_, nao2, nalpha_));

    // Sigma_ij += V_ila ~Y~^TR(la)j
    ZMatrixMap(Sigma[s].get().retptr(tstp), nao_, nao_) -= 
        DMatrixConstMap(Vija_.data(), nao_, naa)
      * ZMatrixMap(X1sija.data(), nao_, naa).transpose();
  }
}
*/

void tti_molGF2SolverSpinDecomp::solve_exch(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  // TV needs to be first since ret,les are just BCs
  solve_exch_tv(tstp, Sigma, G);
  solve_exch_ret(tstp, Sigma, G);
  solve_exch_les(tstp, Sigma, G);
}

void tti_molGF2SolverSpinDecomp::solve(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
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
  
  // Do the contractions
  solve_bubble(tstp, Sigma, G);
  solve_exch(tstp, Sigma, G);
}

void tti_molGF2SolverSpinDecomp::solve_loop(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
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

  // Lesser bubble  
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

                      for(int a=0; a<nalpha_; ++a){
                        for(int b=0; b<nalpha_; ++b){
                          Sigma[s].get().lesptr(-tstp)[i*nao_+j] += Vija_(i,l,a) * Vija_(n,p,a) * Vija_(k,j,b) * Vija_(q,m,b) * G[s].get().lesptr(-tstp)[l*nao_+k] * (G[sp].get().retptr(tstp)[m*nao_+n]-std::conj(G[sp].get().lesptr(-tstp)[n*nao_+m])) * G[sp].get().lesptr(-tstp)[p*nao_+q];
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

  // Retarded Bubble
  for(int s=0; s<ns_; ++s){
    cplx *sigRptr = Sigma[s].get().retptr(tstp);
    cplx *GsLptr = G[s].get().lesptr(-tstp);
    cplx *GsRptr = G[s].get().retptr(tstp);
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        cplx *sigRptrij = sigRptr + i*nao_ + j;

        for(int sp=0; sp<ns_; ++sp){
        cplx *GspLptr = G[sp].get().lesptr(-tstp);
        cplx *GspRptr = G[sp].get().retptr(tstp);

          for(int l=0; l<nao_; ++l){
            for(int n=0; n<nao_; ++n){
              for(int p=0; p<nao_; ++p){
                for(int k=0; k<nao_; ++k){
                  cplx GsRlk = GsRptr[l*nao_+k];
                  cplx GsLkl = std::conj(GsLptr[k*nao_+l]);
                  for(int q=0; q<nao_; ++q){
                    cplx GspRpq = GspRptr[p*nao_+q];
                    cplx GspLqp = std::conj(GspLptr[q*nao_+p]);
                    for(int m=0; m<nao_; ++m){
                      cplx GspLmn = GspLptr[m*nao_+n];
                      cplx GspRnm = std::conj(GspRptr[n*nao_+m]);

                      for(int a=0; a<nalpha_; ++a){
                        for(int b=0; b<nalpha_; ++b){
                          *sigRptrij += Vija_(i,l,a) * Vija_(n,p,a)
                            * Vija_(k,j,b) * Vija_(q,m,b)
                            *( GsRlk * GspLmn * (GspRpq-GspLqp)
                                - GsLkl * (-GspRnm * GspLqp + GspLmn * GspRpq));
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

  // TV bubble
  for(int t=0; t<=ntau; ++t){
    for(int s=0; s<ns_; ++s){
      cplx *sigTVptr = Sigma[s].get().tvptr(tstp,t);
      cplx *GsTVptr = G[s].get().tvptr(tstp,t);
      for(int i=0; i<nao_; ++i){
        for(int j=0; j<nao_; ++j){
          cplx *sigTVptrij = sigTVptr+i*nao_+j;
          
          for(int sp=0; sp<ns_; ++sp){
            cplx *GspTVptr = G[sp].get().tvptr(tstp,t);
            cplx *GspVTptr = G[sp].get().tvptr(tstp,ntau-t);

            for(int l=0; l<nao_; ++l){
              for(int n=0; n<nao_; ++n){
                for(int p=0; p<nao_; ++p){
                  for(int k=0; k<nao_; ++k){
                    cplx GsTVlk = GsTVptr[l*nao_+k];
                    for(int q=0; q<nao_; ++q){
                      cplx GspTVpq = GspTVptr[p*nao_+q];
                      for(int m=0; m<nao_; ++m){
                        cplx GspVTnm = -1*(double)sig*std::conj(GspVTptr[n*nao_+m]);

                        for(int a=0; a<nalpha_; ++a){
                          for(int b=0; b<nalpha_; ++b){
                            *sigTVptrij += Vija_(i,l,a) * Vija_(n,p,a) * Vija_(k,j,b) * Vija_(q,m,b) * GsTVlk * GspVTnm * GspTVpq;
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
  }

  // Lesser exch
  for(int s=0; s<ns_; ++s){
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        
          for(int l=0; l<nao_; ++l){
            for(int n=0; n<nao_; ++n){
              for(int p=0; p<nao_; ++p){
                for(int k=0; k<nao_; ++k){
                  for(int q=0; q<nao_; ++q){
                    for(int m=0; m<nao_; ++m){

                      for(int a=0; a<nalpha_; ++a){
                        for(int b=0; b<nalpha_; ++b){
                          Sigma[s].get().lesptr(-tstp)[i*nao_+j] -= Vija_(i,l,a) * Vija_(n,p,a) * Vija_(q,j,b) * Vija_(k,m,b) * G[s].get().lesptr(-tstp)[l*nao_+k] * (G[s].get().retptr(tstp)[m*nao_+n]-std::conj(G[s].get().lesptr(-tstp)[n*nao_+m])) * G[s].get().lesptr(-tstp)[p*nao_+q];
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

  // TV exch
  for(int t=0; t<=ntau; ++t){
    for(int s=0; s<ns_; ++s){
      cplx *sigTVptr = Sigma[s].get().tvptr(tstp,t);
      cplx *GsTVptr = G[s].get().tvptr(tstp,t);
      cplx *GsVTptr = G[s].get().tvptr(tstp,ntau-t);
      for(int i=0; i<nao_; ++i){
        for(int j=0; j<nao_; ++j){
          cplx *sigTVptrij = sigTVptr+i*nao_+j;
          
            for(int l=0; l<nao_; ++l){
              for(int n=0; n<nao_; ++n){
                for(int p=0; p<nao_; ++p){
                  for(int k=0; k<nao_; ++k){
                    cplx GsTVlk = GsTVptr[l*nao_+k];
                    for(int q=0; q<nao_; ++q){
                      cplx GsTVpq = GsTVptr[p*nao_+q];
                      for(int m=0; m<nao_; ++m){
                        cplx GsVTnm = -1*(double)sig*std::conj(GsVTptr[n*nao_+m]);

                        for(int a=0; a<nalpha_; ++a){
                          for(int b=0; b<nalpha_; ++b){
                            *sigTVptrij -= Vija_(i,l,a) * Vija_(n,p,a) * Vija_(q,j,b) * Vija_(k,m,b) * GsTVlk * GsVTnm * GsTVpq;
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

  // Retarded exch
  for(int s=0; s<ns_; ++s){
    cplx *sigRptr = Sigma[s].get().retptr(tstp);
    cplx *GsLptr = G[s].get().lesptr(-tstp);
    cplx *GsRptr = G[s].get().retptr(tstp);
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        cplx *sigRptrij = sigRptr + i*nao_ + j;

          for(int l=0; l<nao_; ++l){
            for(int n=0; n<nao_; ++n){
              for(int p=0; p<nao_; ++p){
                for(int k=0; k<nao_; ++k){
                  cplx GsRlk = GsRptr[l*nao_+k];
                  cplx GsLkl = std::conj(GsLptr[k*nao_+l]);
                  for(int q=0; q<nao_; ++q){
                    cplx GsRpq = GsRptr[p*nao_+q];
                    cplx GsLqp = std::conj(GsLptr[q*nao_+p]);
                    for(int m=0; m<nao_; ++m){
                      cplx GsLmn = GsLptr[m*nao_+n];
                      cplx GsRnm = std::conj(GsRptr[n*nao_+m]);

                      for(int a=0; a<nalpha_; ++a){
                        for(int b=0; b<nalpha_; ++b){
                          *sigRptrij -= Vija_(i,l,a) * Vija_(n,p,a)
                            * Vija_(q,j,b) * Vija_(k,m,b)
                            *( GsRlk * GsLmn * (GsRpq-GsLqp)
                                - GsLkl * (-GsRnm * GsLqp + GsLmn * GsRpq));
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

