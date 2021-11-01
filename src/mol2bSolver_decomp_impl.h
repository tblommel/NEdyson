//////////////////////////////////////////////////////////////////////////
//                           Decomp Functions                           //
//////////////////////////////////////////////////////////////////////////


void molGF2SolverDecomp::solve_bubble_les(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;
  for(int t=0; t<=tstp; ++t){
    // X^<_naq = V_nap * G^<_pq
    ZMatrixMap(X1ija.data(), naa, nao_).noalias() = DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G.lesptr(t,tstp),nao_,nao_);
    // Y^<_nqa = X^<_naq
    for(int n=0; n<nao_; n++){
      ZMatrixMap(Y1ija.data() +  n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X1ija.data() +  n*naa, nalpha_, nao_).transpose();
    }

    // X^>_nqb(tstp,t) = (G^>)^T(tstp,t) V^T_m(qb)
    //                   = (G^R)^T(tstp,t) + (G^<)^T(tstp,t) ...
    //                   = (G^R)^T(tstp,t) - (G^<)^*(t,tstp) ...
    ZMatrixMap(X1ija.data(), nao_, naa).noalias() =
      (ZMatrixMap(G.retptr(tstp,t),nao_,nao_).transpose()
     - ZMatrixMap(G.lesptr(t,tstp),nao_,nao_).conjugate())
     * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();


    // P^<_ab(t, tstp) = (Y^T)_a(nq) * X_(nq)b
    ZMatrixMap(P1ab.data(), nalpha_, nalpha_).noalias() =
        ZMatrixMap(Y1ija.data(), nao2, nalpha_).transpose()
      * ZMatrixMap(X1ija.data(), nao2, nalpha_);


    // Z^<_kja = V_kjb * (P^<_ab)^T
    ZMatrixMap(X1ija.data(), nao2, nalpha_).noalias() =
        DMatrixConstMap(Vija_.data(), nao2, nalpha_)
      * ZMatrixMap(P1ab.data(), nalpha_, nalpha_).transpose();


    // Z^<_kaj = Z^<_kja
    for(int k=0; k<nao_; ++k){
      ZMatrixMap(X2ija.data() + k*naa, nalpha_, nao_).noalias() =
          ZMatrixMap(X1ija.data() + k*naa, nao_, nalpha_).transpose();
    }

    // Sigma_ij = Y_ika * Z_kaj
    ZMatrixMap(Sigma.lesptr(t, tstp), nao_, nao_).noalias() +=
      2 * ZMatrixMap(Y1ija.data(), nao_, naa)
      * ZMatrixMap(X2ija.data(), naa, nao_);
  }
}

void molGF2SolverDecomp::solve_bubble_tv(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;
  int ntau = G.ntau();
  int sig = G.sig();
  for(int t=0; t<=ntau; ++t){
    // X^\rceil_naq = V_nap * G^\rceil_pq
    ZMatrixMap(X1ija.data(), naa, nao_).noalias() = DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G.tvptr(tstp,t),nao_,nao_);
    // Y^\rceil_nqa = X^\rceil_naq
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y1ija.data() + n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X1ija.data() + n*naa, nalpha_, nao_).transpose();
    }

    // X^\lceil_nqb(tstp,t) = (G^\lceil)^T(t,tstp) V^T_m(qb)
    //                    = -sig(G^\rceil)^*(tstp,ntau-t) ...
    ZMatrixMap(X1ija.data(), nao_, naa).noalias() =
      -(double)sig * ZMatrixMap(G.tvptr(tstp,ntau-t),nao_,nao_).conjugate()
      * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

    // P^\rceil_ab(t, tstp) = (Y^T)_a(nq) * X_(nq)b
    ZMatrixMap(P1ab.data(), nalpha_, nalpha_).noalias() =
        ZMatrixMap(Y1ija.data(), nao2, nalpha_).transpose()
      * ZMatrixMap(X1ija.data(), nao2, nalpha_);

    // Z^\rceil_kja = V_kjb * (P^<_ab)^T
    ZMatrixMap(X1ija.data(), nao2, nalpha_).noalias() =
        DMatrixConstMap(Vija_.data(), nao2, nalpha_)
      * ZMatrixMap(P1ab.data(), nalpha_, nalpha_).transpose();

    // Z^\rceil_kaj = Z^\rceil_kja
    for(int k=0; k<nao_; ++k){
      ZMatrixMap(X2ija.data() + k*naa, nalpha_, nao_).noalias() =
          ZMatrixMap(X1ija.data() + k*naa, nao_, nalpha_).transpose();
    }

    // Sigma^\rceil_ij = Y^\rceil_ika * Z^\rceil_kaj
    ZMatrixMap(Sigma.tvptr(tstp, t), nao_, nao_).noalias() +=
      2 * ZMatrixMap(Y1ija.data(), nao_, naa)
      * ZMatrixMap(X2ija.data(), naa, nao_);
  }
}


void molGF2SolverDecomp::solve_bubble_ret(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;
  for(int t=0; t<=tstp; ++t){
    //  X^<_naq(tstp,t) = V_nap * G^<_pq(tstp,t)
    //                  = V_nap * -(G^<_pq(t,tstp))^\dagger
    //
    //  X^R_naq(tstp,t) = V_nap * G^R_pq(tstp,t)
    ZMatrixMap(X1ija.data(), naa, nao_).noalias() = -1 * DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G.lesptr(t,tstp),nao_,nao_).adjoint();
    ZMatrixMap(X2ija.data(), naa, nao_).noalias() =  DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G.retptr(tstp,t),nao_,nao_);
    // Y^<_nqa = X^<_naq
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y1ija.data() + n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X1ija.data() + n*naa, nalpha_, nao_).transpose();
      ZMatrixMap(Y2ija.data() + n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X2ija.data() + n*naa, nalpha_, nao_).transpose();
    }

    
    // X^A_nqb(t,tstp) = (G^A_mn)^T(t,tstp) * V^T_m(qb)
    //                 = (G^R_mn)^*(tstp,t) * V^T_m(qb)
    //
    // X^<_nqb(t,tstp) = (G^<_mn)^T(t,tstp) * V^T_m(qb)
    //                   
    // X^>_nqb(t,tstp) = (G^>_mn)^T(t,tstp) * V^T_m(qb)
    //                   = (G^R_mn(t,tstp) + G^<_mn(t,tstp))^T ...
    //                   = (-G^R_mn(tstp,t)^* + G^<_mn(t,tstp)^T) ...
    ZMatrixMap(X1ija.data(), nao_, naa).noalias() =
       ZMatrixMap(G.retptr(tstp,t),nao_,nao_).conjugate()
     * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

    ZMatrixMap(X2ija.data(), nao_, naa).noalias() = 
       ZMatrixMap(G.lesptr(t,tstp),nao_,nao_).transpose() 
     * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

    ZMatrixMap(X3ija.data(), nao_, naa).noalias() =
      (ZMatrixMap(G.lesptr(t,tstp),nao_,nao_).transpose()
     - ZMatrixMap(G.retptr(tstp,t),nao_,nao_).conjugate())
     * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();


    // P^<_ab(tstp, t) = (Y^<T)_a(nq) * X^>_(nq)b
    // P^R_ab(tstp, t) = (Y^RT)_a(nq) * X^<_(nq)b + (Y^<T)_a(nq) * X^A_(nq)b
    ZMatrixMap(P1ab.data(), nalpha_, nalpha_).noalias() =
        ZMatrixMap(Y1ija.data(), nao2, nalpha_).transpose()
      * ZMatrixMap(X3ija.data(), nao2, nalpha_);
    ZMatrixMap(P2ab.data(), nalpha_, nalpha_).noalias() = 
        ZMatrixMap(Y2ija.data(), nao2, nalpha_).transpose()
      * ZMatrixMap(X2ija.data(), nao2, nalpha_)
      + ZMatrixMap(Y1ija.data(), nao2, nalpha_).transpose()
      * ZMatrixMap(X1ija.data(), nao2, nalpha_);


    // Z^<_kja = V_kjb * (P^<_ab)^T
    // Z^R_kja = V_kjb * (P^R_ab)^T
    ZMatrixMap(X1ija.data(), nao2, nalpha_).noalias() =
        DMatrixConstMap(Vija_.data(), nao2, nalpha_)
      * ZMatrixMap(P1ab.data(), nalpha_, nalpha_).transpose();
    ZMatrixMap(X2ija.data(), nao2, nalpha_).noalias() =
        DMatrixConstMap(Vija_.data(), nao2, nalpha_)
      * ZMatrixMap(P2ab.data(), nalpha_, nalpha_).transpose();

    // Z^<_kaj = Z^<_kja
    // Z^R_kaj = Z^R_kja
    for(int k=0; k<nao_; ++k){
      ZMatrixMap(X3ija.data() + k*naa, nalpha_, nao_).noalias() =
          ZMatrixMap(X2ija.data() + k*naa, nao_, nalpha_).transpose();
      ZMatrixMap(X2ija.data() + k*naa, nalpha_, nao_).noalias() =
          ZMatrixMap(X1ija.data() + k*naa, nao_, nalpha_).transpose();
    }

    // Sigma^R_ij = Y^R_ika * Z^R_kaj + Y^< * Z^R + Y^R * Z^<
    ZMatrixMap(Sigma.retptr(tstp, t), nao_, nao_).noalias() += 
    2*(ZMatrixMap(Y2ija.data(), nao_, naa)
      * ZMatrixMap(X3ija.data(), naa, nao_)
      + ZMatrixMap(Y1ija.data(), nao_, naa)
      * ZMatrixMap(X3ija.data(), naa, nao_)
      + ZMatrixMap(Y2ija.data(), nao_, naa)
      * ZMatrixMap(X2ija.data(), naa, nao_));
  }
}


void molGF2SolverDecomp::solve_bubble(int tstp, GREEN &Sigma, GREEN &G) const {
  solve_bubble_les(tstp, Sigma, G);
  solve_bubble_ret(tstp, Sigma, G);
  solve_bubble_tv(tstp, Sigma, G);
}


void molGF2SolverDecomp::solve_exch_les(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int naa = nao_ * nalpha_;
  for(int t=0; t<=tstp; ++t){
    // X^>_nkb = (G^>T)_nm (V^T)_m(kb)
    //         = (G^RTnm - G^L*mn) V^T_m(kb)
    ZMatrixMap(X1ija.data(), nao_, naa) = 
        (ZMatrixMap(G.retptr(tstp,t), nao_, nao_).transpose()
      -  ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).conjugate())
      *  DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

    // Transpose last two arguments of X
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(X2ija.data() + n*naa, nalpha_, nao_) =
          ZMatrixMap(X1ija.data() + n*naa, nao_, nalpha_).transpose();
    }

    // X^<_nbl = X^>_nbk * G^<T_kl
    ZMatrixMap(X1ija.data(), naa, nao_) = 
        ZMatrixMap(X2ija.data(), naa, nao_)
      * ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).transpose();
    
    // Z^<_nqjl = V_qjb * X^<_nbl
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Zijkl.data() + n*nao3, nao2, nao_) =
          DMatrixConstMap(Vija_.data(), nao2, nalpha_)
        * ZMatrixMap(X1ija.data() + n*naa, nalpha_, nao_);
    }

    // Y^<_naq = V_nap G^<_pq
    ZMatrixMap(Y1ija.data(), naa, nao_) =
        DMatrixConstMap(Viaj_.data(), naa, nao_)
      * ZMatrixMap(G.lesptr(t,tstp), nao_, nao_);

    // Transpose last two arguments of Y
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y2ija.data() + n*naa, nao_, nalpha_) = 
          ZMatrixMap(Y1ija.data() + n*naa, nalpha_, nao_).transpose();
    }

    // X_jla = Z^T_(jl)(nq) * Y_nqa
    ZMatrixMap(X1ija.data(), nao2, nalpha_) =
        ZMatrixMap(Zijkl.data(), nao2, nao2).transpose()
      * ZMatrixMap(Y2ija.data(), nao2, nalpha_);
    
    // Sigma_ij += V_ila X^T(la)j
    ZMatrixMap(Sigma.lesptr(t,tstp), nao_, nao_) -=
        DMatrixConstMap(Vija_.data(), nao_, naa)
      * ZMatrixMap(X1ija.data(), nao_, naa).transpose();
  }
}


void molGF2SolverDecomp::solve_exch_tv(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int naa = nao_ * nalpha_;
  int ntau = G.ntau();
  int sig = G.sig();

  for(int t=0; t<=ntau; ++t){
    // X^\lceil_nkb = (G^\lceilT)_nm (V^T)_m(kb)
    //              = -sig (G^\rceil*_nm(ntau-t)  V^T_m(kb)
    ZMatrixMap(X1ija.data(), nao_, naa) =
      -(double)sig * ZMatrixMap(G.tvptr(tstp,ntau-t), nao_, nao_).conjugate()
      * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

    // Transpose last two arguments of X
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(X2ija.data() + n*naa, nalpha_, nao_) = 
          ZMatrixMap(X1ija.data() + n*naa, nao_, nalpha_).transpose();
    }

    // X^\rceil_nbl = X^lceil_nbk * G^rceilT_kl
    ZMatrixMap(X1ija.data(), naa, nao_) =
        ZMatrixMap(X2ija.data(), naa, nao_)
      * ZMatrixMap(G.tvptr(tstp,t), nao_, nao_).transpose();
    
    // Z^\rceil_nqjl = V_qjb * X^\rceil_nbl
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Zijkl.data() + n*nao3, nao2, nao_) = 
          DMatrixConstMap(Vija_.data(), nao2, nalpha_)
        * ZMatrixMap(X1ija.data() + n*naa, nalpha_, nao_);
    }

    // Y^\rceil_naq = V_nap G^\rceil_pq
    ZMatrixMap(Y1ija.data(), naa, nao_) = 
        DMatrixConstMap(Viaj_.data(), naa, nao_)
      * ZMatrixMap(G.tvptr(tstp,t), nao_, nao_);

    // Transpose last two arguments of Y
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y2ija.data() + n*naa, nao_, nalpha_) =
          ZMatrixMap(Y1ija.data() + n*naa, nalpha_, nao_).transpose();
    }

    // X_jla = Z^T_(jl)(nq) * Y_nqa
    ZMatrixMap(X1ija.data(), nao2, nalpha_) =
        ZMatrixMap(Zijkl.data(), nao2, nao2).transpose()
      * ZMatrixMap(Y2ija.data(), nao2, nalpha_);
 
    // Sigma_ij += V_ila X^T(la)j
    ZMatrixMap(Sigma.tvptr(tstp,t), nao_, nao_) -= 
        DMatrixConstMap(Vija_.data(), nao_, naa)
      * ZMatrixMap(X1ija.data(), nao_, naa).transpose();
  }
}



void molGF2SolverDecomp::solve_exch_ret(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int naa = nao_ * nalpha_;

  for(int t=0; t<=tstp; ++t){
    // X^>_nkb = (G^>T_s)_nm (V^T)_m(kb)         
    //         = (-G^R*_snm + G^LT_smn) V^T_m(kb)
    // X^<_nkb = (G^<T_s)_nm (V^T)_m(kb)          
    // X^A_nkb = (G^AT_s)_nm (V^T)_m(kb)
    //         = (G^T*_s)_nm (V^T)_m(kb)
    //
    // Transpose last two arguments of X  
    // > -> Y1
    // < -> Y2
    // A -< X3       
    ZMatrixMap(X2ija.data(), nao_, naa) =
        (-ZMatrixMap(G.retptr(tstp,t), nao_, nao_).conjugate()
         +ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).transpose())
      *  DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y1ija.data() + n*naa, nalpha_, nao_) =
          ZMatrixMap(X2ija.data() + n*naa, nao_, nalpha_).transpose();
    }

    ZMatrixMap(X2ija.data(), nao_, naa) = 
        ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).transpose()
      * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y2ija.data() + n*naa, nalpha_, nao_) = 
          ZMatrixMap(X2ija.data() + n*naa, nao_, nalpha_).transpose();
    }

    ZMatrixMap(X2ija.data(), nao_, naa) =
        ZMatrixMap(G.retptr(tstp,t), nao_, nao_).conjugate()
      * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(X3ija.data() + n*naa, nalpha_, nao_) =
          ZMatrixMap(X2ija.data() + n*naa, nao_, nalpha_).transpose();
    }

    // ~X~^<_nbl = X^>_nbk * G^<T_kl
    //           = - X^>_nbk * G^<*_kl                    -> X1
    // ~X~^R_nbl = X^<_nbk * G^RT_kl + X^A_nbk * G^<T_kl
    //           = X^<_nbk * G^RT_kl - X^A_nbk * G^<*_kl  -> X2
    ZMatrixMap(X1ija.data(), naa, nao_) =
      - ZMatrixMap(Y1ija.data(), naa, nao_)
      * ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).conjugate();
    ZMatrixMap(X2ija.data(), naa, nao_) =
        ZMatrixMap(Y2ija.data(), naa, nao_)
      * ZMatrixMap(G.retptr(tstp,t), nao_, nao_).transpose()
      - ZMatrixMap(X3ija.data(), naa, nao_)
      * ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).conjugate();


    // Y^<_naq = V_nap G^<_pq
    //         = - V_nap G^<\dagger_pq 
    // Y^R_naq = V_nap G^R_pq          
    //
    // Transpose last two arguments
    // < -> Y1
    // R -> Y2
    ZMatrixMap(X3ija.data(), naa, nao_) = 
      - DMatrixConstMap(Viaj_.data(), naa, nao_)
      * ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).adjoint();
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y1ija.data() + n*naa, nao_, nalpha_) =
          ZMatrixMap(X3ija.data() + n*naa, nalpha_, nao_).transpose();
    }
  
    ZMatrixMap(X3ija.data(), naa, nao_) =
        DMatrixConstMap(Viaj_.data(), naa, nao_)
      * ZMatrixMap(G.retptr(tstp,t), nao_, nao_);
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Y2ija.data() + n*naa, nao_, nalpha_) =
          ZMatrixMap(X3ija.data() + n*naa, nalpha_, nao_).transpose();
    }

    // Z^<_nqjl = V_qjb * ~X~^<_nbl
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Zijkl.data() + n*nao3, nao2, nao_) =
          DMatrixConstMap(Vija_.data(), nao2, nalpha_)
        * ZMatrixMap(X1ija.data() + n*naa, nalpha_, nao_);
    }

    // ~Y~^R_jla = Z^T^<_(jl)(nq) * Y^R_nqa   -> X3
    ZMatrixMap(X3ija.data(), nao2, nalpha_) =
        ZMatrixMap(Zijkl.data(), nao2, nao2).transpose()
      * ZMatrixMap(Y2ija.data(), nao2, nalpha_);

    // Z^R_nqjl = V_qjb * ~X~^R_nbl
    for(int n=0; n<nao_; ++n){
      ZMatrixMap(Zijkl.data() + n*nao3, nao2, nao_) =
          DMatrixConstMap(Vija_.data(), nao2, nalpha_)
        * ZMatrixMap(X2ija.data() + n*naa, nalpha_, nao_);
    }

    // ~Y~^R_jla += Z^T^R_(jl)(nq) * Y^R_nqa
    //            + Z^T^R_(jl)(nq) * Y^<_nqa  +-> X3
    ZMatrixMap(X3ija.data(), nao2, nalpha_) +=
        ZMatrixMap(Zijkl.data(), nao2, nao2).transpose()
      * ( ZMatrixMap(Y2ija.data(), nao2, nalpha_)
         +ZMatrixMap(Y1ija.data(), nao2, nalpha_));

    // Sigma_ij += V_ila ~Y~^TR(la)j
    ZMatrixMap(Sigma.retptr(tstp,t), nao_, nao_) -=
        DMatrixConstMap(Vija_.data(), nao_, naa)
      * ZMatrixMap(X3ija.data(), nao_, naa).transpose();
  }
}


void molGF2SolverDecomp::solve_exch(int tstp, GREEN &Sigma, GREEN &G) const {
  solve_exch_les(tstp, Sigma, G);
  solve_exch_ret(tstp, Sigma, G);
  solve_exch_tv(tstp, Sigma, G);
}


void molGF2SolverDecomp::solve(int tstp, GREEN &Sigma, GREEN &G) const {
  assert(G.sig() == Sigma.sig());

  assert(G.ntau() == Sigma.ntau());

  assert(tstp <= G.nt());
  assert(tstp <= Sigma.nt());

  assert(G.size1() == nao_);
  assert(Sigma.size1() == nao_);

  // Set self-energy zero
  Sigma.set_tstp_zero(tstp);
  
  // Perform contractions
  solve_bubble(tstp, Sigma, G);
  solve_exch(tstp, Sigma, G);
}

void molGF2SolverDecomp::solve_loop(int tstp, GREEN &Sigma, GREEN &G) const {
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

  // Lesser bubble  
  for(int t=0; t<=tstp; ++t){
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
                            Sigma.lesptr(t,tstp)[i*nao_+j] += 2*Vija_(i,l,a) * Vija_(n,p,a) * Vija_(k,j,b) * Vija_(q,m,b) * G.lesptr(t,tstp)[l*nao_+k] * (G.retptr(tstp,t)[m*nao_+n]-std::conj(G.lesptr(t,tstp)[n*nao_+m])) * G.lesptr(t,tstp)[p*nao_+q];
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

                        for(int a=0; a<nalpha_; ++a){
                          for(int b=0; b<nalpha_; ++b){
                            *sigRptrij += 2*Vija_(i,l,a) * Vija_(n,p,a) 
                              * Vija_(k,j,b) * Vija_(q,m,b)
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

  }

  // TV bubble
  for(int t=0; t<=ntau; ++t){
      cplx *sigTVptr = Sigma.tvptr(tstp,t);
      cplx *GTVptr = G.tvptr(tstp,t);
      cplx *GVTptr = G.tvptr(tstp,ntau-t);
      for(int i=0; i<nao_; ++i){
        for(int j=0; j<nao_; ++j){
          cplx *sigTVptrij = sigTVptr+i*nao_+j;

            for(int l=0; l<nao_; ++l){
              for(int n=0; n<nao_; ++n){
                for(int p=0; p<nao_; ++p){
                  for(int k=0; k<nao_; ++k){
                    cplx GTVlk = GTVptr[l*nao_+k];
                    for(int q=0; q<nao_; ++q){
                      cplx GTVpq = GTVptr[p*nao_+q];
                      for(int m=0; m<nao_; ++m){
                        cplx GVTnm = -1*(double)sig*std::conj(GVTptr[n*nao_+m]);

                        for(int a=0; a<nalpha_; ++a){
                          for(int b=0; b<nalpha_; ++b){
                            *sigTVptrij += 2*Vija_(i,l,a) * Vija_(n,p,a) * Vija_(k,j,b) * Vija_(q,m,b) * GTVlk * GVTnm * GTVpq;
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
  for(int t=0; t<=tstp; ++t){
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
                            Sigma.lesptr(t,tstp)[i*nao_+j] -= Vija_(i,l,a) * Vija_(n,p,a) * Vija_(q,j,b) * Vija_(k,m,b) * G.lesptr(t,tstp)[l*nao_+k] * (G.retptr(tstp,t)[m*nao_+n]-std::conj(G.lesptr(t,tstp)[n*nao_+m])) * G.lesptr(t,tstp)[p*nao_+q];
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
      cplx *sigTVptr = Sigma.tvptr(tstp,t);
      cplx *GTVptr = G.tvptr(tstp,t);
      cplx *GVTptr = G.tvptr(tstp,ntau-t);
      for(int i=0; i<nao_; ++i){
        for(int j=0; j<nao_; ++j){
          cplx *sigTVptrij = sigTVptr+i*nao_+j;
          
            for(int l=0; l<nao_; ++l){
              for(int n=0; n<nao_; ++n){
                for(int p=0; p<nao_; ++p){
                  for(int k=0; k<nao_; ++k){
                    cplx GTVlk = GTVptr[l*nao_+k];
                    for(int q=0; q<nao_; ++q){
                      cplx GTVpq = GTVptr[p*nao_+q];
                      for(int m=0; m<nao_; ++m){
                        cplx GVTnm = -1*(double)sig*std::conj(GVTptr[n*nao_+m]);

                        for(int a=0; a<nalpha_; ++a){
                          for(int b=0; b<nalpha_; ++b){
                            *sigTVptrij -= Vija_(i,l,a) * Vija_(n,p,a) * Vija_(q,j,b) * Vija_(k,m,b) * GTVlk * GVTnm * GTVpq;
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

                      for(int a=0; a<nalpha_; ++a){
                        for(int b=0; b<nalpha_; ++b){
                          *sigRptrij -= Vija_(i,l,a) * Vija_(n,p,a)
                            * Vija_(q,j,b) * Vija_(k,m,b)
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
  }
}


