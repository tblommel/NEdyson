//
// Created by tblommel on 6/2/20
//

#include "mol2bSolver.h"

namespace NEdyson{

//////////////////////////////////////////////////////////////////////////
//                      Spin Decomp Functions                           //
//////////////////////////////////////////////////////////////////////////

void molGF2SolverSpinDecomp::solve_bubble_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;
  for(int t=0; t<=tstp; ++t){
    // X^<_s'naq = V_nap * G^<_s'pq
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, naa, nao_).noalias() = DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_);
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
        (ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_).transpose() 
       - ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).conjugate()) 
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
      ZMatrixMap(Sigma[s].get().lesptr(t, tstp), nao_, nao_).noalias() += 
          ZMatrixMap(Y1sija.data() + s*nao2*nalpha_, nao_, naa)
        * ZMatrixMap(X2sija.data(), naa, nao_);
    }
  }
}

void molGF2SolverSpinDecomp::solve_bubble_tv(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
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


void molGF2SolverSpinDecomp::solve_bubble_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;
  for(int t=0; t<=tstp; ++t){
    //  X^<_s'naq(tstp,t) = V_nap * G^<_s'pq(tstp,t)
    //                    = V_nap * -(G^<_s'pq(t,tstp))^\dagger
    //
    //  X^R_s'naq(tstp,t) = V_nap * G^R_s'pq(tstp,t)
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, naa, nao_).noalias() = -1 * DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).adjoint();
      ZMatrixMap(X2sija.data() + sp*nao2*nalpha_, naa, nao_).noalias() =  DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_);
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
         ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_).conjugate() 
       * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

      ZMatrixMap(X2sija.data() + sp*nao2*nalpha_, nao_, naa).noalias() = 
         ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).transpose() 
       * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

      ZMatrixMap(X3sija.data() + sp*nao2*nalpha_, nao_, naa).noalias() = 
        (ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).transpose() 
       - ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_).conjugate())
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
      ZMatrixMap(Sigma[s].get().retptr(tstp, t), nao_, nao_).noalias() += 
          ZMatrixMap(Y2sija.data() + s*nao2*nalpha_, nao_, naa)
        * ZMatrixMap(X2sija.data() + nao2*nalpha_, naa, nao_)
        + ZMatrixMap(Y1sija.data() + s*nao2*nalpha_, nao_, naa)
        * ZMatrixMap(X2sija.data() + nao2*nalpha_, naa, nao_)
        + ZMatrixMap(Y2sija.data() + s*nao2*nalpha_, nao_, naa)
        * ZMatrixMap(X1sija.data() + nao2*nalpha_, naa, nao_);
    }
  }
}


void molGF2SolverSpinDecomp::solve_bubble(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  solve_bubble_les(tstp, Sigma, G);
  solve_bubble_ret(tstp, Sigma, G);
  solve_bubble_tv(tstp, Sigma, G);
}


void molGF2SolverSpinDecomp::solve_exch_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int naa = nao_ * nalpha_;
  for(int t=0; t<=tstp; ++t){
    for(int s=0; s<ns_; ++s){
      // X^>_nkb = (G^>T_s)_nm (V^T)_m(kb)
      //         = (G^RT_snm - G^L*_smn) V^T_m(kb)
      ZMatrixMap(X1sija.data(), nao_, naa) = 
          (ZMatrixMap(G[s].get().retptr(tstp,t), nao_, nao_).transpose()
        -  ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_).conjugate())
        *  DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

      // Transpose last two arguments of X
      for(int n=0; n<nao_; ++n){
        ZMatrixMap(X2sija.data() + n*naa, nalpha_, nao_) = 
            ZMatrixMap(X1sija.data() + n*naa, nao_, nalpha_).transpose();
      }

      // X^<_nbl = X^>_nbk * G^<T_kl
      ZMatrixMap(X1sija.data(), naa, nao_) = 
          ZMatrixMap(X2sija.data(), naa, nao_)
        * ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_).transpose();
      
      // Z^<_nqjl = V_qjb * X^<_nbl
      for(int n=0; n<nao_; ++n){
        ZMatrixMap(Zijkl.data() + n*nao3, nao2, nao_) = 
            DMatrixConstMap(Vija_.data(), nao2, nalpha_)
          * ZMatrixMap(X1sija.data() + n*naa, nalpha_, nao_);
      }

      // Y^<_naq = V_nap G^<_pq
      ZMatrixMap(Y1sija.data(), naa, nao_) = 
          DMatrixConstMap(Viaj_.data(), naa, nao_)
        * ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_);

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
      ZMatrixMap(Sigma[s].get().lesptr(t,tstp), nao_, nao_) -= 
          DMatrixConstMap(Vija_.data(), nao_, naa)
        * ZMatrixMap(X1sija.data(), nao_, naa).transpose();
    }
  }
}


void molGF2SolverSpinDecomp::solve_exch_tv(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
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



void molGF2SolverSpinDecomp::solve_exch_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int naa = nao_ * nalpha_;

  for(int t=0; t<=tstp; ++t){
    for(int s=0; s<ns_; ++s){
      // X^>_nkb = (G^>T_s)_nm (V^T)_m(kb)         
      //         = (-G^R*_snm + G^LT_smn) V^T_m(kb) -> X10
      // X^<_nkb = (G^<T_s)_nm (V^T)_m(kb)          -> X20
      // X^A_nkb = (G^AT_s)_nm (V^T)_m(kb)
      //         = (G^T*_s)_nm (V^T)_m(kb)          -> X30
      ZMatrixMap(X1sija.data(), nao_, naa) = 
          (-ZMatrixMap(G[s].get().retptr(tstp,t), nao_, nao_).conjugate()
           +ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_).transpose())
        *  DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
      ZMatrixMap(X2sija.data(), nao_, naa) = 
          ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_).transpose()
        * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
      ZMatrixMap(X3sija.data(), nao_, naa) = 
          ZMatrixMap(G[s].get().retptr(tstp,t), nao_, nao_).conjugate()
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
        * ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_).conjugate();
      ZMatrixMap(X2sija.data(), naa, nao_) = 
          ZMatrixMap(X2sija.data() + naa*nao_, naa, nao_)
        * ZMatrixMap(G[s].get().retptr(tstp,t), nao_, nao_).transpose()
        - ZMatrixMap(X3sija.data() + naa*nao_, naa, nao_)
        * ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_).conjugate();


      // Y^<_naq = V_nap G^<_pq
      //         = - V_nap G^<\dagger_pq -> Y10
      // Y^R_naq = V_nap G^R_pq          -> Y20
      ZMatrixMap(Y1sija.data(), naa, nao_) = 
        - DMatrixConstMap(Viaj_.data(), naa, nao_)
        * ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_).adjoint();
      ZMatrixMap(Y2sija.data(), naa, nao_) = 
          DMatrixConstMap(Viaj_.data(), naa, nao_)
        * ZMatrixMap(G[s].get().retptr(tstp,t), nao_, nao_);

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
      ZMatrixMap(Sigma[s].get().retptr(tstp,t), nao_, nao_) -= 
          DMatrixConstMap(Vija_.data(), nao_, naa)
        * ZMatrixMap(X1sija.data(), nao_, naa).transpose();
    }
  }
}


void molGF2SolverSpinDecomp::solve_exch(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  solve_exch_les(tstp, Sigma, G);
  solve_exch_ret(tstp, Sigma, G);
  solve_exch_tv(tstp, Sigma, G);
}

void molGF2SolverSpinDecomp::solve(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
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
  
  solve_bubble(tstp, Sigma, G);
  solve_exch(tstp, Sigma, G);
}

void molGF2SolverSpinDecomp::solve_loop(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
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
  for(int t=0; t<=tstp; ++t){
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
                            Sigma[s].get().lesptr(t,tstp)[i*nao_+j] += Vija_(i,l,a) * Vija_(n,p,a) * Vija_(k,j,b) * Vija_(q,m,b) * G[s].get().lesptr(t,tstp)[l*nao_+k] * (G[sp].get().retptr(tstp,t)[m*nao_+n]-std::conj(G[sp].get().lesptr(t,tstp)[n*nao_+m])) * G[sp].get().lesptr(t,tstp)[p*nao_+q];
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

  // Retarded Bubble
  for(int t=0; t<=tstp; ++t){
    for(int s=0; s<ns_; ++s){
      cplx *sigRptr = Sigma[s].get().retptr(tstp,t);
      cplx *GsLptr = G[s].get().lesptr(t,tstp);
      cplx *GsRptr = G[s].get().retptr(tstp,t);
      for(int i=0; i<nao_; ++i){
        for(int j=0; j<nao_; ++j){
          cplx *sigRptrij = sigRptr + i*nao_ + j;
  
          for(int sp=0; sp<ns_; ++sp){
          cplx *GspLptr = G[sp].get().lesptr(t,tstp);
          cplx *GspRptr = G[sp].get().retptr(tstp,t);

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
  for(int t=0; t<=tstp; ++t){
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
                            Sigma[s].get().lesptr(t,tstp)[i*nao_+j] -= Vija_(i,l,a) * Vija_(n,p,a) * Vija_(q,j,b) * Vija_(k,m,b) * G[s].get().lesptr(t,tstp)[l*nao_+k] * (G[s].get().retptr(tstp,t)[m*nao_+n]-std::conj(G[s].get().lesptr(t,tstp)[n*nao_+m])) * G[s].get().lesptr(t,tstp)[p*nao_+q];
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
  for(int t=0; t<=tstp; ++t){
    for(int s=0; s<ns_; ++s){
      cplx *sigRptr = Sigma[s].get().retptr(tstp,t);
      cplx *GsLptr = G[s].get().lesptr(t,tstp);
      cplx *GsRptr = G[s].get().retptr(tstp,t);
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
}

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
  Sigma.set_tstp_zero(tstp);
  
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


//////////////////////////////////////////////////////////////////////////
//                      Spin Functions                                  //
//////////////////////////////////////////////////////////////////////////

void molGF2SolverSpin::solve_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao_ * nao2;

  auto A1_Aa = ZMatrixMap(A1_aaa.data(), nao2, nao_);
  auto A1_Q = ZColVectorMap(A1_aaa.data(), nao3);

  auto C1_aA = ZMatrixMap(C1_aaa.data(), nao_, nao2);
  auto C1_Aa = ZMatrixMap(C1_aaa.data(), nao2, nao_);
  auto B1_aA = ZMatrixMap(B1_aaa.data(), nao_, nao2);
  auto B1_Aa = ZMatrixMap(B1_aaa.data(), nao2, nao_);

  auto Uexch_i_jkl = DMatrixConstMap(Uijkl_exch_.data(), nao_, nao3);
  auto U_i_jkl = DMatrixConstMap(Uijkl_.data(), nao_, nao3);

  for(int t=0; t<=tstp; ++t){
    for(int s=0; s<ns_; ++s){
      auto gL_s_lk = ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_);

      for(int i=0; i<nao_; ++i){
        cplx *sigptr = Sigma[s].get().lesptr(t,tstp) + i*nao_;
        auto Ui_l_np = DMatrixConstMap(Uijkl_.data() + i*nao3, nao_, nao2);
        C1_aA.noalias() = gL_s_lk.transpose() * Ui_l_np;

        for(int sp=0; sp<ns_; ++sp){
          B1_aA.noalias() = (C1_Aa * ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_)).transpose();
          A1_Aa.noalias() = B1_Aa * 
                              (-ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).conjugate()
                               +ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_).transpose());
          ZColVectorMap(sigptr, nao_).noalias() += Uexch_i_jkl * A1_Q;

          if(s==sp) {
            ZColVectorMap(sigptr, nao_).noalias() -= U_i_jkl * A1_Q;
          }
        }
      }
    }
  }
}


void molGF2SolverSpin::solve_tv(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
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


void molGF2SolverSpin::solve_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
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

  for(int t=0; t<=tstp; ++t){
    for(int s=0; s<ns_; ++s){
      auto gR_s_lk = ZMatrixMap(G[s].get().retptr(tstp,t), nao_, nao_);
      auto gL_s_lk = ZMatrixMap(G[s].get().lesptr(t,tstp), nao_, nao_);

      for(int i=0; i<nao_; ++i){
        cplx *sigptr = Sigma[s].get().retptr(tstp,t) + i*nao_;
        auto Ui_l_np = DMatrixConstMap(Uijkl_.data() + i*nao3, nao_, nao2);
        C1_aA.noalias() = gR_s_lk.transpose() * Ui_l_np;
        C2_aA.noalias() = -gL_s_lk.conjugate() * Ui_l_np;

        for(int sp=0; sp<ns_; ++sp){
          auto gL_sp_lk = ZMatrixMap(G[sp].get().lesptr(t,tstp), nao_, nao_);
          auto gR_sp_lk = ZMatrixMap(G[sp].get().retptr(tstp,t), nao_, nao_);
          
          B2_aA.noalias() = - (C2_Aa * gL_sp_lk.adjoint()).transpose();
          B1_aA.noalias() = (C1_Aa * gR_sp_lk - C1_Aa * gL_sp_lk.adjoint() + C2_Aa * gR_sp_lk).transpose();

          A1_Aa.noalias() = B1_Aa * ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).transpose()
                          + B2_Aa * ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_).conjugate();
          ZColVectorMap(sigptr, nao_).noalias() += Uexch_i_jkl * A1_Q;

          if(s==sp) {
            ZColVectorMap(sigptr, nao_).noalias() -= U_i_jkl * A1_Q;
          }
        }
      }
    }
  }
}


void molGF2SolverSpin::solve(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
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

  solve_les(tstp, Sigma, G);
  solve_ret(tstp, Sigma, G);
  solve_tv(tstp, Sigma, G);
}


void molGF2SolverSpin::solve_loop(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
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
  for(int t=0; t<=tstp; ++t){
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
                        Sigma[s].get().lesptr(t,tstp)[i*nao_+j] += Uijkl_(i,l,n,p) * Uijkl_(j,k,q,m) * G[s].get().lesptr(t,tstp)[l*nao_+k] * (G[sp].get().retptr(tstp,t)[m*nao_+n]-std::conj(G[sp].get().lesptr(t,tstp)[n*nao_+m])) * G[sp].get().lesptr(t,tstp)[p*nao_+q];
                        if(s==sp){
                          Sigma[s].get().lesptr(t,tstp)[i*nao_+j] -= Uijkl_(i,l,n,p) * Uijkl_(j,q,k,m) * G[s].get().lesptr(t,tstp)[l*nao_+k] * (G[sp].get().retptr(tstp,t)[m*nao_+n]-std::conj(G[sp].get().lesptr(t,tstp)[n*nao_+m])) * G[sp].get().lesptr(t,tstp)[p*nao_+q];
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
  for(int t=0; t<=tstp; ++t){
    for(int s=0; s<ns_; ++s){
    cplx *sigRptr = Sigma[s].get().retptr(tstp,t);
    cplx *GsLptr = G[s].get().lesptr(t,tstp);
    cplx *GsRptr = G[s].get().retptr(tstp,t);
      for(int i=0; i<nao_; ++i){
        for(int j=0; j<nao_; ++j){
          cplx *sigRptrij = sigRptr+i*nao_+j;
          for(int sp=0; sp<ns_; ++sp){
            cplx *GspLptr = G[sp].get().lesptr(t,tstp);
            cplx *GspRptr = G[sp].get().retptr(tstp,t);
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
}


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
  Sigma.set_tstp_zero(tstp);
  
  solve_les(tstp, Sigma, G);
  solve_tv(tstp, Sigma, G);
  solve_ret(tstp, Sigma, G);
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


//////////////////////////////////////////////////////////////////////////
//                      TTI Spin Decomp Functions                       //
//////////////////////////////////////////////////////////////////////////

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


void tti_molGF2SolverSpinDecomp::solve_bubble(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  solve_bubble_ret(tstp, Sigma, G);
  solve_bubble_tv(tstp, Sigma, G);
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


void tti_molGF2SolverSpinDecomp::solve_exch(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const {
  solve_exch_ret(tstp, Sigma, G);
  solve_exch_tv(tstp, Sigma, G);
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

  Sigma[0].get().set_tstp_zero(tstp);
  Sigma[1].get().set_tstp_zero(tstp);
  
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


//////////////////////////////////////////////////////////////////////////
//                           TTI Decomp Functions                       //
//////////////////////////////////////////////////////////////////////////

void tti_molGF2SolverDecomp::solve_bubble_les(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;

  // X^<_naq = V_nap * G^<_pq
  ZMatrixMap(X1ija.data(), naa, nao_).noalias() = DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G.lesptr(-tstp),nao_,nao_);
  // Y^<_nqa = X^<_naq
  for(int n=0; n<nao_; n++){
    ZMatrixMap(Y1ija.data() +  n*naa, nao_, nalpha_).noalias() = ZMatrixMap(X1ija.data() +  n*naa, nalpha_, nao_).transpose();
  }

  // X^>_nqb(tstp,t) = (G^>)^T(tstp,t) V^T_m(qb)
  //                   = (G^R)^T(tstp,t) + (G^<)^T(tstp,t) ...
  //                   = (G^R)^T(tstp,t) - (G^<)^*(t,tstp) ...
  ZMatrixMap(X1ija.data(), nao_, naa).noalias() = 
    (ZMatrixMap(G.retptr(tstp),nao_,nao_).transpose() 
   - ZMatrixMap(G.lesptr(-tstp),nao_,nao_).conjugate()) 
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
  ZMatrixMap(Sigma.lesptr(-tstp), nao_, nao_).noalias() += 
    2 * ZMatrixMap(Y1ija.data(), nao_, naa)
    * ZMatrixMap(X2ija.data(), naa, nao_);
}

void tti_molGF2SolverDecomp::solve_bubble_tv(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
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


void tti_molGF2SolverDecomp::solve_bubble_ret(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
  int nao2 = nao_ * nao_;
  int naa = nao_ * nalpha_;
  //  X^<_naq(tstp,t) = V_nap * G^<_pq(tstp,t)
  //                  = V_nap * -(G^<_pq(t,tstp))^\dagger
  //
  //  X^R_naq(tstp,t) = V_nap * G^R_pq(tstp,t)
  ZMatrixMap(X1ija.data(), naa, nao_).noalias() = -1 * DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G.lesptr(-tstp),nao_,nao_).adjoint();
  ZMatrixMap(X2ija.data(), naa, nao_).noalias() =  DMatrixConstMap(Viaj_.data(), naa, nao_) * ZMatrixMap(G.retptr(tstp),nao_,nao_);
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
     ZMatrixMap(G.retptr(tstp),nao_,nao_).conjugate() 
   * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

  ZMatrixMap(X2ija.data(), nao_, naa).noalias() = 
     ZMatrixMap(G.lesptr(-tstp),nao_,nao_).transpose() 
   * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();

  ZMatrixMap(X3ija.data(), nao_, naa).noalias() = 
    (ZMatrixMap(G.lesptr(-tstp),nao_,nao_).transpose() 
   - ZMatrixMap(G.retptr(tstp),nao_,nao_).conjugate())
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
  ZMatrixMap(Sigma.retptr(tstp), nao_, nao_).noalias() += 
  2*(ZMatrixMap(Y2ija.data(), nao_, naa)
    * ZMatrixMap(X3ija.data(), naa, nao_)
    + ZMatrixMap(Y1ija.data(), nao_, naa)
    * ZMatrixMap(X3ija.data(), naa, nao_)
    + ZMatrixMap(Y2ija.data(), nao_, naa)
    * ZMatrixMap(X2ija.data(), naa, nao_));
}


void tti_molGF2SolverDecomp::solve_bubble(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
  solve_bubble_ret(tstp, Sigma, G);
  solve_bubble_tv(tstp, Sigma, G);
}


void tti_molGF2SolverDecomp::solve_exch_les(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
  ZMatrixMap(Sigma.lesptr(-tstp), nao_, nao_) = -ZMatrixMap(Sigma.tvptr(tstp,0), nao_, nao_).adjoint();
}


void tti_molGF2SolverDecomp::solve_exch_tv(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
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



void tti_molGF2SolverDecomp::solve_exch_ret(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;
  int naa = nao_ * nalpha_;

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
      (-ZMatrixMap(G.retptr(tstp), nao_, nao_).conjugate()
       +ZMatrixMap(G.lesptr(-tstp), nao_, nao_).transpose())
    *  DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
  for(int n=0; n<nao_; ++n){
    ZMatrixMap(Y1ija.data() + n*naa, nalpha_, nao_) = 
        ZMatrixMap(X2ija.data() + n*naa, nao_, nalpha_).transpose();
  }

  ZMatrixMap(X2ija.data(), nao_, naa) = 
      ZMatrixMap(G.lesptr(-tstp), nao_, nao_).transpose()
    * DMatrixConstMap(Viaj_.data(), naa, nao_).transpose();
  for(int n=0; n<nao_; ++n){
    ZMatrixMap(Y2ija.data() + n*naa, nalpha_, nao_) = 
        ZMatrixMap(X2ija.data() + n*naa, nao_, nalpha_).transpose();
  }

  ZMatrixMap(X2ija.data(), nao_, naa) = 
      ZMatrixMap(G.retptr(tstp), nao_, nao_).conjugate()
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
    * ZMatrixMap(G.lesptr(-tstp), nao_, nao_).conjugate();
  ZMatrixMap(X2ija.data(), naa, nao_) = 
      ZMatrixMap(Y2ija.data(), naa, nao_)
    * ZMatrixMap(G.retptr(tstp), nao_, nao_).transpose()
    - ZMatrixMap(X3ija.data(), naa, nao_)
    * ZMatrixMap(G.lesptr(-tstp), nao_, nao_).conjugate();


  // Y^<_naq = V_nap G^<_pq
  //         = - V_nap G^<\dagger_pq 
  // Y^R_naq = V_nap G^R_pq          
  //
  // Transpose last two arguments
  // < -> Y1
  // R -> Y2
  ZMatrixMap(X3ija.data(), naa, nao_) = 
    - DMatrixConstMap(Viaj_.data(), naa, nao_)
    * ZMatrixMap(G.lesptr(-tstp), nao_, nao_).adjoint();
  for(int n=0; n<nao_; ++n){
    ZMatrixMap(Y1ija.data() + n*naa, nao_, nalpha_) = 
        ZMatrixMap(X3ija.data() + n*naa, nalpha_, nao_).transpose();
  }

  ZMatrixMap(X3ija.data(), naa, nao_) = 
      DMatrixConstMap(Viaj_.data(), naa, nao_)
    * ZMatrixMap(G.retptr(tstp), nao_, nao_);
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
  ZMatrixMap(Sigma.retptr(tstp), nao_, nao_) -= 
      DMatrixConstMap(Vija_.data(), nao_, naa)
    * ZMatrixMap(X3ija.data(), nao_, naa).transpose();
}


void tti_molGF2SolverDecomp::solve_exch(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
  solve_exch_ret(tstp, Sigma, G);
  solve_exch_tv(tstp, Sigma, G);
  solve_exch_les(tstp, Sigma, G);
}


void tti_molGF2SolverDecomp::solve(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
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
  Sigma.set_tstp_zero(tstp);
  
  solve_bubble(tstp, Sigma, G);
  solve_exch(tstp, Sigma, G);
}

void tti_molGF2SolverDecomp::solve_loop(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
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
  Sigma.set_tstp_zero(tstp);

  int sig = G.sig();
  int ntau= G.ntau();

  // Lesser bubble  
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
                        Sigma.lesptr(-tstp)[i*nao_+j] += 2*Vija_(i,l,a) * Vija_(n,p,a) * Vija_(k,j,b) * Vija_(q,m,b) * G.lesptr(-tstp)[l*nao_+k] * (G.retptr(tstp)[m*nao_+n]-std::conj(G.lesptr(-tstp)[n*nao_+m])) * G.lesptr(-tstp)[p*nao_+q];
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
  cplx *sigRptr = Sigma.retptr(tstp);
  cplx *GLptr = G.lesptr(-tstp);
  cplx *GRptr = G.retptr(tstp);
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
                        Sigma.lesptr(-tstp)[i*nao_+j] -= Vija_(i,l,a) * Vija_(n,p,a) * Vija_(q,j,b) * Vija_(k,m,b) * G.lesptr(-tstp)[l*nao_+k] * (G.retptr(tstp)[m*nao_+n]-std::conj(G.lesptr(-tstp)[n*nao_+m])) * G.lesptr(-tstp)[p*nao_+q];
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
  sigRptr = Sigma.retptr(tstp);
  GLptr = G.lesptr(-tstp);
  GRptr = G.retptr(tstp);
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

  Sigma[0].get().set_tstp_zero(tstp);
  Sigma[1].get().set_tstp_zero(tstp);

  solve_ret(tstp, Sigma, G);
  solve_tv(tstp, Sigma, G);
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


//////////////////////////////////////////////////////////////////////////
//                         TTI Full Functions                           //
//////////////////////////////////////////////////////////////////////////

void tti_molGF2Solver::solve_les(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
  ZMatrixMap(Sigma.lesptr(-tstp), nao_, nao_) = -ZMatrixMap(Sigma.tvptr(tstp,0), nao_, nao_).adjoint();
}


void tti_molGF2Solver::solve_tv(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
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


void tti_molGF2Solver::solve_ret(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
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

  auto gles = ZMatrixMap(G.lesptr(-tstp), nao_, nao_);
  auto gret = ZMatrixMap(G.retptr(tstp), nao_, nao_);

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
    ZRowVectorMap(Sigma.retptr(tstp) + i * nao_, nao_) += Uexch_j_qkm * A1_Q;
  }
}


void tti_molGF2Solver::solve(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
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
  Sigma.set_tstp_zero(tstp);
  
  solve_tv(tstp, Sigma, G);
  solve_les(tstp, Sigma, G);
  solve_ret(tstp, Sigma, G);
}

void tti_molGF2Solver::solve_loop(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const {
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
  Sigma.set_tstp_zero(tstp);

  int sig = G.sig();
  int ntau= G.ntau();

  // Lesser
  for(int i=0; i<nao_; ++i){
    for(int j=0; j<nao_; ++j){          
      for(int l=0; l<nao_; ++l){
        for(int n=0; n<nao_; ++n){
          for(int p=0; p<nao_; ++p){
            for(int k=0; k<nao_; ++k){
              for(int q=0; q<nao_; ++q){
                for(int m=0; m<nao_; ++m){
                  Sigma.lesptr(-tstp)[i*nao_+j] += Uijkl_(i,l,n,p) * (2*Uijkl_(j,k,q,m)-Uijkl_(j,q,k,m)) * G.lesptr(-tstp)[l*nao_+k] * (G.retptr(tstp)[m*nao_+n]-std::conj(G.lesptr(-tstp)[n*nao_+m])) * G.lesptr(-tstp)[p*nao_+q];
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
  cplx *sigRptr = Sigma.retptr(tstp);
  cplx *GLptr = G.lesptr(-tstp);
  cplx *GRptr = G.retptr(tstp);
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




}// namespace NEdyson
