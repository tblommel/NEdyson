//
// Created by tblommel on 6/2/20
//

#include "mol2bSolver.h"

namespace NEdyson{

//////////////////////////////////////////////////////////////////////////
//                      Spin Decomp Functions                           //
//////////////////////////////////////////////////////////////////////////

void molGF2SolverSpinDecomp::TransposeV()
{ 
  for(int i=0; i<nao_; i++){
    for(int j=0; j<nao_; j++){
      for(int a=0; a<nalpha_; a++){
        Viaj_(i,a,j) = Vija_(i,j,a);
      }
    }
  }
}


void molGF2SolverSpinDecomp::solve_bubble_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  for(int t=0; t<=tstp; ++t){
    // X^<_s'naq = V_nap * G^<_s'pq
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_*nalpha_, nao_).noalias() = DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_) * ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_);
    }
    // Y^<_s'nqa = X^<_s'naq
    for(int sp=0; sp<ns_; ++sp){
      for(int n=0; n<nao_; n++){
        ZMatrixMap(Y1sija.data() + sp*nao2*nalpha_ + n*nalpha_*nao_, nao_, nalpha_).noalias() = ZMatrixMap(X1sija.data() + sp*nao2*nalpha_ + n*nalpha_*nao_, nalpha_, nao_).transpose();
      }
    }

    // X^>_s'nqb(tstp,t) = (G^>_s')^T(tstp,t) V^T_m(qb)
    //                   = (G^R_s')^T(tstp,t) + (G^<_s')^T(tstp,t) ...
    //                   = (G^R_s')^T(tstp,t) - (G^<_s')^*(t,tstp) ...
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_, nao_*nalpha_).noalias() = 
        (ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_).transpose() 
       - ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).conjugate()) 
       * DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_).transpose();
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
      ZMatrixMap(X2sija.data() + k*nao_*nalpha_, nalpha_, nao_).noalias() =
          ZMatrixMap(X1sija.data() + k*nao_*nalpha_, nao_, nalpha_).transpose();
    }

    // Sigma_sij = Y_sika * Z_kaj
    for(int s=0; s<ns_; ++s){
      ZMatrixMap(Sigma[s].get().lesptr(t, tstp), nao_, nao_).noalias() += 
          ZMatrixMap(Y1sija.data() + s*nao2*nalpha_, nao_, nao_*nalpha_)
        * ZMatrixMap(X2sija.data(), nao_*nalpha_, nao_);
    }
  }
}

void molGF2SolverSpinDecomp::solve_bubble_tv(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  int ntau = G[0].get().ntau();
  int sig = G[0].get().sig();
  for(int t=0; t<=ntau; ++t){
    // X^\rceil_s'naq = V_nap * G^\rceil_s'pq
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_*nalpha_, nao_).noalias() = DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_) * ZMatrixMap(G[sp].get().tvptr(tstp,t),nao_,nao_);
    }
    // Y^\rceil_s'nqa = X^\rceil_s'naq
    for(int sp=0; sp<ns_; ++sp){
      for(int n=0; n<nao_; ++n){
        ZMatrixMap(Y1sija.data() + sp*nao2*nalpha_ + n*nalpha_*nao_, nao_, nalpha_).noalias() = ZMatrixMap(X1sija.data() + sp*nao2*nalpha_ + n*nalpha_*nao_, nalpha_, nao_).transpose();
      }
    }

    // X^\lceil_s'nqb(tstp,t) = (G^\lceil_s')^T(t,tstp) V^T_m(qb)
    //                    = -sig(G^\rceil_s')^*(tstp,ntau-t) ...
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_, nao_*nalpha_).noalias() = 
       -(double)sig * ZMatrixMap(G[sp].get().tvptr(tstp,ntau-t),nao_,nao_).conjugate()
       * DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_).transpose();
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
      ZMatrixMap(X2sija.data() + k*nao_*nalpha_, nalpha_, nao_).noalias() =
          ZMatrixMap(X1sija.data() + k*nao_*nalpha_, nao_, nalpha_).transpose();
    }

    // Sigma^\rceil_sij = Y^\rceil_sika * Z^\rceil_kaj
    for(int s=0; s<ns_; ++s){
      ZMatrixMap(Sigma[s].get().tvptr(tstp, t), nao_, nao_).noalias() += 
          ZMatrixMap(Y1sija.data() + s*nao2*nalpha_, nao_, nao_*nalpha_)
        * ZMatrixMap(X2sija.data(), nao_*nalpha_, nao_);
    }
  }
}


void molGF2SolverSpinDecomp::solve_bubble_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  int nao2 = nao_ * nao_;
  for(int t=0; t<=tstp; ++t){
    //  X^<_s'naq(tstp,t) = V_nap * G^<_s'pq(tstp,t)
    //                    = V_nap * -(G^<_s'pq(t,tstp))^\dagger
    //
    //  X^R_s'naq(tstp,t) = V_nap * G^R_s'pq(tstp,t)
    for(int sp = 0; sp<ns_; ++sp){
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_*nalpha_, nao_).noalias() = -1 * DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_) * ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).adjoint();
      ZMatrixMap(X2sija.data() + sp*nao2*nalpha_, nao_*nalpha_, nao_).noalias() =  DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_) * ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_);
    }
    // Y^<_s'nqa = X^<_s'naq
    for(int sp=0; sp<ns_; ++sp){
      for(int n=0; n<nao_; ++n){
        ZMatrixMap(Y1sija.data() + sp*nao2*nalpha_ + n*nalpha_*nao_, nao_, nalpha_).noalias() = ZMatrixMap(X1sija.data() + sp*nao2*nalpha_ + n*nalpha_*nao_, nalpha_, nao_).transpose();
        ZMatrixMap(Y2sija.data() + sp*nao2*nalpha_ + n*nalpha_*nao_, nao_, nalpha_).noalias() = ZMatrixMap(X2sija.data() + sp*nao2*nalpha_ + n*nalpha_*nao_, nalpha_, nao_).transpose();
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
      ZMatrixMap(X1sija.data() + sp*nao2*nalpha_, nao_, nao_*nalpha_).noalias() = 
         ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_).conjugate() 
       * DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_).transpose();

      ZMatrixMap(X2sija.data() + sp*nao2*nalpha_, nao_, nao_*nalpha_).noalias() = 
         ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).transpose() 
       * DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_).transpose();

      ZMatrixMap(X3sija.data() + sp*nao2*nalpha_, nao_, nao_*nalpha_).noalias() = 
        (ZMatrixMap(G[sp].get().lesptr(t,tstp),nao_,nao_).transpose() 
       - ZMatrixMap(G[sp].get().retptr(tstp,t),nao_,nao_).conjugate())
       * DMatrixConstMap(Viaj_.data(), nao_*nalpha_, nao_).transpose();
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
      ZMatrixMap(X1sija.data() + nao2*nalpha_ + k*nao_*nalpha_, nalpha_, nao_).noalias() =
          ZMatrixMap(X1sija.data() + k*nao_*nalpha_, nao_, nalpha_).transpose();
      ZMatrixMap(X2sija.data() + nao2*nalpha_ + k*nao_*nalpha_, nalpha_, nao_).noalias() =
          ZMatrixMap(X2sija.data() + k*nao_*nalpha_, nao_, nalpha_).transpose();
    }

    // Sigma^R_sij = Y^R_sika * Z^R_kaj + Y^< * Z^R + Y^R * Z^<
    for(int s=0; s<ns_; ++s){
      ZMatrixMap(Sigma[s].get().retptr(tstp, t), nao_, nao_).noalias() += 
          ZMatrixMap(Y2sija.data() + s*nao2*nalpha_, nao_, nao_*nalpha_)
        * ZMatrixMap(X2sija.data() + nao2*nalpha_, nao_*nalpha_, nao_)
        + ZMatrixMap(Y1sija.data() + s*nao2*nalpha_, nao_, nao_*nalpha_)
        * ZMatrixMap(X2sija.data() + nao2*nalpha_, nao_*nalpha_, nao_)
        + ZMatrixMap(Y2sija.data() + s*nao2*nalpha_, nao_, nao_*nalpha_)
        * ZMatrixMap(X1sija.data() + nao2*nalpha_, nao_*nalpha_, nao_);
    }
  }
}


void molGF2SolverSpinDecomp::solve_bubble(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
  solve_bubble_les(tstp, Sigma, G);
  solve_bubble_ret(tstp, Sigma, G);
  solve_bubble_tv(tstp, Sigma, G);
}

void molGF2SolverSpinDecomp::solve_exch(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const {
//  solve_exch_les(tstp, Sigma, G);
//  solve_exch_ret(tstp, Sigma, G);
//  solve_exch_tv(tstp, Sigma, G);
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
  //solve_exch(tstp, Sigma, G);
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
}


}// namespace NEdyson
