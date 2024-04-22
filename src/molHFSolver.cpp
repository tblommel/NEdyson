//
// Created by Thomas Blommel on 5/21/20
// This implements teh functions that solve for the HF self energy using the non-decomp interaction tensor
//

#ifndef HF_SIGMA_IMPL
#define HF_SIGMA_IMPL

#include "molHFSolver.h"

namespace NEdyson {

void molHFSolver::solve_HF_loop(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const {
//  assert(tstp < hmf.shape()[0]);
//  assert(tstp >= 0);
//  assert(hmf.shape()[1] == nao_);
//  assert(hmf.shape()[2] == nao_);
//  assert(rho.shape()[1] == nao_);
//  assert(rho.shape()[2] == nao_);
//
//  ZMatrixMap hmfMap = ZMatrixMap(hmf.data()+tstp*nao_*nao_,nao_,nao_);
//
//  for(int i=0; i<nao_; i++){
//    for(int j=0; j<nao_; j++){
//      for(int k=0; k<nao_; k++){
//        for(int l=0; l<nao_; l++){
//          hmfMap(i,j) += rho(k,l) * (2*Uijkl_(i,j,k,l) - Uijkl_(i,l,k,j));
//        }
//      }
//    }
//  }
}

  
void molHFSolver::solve_HF(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const {
  assert(tstp < hmf.shape()[0]);
  assert(tstp >= 0);
  assert(hmf.shape()[1] == nao_);
  assert(hmf.shape()[2] == nao_);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);

  int nao2 = nao_*nao_;
  int nao3 = nao2*nao_;
  
//  rho_T = ZMatrixConstMap(rho.data(),nao_,nao_).transpose();
//  ZRowVectorMap rho_T_flat(rho_T.data(), nao2);
//
//  for(int i=0; i<nao_; ++i){
//    DMatrixConstMap U_i_lk_j = DMatrixConstMap(Uijkl_.data() + i*nao3, nao2, nao_);
//    ZColVectorMap(hmf.data()+tstp*nao2+i*nao_, nao_).noalias() -= rho_T_flat * U_i_lk_j;
//  }
  
//  ZColVectorMap(hmf.data()+tstp*nao2, nao2).noalias() += 2*DMatrixConstMap(Uijkl_.data(),nao2,nao2) * ZColVectorConstMap(rho.data(), nao2);

//  ZMatrixMap(hmf.data() + tstp*nao2, nao2).diagonal() += U_ * ZMatrixMap(rho.data(), nao_, nao_).diagonal();

}

void molHFSolver::solve_HF(int tstp, ZMatrix &hmf, ZMatrix &rho) const {
  assert(tstp >= 0);

  int nao2 = nao_*nao_;
  int nao3 = nao2*nao_;

//  ZMatrixMap(hmf.data() + tstp*nao2, nao2).diagonal() += U_ * ZMatrixMap(rho.data(), nao_, nao_).diagonal();

}



void molHFSolverDecomp::solve_HF_loop(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const {
  assert(tstp < hmf.shape()[0]);
  assert(tstp >= 0);
  assert(hmf.shape()[1] == nao_);
  assert(hmf.shape()[2] == nao_);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);
  
  ZMatrixMap hmfMap = ZMatrixMap(hmf.data()+tstp*nao_*nao_,nao_,nao_);

  for(int i=0;i<nao_;i++){
    for(int j=0;j<nao_;j++){
      for(int k=0;k<nao_;k++){
        for(int l=0;l<nao_;l++){
          for(int a=0;a<nalpha_;a++){
            hmfMap(i,j) += rho(k,l) * (2 * Vija_(i,j,a) * Vija_(k,l,a) - Vija_(i,l,a) * Vija_(k,j,a));
          }
        }
      }
    }
  }
}

void molHFSolverDecomp::solve_HF(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const {
  assert(tstp < hmf.shape()[0]);
  assert(tstp >= 0);
  assert(hmf.shape()[1] == nao_);
  assert(hmf.shape()[2] == nao_);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);
  
  int nao2 = nao_ * nao_;

  // Hartree Diagram
  auto VmapIa = DMatrixConstMap(Vija_.data(), nao2, nalpha_);
  auto RhomapI = ZRowVectorConstMap(rho.data(), nao2);
  Xa_ = RhomapI * VmapIa;
  auto hmfmapI = ZRowVectorMap(hmf.data()+tstp*nao2, nao2);
  hmfmapI.noalias() += 2 * Xa_ * VmapIa.transpose();

  // Fock Diagram
  auto tmp_12_3 = ZMatrixMap(tmp_.data(), nao2, nao_);
  auto tmp_1_23 = ZMatrixMap(tmp_.data(), nao_, nao2);
  auto V_kj_a_T = DMatrixConstMap(Vija_.data(), nao2, nalpha_).transpose();
  rho_T = ZMatrixConstMap(rho.data(), nao_, nao_).transpose();
  ZRowVectorMap rho_T_flat(rho_T.data(), nao2);
  
  for(int i=0;i<nao_;i++){
    auto V_i_l_a = DMatrixConstMap(Vija_.data() + i*nao_*nalpha_, nao_, nalpha_);
    tmp_1_23 = V_i_l_a * V_kj_a_T;
    ZRowVectorMap(hmf.data() + tstp*nao2 + i*nao_, nao_).noalias() -= rho_T_flat * tmp_12_3;
  }
}

  
void molHFSolverSpin::solve_HF_loop(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const {
  assert(tstp >= 0);
  assert(hmf.shape()[1] > tstp);
  assert(hmf.shape()[0] == 2);
  assert(hmf.shape()[2] == nao_);
  assert(hmf.shape()[3] == nao_);
  assert(rho.shape()[0] == 2);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);
  
  for(int s=0; s<ns_; ++s){
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        for(int k=0; k<nao_; ++k){
          for(int l=0; l<nao_; ++l){
            hmf(s,tstp,i,j) -= Uijkl_(i,l,k,j) * rho(s, k, l);
            for(int sp=0; sp<ns_; sp++){
              hmf(s,tstp,i,j) += rho(sp,k,l) * Uijkl_(i,j,k,l);
            }
          }
        }
      }
    }
  }
}


void molHFSolverSpin::solve_HF(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const {
  assert(tstp >= 0);
  assert(hmf.shape()[1] > tstp);
  assert(hmf.shape()[0] == 2);
  assert(hmf.shape()[2] == nao_);
  assert(hmf.shape()[3] == nao_);
  assert(rho.shape()[0] == 2);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);

  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;

  int nt = hmf.shape()[1];
  
  for(int s=0; s<ns_; ++s){
    rho_T = ZMatrixConstMap(rho.data() + s*nao2, nao_, nao_).transpose();
    ZRowVectorMap rho_T_flat(rho_T.data(), nao2);

    // Fock
    for(int i=0; i<nao_; ++i){
      DMatrixConstMap Ui_lk_j(Uijkl_.data() + i*nao3, nao2, nao_);
      ZColVectorMap(hmf.data() + s*nt*nao2 + tstp*nao2 + i*nao_, nao_).noalias() -= rho_T_flat * Ui_lk_j;
    }

    // Hartree
    for(int sp=0; sp<ns_; sp++){
      ZColVectorMap(hmf.data() + s*nt*nao2 + tstp*nao2, nao2).noalias() += DMatrixConstMap(Uijkl_.data(),nao2,nao2) * ZColVectorConstMap(rho.data() + sp*nao2 ,nao2);
    }
  }
}


void molHFSolverSpinDecomp::solve_HF_loop(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const {
  assert(tstp >= 0);
  assert(hmf.shape()[1] > tstp);
  assert(hmf.shape()[0] == 2);
  assert(hmf.shape()[2] == nao_);
  assert(hmf.shape()[3] == nao_);
  assert(rho.shape()[0] == 2);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);
  
  for(int s=0; s<ns_; ++s){
    for(int i=0; i<nao_; ++i){
      for(int j=0; j<nao_; ++j){
        for(int a=0; a<nalpha_; ++a){
          for(int sp=0; sp<ns_; ++sp){
            for(int k=0; k<nao_; ++k){
              for(int l=0; l<nao_; ++l){
                hmf(s,tstp,i,j) += Vija_(i,j,a) * Vija_(k,l,a) * rho(sp, k, l);
                if(s==sp)
                  hmf(s,tstp,i,j) -= Vija_(i,l,a) * Vija_(k,j,a) * rho(s,k,l);
              }
            }
          }
        }
      }
    }
  }
}


void molHFSolverSpinDecomp::solve_HF(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const {
  assert(tstp >= 0);
  assert(hmf.shape()[1] > tstp);
  assert(hmf.shape()[0] == 2);
  assert(hmf.shape()[2] == nao_);
  assert(hmf.shape()[3] == nao_);
  assert(rho.shape()[0] == 2);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);

  int nao2 = nao_*nao_;
  int nt = hmf.shape()[1];

  rho_T = ZMatrixConstMap(rho.data(), nao_, nao_) + ZMatrixConstMap(rho.data() + nao2, nao_, nao_);

  // Hartree
  auto VmapIa = DMatrixConstMap(Vija_.data(),nao2, nalpha_);
  auto RhomapI= ZRowVectorConstMap(rho_T.data(), nao2);
  Xa_ = RhomapI * VmapIa;
  ZRowVectorMap(hmf.data() + tstp*nao2, nao2).noalias() += Xa_ * VmapIa.transpose();
  ZRowVectorMap(hmf.data() + nt*nao2 + tstp*nao2, nao2).noalias() += Xa_ * VmapIa.transpose();

  // Fock
  auto tmp_12_3 = ZMatrixMap(tmp_.data(), nao2, nao_);
  auto tmp_1_23 = ZMatrixMap(tmp_.data(), nao_, nao2);
  auto V_kj_a_T = DMatrixConstMap(Vija_.data(), nao2, nalpha_).transpose();
  ZRowVectorMap rho_T_flat(rho_T.data(), nao2);
  
  for(int i=0; i<nao_; i++){
    auto V_i_l_a = DMatrixConstMap(Vija_.data() + i*nalpha_*nao_, nao_, nalpha_);
    tmp_1_23 = V_i_l_a * V_kj_a_T;
    for(int s=0; s<ns_; s++){
      rho_T = ZMatrixConstMap(rho.data() + s*nao2, nao_,nao_).transpose();
      ZRowVectorMap(hmf.data() + s*nt*nao2 + tstp*nao2 + i*nao_, nao_).noalias() -= rho_T_flat * tmp_12_3;
    }    
  }
}

} // namespace
#endif
