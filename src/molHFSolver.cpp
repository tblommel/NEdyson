//
// Created by Thomas Blommel on 5/21/20
// This implements teh functions that solve for the HF self energy using the non-decomp interaction tensor
//

#ifndef HF_SIGMA_IMPL
#define HF_SIGMA_IMPL

#include "molHFSolver.h"

namespace NEdyson {

void molHFSolver::solve_HF_loop(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const {
  assert(tstp <= hmf.shape()[0]);
  assert(tstp >= 0);
  assert(hmf.shape()[1] == nao_);
  assert(hmf.shape()[2] == nao_);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);

  ZMatrixMap hmfMap = ZMatrixMap(hmf.data()+tstp*nao_*nao_,nao_,nao_);

  for(int i=0; i<nao_; i++){
    for(int j=0; j<nao_; j++){
      for(int k=0; k<nao_; k++){
        for(int l=0; l<nao_; l++){
          hmfMap(i,j) += rho(k,l) * (2*Vijkl_(i,j,k,l) - Vijkl_(i,l,k,j));
        }
      }
    }
  }
}

  
void molHFSolver::solve_HF(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const {
  assert(tstp <= hmf.shape()[0]);
  assert(tstp >= 0);
  assert(hmf.shape()[1] == nao_);
  assert(hmf.shape()[2] == nao_);
  assert(rho.shape()[1] == nao_);
  assert(rho.shape()[2] == nao_);

  int nao2 = nao_*nao_;
  int nao3 = nao2*nao_;
  
  rho_T = ZMatrixConstMap(rho.data(),nao_,nao_).transpose();
  ZRowVectorMap rho_T_flat(rho_T.data(), nao2);

  for(int i=0; i<nao_; ++i){
    DMatrixConstMap V_i_lk_j = DMatrixConstMap(Vijkl_.data() + i*nao3, nao2, nao_);
    ZColVectorMap(hmf.data()+tstp*nao2+i*nao_, nao_).noalias() -= rho_T_flat * V_i_lk_j;
  }
  
  ZColVectorMap(hmf.data()+tstp*nao2, nao2).noalias() += 2*DMatrixConstMap(Vijkl_.data(),nao2,nao2) * ZColVectorConstMap(rho.data(), nao2);
}


void molHFSolverDecomp::solve_HF_loop(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const {
  assert(tstp <= hmf.shape()[0]);
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
  assert(tstp <= hmf.shape()[0]);
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

  
void molHFSolverSpin::solve_HF(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const {





}






} // namespace
#endif
