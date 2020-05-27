//
// Created by Thomas Blommel on 5/21/20
// This implements teh functions that solve for the HF self energy using the non-decomp interaction tensor
//

#ifndef HF_SIGMA_IMPL
#define HF_SIGMA_IMPL

#include "molHFSolver.h"

namespace NEdyson {

void molHFSolver::solve_HF_loop(int tstp, ZMatrix &rho, CFUNC &hmf) const {
  assert(tstp <= hmf.nt());
  assert(tstp >= 0);
  assert(hmf.size1() == size1_);
  assert(rho.rows() == size1_);
  assert(rho.cols() == size1_);

  ZMatrixMap hmfMap = ZMatrixMap(hmf.ptr(tstp),size1_,size1_);

  for(int i=0; i<size1_; i++){
    for(int j=0; j<size1_; j++){
      for(int k=0; k<size1_; k++){
        for(int l=0; l<size1_; l++){
          hmfMap(i,j) += rho(k,l) * (2*Vijkl_(i,j,k,l) - Vijkl_(i,l,k,j));
        }
      }
    }
  }
}

  
void molHFSolver::solve_HF(int tstp, ZMatrix &rho, CFUNC &hmf) const {
  assert(tstp <= hmf.nt());
  assert(tstp >= 0);
  assert(hmf.size1() == size1_);
  assert(rho.rows() == size1_);
  assert(rho.cols() == size1_);

  int size12 = size1_*size1_;
  int size13 = size12*size1_;
  
  rho_T = rho.transpose();
  ZRowVectorMap rho_T_flat(rho_T.data(), size12);

  for(int i=0; i<size1_; ++i){
    DMatrixConstMap V_i_lk_j = DMatrixConstMap(Vijkl_.data() + i*size13, size12, size1_);
    ZColVectorMap(hmf.ptr(tstp)+i*size1_, size1_).noalias() -= rho_T_flat * V_i_lk_j;
  }
  
  ZColVectorMap(hmf.ptr(tstp), size12).noalias() += 2*DMatrixConstMap(Vijkl_.data(),size12,size12) * ZColVectorMap(rho.data(), size12);
}


/*
void molHFSolverSpin::solve_HF_loop(int tstp, ZMatrix &rhoU, CFUNC &hmfU, ZMatrix &rhoD, CFUNC &hmfD) const {
  assert(tstp <= hmfU.nt());
  assert(tstp <= hmfD.nt());
  assert(tstp >= 0);
  assert(hmfU.size1() == size1_);
  assert(rhoU.rows() == size1_);
  assert(rhoU.cols() == size1_);
  assert(hmfD.size1() == size1_);
  assert(rhoD.rows() == size1_);
  assert(rhoD.cols() == size1_);

  ZMatrixMap hmfMap = ZMatrixMap(hmf.ptr(tstp),size1_,size1_);

  for(int i=0; i<size1_; i++){
    for(int j=0; j<size1_; j++){
      for(int k=0; k<size1_; k++){
        for(int l=0; l<size1_; l++){
          hmfMap(i,j) += rho(k,l) * (2*Vijkl_(i,j,k,l) - Vijkl_(i,l,k,j));
        }
      }
    }
  }
}

  
void molHFSolverSpin::solve_HF(int tstp, ZMatrix &rhoU, CFUNC &hmfU, ZMatrix &rhoD, CFUNC &hmfD) const 
{
  assert(tstp <= hmfU.nt());
  assert(tstp <= hmfD.nt());
  assert(tstp >= 0);
  assert(hmfU.size1() == size1_);
  assert(rhoU.rows() == size1_);
  assert(rhoU.cols() == size1_);
  assert(hmfD.size1() == size1_);
  assert(rhoD.rows() == size1_);
  assert(rhoD.cols() == size1_);

  int size12 = size1_*size1_;
  int size13 = size12*size1_;
  
  rho_T = rho.transpose();
  ZRowVectorMap rho_T_flat(rho_T.data(), size12);

  for(int i=0; i<size1_; ++i){
    DMatrixConstMap V_i_lk_j = DMatrixConstMap(Vijkl_.data() + i*size13, size12, size1_);
    ZColVectorMap(hmf.ptr(tstp)+i*size1_, size1_).noalias() -= rho_T_flat * V_i_lk_j;
  }
  
  ZColVectorMap(hmf.ptr(tstp), size12).noalias() += 2*DMatrixConstMap(Vijkl_.data(),size12,size12) * ZColVectorMap(rho.data(), size12);
}

*/
} // namespace
#endif
