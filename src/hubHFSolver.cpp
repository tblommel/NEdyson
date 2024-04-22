//
// Created by Thomas Blommel on 5/21/20
// This implements teh functions that solve for the HF self energy using the non-decomp interaction tensor
//

#ifndef hubHF_SIGMA_IMPL
#define hubHF_SIGMA_IMPL

#include "hubHFSolver.h"

namespace NEdyson {


void hubHFSolver::solve_HF(int tstp, cplx *hmf, ZMatrix &rho) const {
  assert(tstp >= 0);

  int nao2 = nao_*nao_;
  int nao3 = nao2*nao_;

  ZMatrixMap(hmf + tstp*nao2, nao_, nao_).diagonal() += U_ * ZMatrixMap(rho.data(), nao_, nao_).diagonal();

}

} // namespace
#endif
