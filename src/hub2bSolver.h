//
// Created by tblommel on 6/1/20
//

#ifndef hubNEDYSON_SIGMA_GF2_DECL
#define hubNEDYSON_SIGMA_GF2_DECL

#include "greens.h"
#include "utils.h"
#include "hubHFSolver.h"
#include <chrono>

namespace NEdyson{

/////////////////////////////////////////////////////////////////////////////
////                               FULL                                  ////
/////////////////////////////////////////////////////////////////////////////
class hubGF2Solver : public hubHFSolver {
private:
  mutable ZTensor<3> A1_aaa;
  mutable ZTensor<3> B1_aaa;
  mutable ZTensor<3> B2_aaa;

  void solve_les(int tstp, GREEN &Sigma, GREEN &G) const;
  void solve_ret(int tstp, GREEN &Sigma, GREEN &G) const;

public:
  hubGF2Solver(double U, int nao)
      : hubHFSolver(U, nao),
        B1_aaa(nao_, nao_, nao_),
        B2_aaa(nao_, nao_, nao_),
        A1_aaa(nao_, nao_, nao_) {};

  
  void solve(int tstp, GREEN &Sigma, GREEN &G) const;

}; // class molGF2SolverSpinDecomp

} // namespace NEdyson
#endif
