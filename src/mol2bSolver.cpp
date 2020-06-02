//
// Created by tblommel on 6/2/20
//

#include "mol2bSolver.h"

namespace NEdyson{

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

}


}// namespace NEdyson
