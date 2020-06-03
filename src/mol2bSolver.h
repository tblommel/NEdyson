//
// Created by tblommel on 6/1/20
//

#ifndef NEDYSON_SIGMA_GF2_DECL
#define NEDYSON_SIGMA_GF2_DECL

#include "greens.h"
#include "utils.h"
#include "molHFSolver.h"

namespace NEdyson{

class molGF2SolverSpinDecomp : public molHFSolverSpinDecomp {
private:
  mutable ZTensor<4> X1sija;
  mutable ZTensor<4> X2sija;
  mutable ZTensor<4> X3sija;
  mutable ZTensor<4> Y1sija;
  mutable ZTensor<4> Y2sija;
  mutable ZTensor<2> P1ab;
  mutable ZTensor<2> P2ab;
  mutable ZTensor<4> Zijkl;
  DTensor<3> Viaj_;

  void solve_bubble(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_bubble_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_bubble_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_bubble_tv(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

  void solve_exch(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_exch_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_exch_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_exch_tv(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

  void TransposeV();

public:
  molGF2SolverSpinDecomp(const DTensor<3> &Vija)
      : molHFSolverSpinDecomp(Vija),
        X1sija(2,nao_, nao_, nalpha_),
        X2sija(2,nao_, nao_, nalpha_),
        X3sija(2,nao_, nao_, nalpha_),
        Y1sija(2,nao_, nao_, nalpha_),
        Y2sija(2,nao_, nao_, nalpha_),
        P1ab(nalpha_,nalpha_),
        P2ab(nalpha_,nalpha_),
        Viaj_(nao_, nalpha_, nao_),
        Zijkl(nao_, nao_, nao_, nao_) { TransposeV(); };

  
  void solve(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

  void solve_loop(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

}; // class molGF2SolverSpinDecomp


} // namespace NEdyson
#endif
