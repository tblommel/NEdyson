//
// Created by tblommel on 6/1/20
//

#ifndef NEDYSON_SIGMA_GF2_DECL
#define NEDYSON_SIGMA_GF2_DECL

#include "greens.h"
#include "utils.h"
#include "molHFSolver.h"

namespace NEdyson{

/////////////////////////////////////////////////////////////////////////////
////                      SPIN DECOMP                                    ////
/////////////////////////////////////////////////////////////////////////////

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



/////////////////////////////////////////////////////////////////////////////
////                           DECOMP                                    ////
/////////////////////////////////////////////////////////////////////////////
class molGF2SolverDecomp : public molHFSolverDecomp {
private:
  mutable ZTensor<3> X1ija;
  mutable ZTensor<3> X2ija;
  mutable ZTensor<3> X3ija;
  mutable ZTensor<3> Y1ija;
  mutable ZTensor<3> Y2ija;
  mutable ZTensor<2> P1ab;
  mutable ZTensor<2> P2ab;
  mutable ZTensor<4> Zijkl;
  DTensor<3> Viaj_;

  void solve_bubble(int tstp,     GREEN &Sigma, GREEN &G) const;
  void solve_bubble_les(int tstp, GREEN &Sigma, GREEN &G) const;
  void solve_bubble_ret(int tstp, GREEN &Sigma, GREEN &G) const;
  void solve_bubble_tv(int tstp,  GREEN &Sigma, GREEN &G) const;

  void solve_exch(int tstp,       GREEN &Sigma, GREEN &G) const;
  void solve_exch_les(int tstp,   GREEN &Sigma, GREEN &G) const;
  void solve_exch_ret(int tstp,   GREEN &Sigma, GREEN &G) const;
  void solve_exch_tv(int tstp,    GREEN &Sigma, GREEN &G) const;

  void TransposeV();

public:
  molGF2SolverDecomp(const DTensor<3> &Vija)
      : molHFSolverDecomp(Vija),
        X1ija(nao_, nao_, nalpha_),
        X2ija(nao_, nao_, nalpha_),
        X3ija(nao_, nao_, nalpha_),
        Y1ija(nao_, nao_, nalpha_),
        Y2ija(nao_, nao_, nalpha_),
        P1ab(nalpha_,nalpha_),
        P2ab(nalpha_,nalpha_),
        Viaj_(nao_, nalpha_, nao_),
        Zijkl(nao_, nao_, nao_, nao_) { TransposeV(); };

  
  void solve(int tstp, GREEN &Sigma, GREEN &G) const;

  void solve_loop(int tstp, GREEN &Sigma, GREEN &G) const;

}; // class molGF2SolverSpinDecomp



/////////////////////////////////////////////////////////////////////////////
////                           FULL SPIN                                 ////
/////////////////////////////////////////////////////////////////////////////
class molGF2SolverSpin : public molHFSolverSpin {
private:
  mutable ZTensor<3> C1_aaa;
  mutable ZTensor<3> C2_aaa;
  mutable ZTensor<3> C3_aaa;
  mutable ZTensor<3> A1_aaa;
  mutable ZTensor<3> A2_aaa;
  mutable ZTensor<3> B1_aaa;
  mutable ZTensor<3> B2_aaa;
  DTensor<4> Uijkl_exch_;

  void solve_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_tv(int tstp,  std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

  void make_U_exch();

public:
  molGF2SolverSpin(const DTensor<4> &U_int)
      : molHFSolverSpin(U_int),
        C1_aaa(nao_, nao_, nao_),
        C2_aaa(nao_, nao_, nao_),
        C3_aaa(nao_, nao_, nao_),
        B1_aaa(nao_, nao_, nao_),
        B2_aaa(nao_, nao_, nao_),
        A1_aaa(nao_, nao_, nao_),
        A2_aaa(nao_, nao_, nao_),
        Uijkl_exch_(nao_, nao_, nao_, nao_) { make_U_exch(); };

  
  void solve(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

  void solve_loop(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

}; // class molGF2SolverSpinDecomp




} // namespace NEdyson
#endif
