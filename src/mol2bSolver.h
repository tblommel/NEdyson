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
  const DTensor<3> &Viaj_;

  void solve_bubble(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_bubble_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_bubble_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_bubble_tv(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

  void solve_exch(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_exch_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_exch_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_exch_tv(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

public:
  molGF2SolverSpinDecomp(const DTensor<3> &Vija, const DTensor<3> &Viaj)
      : molHFSolverSpinDecomp(Vija),
        X1sija(2,nao_, nao_, nalpha_),
        X2sija(2,nao_, nao_, nalpha_),
        X3sija(2,nao_, nao_, nalpha_),
        Y1sija(2,nao_, nao_, nalpha_),
        Y2sija(2,nao_, nao_, nalpha_),
        P1ab(nalpha_,nalpha_),
        P2ab(nalpha_,nalpha_),
        Viaj_(Viaj),
        Zijkl(nao_, nao_, nao_, nao_) {};

  
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
  const DTensor<3> &Viaj_;

  void solve_bubble(int tstp,     GREEN &Sigma, GREEN &G) const;
  void solve_bubble_les(int tstp, GREEN &Sigma, GREEN &G) const;
  void solve_bubble_ret(int tstp, GREEN &Sigma, GREEN &G) const;
  void solve_bubble_tv(int tstp,  GREEN &Sigma, GREEN &G) const;

  void solve_exch(int tstp,       GREEN &Sigma, GREEN &G) const;
  void solve_exch_les(int tstp,   GREEN &Sigma, GREEN &G) const;
  void solve_exch_ret(int tstp,   GREEN &Sigma, GREEN &G) const;
  void solve_exch_tv(int tstp,    GREEN &Sigma, GREEN &G) const;

public:
  molGF2SolverDecomp(const DTensor<3> &Vija, const DTensor<3> &Viaj)
      : molHFSolverDecomp(Vija),
        X1ija(nao_, nao_, nalpha_),
        X2ija(nao_, nao_, nalpha_),
        X3ija(nao_, nao_, nalpha_),
        Y1ija(nao_, nao_, nalpha_),
        Y2ija(nao_, nao_, nalpha_),
        P1ab(nalpha_,nalpha_),
        P2ab(nalpha_,nalpha_),
        Viaj_(Viaj),
        Zijkl(nao_, nao_, nao_, nao_) {};

  
  void solve(int tstp, GREEN &Sigma, GREEN &G) const;

  void solve_loop(int tstp, GREEN &Sigma, GREEN &G) const;

}; // class molGF2SolverDecomp



/////////////////////////////////////////////////////////////////////////////
////                           FULL SPIN                                 ////
/////////////////////////////////////////////////////////////////////////////
class molGF2SolverSpin : public molHFSolverSpin {
private:
  mutable ZTensor<3> C1_aaa;
  mutable ZTensor<3> C2_aaa;
  mutable ZTensor<3> A1_aaa;
  mutable ZTensor<3> B1_aaa;
  mutable ZTensor<3> B2_aaa;
  const DTensor<4> &Uijkl_exch_;

  void solve_les(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_ret(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;
  void solve_tv(int tstp,  std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

public:
  molGF2SolverSpin(const DTensor<4> &U_int, const DTensor<4> &U_exch)
      : molHFSolverSpin(U_int),
        C1_aaa(nao_, nao_, nao_),
        C2_aaa(nao_, nao_, nao_),
        B1_aaa(nao_, nao_, nao_),
        B2_aaa(nao_, nao_, nao_),
        A1_aaa(nao_, nao_, nao_),
        Uijkl_exch_(U_exch) {};

  
  void solve(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

  void solve_loop(int tstp, std::vector<std::reference_wrapper<GREEN>> &Sigma, std::vector<std::reference_wrapper<GREEN>> &G) const;

}; // class molGF2SolverSpinDecomp


/////////////////////////////////////////////////////////////////////////////
////                               FULL                                  ////
/////////////////////////////////////////////////////////////////////////////
class molGF2Solver : public molHFSolver {
private:
  mutable ZTensor<3> A1_aaa;
  mutable ZTensor<3> B1_aaa;
  mutable ZTensor<3> B2_aaa;
  const DTensor<4> &Uijkl_exch_;

  void solve_les(int tstp, GREEN &Sigma, GREEN &G) const;
  void solve_ret(int tstp, GREEN &Sigma, GREEN &G) const;
  void solve_tv(int tstp,  GREEN &Sigma, GREEN &G) const;

public:
  molGF2Solver(const DTensor<4> &U_int, const DTensor<4> &U_exch)
      : molHFSolver(U_int),
        B1_aaa(nao_, nao_, nao_),
        B2_aaa(nao_, nao_, nao_),
        A1_aaa(nao_, nao_, nao_),
        Uijkl_exch_(U_exch) {};

  
  void solve(int tstp, GREEN &Sigma, GREEN &G) const;

  void solve_loop(int tstp, GREEN &Sigma, GREEN &G) const;

}; // class molGF2SolverSpinDecomp


/////////////////////////////////////////////////////////////////////////////
////                      TTI SPIN DECOMP                                ////
/////////////////////////////////////////////////////////////////////////////

class tti_molGF2SolverSpinDecomp {
private:
  const int nao_;
  const int nalpha_;
  static constexpr int ns_ = 2;
    
  const DTensor<3> &Vija_;

  mutable ZTensor<4> X1sija;
  mutable ZTensor<4> X2sija;
  mutable ZTensor<4> X3sija;
  mutable ZTensor<4> Y1sija;
  mutable ZTensor<4> Y2sija;
  mutable ZTensor<2> P1ab;
  mutable ZTensor<2> P2ab;
  mutable ZTensor<4> Zijkl;
  const DTensor<3> &Viaj_;

  void solve_bubble(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;
  void solve_bubble_les(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;
  void solve_bubble_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;
  void solve_bubble_tv(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;

  void solve_exch(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;
  void solve_exch_les(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;
  void solve_exch_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;
  void solve_exch_tv(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;

public:
  tti_molGF2SolverSpinDecomp(const DTensor<3> &Vija, const DTensor<3> &Viaj)
      : nao_(Vija.shape()[0]),
        nalpha_(Vija.shape()[2]),
        Vija_(Vija),
        X1sija(2,nao_, nao_, nalpha_),
        X2sija(2,nao_, nao_, nalpha_),
        X3sija(2,nao_, nao_, nalpha_),
        Y1sija(2,nao_, nao_, nalpha_),
        Y2sija(2,nao_, nao_, nalpha_),
        P1ab(nalpha_,nalpha_),
        P2ab(nalpha_,nalpha_),
        Viaj_(Viaj),
        Zijkl(nao_, nao_, nao_, nao_) {};

  
  void solve(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;

  void solve_loop(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;

  int nao() const { return nao_; }
  int nalpha() const { return nalpha_; }
  const DTensor<3> &Vija() const { return Vija_; }
}; // class tti_molGF2SolverSpinDecomp

/////////////////////////////////////////////////////////////////////////////
////                         TTI DECOMP                                  ////
/////////////////////////////////////////////////////////////////////////////
class tti_molGF2SolverDecomp {
private:
  const int nao_;
  const int nalpha_;
    
  const DTensor<3> &Vija_;

  mutable ZTensor<3> X1ija;
  mutable ZTensor<3> X2ija;
  mutable ZTensor<3> X3ija;
  mutable ZTensor<3> Y1ija;
  mutable ZTensor<3> Y2ija;
  mutable ZTensor<2> P1ab;
  mutable ZTensor<2> P2ab;
  mutable ZTensor<4> Zijkl;
  const DTensor<3> &Viaj_;

  void solve_bubble(int tstp,     TTI_GREEN &Sigma, TTI_GREEN &G) const;
  void solve_bubble_les(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const;
  void solve_bubble_ret(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const;
  void solve_bubble_tv(int tstp,  TTI_GREEN &Sigma, TTI_GREEN &G) const;

  void solve_exch(int tstp,       TTI_GREEN &Sigma, TTI_GREEN &G) const;
  void solve_exch_les(int tstp,   TTI_GREEN &Sigma, TTI_GREEN &G) const;
  void solve_exch_ret(int tstp,   TTI_GREEN &Sigma, TTI_GREEN &G) const;
  void solve_exch_tv(int tstp,    TTI_GREEN &Sigma, TTI_GREEN &G) const;

public:
  tti_molGF2SolverDecomp(const DTensor<3> &Vija, const DTensor<3> &Viaj)
      : nao_(Vija.shape()[0]),
        nalpha_(Vija.shape()[2]),
        Vija_(Vija),
        X1ija(nao_, nao_, nalpha_),
        X2ija(nao_, nao_, nalpha_),
        X3ija(nao_, nao_, nalpha_),
        Y1ija(nao_, nao_, nalpha_),
        Y2ija(nao_, nao_, nalpha_),
        P1ab(nalpha_,nalpha_),
        P2ab(nalpha_,nalpha_),
        Viaj_(Viaj),
        Zijkl(nao_, nao_, nao_, nao_) {};

  
  void solve(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const;

  void solve_loop(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const;

  int nao() const { return nao_; }
  int nalpha() const { return nalpha_; }
  const DTensor<3> &Vija() const { return Vija_; }

}; // class tti_molGF2SolverDecomp

/////////////////////////////////////////////////////////////////////////////
////                           TTI FULL SPIN                             ////
/////////////////////////////////////////////////////////////////////////////
class tti_molGF2SolverSpin {
private:
  const int nao_;
  static constexpr int ns_ = 2;
  
  const DTensor<4> &Uijkl_;

  mutable ZTensor<3> C1_aaa;
  mutable ZTensor<3> C2_aaa;
  mutable ZTensor<3> A1_aaa;
  mutable ZTensor<3> B1_aaa;
  mutable ZTensor<3> B2_aaa;
  const DTensor<4> &Uijkl_exch_;

  void solve_les(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;
  void solve_ret(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;
  void solve_tv(int tstp,  std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;

public:
  tti_molGF2SolverSpin(const DTensor<4> &U_int, const DTensor<4> &Uijkl_exch)
      : nao_(U_int.shape()[0]),
        Uijkl_(U_int),
        C1_aaa(nao_, nao_, nao_),
        C2_aaa(nao_, nao_, nao_),
        B1_aaa(nao_, nao_, nao_),
        B2_aaa(nao_, nao_, nao_),
        A1_aaa(nao_, nao_, nao_),
        Uijkl_exch_(Uijkl_exch) {};

  
  void solve(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;

  void solve_loop(int tstp, std::vector<std::reference_wrapper<TTI_GREEN>> &Sigma, std::vector<std::reference_wrapper<TTI_GREEN>> &G) const;

  int nao() const { return nao_; }
  const DTensor<4> &Uijkl() const { return Uijkl_; }
}; // classt tti_molGF2SolverSpin


/////////////////////////////////////////////////////////////////////////////
////                            TTI FULL                                 ////
/////////////////////////////////////////////////////////////////////////////
class tti_molGF2Solver {
private:
  const int nao_;
  
  const DTensor<4> &Uijkl_;

  mutable ZTensor<3> A1_aaa;
  mutable ZTensor<3> B1_aaa;
  mutable ZTensor<3> B2_aaa;
  const DTensor<4> &Uijkl_exch_;

  void solve_les(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const;
  void solve_ret(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const;
  void solve_tv(int tstp,  TTI_GREEN &Sigma, TTI_GREEN &G) const;

public:
  tti_molGF2Solver(const DTensor<4> &U_int, const DTensor<4> &u_exch)
      : nao_(U_int.shape()[0]),
        Uijkl_(U_int),
        B1_aaa(nao_, nao_, nao_),
        B2_aaa(nao_, nao_, nao_),
        A1_aaa(nao_, nao_, nao_),
        Uijkl_exch_(u_exch) {};

  
  void solve(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const;

  void solve_loop(int tstp, TTI_GREEN &Sigma, TTI_GREEN &G) const;

  int nao() const { return nao_; }
  const DTensor<4> &Uijkl() const { return Uijkl_; }

}; // class tti_molGF2Solver
} // namespace NEdyson
#endif
