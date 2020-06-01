//
// Created by Thomas Blommel on 5/21/20
// This defines teh functions that solve for the HF self energy using the non-decomp interaction tensor
//

#ifndef HF_SIGMA_DEFN
#define HF_SIGMA_DEFN

#include "function.h"
#include "greens.h"
#include "utils.h"

namespace NEdyson{

class molHFSolver{
protected:
  const int nao_;
  const DTensor<4> &Vijkl_;

private:
  mutable ZMatrix rho_T;

public:
  //Constructor
  molHFSolver(const DTensor<4> &Vijkl) : nao_(Vijkl.shape()[0]), Vijkl_(Vijkl), rho_T(nao_,nao_) {};

  int nao() const { return nao_; }
  
  void solve_HF_loop(int tstp, const ZTensor<2> &rho, ZTensor<3> &hmf) const;
  
  void solve_HF(int tstp, const ZTensor<2> &rho, ZTensor<3> &hmf) const;

  virtual ~molHFSolver(){};

};


/*
class molHFSolverSpin{
protected:
  const int size1_;
  const DTensor<4> &Vijkl_;
  static constexpr int ns_ = 2;

private:
  mutable ZMatrix rho_T;

public:
  //Constructor
  molHFSolverSpin(const DTensor<4> &Vijkl) : size1_(Vijkl.shape()[0]), Vijkl_(Vijkl), rho_T(size1_,size1_) {};

  int nao() const { return size1_; }
  
  void solve_HF(int tstp, ZMatrix &rhoU, CFUNC &hmfU, ZMatrix &rhoD, CFUNC &hmfD) const ;
  void solve_HF_loop(int tstp, ZMatrix &rhoU, CFUNC &hmfU, ZMatrix &rhoD, CFUNC &hmfD) const ;
  

  virtual ~molHFSolver(){};

};*/

} // namespace

#endif // HF_SIGMA_DEFN
