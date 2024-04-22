//
// Created by Thomas Blommel on 5/21/20
// This defines teh functions that solve for the HF self energy using the non-decomp interaction tensor
//

#ifndef hubHF_SIGMA_DEFN
#define hubHF_SIGMA_DEFN

#include "greens.h"
#include "utils.h"

namespace NEdyson{

class hubHFSolver{
protected:
  int nao_;
  double U_;
  mutable ZMatrix rho_T;

public:
  //Constructor
  hubHFSolver(double U, int nao) : nao_(nao), 
                                   rho_T(nao,nao) { U_ = U; }

  int nao() const { return nao_; }
  
  void solve_HF(int tstp, cplx *hmf, ZMatrix &rho) const;

  virtual ~hubHFSolver(){};

};


} // namespace

#endif // HF_SIGMA_DEFN
