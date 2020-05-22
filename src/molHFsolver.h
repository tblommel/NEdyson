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

class HFSolver{
protected:
  const int size1_;
  const DTensor<4> &Vijkl_;

private:
  mutable DMatrix rho_lk;

public:
  //Constructor
  HFSolver(const DTensor<4> &Vijkl) : size1_(Vijkl.shape()[0]), Vijkl_(Vijkl), rho_lk(size1_,size1_) {};











};

} // namespace

#endif // HF_SIGMA_DEFN
