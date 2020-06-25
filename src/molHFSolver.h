//
// Created by Thomas Blommel on 5/21/20
// This defines teh functions that solve for the HF self energy using the non-decomp interaction tensor
//

#ifndef HF_SIGMA_DEFN
#define HF_SIGMA_DEFN

#include "greens.h"
#include "utils.h"

namespace NEdyson{

class molHFSolver{
protected:
  const int nao_;
  const DTensor<4> &Uijkl_;

private:
  mutable ZMatrix rho_T;

public:
  //Constructor
  molHFSolver(const DTensor<4> &Uijkl) : nao_(Uijkl.shape()[0]), 
                                         Uijkl_(Uijkl),
                                         rho_T(nao_,nao_) {}

  int nao() const { return nao_; }
  
  void solve_HF_loop(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const;
  
  void solve_HF(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const;

  virtual ~molHFSolver(){};

  const DTensor<4> &Uijkl() { return Uijkl_; }

};



class molHFSolverDecomp{
protected:
  const int nao_;
  const int nalpha_;
  const DTensor<3> &Vija_;

  mutable ZTensor<3> tmp_;
  mutable ZRowVector Xa_;
  mutable ZMatrix rho_T;

public:
  molHFSolverDecomp(const DTensor<3> &Vija) : nao_(Vija.shape()[0]),
                                              nalpha_(Vija.shape()[2]),
                                              Vija_(Vija),
                                              tmp_(nao_,nao_,nao_),
                                              Xa_(nalpha_),
                                              rho_T(nao_,nao_) {}

  int nao() const {return nao_;}
  
  int nalpha() const {return nalpha_;};

  void solve_HF(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const;

  void solve_HF_loop(int tstp, ZTensor<3> &hmf, const ZTensor<2> &rho) const;

};


class molHFSolverSpin{
protected:
  const int nao_;
  static constexpr int ns_ = 2;
    
  const DTensor<4> &Uijkl_;

private:
  mutable ZMatrix rho_T;

public:
  molHFSolverSpin(const DTensor<4> &U_int) : Uijkl_(U_int),
                                             nao_(U_int.shape()[0]),
                                             rho_T(nao_,nao_) {}

  int nao() const {return nao_;}

  // hmf must have size (nt,2,nao,nao)
  // rho must have size (2,nao,nao)
  void solve_HF(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const;

  void solve_HF_loop(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const;
  
  const DTensor<4> &Uijkl() { return Uijkl_; }
};


class molHFSolverSpinDecomp{
protected:
  const int nao_;
  const int nalpha_;
  static constexpr int ns_ = 2;

  const DTensor<3> &Vija_;

private:
  mutable ZMatrix rho_T;
  mutable ZTensor<3> tmp_;
  mutable ZRowVector Xa_;

public:
  molHFSolverSpinDecomp(const DTensor<3> &Vija) : Vija_(Vija),
                                             nao_(Vija.shape()[0]),
                                             nalpha_(Vija.shape()[2]),
                                             rho_T(nao_,nao_),
                                             tmp_(nao_,nao_,nao_),
                                             Xa_(nalpha_) {}

  int nao() const {return nao_;}
  
  int nalpha() const {return nalpha_;}

  // hmf must have size (nt,2,nao,nao)
  // rho must have size (2,nao,nao)
  void solve_HF(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const;

  void solve_HF_loop(int tstp, ZTensor<4> &hmf, const ZTensor<3> &rho) const;
};




class tti_molHFSolver{
protected:
  const int nao_;
  const DTensor<4> &Uijkl_;

private:
  mutable ZMatrix rho_T;

public:
  //Constructor
  tti_molHFSolver(const DTensor<4> &Uijkl) : nao_(Uijkl.shape()[0]), 
                                         Uijkl_(Uijkl),
                                         rho_T(nao_,nao_) {}

  int nao() const { return nao_; }
  
  void solve_HF_loop(ZTensor<2> &hmf, const ZTensor<2> &rho) const;
  
  void solve_HF(ZTensor<2> &hmf, const ZTensor<2> &rho) const;

  virtual ~tti_molHFSolver(){};

  const DTensor<4> &Uijkl() { return Uijkl_; }

};

class tti_molHFSolverDecomp{
protected:
  const int nao_;
  const int nalpha_;
  const DTensor<3> &Vija_;

  mutable ZTensor<3> tmp_;
  mutable ZRowVector Xa_;
  mutable ZMatrix rho_T;

public:
  tti_molHFSolverDecomp(const DTensor<3> &Vija) : nao_(Vija.shape()[0]),
                                              nalpha_(Vija.shape()[2]),
                                              Vija_(Vija),
                                              tmp_(nao_,nao_,nao_),
                                              Xa_(nalpha_),
                                              rho_T(nao_,nao_) {}

  int nao() const {return nao_;}
  
  int nalpha() const {return nalpha_;};

  void solve_HF(ZTensor<2> &hmf, const ZTensor<2> &rho) const;

  void solve_HF_loop(ZTensor<2> &hmf, const ZTensor<2> &rho) const;

};

class tti_molHFSolverSpin{
protected:
  const int nao_;
  static constexpr int ns_ = 2;
    
  const DTensor<4> &Uijkl_;

private:
  mutable ZMatrix rho_T;

public:
  tti_molHFSolverSpin(const DTensor<4> &U_int) : Uijkl_(U_int),
                                             nao_(U_int.shape()[0]),
                                             rho_T(nao_,nao_) {}

  int nao() const {return nao_;}

  // hmf must have size (2,nao,nao)
  // rho must have size (2,nao,nao)
  void solve_HF(ZTensor<3> &hmf, const ZTensor<3> &rho) const;

  void solve_HF_loop(ZTensor<3> &hmf, const ZTensor<3> &rho) const;
  
  const DTensor<4> &Uijkl() { return Uijkl_; }
};

class tti_molHFSolverSpinDecomp{
protected:
  const int nao_;
  const int nalpha_;
  static constexpr int ns_ = 2;

  const DTensor<3> &Vija_;

private:
  mutable ZMatrix rho_T;
  mutable ZTensor<3> tmp_;
  mutable ZRowVector Xa_;

public:
  tti_molHFSolverSpinDecomp(const DTensor<3> &Vija) : Vija_(Vija),
                                             nao_(Vija.shape()[0]),
                                             nalpha_(Vija.shape()[2]),
                                             rho_T(nao_,nao_),
                                             tmp_(nao_,nao_,nao_),
                                             Xa_(nalpha_) {}

  int nao() const {return nao_;}
  
  int nalpha() const {return nalpha_;}

  // hmf must have size (2,nao,nao)
  // rho must have size (2,nao,nao)
  void solve_HF(ZTensor<3> &hmf, const ZTensor<3> &rho) const;

  void solve_HF_loop(ZTensor<3> &hmf, const ZTensor<3> &rho) const;
};








} // namespace

#endif // HF_SIGMA_DEFN
