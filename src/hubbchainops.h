#ifndef HUBB_OPS_DECL
#define HUBB_OPS_DECL

#include <Eigen/Eigen>
#include <complex>

#include "greens.h"
#include "function.h"
#include "greens_tstp.h"
#include "integration.h"
#include "vie2.h"
#include "utils.h"


namespace Hubb{
  void Ham_MF(int tstp, const GREEN &G, const CFUNC &U, const cdmatrix &h0, CFUNC &hmf);

  void Sigma_2B(int tstp, const GREEN &G, const CFUNC &Ut, GREEN &Sigma);

  void Polarization(int tstp, GREEN &G, GREEN_TSTP &Pol);

  void Bubble2(int tstp, GREEN &Sigma, int s1, int s2, const GREEN &G, int g1, int g2, const GREEN_TSTP &Pol, int p1, int p2);

  void Bubble2(int tstp, GREEN &Sigma, int s1, int s2, const GREEN &G, int g1, int g2, const GREEN &Pol, int p1, int p2);

  void Bubble1(int tstp, GREEN_TSTP &Pol, int p1, int p2, const GREEN &A, int a1, int a2, const GREEN &B, int b1, int b2);

  void Bubble1(int tstp, GREEN &Pol, int p1, int p2, const GREEN_TSTP &A, int a1, int a2, const GREEN &B, int b1, int b2);

  void GenTPP(int tstp, double dt, double beta, const GREEN &G, GREEN &Phi, const CFUNC &Ut, GREEN &UxPhi, GREEN &PhixU, GREEN &TPP, const INTEG &I);
  
  void GenTPP(double dt, double beta, const GREEN &G, GREEN &Phi, const CFUNC &Ut, GREEN &PhixU, GREEN &UxPhi, GREEN &TPP, const INTEG &I);

  void Sigma_TPP(int tstp, const GREEN &G, const CFUNC &Ut, const GREEN &TPP, GREEN &Sigma);
}

#endif
