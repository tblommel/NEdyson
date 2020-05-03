#ifndef DYSON_DECL
#define DYSON_DECL

#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>

#include "elementops.h"
#include "function.h"
#include "greens.h"
#include "integration.h"
#include "integrals.h"
#include "mat_utils.h"
#include "utils.h"

namespace NEdyson{

 
void G0_from_h0(GREEN &G, double mu, const cdmatrix &h0, double beta, double dt);


void G0_from_h0(GREEN &G, const INTEG &I, double mu, const function &eps, double beta, double dt);

 
double mat_fourier(GREEN &G, const GREEN &Sigma, double mu, const cplx *hmf, double beta);


void Extrapolate(const INTEG &I, GREEN &G, int n);


double dyson_start(const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt);


void dyson_step(int n, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt);
}//namespace
#endif
