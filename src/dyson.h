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

// NOT TTI 
void G0_from_h0(GREEN &G, double mu, const ZMatrix &h0, double beta, double dt);
void G0_from_h0(GREEN &G, const INTEG &I, double mu, const function &eps, double beta, double dt);
void G0_from_h0(GREEN &G, double mu, const DTensor<2> &H0, double beta, double h);
void G0_from_h0(GREEN &G, double mu, const double *H0, double beta, double h);

 
double mat_fourier(GREEN &G, const GREEN &Sigma, double mu, const cplx *hmf, double beta);


void Extrapolate(const INTEG &I, GREEN &G, int n);


double dyson_start(const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt);
double dyson_start(const INTEG &I, GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt);
double dyson_start(const INTEG &I, GREEN &G, const GREEN &Sig, const cplx* hmf, double mu, double beta, double dt);


double dyson_step_les(int n, const INTEG &I, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt);

void dyson_step(int n, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt);
void dyson_step(int n, const INTEG &I, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt);
void dyson_step(int n, const INTEG &I, GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt);




// TTI FUNCTIONS
void G0_from_h0(TTI_GREEN &G, double mu, const DTensor<2> &H0, double beta, double h);
void G0_from_h0(TTI_GREEN &G, double mu, const double *H0, double beta, double h);

 
void Extrapolate(const INTEG &I, TTI_GREEN &G, int n);


double dyson_start(const INTEG &I, TTI_GREEN &G, const TTI_GREEN &Sig, const DTensor<2> &hmf, double mu, double beta, double dt);
double dyson_start(const INTEG &I, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt);


void dyson_step(int n, const INTEG &I, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt);
void dyson_step(int n, const INTEG &I, TTI_GREEN &G, const TTI_GREEN &Sig, const DTensor<2> &hmf, double mu, double beta, double dt);






}//namespace
#endif
