#ifndef DYSON_DECL
#define DYSON_DECL

#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>

#include "elementops.hpp"
#include "function.hpp"
#include "greens.hpp"
#include "integration.hpp"
#include "integrals.hpp"

#define cdmatrix Eigen::MatrixXcd
#define cplx std::complex<double>
#define INTEG integration::Integrator

namespace NEdyson{

 
void G0_from_h0(green_func  &G, double mu, cdmatrix &h0, double beta, double dt);


void G0_from_h0(green_func &G, INTEG &I, double mu, function &eps, double beta, double dt);


 
void mat_fourier(green_func  &G, green_func  &Sigma, double mu, cplx *hmf, double beta);


void CTV1(INTEG &I, green_func &A, green_func &B, int n, double dt);

 
void CTV2(INTEG &I, green_func  &A, green_func  &B, int n, int m, double beta, cplx *res);

 
void CTV3(INTEG &I, green_func  &A, green_func  &B, int n, int m, double beta, cplx *res);

 
void Cles2(INTEG &I, green_func  &A, green_func  &B, int n, int m, double dt, cplx *res);

 
void Cles3(INTEG &I, green_func  &A, green_func  &B, int n, int m, double beta, cplx *res);

 
void Extrapolate(INTEG &I, green_func  &G, int n);


void GRstart(INTEG &I, green_func &G, green_func &Sig, function &hmf, double mu, double dt);


void GRstep(int tstp, INTEG &I, green_func &G, green_func &Sig, function &hmf, double mu, double dt);


void GTVstart(INTEG &I, green_func &G, green_func &Sig, function &hmf, double mu, double beta, double dt);


void GTVstep(int tstp, INTEG &I, green_func &G, green_func &Sig, function &hmf, double mu, double beta, double dt);


void GLstep(int n, INTEG &I, green_func &G, green_func &Sig, function &hmf, double mu, double beta, double dt);


void GLstart(INTEG &I, green_func &G, green_func &Sig, function &hmf, double mu, double beta, double dt);


void dyson_start(INTEG &I, green_func &G, green_func &Sig, function &hmf, double mu, double beta, double dt);


void dyson_step(int n, INTEG &I, green_func &G, green_func &Sig, function &hmf, double mu, double beta, double dt);






}//namespace
#endif
