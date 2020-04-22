#ifndef INTEGRALS_DECL
#define INTEGRALS_DECL

#include "elementops.hpp"
#include "greens.hpp"
#include "integration.hpp"
#include <complex>

namespace NEdyson{

void CTV1(integration::Integrator &I, green_func &A, green_func &Acc, green_func &B, int n, double dt);

void CTV2(integration::Integrator &I, green_func &A, green_func &B, int n, int m, double beta, cplx *res);

void CTV3(integration::Integrator &I, green_func &A, green_func &B, int n, int m, double beta, cplx *res);

void Cles2_tstp(integration::Integrator &I, green_func &A, green_func &Acc,green_func &B, int m, double dt, cplx *Q);

void Cles3_tstp(integration::Integrator &I,green_func &A, green_func &B, int m, double beta, cplx *Q);


}//namespace

#endif
