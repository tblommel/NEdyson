#ifndef VIE2_decl
#define VIE2_decl

#include "function.hpp"
#include "greens.hpp"
#include <complex>
#include "elementops.hpp"
#include "integration.hpp"
#include "integrals.hpp"

#define cplx std::complex<double>
#define GREENS NEdyson::green_func
#define CFUNC NEdyson::function
#define INTEG integration::Integrator

namespace NEdyson{








void vie2_mat_fourier(GREENS &G, GREENS &F, GREENS &Fcc, GREENS &Q, double beta, int pcf=20);

void vie2_start(GREENS &TPP, GREENS &PhixU, GREENS &UxPhi, GREENS &Phi, double beta, double dt, INTEG &I);

void vie2_timestep(int tstp, GREENS &TPP, GREENS &PhixU, GREENS &UxPhi, GREENS &Phi, double beta, double dt, INTEG &I); 










} // namespace
#endif
