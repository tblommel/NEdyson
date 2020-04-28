#ifndef VIE2_decl
#define VIE2_decl

#include "function.h"
#include "greens.h"
#include <complex>
#include "elementops.h"
#include "integration.h"
#include "integrals.h"
#include "utils.h"


namespace NEdyson{

void vie2_mat_fourier(GREEN &G, const GREEN &F, const GREEN &Fcc, const GREEN &Q, double beta, int pcf=20);

void vie2_start(GREEN &TPP, const GREEN &PhixU, const GREEN &UxPhi, const GREEN &Phi, double beta, double dt, const INTEG &I);

void vie2_timestep(int tstp, GREEN &TPP, const GREEN &PhixU, const GREEN &UxPhi, const GREEN &Phi, double beta, double dt, const INTEG &I); 

} // namespace
#endif
