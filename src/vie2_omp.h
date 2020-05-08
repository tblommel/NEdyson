#ifndef VIE2_OMP_DECL
#define VIE2_OMP_DECL

#if USE_OMP == 1
#include <omp.h>
#endif

#include "vie2.h"
#include "function.h"
#include "greens.h"
#include <complex>
#include "elementops.h"
#include "integration.h"
#include "integrals.h"
#include "utils.h"


namespace NEdyson{

#if USE_OMP == 1
void vie2_timestep_omp(int threads, int tstp, GREEN &TPP, const GREEN &PhixU, const GREEN &UxPhi, const GREEN &Phi, double beta, double dt, const INTEG &I); 
#endif

} // namespace
#endif
