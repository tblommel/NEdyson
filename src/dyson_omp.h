#ifndef DYSON_OMP_DECL
#define DYSON_OMP_DECL

#include "elementops.h"
#include "function.h"
#include "greens.h"
#include "integration.h"
#include "integrals.h"
#include "utils.h"

namespace NEdyson{

#if USE_OMP == 1

void dyson_step_omp(int omp_threads, int n, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt);

#endif // USE_OMP

} // namespace

#endif // DYSON_OMP_DECL
