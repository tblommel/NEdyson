#ifndef VIE2_OMP_IMPL
#define VIE2_OMP_IMPL

#include "vie2_omp.h"

namespace NEdyson{


void vie2_timestep_omp(int threads, int tstp, GREEN &TPP, const GREEN &PhixU, const GREEN &UxPhi, const GREEN &Phi, double beta, double dt, const INTEG &I){
  int k=I.k(), ntau = TPP.ntau(), size1=TPP.size1(), es=size1*size1, 
  GREEN_TSTP Qtstp(tstp, ntau, size1);
  Q.get_tstp(tstp, Qtstp);





}


} // namespace
#endif
