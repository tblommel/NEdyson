//
// Created by tblommel on 06/09/2020
//

#include "molNEsim.h"

namespace NEdyson {

SimulationBase::SimulationBase(const gfmol::HartreeFock &hf, 
                               int nt, int ntau, int k, double dt, int nw, double wmax,
                               int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps) 
                               : enuc_(hf.enuc()), 
                                 nao_(hf.nao()), 
                                 eKin_(nt+1), 
                                 ePot_(nt+1), 
                                 I(k) { 
  nw_ = nw;
  wmax_ = wmax;
  nt_ = nt;
  ntau_ = ntau;
  dt_ = dt;
  k_ = k;
  bootstrap_converged = false;

  MatMax_ = MatMax;
  MatTol_ = MatTol;
  BootMax_ = BootMax;
  BootTol_ = BootTol;
  CorrSteps_ = CorrSteps;
}

void SimulationBase::run(){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
  
  // Calculate the free GF
  free_gf();

  // Run the gfmol Matsubara solver
  do_mat();

  // Take the coefficients from gfmol and transform them into equidistant mesh
  L_to_Tau();
  
  // Do the bootstrapping routines
  start = std::chrono::system_clock::now();
  do_boot();
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << "Time [Bootstrapping] = "<<elapsed_seconds.count() << "s\n";

  // Do the timestepping routines
  if (bootstrap_converged == true) {
    for(int tstp = k_+1; tstp <= nt_; tstp++){
      std::cout<<tstp<<std::endl;
      do_tstp(tstp);
    }

    // Calculate the spectral function
    do_spectral();
  }
  else {
    std::cout << "Bootstrapping did not converge after " << BootMax_ << " iterations!" << std::endl;
  }
}


template class Simulation<gfmol::ChebyshevRepr>;
template class Simulation<gfmol::IntermediateRepr>;
template class tti_Simulation<gfmol::ChebyshevRepr>;
template class tti_Simulation<gfmol::IntermediateRepr>;

template class DecompSimulation<gfmol::ChebyshevRepr>;
template class DecompSimulation<gfmol::IntermediateRepr>;
template class tti_DecompSimulation<gfmol::ChebyshevRepr>;
template class tti_DecompSimulation<gfmol::IntermediateRepr>;

template class SpinSimulation<gfmol::ChebyshevRepr>;
template class SpinSimulation<gfmol::IntermediateRepr>;
template class tti_SpinSimulation<gfmol::ChebyshevRepr>;
template class tti_SpinSimulation<gfmol::IntermediateRepr>;

template class DecompSpinSimulation<gfmol::ChebyshevRepr>;
template class DecompSpinSimulation<gfmol::IntermediateRepr>;
template class tti_DecompSpinSimulation<gfmol::ChebyshevRepr>;
template class tti_DecompSpinSimulation<gfmol::IntermediateRepr>;

} // namespace NEdyson
