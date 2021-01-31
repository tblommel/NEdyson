//
// Created by tblommel on 06/09/2020
//

#include "hf_efield_molNEsim.h"

namespace NEdyson {

hfefmolSimulationBase::hfefmolSimulationBase(const gfmol::HartreeFock &hf, 
                               int nt, int ntau, int k, double dt,
                               int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps) 
                               : enuc_(hf.enuc()), 
                                 nao_(hf.nao()), 
                                 eKin_(nt+1), 
                                 ePot_(nt+1),
                                 Dyson_hf(nt, ntau, nao_, k) { 
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

void hfefmolSimulationBase::save_base(h5::File &file, const std::string &path) const {
  h5e::dump(file, path + "/energy/nuclear", enuc_);
  h5e::dump(file, path + "/energy/kinetic", eKin_);
  h5e::dump(file, path + "/energy/potential", ePot_);
  
  h5e::dump(file, path + "/solve/params/nao", nao_);
  h5e::dump(file, path + "/solve/params/dt", dt_);
  h5e::dump(file, path + "/solve/params/nt", nt_);
  h5e::dump(file, path + "/solve/params/ntau", ntau_);
  h5e::dump(file, path + "/solve/params/k", k_);

  h5e::dump(file, path + "/solve/params/MatMax", MatMax_);
  h5e::dump(file, path + "/solve/params/MatTol", MatTol_);
  h5e::dump(file, path + "/solve/params/BootMax", BootMax_);
  h5e::dump(file, path + "/solve/params/BootTol", BootTol_);
  h5e::dump(file, path + "/solve/params/CorrSteps", CorrSteps_);
  
  h5e::dump(file, path + "/solve/params/boot_conv", bootstrap_converged);
}


void hfefmolSimulationBase::run(){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
  
  // Run the gfmol Matsubara solver
  do_mat();

  // Calculate the free GF
  free_gf();

  // Take the coefficients from gfmol and transform them into Legendre mesh
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

    do_energy();
  }
  else {
    std::cout << "Bootstrapping did not converge after " << BootMax_ << " iterations!" << std::endl;
  }
}


template class hfefmolSimulation<gfmol::ChebyshevRepr>;

} // namespace NEdyson
