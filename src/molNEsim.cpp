//
// Created by tblommel on 06/09/2020
//

#include "molNEsim.h"

namespace NEdyson {

SimulationBase::SimulationBase(const gfmol::HartreeFock &hf, 
                               const NEdyson::Params &p) 
                               : enuc_(hf.enuc()), 
                                 nao_(hf.nao()), 
                                 eKin_(p.nt+1), 
                                 ePot_(p.nt+1),
                                 Dyson(p.nt, p.ntau, nao_, p.k, p.gfmolmode),
                                 efield_(3, p.nt+1),
                                 Efield_(3, p.nt+1),
                                 dfield_(3, p.nt+1),
                                 dipole_(3, nao_, nao_) {
  nt_ = p.nt;
  ntau_ = p.ntau;
  dt_ = p.dt;
  k_ = p.k;
  bootstrap_converged = false;
  boolPumpProbe_ = p.boolPumpProbe;
  lPumpProbe_ = p.lPumpProbe;
  nPumpProbe_ = p.nPumpProbe;
  mode_ = p.gfmolmode;

  if(boolPumpProbe_) {
    h5e::File ppinp(p.PumpProbe_file);
    Efield_ = h5e::load<DTensor<2>>(ppinp, "/E");
    efield_ = h5e::load<DTensor<2>>(ppinp, "/e");

    h5e::File molinp_file(p.hf_input);
    dipole_ = h5e::load<DTensor<3>>(molinp_file, "/dipole");

    // Bootstrapping is tti so we need to make sure e and E are zero
    // for first k+1 timesteps
    for(int d = 0; d < 3; d++) {
      int erows = efield_.shape()[1];
      int Erows = Efield_.shape()[2];
      double enorm = DRowVectorMap(efield_.data() + d * erows, k_+1).norm();
      double Enorm = DRowVectorMap(Efield_.data() + d * Erows, k_+1).norm();
      if(enorm > 1e-15 || Enorm > 1e-15) {
        std::cout << "applied field is nonzero in the first k timesteps" << std::endl;
        std::exit(0);
      }
    }
  }

  MatMax_ = p.maxiter;
  MatTol_ = p.etol;
  BootMax_ =p.BootMaxIter;
  BootTol_ = p.BootMaxErr;
  CorrSteps_ = p.CorrSteps;
}

void SimulationBase::save_base(h5::File &file, const std::string &path) const {
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
  
  h5e::dump(file, path + "/solve/params/boot_conv", (int)bootstrap_converged);
  if(mode_ == gfmol::Mode::GF2) {
    h5e::dump(file, path + "/solve/params/mode", std::string("GF2"));
  }
  if(mode_ == gfmol::Mode::HF) {
    h5e::dump(file, path + "/solve/params/mode", std::string("HF"));
  }
}

void SimulationBase::save_PP(h5::File &file, const std::string &path) const {
  h5e::dump(file, path + "/l", lPumpProbe_);
  h5e::dump(file, path + "/n", nPumpProbe_);
  h5e::dump(file, path + "/dfield", dfield_);
  h5e::dump(file, path + "/efield", efield_);
  h5e::dump(file, path + "/dt", dt_);
}


void SimulationBase::run(){
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
    std::exit(0);
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
