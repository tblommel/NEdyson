//
// Created by tblommel on 06/08/2020
//

#ifndef HFEF_NE_SIMULATION_BASE_H
#define HFEF_NE_SIMULATION_BASE_H

#include "gfmol/sim.h"
#include "dyson_hf.h"
#include "greens.h"
#include "integration.h"
#include "molHFSolver.h"
#include "utils.h"
#include <chrono>

namespace NEdyson{

// base class
class hfefmolSimulationBase {
protected:
  double enuc_;
  int nao_;

  double dt_;
  int nt_;
  int ntau_;
  int k_;

  int MatMax_;
  double MatTol_;
  int BootMax_;
  double BootTol_;
  int CorrSteps_;

  bool bootstrap_converged;

  DTensor<1> eKin_;
  DTensor<1> ePot_;
  DTensor<1> tot_time;
  DTensor<1> dys_time;
  DTensor<1> gf2_time;

  dyson_hf Dyson_hf;

  void save_base(h5::File &file, const std::string &path) const;
  
public:
  // Base class constructor
  explicit hfefmolSimulationBase(const gfmol::HartreeFock &hf, int nt, int ntau, int k, double dt, int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps);
  
  // Calculate free GF
  virtual void free_gf() = 0;
  
  // Perform a timestep iteration
  virtual void do_tstp(int tstp) = 0;
  
  // Perform a bootstrap iteration
  virtual void do_boot() = 0;

  // Run the Matsubara solver
  virtual void do_mat() = 0;

  // Coefficients to equidistant mesh
  virtual inline void L_to_Tau() = 0;

  // Calculate energy
  virtual void do_energy() = 0;

  // Output Results to HDF5
  virtual void save(h5::File &file, const std::string &path) = 0;

  // Load from HDF5
  virtual void load(const h5::File &file, const std::string &path) = 0;

  // Run a simulation
  void run();

}; // class SimulationBase


// Vanilla Simulation without decomp and spin
template <typename Repr>
class hfefmolSimulation : public hfefmolSimulationBase {
public:
  std::unique_ptr<gfmol::Simulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<molGF2Solver> p_NEgf2_;

  ZTensor<3> hmf;
  const DTensor<2> &h0;
  ZTensor<2> rho;
  GREEN G; 

  // Construct sim
  hfefmolSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  inline void L_to_Tau() override;
  
  void do_energy() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

}; // class Simulation

} // namespace NEdyson

#include "hf_efield_molNEsim_impl.h"

namespace NEdyson{

extern template class hfefmolSimulation<gfmol::ChebyshevRepr>;
}


#endif // header guard
