//
// Created by tblommel on 06/08/2020
//

#ifndef NE_SIMULATION_BASE_H
#define NE_SIMULATION_BASE_H

#include "gfmol/sim.h"
#include "dyson.h"
#include "greens.h"
#include "spectral.h"
#include "integration.h"
#include "molHFSolver.h"
#include "mol2bSolver.h"
#include "utils.h"
#include <chrono>

namespace NEdyson{

// base class
class SimulationBase {
protected:
  double enuc_;
  int nao_;

  double dt_;
  int nt_;
  int ntau_;
  int k_;
  int nw_;
  double wmax_;

  int MatMax_;
  double MatTol_;
  int BootMax_;
  double BootTol_;
  int CorrSteps_;

  bool bootstrap_converged;

  INTEG I;

  DTensor<1> eKin_;
  DTensor<1> ePot_;
  
public:
  // Base class constructor
  explicit SimulationBase(const gfmol::HartreeFock &hf, int nt, int ntau, int k, double dt, int nw, double wmax, int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps);
  
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

  // Calculate the spectral function
  virtual void do_spectral() = 0;

  // Output Results to HDF5
  virtual void save(h5::File &file, const std::string &path) = 0;

  // Load from HDF5
  virtual void load(const h5::File &file, const std::string &path) = 0;

  // Run a simulation
  void run();

}; // class SimulationBase


// Vanilla Simulation without decomp and spin
template <typename Repr>
class Simulation : public SimulationBase {
  std::unique_ptr<gfmol::Simulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<molGF2Solver> p_NEgf2_;

  ZTensor<3> hmf;
  const DTensor<2> &h0;
  ZTensor<2> rho;
  GREEN Sigma;
  GREEN G; 
  SPECT A;

public:
  // Construct sim
  Simulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt, int nw, double wmax,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_spectral() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

}; // class Simulation


// Spin Restricted with Decomp
template <typename Repr>
class DecompSimulation : public SimulationBase {
  std::unique_ptr<gfmol::DecompSimulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<molGF2SolverDecomp> p_NEgf2_;

  ZTensor<3> hmf;
  const DTensor<2> &h0;
  ZTensor<2> rho;
  GREEN Sigma;
  GREEN G; 
  SPECT A;

public:
  // Construct sim
  DecompSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt, int nw, double wmax,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0,
             double decomp_prec = 1e-8);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_spectral() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

}; // class DecompSimulation


// Simulation without decomp and with spin
template <typename Repr>
class SpinSimulation : public SimulationBase {
  std::unique_ptr<gfmol::SpinSimulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<molGF2SolverSpin> p_NEgf2_;

  ZTensor<4> hmf;
  const DTensor<2> &h0;
  ZTensor<3> rho;

  GREEN Sup;
  GREEN Sdown;
  GREEN Gup;
  GREEN Gdown;
  
  std::vector<std::reference_wrapper<GREEN>> G;
  std::vector<std::reference_wrapper<GREEN>> Sigma;
  
  SPECT A;

public:
  // Construct sim
  SpinSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt, int nw, double wmax,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_spectral() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

}; // class SpinSimulation


// Simulation with decomp and with spin
template <typename Repr>
class DecompSpinSimulation : public SimulationBase {
  std::unique_ptr<gfmol::DecompSpinSimulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<molGF2SolverSpinDecomp> p_NEgf2_;

  ZTensor<4> hmf;
  const DTensor<2> &h0;
  ZTensor<3> rho;

  GREEN Sup;
  GREEN Sdown;
  GREEN Gup;
  GREEN Gdown;
  
  std::vector<std::reference_wrapper<GREEN>> G;
  std::vector<std::reference_wrapper<GREEN>> Sigma;
  
  SPECT A;

public:
  // Construct sim
  DecompSpinSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt, int nw, double wmax,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0, double decomp_prec = 1e-8);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_spectral() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

}; // class DecompSpinSimulation

// Factory for simulations
template <typename Repr>
std::unique_ptr<SimulationBase> make_simulation(bool unrestricted,
                                                bool decomposed,
                                                const gfmol::HartreeFock &hf,
                                                const gfmol::RepresentationBase<Repr> &frepr,
                                                const gfmol::RepresentationBase<Repr> &brepr,
                                                gfmol::Mode mode,
                                                int nt, int ntau, int k, double dt, int nw, double wmax,
                                                int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                                                double damping = 0, 
                                                double decomp_prec = 1e-7)
{
  if(unrestricted) {
    if(decomposed) {
      return std::unique_ptr<DecompSpinSimulation<Repr>>(
        new DecompSpinSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, decomp_prec));
    }
    else{
      return std::unique_ptr<SpinSimulation<Repr>>(
        new SpinSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping));
    }
  }
  else {
    if(decomposed) {
      return std::unique_ptr<DecompSimulation<Repr>>(
        new DecompSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, decomp_prec));
    }
    else {
      return std::unique_ptr<Simulation<Repr>>(
        new Simulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping));
    }
  }
  return nullptr;
}

} // namespace NEdyson

#include "molNEsim_impl.h"
#include "molNEsim_decomp_impl.h"
#include "molNEsim_spin_impl.h"
#include "molNEsim_decomp_spin_impl.h"

namespace NEdyson{

extern template class Simulation<gfmol::ChebyshevRepr>;
extern template class Simulation<gfmol::IntermediateRepr>;
extern template class DecompSimulation<gfmol::ChebyshevRepr>;
extern template class DecompSimulation<gfmol::IntermediateRepr>;
extern template class SpinSimulation<gfmol::ChebyshevRepr>;
extern template class SpinSimulation<gfmol::IntermediateRepr>;
extern template class DecompSpinSimulation<gfmol::ChebyshevRepr>;
extern template class DecompSpinSimulation<gfmol::IntermediateRepr>;

}

////////////////////////////////////////////////////////////////////////////////
///////////////////  TTI        ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




template <typename Repr>
class tti_Simulation : public SimulationBase {
  std::unique_ptr<gfmol::Simulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<tti_molGF2Solver> p_NEgf2_;

  ZTensor<2> hmf;
  const DTensor<2> &h0;
  ZTensor<2> rho;
  TTI_GREEN Sigma;
  TTI_GREEN G; 
  TTI_SPECT A;

public:
  // Construct sim
  tti_Simulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt, int nw, double wmax,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_spectral() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

}; // class Simulation


#endif // header guard
