//
// Created by tblommel on 06/08/2020
//

#ifndef NE_SIMULATION_BASE_H
#define NE_SIMULATION_BASE_H

#include "gfmol/sim.h"
#include "dyson.h"
#include "greens.h"
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

  int MatMax_;
  double MatTol_;
  int BootMax_;
  double BootTol_;
  int CorrSteps_;

  bool bootstrap_converged;
  bool hfbool_;

  bool boolPumpProbe_;
  double lPumpProbe_;
  double nPumpProbe_;
  DTensor<3> dipole_;
  DTensor<2> Efield_;
  DTensor<2> efield_;

  DTensor<1> eKin_;
  DTensor<1> ePot_;
  DTensor<1> tot_time;
  DTensor<1> dys_time;
  DTensor<1> gf2_time;

  dyson Dyson;

  void save_base(h5::File &file, const std::string &path) const;
  
public:
  // Base class constructor
  explicit SimulationBase(const gfmol::HartreeFock &hf, int nt, int ntau, int k, double dt, int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps, bool hfbool, bool boolPumpProbe, std::string PumpProbeInput, std::string MolInput, double lPumpProbe, double nPumpProbe);
  
  // Calculate free GF
  virtual void free_gf() = 0;
  
  // Perform a timestep iteration
  virtual void do_tstp(int tstp) = 0;
  
  // Perform a bootstrap iteration
  virtual void do_boot() = 0;

  // Run the Matsubara solver
  virtual void do_mat() = 0;

  // Coefficients to Legendre mesh
  virtual inline void L_to_Tau() = 0;

  // Calculate energy
  virtual void do_energy() = 0;

  // Output Results to HDF5
  virtual void save(h5::File &file, const std::string &path) = 0;

  // Load from HDF5
  virtual void load(const h5::File &file, const std::string &path) = 0;

  // Do Electric-field, dipole contractions
  virtual void Ed_contractions(int tstp) = 0;

  // Run a simulation
  void run();

}; // class SimulationBase


// Vanilla Simulation without decomp and spin
template <typename Repr>
class Simulation : public SimulationBase {
public:
  std::unique_ptr<gfmol::Simulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<molGF2Solver> p_NEgf2_;

  ZTensor<3> hmf;
  const DTensor<2> &h0;
  ZTensor<2> rho;
  GREEN Sigma;
  GREEN G; 
  ZTensor<2> dfield_;

  // Construct sim
  Simulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0, bool hfbool = false, bool boolPumpProbe = false, 
             std::string PumpProbeInp = "", std::string MolInp = "", 
             double lPumpProbe = 0, double nPumpProbe = 0);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  inline void L_to_Tau() override;

  void do_energy() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

  void Ed_contractions(int tstp) override;

}; // class Simulation


// Spin Restricted with Decomp
template <typename Repr>
class DecompSimulation : public SimulationBase {
public:
  std::unique_ptr<gfmol::DecompSimulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<molGF2SolverDecomp> p_NEgf2_;

  ZTensor<3> hmf;
  const DTensor<2> &h0;
  ZTensor<2> rho;
  GREEN Sigma;
  GREEN G; 
  ZTensor<2> dfield_;

  // Construct sim
  DecompSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0,
             double decomp_prec = 1e-8, bool hfbool = false, bool boolPumpProbe = false, 
             std::string PumpProbeInp = "", std::string MolInp = "",
             double lPumpProbe = 0, double nPumpProbe = 0);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_energy() override;
  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

  void Ed_contractions(int tstp) override;
}; // class DecompSimulation


// Simulation without decomp and with spin
template <typename Repr>
class SpinSimulation : public SimulationBase {
public:
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
  ZTensor<2> dfieldu_;
  ZTensor<2> dfieldd_;
  
  std::vector<std::reference_wrapper<GREEN>> G;
  std::vector<std::reference_wrapper<GREEN>> Sigma;
  
  // Construct sim
  SpinSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0, bool hfbool = false, bool boolPumpProbe = false, 
             std::string PumpProbeInp = "", std::string MolInp = "",
             double lPumpProbe = 0, double nPumpProbe = 0);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_energy() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

  void Ed_contractions(int tstp) override;
}; // class SpinSimulation


// Simulation with decomp and with spin
template <typename Repr>
class DecompSpinSimulation : public SimulationBase {
public:
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
  ZTensor<2> dfieldu_;
  ZTensor<2> dfieldd_;
  
  std::vector<std::reference_wrapper<GREEN>> G;
  std::vector<std::reference_wrapper<GREEN>> Sigma;
  
  // Construct sim
  DecompSpinSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0, 
             double decomp_prec = 1e-8, bool hfbool = false, bool boolPumpProbe = false, 
             std::string PumpProbeInp = "", std::string MolInp = "",
             double lPumpProbe = 0, double nPumpProbe = 0);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_energy() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

  void Ed_contractions(int tstp) override;
}; // class DecompSpinSimulation

////////////////////////////////////////////////////////////////////////////////
///////////////////  TTI        ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




template <typename Repr>
class tti_Simulation : public SimulationBase {
public:
  std::unique_ptr<gfmol::Simulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<tti_molGF2Solver> p_NEgf2_;

  const DTensor<2> &h0;
  TTI_GREEN Sigma;
  TTI_GREEN G; 
  // Construct sim
  tti_Simulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0, bool hfbool = false);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_energy() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

  void Ed_contractions(int tstp) override;
}; // class Simulation

// Spin Restricted with Decomp
template <typename Repr>
class tti_DecompSimulation : public SimulationBase {
public:
  std::unique_ptr<gfmol::DecompSimulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<tti_molGF2SolverDecomp> p_NEgf2_;

  const DTensor<2> &h0;
  TTI_GREEN Sigma;
  TTI_GREEN G; 

  // Construct sim
  tti_DecompSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0,
             double decomp_prec = 1e-8, bool hfbool = false);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_energy() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

  void Ed_contractions(int tstp) override;
}; // class DecompSimulation

// Simulation without decomp and with spin
template <typename Repr>
class tti_SpinSimulation : public SimulationBase {
public:
  std::unique_ptr<gfmol::SpinSimulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<tti_molGF2SolverSpin> p_NEgf2_;

  const DTensor<2> &h0;

  TTI_GREEN Sup;
  TTI_GREEN Sdown;
  TTI_GREEN Gup;
  TTI_GREEN Gdown;
  
  std::vector<std::reference_wrapper<TTI_GREEN>> G;
  std::vector<std::reference_wrapper<TTI_GREEN>> Sigma;
  
  // Construct sim
  tti_SpinSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0, bool hfbool = false);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_energy() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

  void Ed_contractions(int tstp) override;
}; // class SpinSimulation


// Simulation with decomp and with spin
template <typename Repr>
class tti_DecompSpinSimulation : public SimulationBase {
public:
  std::unique_ptr<gfmol::DecompSpinSimulation<Repr>> p_MatSim_;
  double beta_;
  double dtau_;
  std::unique_ptr<tti_molGF2SolverSpinDecomp> p_NEgf2_;

  const DTensor<2> &h0;

  TTI_GREEN Sup;
  TTI_GREEN Sdown;
  TTI_GREEN Gup;
  TTI_GREEN Gdown;
  
  std::vector<std::reference_wrapper<TTI_GREEN>> G;
  std::vector<std::reference_wrapper<TTI_GREEN>> Sigma;
  
  // Construct sim
  tti_DecompSpinSimulation(const gfmol::HartreeFock &hf,
             const gfmol::RepresentationBase<Repr> &frepr,
             const gfmol::RepresentationBase<Repr> &brepr,
             int nt, int ntau, int k, double dt,
             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
             gfmol::Mode mode = gfmol::Mode::GF2,
             double damping = 0, double decomp_prec = 1e-8, bool hfbool = false);

  void free_gf() override;

  void do_tstp(int tstp) override;

  void do_boot() override;

  void do_mat() override;

  void do_energy() override;

  inline void L_to_Tau() override;

  void save(h5::File &file, const std::string &path) override;

  void load(const h5::File &file, const std::string &path) override;

  void Ed_contractions(int tstp) override;
}; // class DecompSpinSimulation

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


// Factory for simulations
template <typename Repr>
std::unique_ptr<SimulationBase> make_simulation(bool tti,
                                                bool unrestricted,
                                                bool decomposed,
                                                const gfmol::HartreeFock &hf,
                                                const gfmol::RepresentationBase<Repr> &frepr,
                                                const gfmol::RepresentationBase<Repr> &brepr,
                                                gfmol::Mode mode,
                                                int nt, int ntau, int k, double dt, 
                                                int MatMax, double MatTol, int BootMax, 
                                                double BootTol, int CorrSteps, double damping = 0, 
                                                double decomp_prec = 1e-7, bool hfbool = false, 
                                                bool boolPumpProbe = false, std::string PumpProbeInp= "", 
                                                std::string MolInp = "", double lPumpProbe = 0., double nPumpProbe = 0.)
{
  if(tti) {
    if(unrestricted) {
      if(decomposed) {
        return std::unique_ptr<tti_DecompSpinSimulation<Repr>>(
          new tti_DecompSpinSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, decomp_prec, hfbool));
      }
      else{
        return std::unique_ptr<tti_SpinSimulation<Repr>>(
          new tti_SpinSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, hfbool));
      }
    }
    else {
      if(decomposed) {
        return std::unique_ptr<tti_DecompSimulation<Repr>>(
          new tti_DecompSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, decomp_prec, hfbool));
      }
      else {
        return std::unique_ptr<tti_Simulation<Repr>>(
          new tti_Simulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, hfbool));
      }
    }
  }
  else {
    if(unrestricted) {
      if(decomposed) {
        return std::unique_ptr<DecompSpinSimulation<Repr>>(
          new DecompSpinSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, decomp_prec, hfbool, boolPumpProbe, PumpProbeInp, MolInp, lPumpProbe, nPumpProbe));
      }
      else{
        return std::unique_ptr<SpinSimulation<Repr>>(
          new SpinSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, hfbool, boolPumpProbe, PumpProbeInp, MolInp, lPumpProbe, nPumpProbe));
      }
    }
    else {
      if(decomposed) {
        return std::unique_ptr<DecompSimulation<Repr>>(
          new DecompSimulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, decomp_prec, hfbool, boolPumpProbe, PumpProbeInp, MolInp, lPumpProbe, nPumpProbe));
      }
      else {
        return std::unique_ptr<Simulation<Repr>>(
          new Simulation<Repr>(hf, frepr, brepr, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, mode, damping, hfbool, boolPumpProbe, PumpProbeInp, MolInp, lPumpProbe, nPumpProbe));
      }
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
extern template class tti_Simulation<gfmol::ChebyshevRepr>;
extern template class tti_Simulation<gfmol::IntermediateRepr>;

extern template class DecompSimulation<gfmol::ChebyshevRepr>;
extern template class DecompSimulation<gfmol::IntermediateRepr>;
extern template class tti_DecompSimulation<gfmol::ChebyshevRepr>;
extern template class tti_DecompSimulation<gfmol::IntermediateRepr>;

extern template class SpinSimulation<gfmol::ChebyshevRepr>;
extern template class SpinSimulation<gfmol::IntermediateRepr>;
extern template class tti_SpinSimulation<gfmol::ChebyshevRepr>;
extern template class tti_SpinSimulation<gfmol::IntermediateRepr>;

extern template class DecompSpinSimulation<gfmol::ChebyshevRepr>;
extern template class DecompSpinSimulation<gfmol::IntermediateRepr>;
extern template class tti_DecompSpinSimulation<gfmol::ChebyshevRepr>;
extern template class tti_DecompSpinSimulation<gfmol::IntermediateRepr>;

}


#endif // header guard
