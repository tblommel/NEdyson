//
// Created by tblommel on 6/8/20
//

#ifndef NE_SIMULATION_DECOMP_IMPL_H
#define NE_SIMULATION_DECOMP_IMPL_H

#include <fstream>
#include <boost/math/special_functions/chebyshev.hpp>

namespace NEdyson {

template <typename Repr>
DecompSimulation<Repr>::DecompSimulation(const gfmol::HartreeFock &hf,
                             const gfmol::RepresentationBase<Repr> &frepr,
                             const gfmol::RepresentationBase<Repr> &brepr,
                             int nt, int ntau, int k, double dt,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping,
                             double decomp_prec, bool hfbool, bool boolPumpProbe, 
                             std::string PumpProbeInp, std::string MolInp,
                             double lPumpProbe, double nPumpProbe) : 
                                 SimulationBase(hf, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, hfbool, boolPumpProbe, PumpProbeInp, MolInp, lPumpProbe, nPumpProbe),
                                 hmf(nt+1, nao_, nao_), 
                                 h0(hf.hcore()), 
                                 rho(nao_,nao_)
{
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::DecompSimulation<Repr> >(new gfmol::DecompSimulation<Repr>(hf, frepr, brepr, mode, 0., decomp_prec, hfbool));
      beta_ = p_MatSim_->frepr().beta();
      p_NEgf2_ = std::unique_ptr<molGF2SolverDecomp>(new molGF2SolverDecomp(p_MatSim_->Vija(), p_MatSim_->Viaj()));
  }

  Sigma = GREEN(nt, ntau, nao_, -1);
  G = GREEN(nt, ntau, nao_, -1);
}


template <typename Repr>
void DecompSimulation<Repr>::free_gf() {
  Dyson.G0_from_h0(G, p_MatSim_->mu(), p_MatSim_->fock(), p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void DecompSimulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}

template <typename Repr>
void DecompSimulation<Repr>::Ed_contractions(int tstp) {
  int nao = hmf.shape()[2];
  for(int d = 0; d < 3; d++) {
    ZMatrixMap(hmf.data() + tstp*nao*nao, nao, nao) += (Efield_(d, tstp) + efield_(d, tstp) + dfield_(d, tstp)) * DMatrixMap(dipole_.data() + d*nao*nao, nao, nao);
  }
}

template <typename Repr>
void DecompSimulation<Repr>::do_boot() {
  for(int iter = 0; iter <= BootMax_; iter++){
    double err = 0;

    // Update mean field & self energy
    for(int tstp = 0; tstp <= k_; tstp++){
      G.get_dm(tstp, rho);
      ZMatrixMap(hmf.data() + tstp*nao_*nao_, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
      p_NEgf2_->solve_HF(tstp, hmf, rho);
      if (!hfbool_) p_NEgf2_->solve(tstp, Sigma, G);
      if(boolPumpProbe_) {
        Dyson.dipole_field(tstp, dfield_, G, G, dipole_, lPumpProbe_, nPumpProbe_, dt_);
        Ed_contractions(tstp);
      }
    }

    // Solve G Equation of Motion
    err = Dyson.dyson_start(G, Sigma, hmf, p_MatSim_->mu(), beta_, dt_);

    std::cout<<"Bootstrapping iteration : "<<iter<<" | Error = "<<err<<std::endl;
    if(err<BootTol_){
      bootstrap_converged = true;
      break;
    }
  }
}


template <typename Repr>
void DecompSimulation<Repr>::do_tstp(int tstp) {
  // Predictor
  Dyson.Extrapolate(tstp, G);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    G.get_dm(tstp, rho);
    ZMatrixMap(hmf.data() + tstp*nao_*nao_, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
    p_NEgf2_->solve_HF(tstp, hmf, rho);
    if(!hfbool_) p_NEgf2_->solve(tstp, Sigma, G);

    if(boolPumpProbe_) {
      Dyson.dipole_field(tstp, dfield_, G, G, dipole_, lPumpProbe_, nPumpProbe_, dt_);
      Ed_contractions(tstp);
    }

    Dyson.dyson_step(tstp, G, Sigma, hmf, p_MatSim_->mu(), beta_, dt_);
  }
}


template <typename Repr>
void DecompSimulation<Repr>::save(h5::File &file, const std::string &path) {
  save_base(file, path);

  G.print_to_file(file, path + "/G");
  Sigma.print_to_file(file, path + "/Sigma");
  
  h5e::dump(file, path + "/params/beta", beta_);
  h5e::dump(file, path + "/mu", p_MatSim_->mu());
  h5e::dump(file, path + "/hmf", hmf);
  h5e::dump(file, path + "/params/h0", h0);
  h5e::dump(file, path + "/params/filling", p_MatSim_->filling());
  h5e::dump(file, path + "/tti", 0);
  h5e::dump(file, path + "/energy/EkinM", p_MatSim_->ehf() + p_MatSim_->ekin());
  h5e::dump(file, path + "/energy/EpotM", p_MatSim_->epot());
}


template <typename Repr>
void DecompSimulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}


template <typename Repr>
void DecompSimulation<Repr>::do_energy() {
  int nao2 = nao_*nao_;
  for(int t=0; t<=nt_; t++) {
    G.get_dm(t,rho);

    auto rhomat = ZMatrixMap(rho.data(), nao_, nao_);
    auto hmfmat = ZMatrixMap(hmf.data() + t*nao2, nao_, nao_);
    auto h0mat = DMatrixConstMap(h0.data(),nao_,nao_);

    eKin_(t) = rhomat.cwiseProduct((hmfmat+h0mat).transpose()).sum().real();
    ePot_(t) = 2*Dyson.energy_conv(t, Sigma, G, beta_, dt_);
  }
}

template <>
inline void DecompSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
  int ntau = ntau_;
  int nL = p_MatSim_->frepr().nl();
  DMatrix Trans(ntau+1, nL);
  for(int t=0; t<=ntau; t++){
    double x = Dyson.Convolution().collocation().x_i()(t);
    for(int l=0; l<nL; l++){
      Trans(t,l) = boost::math::chebyshev_t(l,x);
    }
  }
  int nao2 = nao_*nao_;
  for(int i=0; i<nao2; i++){
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(G.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data()+i, nL, Eigen::InnerStride<>(nao2));

    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sigma.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->sigmal().data()+i, nL, Eigen::InnerStride<>(nao2));
  }
}


template <>
inline void DecompSimulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}

///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////


template <typename Repr>
tti_DecompSimulation<Repr>::tti_DecompSimulation(const gfmol::HartreeFock &hf,
                             const gfmol::RepresentationBase<Repr> &frepr,
                             const gfmol::RepresentationBase<Repr> &brepr,
                             int nt, int ntau, int k, double dt,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping,
                             double decomp_prec, bool hfbool) : 
                                 SimulationBase(hf, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, hfbool, false, "", "", 0, 0),
                                 h0(hf.hcore())
{
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::DecompSimulation<Repr> >(new gfmol::DecompSimulation<Repr>(hf, frepr, brepr, mode, 0., decomp_prec, hfbool));
      beta_ = p_MatSim_->frepr().beta();
      p_NEgf2_ = std::unique_ptr<tti_molGF2SolverDecomp>(new tti_molGF2SolverDecomp(p_MatSim_->Vija(), p_MatSim_->Viaj()));
  }

  Sigma = TTI_GREEN(nt, ntau, nao_, -1);
  G = TTI_GREEN(nt, ntau, nao_, -1);
}


template <typename Repr>
void tti_DecompSimulation<Repr>::free_gf() {
  Dyson.G0_from_h0(G, p_MatSim_->mu(), p_MatSim_->fock(), p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void tti_DecompSimulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}

template <typename Repr>
void tti_DecompSimulation<Repr>::Ed_contractions(int tstp) {
  int a = 0;
}



template <typename Repr>
void tti_DecompSimulation<Repr>::do_boot() {
  for(int iter = 0; iter <= BootMax_; iter++){
    double err = 0;

    // Update mean field & self energy
    for(int tstp = 0; tstp <= k_; tstp++){
      if(!hfbool_) p_NEgf2_->solve(tstp, Sigma, G);
    }

    // Solve G Equation of Motion
    err = Dyson.dyson_start(G, Sigma, p_MatSim_->fock(), p_MatSim_->mu(), beta_, dt_);

    std::cout<<"Bootstrapping iteration : "<<iter<<" | Error = "<<err<<std::endl;
    if(err<BootTol_){
      bootstrap_converged = true;
      break;
    }
  }
}

template <typename Repr>
void tti_DecompSimulation<Repr>::do_energy() {
  int nao2 = nao_*nao_;
  ZTensor<2> rho(nao_, nao_);

  for(int t=0; t<=nt_; t++) {
    G.get_dm(t,rho);

    auto rhomat = ZMatrixMap(rho.data(), nao_, nao_);
    auto hmfmat = DMatrixConstMap(p_MatSim_->fock().data(), nao_, nao_);
    auto h0mat = DMatrixConstMap(h0.data(),nao_,nao_);

    eKin_(t) = rhomat.cwiseProduct((hmfmat+h0mat).transpose()).sum().real();
    ePot_(t) = 2*Dyson.energy_conv(t, Sigma, G, beta_, dt_);
  }
}

template <typename Repr>
void tti_DecompSimulation<Repr>::do_tstp(int tstp) {
  // Predictor
  Dyson.Extrapolate(tstp, G);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    if(!hfbool_) p_NEgf2_->solve(tstp, Sigma, G);
    Dyson.dyson_step(tstp, G, Sigma, p_MatSim_->fock(), p_MatSim_->mu(), beta_, dt_);
  }
}


template <typename Repr>
void tti_DecompSimulation<Repr>::save(h5::File &file, const std::string &path) {
  save_base(file, path);

  G.print_to_file(file, path + "/G");
  Sigma.print_to_file(file, path + "/Sigma");
  
  h5e::dump(file, path + "/params/beta", beta_);
  h5e::dump(file, path + "/mu", p_MatSim_->mu());
  h5e::dump(file, path + "/hmf", p_MatSim_->fock());
  h5e::dump(file, path + "/params/h0", h0);
  h5e::dump(file, path + "/params/filling", p_MatSim_->filling());
  h5e::dump(file, path + "/tti", 1);
  h5e::dump(file, path + "/energy/EkinM", p_MatSim_->ehf() + p_MatSim_->ekin());
  h5e::dump(file, path + "/energy/EpotM", p_MatSim_->epot());
}


template <typename Repr>
void tti_DecompSimulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}


template <>
inline void tti_DecompSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
  int ntau = ntau_;
  int nL = p_MatSim_->frepr().nl();
  DMatrix Trans(ntau+1, nL);
  for(int t=0; t<=ntau; t++){
    double x = Dyson.Convolution().collocation().x_i()(t);
    for(int l=0; l<nL; l++){
      Trans(t,l) = boost::math::chebyshev_t(l,x);
    }
  }
  int nao2 = nao_*nao_;
  for(int i=0; i<nao2; i++){
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(G.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data()+i, nL, Eigen::InnerStride<>(nao2));

    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sigma.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->sigmal().data()+i, nL, Eigen::InnerStride<>(nao2));
  }
}


template <>
inline void tti_DecompSimulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}

} // namespace NEdyson
#endif // header guard
