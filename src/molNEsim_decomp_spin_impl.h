//
// Created by tblommel on 6/8/20
//

#ifndef NE_SIMULATION_DECOMP_SPIN_IMPL_H
#define NE_SIMULATION_DECOMP_SPIN_IMPL_H

#include <fstream>
#include <boost/math/special_functions/chebyshev.hpp>

namespace NEdyson {

template <typename Repr>
DecompSpinSimulation<Repr>::DecompSpinSimulation(const gfmol::HartreeFock &hf,
                             const gfmol::RepresentationBase<Repr> &frepr,
                             const gfmol::RepresentationBase<Repr> &brepr,
                             const Params &p) : 
                                 SimulationBase(hf, p),
                                 hmf(2, p.nt + 1, nao_, nao_), 
                                 h0(hf.hcore()), 
                                 rho(2, nao_, nao_)
{
  p_MatSim_ = std::unique_ptr<gfmol::DecompSpinSimulation<Repr> >(new gfmol::DecompSpinSimulation<Repr>(hf, frepr, brepr, p.gfmolmode, 0.));
  beta_ = p_MatSim_->frepr().beta();

  switch (p.gfmolmode) {
    case gfmol::Mode::HF:
      p_NEgf2_ = std::unique_ptr<molGF2SolverSpinDecomp>(new molGF2SolverSpinDecomp(p_MatSim_->Vija() ));
      Sup = GREEN(0, 0, 1, -1);
      Sdown = GREEN(0, 0, 1, -1);
      break;
    case gfmol::Mode::GF2:
      p_NEgf2_ = std::unique_ptr<molGF2SolverSpinDecomp>(new molGF2SolverSpinDecomp(p_MatSim_->Vija(), dynamic_cast<gfmol::DecompGF2SolverSpin *>(p_MatSim_->p_sigma().get())->Viaj() ));
      Sup = GREEN(p.nt, p.ntau, nao_, -1);
      Sdown = GREEN(p.nt, p.ntau, nao_, -1);
      break;
    default:
      throw std::runtime_error("[Simulation] Unknown mode.");
  }

  Gup = GREEN(p.nt, p.ntau, nao_, -1);
  Gdown = GREEN(p.nt, p.ntau, nao_, -1);

  G = {Gup, Gdown};
  Sigma = {Sup, Sdown};
}


template <typename Repr>
void DecompSpinSimulation<Repr>::free_gf() {
  Dyson.G0_from_h0(Gup, p_MatSim_->mu()[0], p_MatSim_->fock().data(), p_MatSim_->frepr().beta(), dt_);
  Dyson.G0_from_h0(Gdown, p_MatSim_->mu()[1], p_MatSim_->fock().data() + nao_*nao_, p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void DecompSpinSimulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}

template <typename Repr>
void DecompSpinSimulation<Repr>::Ed_contractions(int tstp) {
  int nao = hmf.shape()[3];
  for(int d = 0; d < 3; d++) {
    ZMatrixMap(hmf.data() + tstp*nao*nao, nao, nao) += (Efield_(d, tstp) + efield_(d, tstp) + dfield_(d, tstp)) * DMatrixMap(dipole_.data() + d*nao*nao, nao, nao);
    ZMatrixMap(hmf.data() + (nt_+1)*nao_*nao_ + tstp*nao*nao, nao, nao) += (Efield_(d, tstp) + efield_(d, tstp) + dfield_(d, tstp)) * DMatrixMap(dipole_.data() + d*nao*nao, nao, nao);
  }
}

template <typename Repr>
void DecompSpinSimulation<Repr>::do_boot() {
  int nao2 = nao_*nao_;
  for(int iter = 0; iter <= BootMax_; iter++){
    double err = 0;

    // Update mean field & self energy
    for(int tstp = 0; tstp <= k_; tstp++){
      ZMatrixMap(hmf.data() + tstp*nao2, nao_, nao_) = DMatrixConstMap(p_MatSim_->fock().data(),nao_,nao_);
      ZMatrixMap(hmf.data() + (nt_+1)*nao2 + tstp*nao2, nao_, nao_) = DMatrixConstMap(p_MatSim_->fock().data() + nao2, nao_, nao_);

      if(mode_ == gfmol::Mode::GF2) p_NEgf2_->solve(tstp, Sigma, G);
    }

    // Solve G Equation of Motion
    err = Dyson.dyson_start(Gup, Sup, hmf.data(), p_MatSim_->mu()[0], beta_, dt_);
    err += Dyson.dyson_start(Gdown, Sdown, hmf.data() + (nt_+1)*nao2, p_MatSim_->mu()[1], beta_, dt_);

    std::cout<<"Bootstrapping iteration : "<<iter<<" | Error = "<<err<<std::endl;
    if(err<BootTol_){
      bootstrap_converged = true;
      break;
    }
  }
}


template <typename Repr>
void DecompSpinSimulation<Repr>::do_tstp(int tstp) {
  int nao2 = nao_ * nao_;
  // Predictor
  Dyson.Extrapolate(tstp, Gup);
  Dyson.Extrapolate(tstp, Gdown);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    Gup.get_dm(tstp, rho.data());
    Gdown.get_dm(tstp, rho.data() + nao2);

    ZMatrixMap(hmf.data() + tstp*nao2, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
    ZMatrixMap(hmf.data() + (nt_+1)*nao2 + tstp*nao2, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);

    p_NEgf2_->solve_HF(tstp, hmf, rho);
    if(mode_ == gfmol::Mode::GF2) p_NEgf2_->solve(tstp, Sigma, G);
    if(boolPumpProbe_) {
      Dyson.dipole_field(tstp, dfield_, Gup, Gdown, dipole_, lPumpProbe_, nPumpProbe_, dt_);
      Ed_contractions(tstp);
    }

    Dyson.dyson_step(tstp, Gup, Sup, hmf.data(), p_MatSim_->mu()[0], beta_, dt_);
    Dyson.dyson_step(tstp, Gdown, Sdown, hmf.data() + (nt_+1)*nao2, p_MatSim_->mu()[1], beta_, dt_);
  }
}


template <typename Repr>
void DecompSpinSimulation<Repr>::save(h5::File &file, const std::string &path) {
  save_base(file, path);

  Gup.print_to_file(file, path + "/G/up");
  Sup.print_to_file(file, path + "/Sigma/up");
  Gdown.print_to_file(file, path + "/G/down");
  Sdown.print_to_file(file, path + "/Sigma/down");
  
  h5e::dump(file, path + "/params/beta", beta_);
  h5e::dump(file, path + "/mu", std::vector<double>{p_MatSim_->mu()[0], p_MatSim_->mu()[1]});
  h5e::dump(file, path + "/hmf", hmf);
  h5e::dump(file, path + "/params/h0", h0);
  h5e::dump(file, path + "/params/filling", std::vector<double>{p_MatSim_->filling()[0], p_MatSim_->filling()[1]});
  h5e::dump(file, path + "/tti", 0);
  h5e::dump(file, path + "/energy/EkinM", p_MatSim_->ehf() + p_MatSim_->ekin());
  h5e::dump(file, path + "/energy/EpotM", p_MatSim_->epot());
}

template <typename Repr>
void DecompSpinSimulation<Repr>::do_energy() {
  int nao2 = nao_*nao_;
  for(int t=0; t<=nt_; t++) {
    eKin_(t) = 0;
    ePot_(t) = 0;
    for(int s=0; s<2; s++) {
      G[s].get().get_dm(t,rho.data());

      auto rhomat = ZMatrixMap(rho.data(), nao_, nao_);
      auto hmfmat = ZMatrixMap(hmf.data() + s*(nt_+1)*nao2 + t*nao2, nao_, nao_);
      auto h0mat = DMatrixConstMap(h0.data(),nao_,nao_);

      eKin_(t) += 0.5*rhomat.cwiseProduct((hmfmat+h0mat).transpose()).sum().real();
      ePot_(t) += Dyson.energy_conv(t, Sigma[s], G[s], beta_, dt_);
    }
  }
}

template <typename Repr>
void DecompSpinSimulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}


template <>
inline void DecompSpinSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
  int nL = p_MatSim_->frepr().nl();
  DMatrix Trans(ntau_ + 1, nL);
  for(int t=0; t<=ntau_; t++){
    double x = Dyson.Convolution().collocation().x_i()(t);
    for(int l=0; l<nL; l++){
      Trans(t,l) = boost::math::chebyshev_t(l,x);
    }
  }
  int nao2 = nao_*nao_;
  for(int i=0; i<nao2; i++){
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Gup.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data() + i, nL, Eigen::InnerStride<>(nao2));

    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Gdown.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data() + nL*nao_*nao_ + i, nL, Eigen::InnerStride<>(nao2));

    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sup.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->sigmal().data()+i, nL, Eigen::InnerStride<>(nao2));

    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sdown.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->sigmal().data() + nL*nao_*nao_ + i, nL, Eigen::InnerStride<>(nao2));
  }
}


template <>
inline void DecompSpinSimulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////


template <typename Repr>
tti_DecompSpinSimulation<Repr>::tti_DecompSpinSimulation(const gfmol::HartreeFock &hf,
                             const gfmol::RepresentationBase<Repr> &frepr,
                             const gfmol::RepresentationBase<Repr> &brepr,
                             const Params &p) : 
                                 SimulationBase(hf, p),
                                 h0(hf.hcore())
{
  p_MatSim_ = std::unique_ptr<gfmol::DecompSpinSimulation<Repr> >(new gfmol::DecompSpinSimulation<Repr>(hf, frepr, brepr, p.gfmolmode, 0.));
  beta_ = p_MatSim_->frepr().beta();

  switch (p.gfmolmode) {
    case gfmol::Mode::HF:
      p_NEgf2_ = std::unique_ptr<tti_molGF2SolverSpinDecomp>(new tti_molGF2SolverSpinDecomp(p_MatSim_->Vija() ));
      Sup = TTI_GREEN(0, 0, 1, -1);
      Sdown = TTI_GREEN(0, 0, 1, -1);
      break;
    case gfmol::Mode::GF2:
      p_NEgf2_ = std::unique_ptr<tti_molGF2SolverSpinDecomp>(new tti_molGF2SolverSpinDecomp(p_MatSim_->Vija(), dynamic_cast<gfmol::DecompGF2SolverSpin *>(p_MatSim_->p_sigma().get())->Viaj() ));
      Sup = TTI_GREEN(p.nt, p.ntau, nao_, -1);
      Sdown = TTI_GREEN(p.nt, p.ntau, nao_, -1);
      break;
    default:
      throw std::runtime_error("[Simulation] Unknown mode.");
  }

  Gup = TTI_GREEN(p.nt, p.ntau, nao_, -1);
  Gdown = TTI_GREEN(p.nt, p.ntau, nao_, -1);

  G = {Gup, Gdown};
  Sigma = {Sup, Sdown};
}


template <typename Repr>
void tti_DecompSpinSimulation<Repr>::free_gf() {
  Dyson.G0_from_h0(Gup, p_MatSim_->mu()[0], p_MatSim_->fock().data(), p_MatSim_->frepr().beta(), dt_);
  Dyson.G0_from_h0(Gdown, p_MatSim_->mu()[1], p_MatSim_->fock().data() + nao_*nao_, p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void tti_DecompSpinSimulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}


template <typename Repr>
void tti_DecompSpinSimulation<Repr>::Ed_contractions(int tstp) {
  int a = 0;
}



template <typename Repr>
void tti_DecompSpinSimulation<Repr>::do_boot() {
  int nao2 = nao_*nao_;
  for(int iter = 0; iter <= BootMax_; iter++){
    double err = 0;

    // Update mean field & self energy
    for(int tstp = 0; tstp <= k_; tstp++){
      if(mode_ == gfmol::Mode::GF2) p_NEgf2_->solve(tstp, Sigma, G);
    }

    // Solve G Equation of Motion
    err = Dyson.dyson_start(Gup, Sup, p_MatSim_->fock().data(), p_MatSim_->mu()[0], beta_, dt_);
    err += Dyson.dyson_start(Gdown, Sdown, p_MatSim_->fock().data() + nao2, p_MatSim_->mu()[1], beta_, dt_);

    std::cout<<"Bootstrapping iteration : "<<iter<<" | Error = "<<err<<std::endl;
    if(err<BootTol_){
      bootstrap_converged = true;
      break;
    }
  }
}


template <typename Repr>
void tti_DecompSpinSimulation<Repr>::do_tstp(int tstp) {
  int nao2 = nao_ * nao_;
  // Predictor
  Dyson.Extrapolate(tstp, Gup);
  Dyson.Extrapolate(tstp, Gdown);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    if(mode_ == gfmol::Mode::GF2) p_NEgf2_->solve(tstp, Sigma, G);

    Dyson.dyson_step(tstp, Gup, Sup, p_MatSim_->fock().data(), p_MatSim_->mu()[0], beta_, dt_);
    Dyson.dyson_step(tstp, Gdown, Sdown, p_MatSim_->fock().data() + nao2, p_MatSim_->mu()[1], beta_, dt_);
  }
}


template <typename Repr>
void tti_DecompSpinSimulation<Repr>::save(h5::File &file, const std::string &path) {
  save_base(file, path);

  Gup.print_to_file(file, path + "/G/up");
  Sup.print_to_file(file, path + "/Sigma/up");
  Gdown.print_to_file(file, path + "/G/down");
  Sdown.print_to_file(file, path + "/Sigma/down");
  
  h5e::dump(file, path + "/params/beta", beta_);
  h5e::dump(file, path + "/mu", std::vector<double>{p_MatSim_->mu()[0], p_MatSim_->mu()[1]});
  h5e::dump(file, path + "/hmf", p_MatSim_->fock());
  h5e::dump(file, path + "/params/h0", h0);
  h5e::dump(file, path + "/params/filling", std::vector<double>{p_MatSim_->filling()[0], p_MatSim_->filling()[1]});
  h5e::dump(file, path + "/tti", 1);
  h5e::dump(file, path + "/energy/EkinM", p_MatSim_->ehf() + p_MatSim_->ekin());
  h5e::dump(file, path + "/energy/EpotM", p_MatSim_->epot());
}

template <typename Repr>
void tti_DecompSpinSimulation<Repr>::do_energy() {
  int nao2 = nao_*nao_;
  ZTensor<2> rho(nao_, nao_);

  for(int t=0; t<=nt_; t++) {
    eKin_(t) = 0;
    ePot_(t) = 0;
    for(int s=0; s<2; s++) {
      G[s].get().get_dm(t,rho.data());

      auto rhomat = ZMatrixMap(rho.data(), nao_, nao_);
      auto hmfmat = DMatrixConstMap(p_MatSim_->fock().data()+s*nao2, nao_, nao_);
      auto h0mat = DMatrixConstMap(h0.data(),nao_,nao_);

      eKin_(t) += 0.5*rhomat.cwiseProduct((hmfmat+h0mat).transpose()).sum().real();
      ePot_(t) += Dyson.energy_conv(t, Sigma[s], G[s], beta_, dt_);
    }
  }
}

template <typename Repr>
void tti_DecompSpinSimulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}


template <>
inline void tti_DecompSpinSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
  int nL = p_MatSim_->frepr().nl();
  DMatrix Trans(ntau_+1, nL);
  for(int t=0; t<=ntau_; t++){
    double x = Dyson.Convolution().collocation().x_i()(t);
    for(int l=0; l<nL; l++){
      Trans(t,l) = boost::math::chebyshev_t(l,x);
    }
  }
  int nao2 = nao_*nao_;
  for(int i=0; i<nao2; i++){
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Gup.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data() + i, nL, Eigen::InnerStride<>(nao2));

    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Gdown.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data() + nL*nao_*nao_ + i, nL, Eigen::InnerStride<>(nao2));

    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sup.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->sigmal().data()+i, nL, Eigen::InnerStride<>(nao2));

    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sdown.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->sigmal().data() + nL*nao_*nao_ + i, nL, Eigen::InnerStride<>(nao2));
  }
}


template <>
inline void tti_DecompSpinSimulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}
} // namespace NEdyson
#endif // header guard
