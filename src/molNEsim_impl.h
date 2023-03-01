//
// Created by tblommel on 6/8/20
//

#ifndef NE_SIMULATION_IMPL_H
#define NE_SIMULATION_IMPL_H

#include <fstream>
#include <boost/math/special_functions/chebyshev.hpp>

namespace NEdyson {

template <typename Repr>
Simulation<Repr>::Simulation(const gfmol::HartreeFock &hf,
                             const gfmol::RepresentationBase<Repr> &frepr,
                             const gfmol::RepresentationBase<Repr> &brepr,
                             const Params &p) : 
                                 SimulationBase(hf, p),
                                 hmf(p.nt + 1, nao_, nao_), 
                                 h0(hf.hcore()), 
                                 rho(nao_, nao_)
{
  switch (p.gfmolmode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::Simulation<Repr> >(new gfmol::Simulation<Repr>(hf, frepr, brepr, p.gfmolmode, 0., p.hfbool));
      beta_ = p_MatSim_->frepr().beta();
      p_NEgf2_ = std::unique_ptr<molGF2Solver>(new molGF2Solver(hf.uchem(), dynamic_cast<gfmol::GF2Solver *>(p_MatSim_->p_sigma().get())->Vijkl_exch() ));
  }

  Sigma = GREEN(p.nt, p.ntau, nao_, -1);
  G = GREEN(p.nt, p.ntau, nao_, -1);
}


template <typename Repr>
void Simulation<Repr>::free_gf() {
  Dyson.G0_from_h0(G, p_MatSim_->mu(), p_MatSim_->fock(), p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void Simulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}

template <typename Repr>
void Simulation<Repr>::Ed_contractions(int tstp) {
  int nao = hmf.shape()[2];
  for(int d = 0; d < 3; d++) {
    ZMatrixMap(hmf.data() + tstp*nao*nao, nao, nao) += (Efield_(d, tstp) + efield_(d, tstp) + dfield_(d, tstp)) * DMatrixMap(dipole_.data() + d*nao*nao, nao, nao);
  }
}


template <typename Repr>
void Simulation<Repr>::do_boot() {
  for(int iter = 0; iter <= BootMax_; iter++){
    double err = 0;

    // Update mean field & self energy
    for(int tstp = 0; tstp <= k_; tstp++){
      ZMatrixMap(hmf.data() + tstp*nao_*nao_, nao_, nao_) = DMatrixConstMap(p_MatSim_->fock().data(),nao_,nao_);
      if(!hfbool_)  p_NEgf2_->solve(tstp, Sigma, G);
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
void Simulation<Repr>::do_tstp(int tstp) {
  // Predictor
  Dyson.Extrapolate(tstp, G);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    G.get_dm(tstp, rho);
    ZMatrixMap(hmf.data() + tstp*nao_*nao_, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
    
    // DEBUG  RAMP
    ZMatrix A = ZMatrix::Zero(2,2);
    A(0,0) = 1;
//    ZMatrixMap(hmf.data() + tstp*nao_*nao_, nao_, nao_) += 1 * (tanh(10*(tstp*dt_-1))+1) * A + 1 * (tanh(10*(tstp*dt_-3))+1) * A;
    // END DEBUG

    // HF Contractions
    p_NEgf2_->solve_HF(tstp, hmf, rho);

    // 2B Contractions
    if(!hfbool_) p_NEgf2_->solve(tstp, Sigma, G);

    // Efield contractions
    if(boolPumpProbe_) {
      Dyson.dipole_field(tstp, dfield_, G, G, dipole_, lPumpProbe_, nPumpProbe_, dt_);
      Ed_contractions(tstp);
    }

    Dyson.dyson_step(tstp, G, Sigma, hmf, p_MatSim_->mu(), beta_, dt_);
  }
}


template <typename Repr>
void Simulation<Repr>::save(h5::File &file, const std::string &path) {
  save_base(file, path);

  G.print_to_file(file, path + "/G");
  Sigma.print_to_file(file, path + "/Sigma");

  h5e::dump(file, path + "/params/beta", beta_);
  h5e::dump(file, path + "/mu", p_MatSim_->mu());
  h5e::dump(file, path + "/hmf", hmf);
  h5e::dump(file, path + "/params/h0", h0);
  h5e::dump(file, path + "/h0minusMu", DMatrixConstMap(h0.data(), nao_, nao_) - p_MatSim_->mu()*DMatrix::Identity(nao_, nao_));
  h5e::dump(file, path + "/GL00", -ZMatrixMap(G.lesptr(0,0), nao_, nao_).adjoint());
  h5e::dump(file, path + "/GG00",  ZMatrixMap(G.retptr(0,0), nao_, nao_) - ZMatrixMap(G.lesptr(0,0), nao_, nao_).adjoint());
  h5e::dump(file, path + "/SL00", -ZMatrixMap(Sigma.lesptr(0,0), nao_, nao_).adjoint());
  h5e::dump(file, path + "/SG00",  ZMatrixMap(Sigma.retptr(0,0), nao_, nao_) - ZMatrixMap(Sigma.lesptr(0,0), nao_, nao_).adjoint());
  h5e::dump(file, path + "/params/filling", p_MatSim_->filling());
  h5e::dump(file, path + "/tti", 0);
  
  h5e::dump(file, path + "/energy/EkinM", p_MatSim_->ehf() + p_MatSim_->ekin());
  h5e::dump(file, path + "/energy/EpotM", p_MatSim_->epot());

  ZTensor<3> densm(nt_+1, nao_, nao_);
  for(int i = 0; i<nao_; i++){
    for(int j = 0; j<nao_; j++){
      for(int t = 0; t<=nt_; t++){
        densm(t,i,j) = -cplx(0.,1.)*G.lesptr(t,t)[i*nao_+j];
      }
    }
  }
  h5e::dump(file, path + "/rho", densm);

  h5e::dump(file, path + "/rhoM", p_MatSim_->rho());
  h5e::dump(file, path + "/hmfM", p_MatSim_->fock());

  h5e::dump(file, path + "/U", p_NEgf2_->Uijkl());
  h5e::dump(file, path + "/U_flat", DColVectorConstMap(p_NEgf2_->Uijkl().data(), nao_*nao_*nao_*nao_));
  h5e::dump(file, path + "/U_ex_flat", DColVectorConstMap(dynamic_cast<gfmol::GF2Solver *>(p_MatSim_->p_sigma().get())->Vijkl_exch().data(), nao_*nao_*nao_*nao_));

//  ZTensor<3> coeff(ntau_+1, nao_, nao_);
//  Dyson.Convolution().collocation().to_spectral(coeff, ZTensorView<3>(G.matptr(0), ntau_+1, nao_, nao_));
//  h5e::dump(file, path + "/Gcoeff", coeff);
//  Dyson.Convolution().collocation().to_spectral(coeff, ZTensorView<3>(Sigma.matptr(0), ntau_+1, nao_, nao_));
//  h5e::dump(file, path + "/Scoeff", coeff);
}


template <typename Repr>
void Simulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}


template <typename Repr>
void Simulation<Repr>::do_energy() {
  int nao2 = nao_*nao_;
  for(int t=0; t<=nt_; t++) {
    G.get_dm(t, rho);

    auto rhomat = ZMatrixMap(rho.data(), nao_, nao_);
    auto hmfmat = ZMatrixMap(hmf.data() + t*nao2, nao_, nao_);
    auto h0mat = DMatrixConstMap(h0.data(),nao_,nao_);

    eKin_(t) = rhomat.cwiseProduct((hmfmat+h0mat).transpose()).sum().real();
    if(boolPumpProbe_) {
      for(int i=0; i<3; i++) {
        eKin_(t) += rhomat.cwiseProduct((Efield_(i,t) + efield_(i,t) + dfield_(i,t))
                                        * DMatrixMap(dipole_.data() + i*nao2, nao_, nao_).transpose()).sum().real();
      }
    }
    ePot_(t) = 2*Dyson.energy_conv(t, Sigma, G, beta_, dt_);
  }
}


template <>
inline void Simulation<gfmol::ChebyshevRepr>::L_to_Tau(){
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
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(G.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data()+i, nL, Eigen::InnerStride<>(nao2));
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sigma.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->sigmal().data()+i, nL, Eigen::InnerStride<>(nao2));
  }
}


template <>
inline void Simulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}


//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////








template <typename Repr>
tti_Simulation<Repr>::tti_Simulation(const gfmol::HartreeFock &hf,
                             const gfmol::RepresentationBase<Repr> &frepr,
                             const gfmol::RepresentationBase<Repr> &brepr,
                             const Params &p) : 
                                 SimulationBase(hf, p),
                                 h0(hf.hcore()) 
{
  switch (p.gfmolmode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::Simulation<Repr> >(new gfmol::Simulation<Repr>(hf, frepr, brepr, p.gfmolmode, 0., p.hfbool));
      beta_ = p_MatSim_->frepr().beta();
      p_NEgf2_ = std::unique_ptr<tti_molGF2Solver>(new tti_molGF2Solver(hf.uchem(), dynamic_cast<gfmol::GF2Solver *>(p_MatSim_->p_sigma().get())->Vijkl_exch() ));
  }

  Sigma = TTI_GREEN(p.nt, p.ntau, nao_, -1);
  G = TTI_GREEN(p.nt, p.ntau, nao_, -1);
}


template <typename Repr>
void tti_Simulation<Repr>::free_gf() {
  Dyson.G0_from_h0(G, p_MatSim_->mu(), p_MatSim_->fock(), p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void tti_Simulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}

template <typename Repr>
void tti_Simulation<Repr>::Ed_contractions(int tstp) {
  int a = 0;
}

template <typename Repr>
void tti_Simulation<Repr>::do_boot() {
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
void tti_Simulation<Repr>::do_tstp(int tstp) {
  // Predictor
  Dyson.Extrapolate(tstp, G);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {

    if(!hfbool_) p_NEgf2_->solve(tstp, Sigma, G);
    Dyson.dyson_step(tstp, G, Sigma, p_MatSim_->fock(), p_MatSim_->mu(), beta_, dt_);
  }
}

template <typename Repr>
double CoeffMax(ZTensorView<1> GetMaxOf) {
  int length = GetMaxOf.shape()[0];
  double ret = std::abs(GetMaxOf(0));
  for(int i = 1; i < length; i++) {
    if(std::abs(GetMaxOf(i)) > ret) ret = std::abs(GetMaxOf(i));
  }
  return ret;
}

template <typename Repr>
void tti_Simulation<Repr>::save(h5::File &file, const std::string &path) {
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
void tti_Simulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}

template <typename Repr>
void tti_Simulation<Repr>::do_energy() {
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

template <>
inline void tti_Simulation<gfmol::ChebyshevRepr>::L_to_Tau(){
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
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(G.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data()+i, nL, Eigen::InnerStride<>(nao2));
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sigma.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->sigmal().data()+i, nL, Eigen::InnerStride<>(nao2));
  }
}


template <>
inline void tti_Simulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}

} // namespace NEdyson
#endif // header guard
