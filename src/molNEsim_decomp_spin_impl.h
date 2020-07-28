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
                             int nt, int ntau, int k, double dt, int nw, double wmax,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping,
                             double decomp_prec) : 
                                 SimulationBase(hf, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps),
                                 hmf(2, nt+1, nao_, nao_), 
                                 h0(hf.hcore()), 
                                 rho(2,nao_,nao_)
{
  int nl = frepr.nl();
  const size_t size_MB = 1024*1024;
  size_t mem = 0;
  // Members of gfmol::sim
  mem += 2*5*nao_*nao_*sizeof(double); 
  mem += 2*5*nao_*nao_*nl*sizeof(double);
  mem += 2*3*nao_*nao_*nl*sizeof(cplx);
  // gfmol::SESolver
  mem += 2*4*nao_*nao_*nao_*sizeof(double); // V
  mem += 2*(2*4+1)*nao_*nao_*nao_*sizeof(double); // X,Y,tmp
  mem += (4*4+1)*nao_*nao_*sizeof(double);      // P, rho
  mem += nao_*nao_*nao_*nao_*sizeof(double);// Z
  // gfmol::repn
  mem += 2*2*nl*sizeof(double);
  mem += 2*2*nl*sizeof(int);
  mem += 2*3*nl*nl*sizeof(double);
  mem += 2*2*nl*nl*sizeof(cplx);
  // Members of NEdyson::sim
  mem += 2*(nt+1)*(nao_*nao_+2)*sizeof(double);     // hmf, energy
  mem += 2*nao_*nao_*sizeof(double);                // rho
  mem += 4*(nt+1)*(nt+2)/2*nao_*nao_*sizeof(cplx);// G,S, R,<
  mem += 4*(ntau+1)*(nt+1)*nao_*nao_*sizeof(cplx);// G,S, tv
  mem += 4*(ntau+1)*nao_*nao_*sizeof(cplx);       // G,S, M
  mem += nw*(nt+1)*nao_*nao_*sizeof(cplx);        // A
  // molNEgf2
  mem += nao_*nao_*nao_*sizeof(cplx);   // tmp
  mem += 10*nao_*nao_*4*nao_*sizeof(cplx); // X, Y
  mem += 2*4*4*nao_*nao_*sizeof(cplx); // P
  mem += nao_*nao_*nao_*nao_*sizeof(cplx); // Z
  // dyson
  mem += k*nao_*nao_*(k+2)*sizeof(cplx);
  mem += (2*nt+ntau)*nao_*nao_*sizeof(cplx); // temporary integral storage
  
  std::cout<< " Approximate memory needed for simulation : " << std::ceil(mem / (double)size_MB) << " MB"<<std::endl;
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::DecompSpinSimulation<Repr> >(new gfmol::DecompSpinSimulation<Repr>(hf, frepr, brepr, mode, 0.));
      beta_ = p_MatSim_->frepr().beta();
      dtau_ = beta_/ntau;
      p_NEgf2_ = std::unique_ptr<molGF2SolverSpinDecomp>(new molGF2SolverSpinDecomp(p_MatSim_->Vija(), p_MatSim_->Viaj()));
  }

  Sup = GREEN(nt, ntau, nao_, -1);
  Sdown = GREEN(nt, ntau, nao_, -1);
  Gup = GREEN(nt, ntau, nao_, -1);
  Gdown = GREEN(nt, ntau, nao_, -1);

  G = {Gup, Gdown};
  Sigma = {Sup, Sdown};

  A = SPECT();
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
void DecompSpinSimulation<Repr>::do_boot() {
  int nao2 = nao_*nao_;
  for(int iter = 0; iter <= BootMax_; iter++){
    double err = 0;

    // Update mean field & self energy
    for(int tstp = 0; tstp <= k_; tstp++){
      Gup.get_dm(tstp, rho.data());
      Gdown.get_dm(tstp, rho.data() + nao2);
        
      ZMatrixMap(hmf.data() + tstp*nao2, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
      ZMatrixMap(hmf.data() + (nt_+1)*nao2 + tstp*nao2, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);

      p_NEgf2_->solve_HF(tstp, hmf, rho);

      p_NEgf2_->solve(tstp, Sigma, G);
    }

    // Solve G Equation of Motion
    err = Dyson.dyson_start(Gup, Sup, hmf.data(), p_MatSim_->mu()[0], beta_, dt_);
    err = Dyson.dyson_start(Gdown, Sdown, hmf.data() + (nt_+1)*nao2, p_MatSim_->mu()[1], beta_, dt_);

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
    p_NEgf2_->solve(tstp, Sigma, G);

    Dyson.dyson_step(tstp, Gup, Sup, hmf.data(), p_MatSim_->mu()[0], beta_, dt_);
    Dyson.dyson_step(tstp, Gdown, Sdown, hmf.data() + (nt_+1)*nao2, p_MatSim_->mu()[1], beta_, dt_);
  }
}


template <typename Repr>
void DecompSpinSimulation<Repr>::save(h5::File &file, const std::string &path) {
  save_base(file, path);

  Gup.print_to_file(file, path + "/G/up");
  Sup.print_to_file(file, path + "/Sigma/up");
  Gup.print_to_file(file, path + "/G/up");
  Sup.print_to_file(file, path + "/Sigma/up");
  A.print_to_file(file, path + "/A/up");
  A.AfromG(Gdown,nw_,wmax_,dt_);
  A.print_to_file(file, path + "/A/down");
  
  h5e::dump(file, path + "/params/beta", beta_);
  h5e::dump(file, path + "/params/dtau", dtau_);
  h5e::dump(file, path + "/mu", std::vector<double>{p_MatSim_->mu()[0], p_MatSim_->mu()[1]});
  h5e::dump(file, path + "/hmf", hmf);
  h5e::dump(file, path + "/params/h0", h0);
  h5e::dump(file, path + "/params/filling", std::vector<double>{p_MatSim_->filling()[0], p_MatSim_->filling()[1]});
  h5e::dump(file, path + "/tti", 0);
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


template <typename Repr>
void DecompSpinSimulation<Repr>::do_spectral() {
  A.AfromG(Gup,nw_,wmax_,dt_);
}


template <>
inline void DecompSpinSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
  int ntau = ntau_;
  int nL = p_MatSim_->frepr().nl();
  DMatrix Trans(ntau+1, nL);
  for(int t=0; t<=ntau; t++){
    double x = (2.*t-ntau)/ntau;
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
                             int nt, int ntau, int k, double dt, int nw, double wmax,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping,
                             double decomp_prec) : 
                                 SimulationBase(hf, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps),
                                 h0(hf.hcore())
{
  int nl = frepr.nl();
  const size_t size_MB = 1024*1024;
  size_t mem = 0;
  // Members of gfmol::sim
  mem += 2*5*nao_*nao_*sizeof(double); 
  mem += 2*5*nao_*nao_*nl*sizeof(double);
  mem += 2*3*nao_*nao_*nl*sizeof(cplx);
  // gfmol::SESolver
  mem += 2*4*nao_*nao_*nao_*sizeof(double); // V
  mem += 2*(2*4+1)*nao_*nao_*nao_*sizeof(double); // X,Y,tmp
  mem += (4*4+1)*nao_*nao_*sizeof(double);      // P, rho
  mem += nao_*nao_*nao_*nao_*sizeof(double);// Z
  // gfmol::repn
  mem += 2*2*nl*sizeof(double);
  mem += 2*2*nl*sizeof(int);
  mem += 2*3*nl*nl*sizeof(double);
  mem += 2*2*nl*nl*sizeof(cplx);
  // Members of NEdyson::sim
  mem += 2*(nt+1)*(nao_*nao_+2)*sizeof(double);     // hmf, energy
  mem += 2*nao_*nao_*sizeof(double);                // rho
  mem += 4*(nt+1)*nao_*nao_*sizeof(cplx);// G,S, R,<
  mem += 4*(ntau+1)*(nt+1)*nao_*nao_*sizeof(cplx);// G,S, tv
  mem += 4*(ntau+1)*nao_*nao_*sizeof(cplx);       // G,S, M
  mem += nw*(nt+1)*nao_*nao_*sizeof(cplx);        // A
  // molNEgf2
  mem += nao_*nao_*nao_*sizeof(cplx);   // tmp
  mem += 10*nao_*nao_*4*nao_*sizeof(cplx); // X, Y
  mem += 2*4*4*nao_*nao_*sizeof(cplx); // P
  mem += nao_*nao_*nao_*nao_*sizeof(cplx); // Z
  // dyson
  mem += k*nao_*nao_*(k+2)*sizeof(cplx);
  mem += (2*nt+ntau)*nao_*nao_*sizeof(cplx); // temporary integral storage
  
  std::cout<< " Approximate memory needed for simulation : " << std::ceil(mem / (double)size_MB) << " MB"<<std::endl;
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::DecompSpinSimulation<Repr> >(new gfmol::DecompSpinSimulation<Repr>(hf, frepr, brepr, mode, 0.));
      beta_ = p_MatSim_->frepr().beta();
      dtau_ = beta_/ntau;
      p_NEgf2_ = std::unique_ptr<tti_molGF2SolverSpinDecomp>(new tti_molGF2SolverSpinDecomp(p_MatSim_->Vija(), p_MatSim_->Viaj()));
  }

  Sup = TTI_GREEN(nt, ntau, nao_, -1);
  Sdown = TTI_GREEN(nt, ntau, nao_, -1);
  Gup = TTI_GREEN(nt, ntau, nao_, -1);
  Gdown = TTI_GREEN(nt, ntau, nao_, -1);

  G = {Gup, Gdown};
  Sigma = {Sup, Sdown};

  A = TTI_SPECT();
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
void tti_DecompSpinSimulation<Repr>::do_boot() {
  int nao2 = nao_*nao_;
  for(int iter = 0; iter <= BootMax_; iter++){
    double err = 0;

    // Update mean field & self energy
    for(int tstp = 0; tstp <= k_; tstp++){
      p_NEgf2_->solve(tstp, Sigma, G);
    }

    // Solve G Equation of Motion
    err = Dyson.dyson_start(Gup, Sup, p_MatSim_->fock().data(), p_MatSim_->mu()[0], beta_, dt_);
    err = Dyson.dyson_start(Gdown, Sdown, p_MatSim_->fock().data() + nao2, p_MatSim_->mu()[1], beta_, dt_);

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
    p_NEgf2_->solve(tstp, Sigma, G);

    Dyson.dyson_step(tstp, Gup, Sup, p_MatSim_->fock().data(), p_MatSim_->mu()[0], beta_, dt_);
    Dyson.dyson_step(tstp, Gdown, Sdown, p_MatSim_->fock().data() + nao2, p_MatSim_->mu()[1], beta_, dt_);
  }
}


template <typename Repr>
void tti_DecompSpinSimulation<Repr>::save(h5::File &file, const std::string &path) {
  save_base(file, path);

  Gup.print_to_file(file, path + "/G/up");
  Sup.print_to_file(file, path + "/Sigma/up");
  Gup.print_to_file(file, path + "/G/up");
  Sup.print_to_file(file, path + "/Sigma/up");
  A.print_to_file(file, path + "/A/up");
  A.AfromG(Gdown,nw_,wmax_,dt_);
  A.print_to_file(file, path + "/A/down");
  
  h5e::dump(file, path + "/params/beta", beta_);
  h5e::dump(file, path + "/params/dtau", dtau_);
  h5e::dump(file, path + "/mu", std::vector<double>{p_MatSim_->mu()[0], p_MatSim_->mu()[1]});
  h5e::dump(file, path + "/hmf", p_MatSim_->fock());
  h5e::dump(file, path + "/params/h0", h0);
  h5e::dump(file, path + "/params/filling", std::vector<double>{p_MatSim_->filling()[0], p_MatSim_->filling()[1]});
  h5e::dump(file, path + "/tti", 1);
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


template <typename Repr>
void tti_DecompSpinSimulation<Repr>::do_spectral() {
  A.AfromG(Gup,nw_,wmax_,dt_);
}


template <>
inline void tti_DecompSpinSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
  int ntau = ntau_;
  int nL = p_MatSim_->frepr().nl();
  DMatrix Trans(ntau+1, nL);
  for(int t=0; t<=ntau; t++){
    double x = (2.*t-ntau)/ntau;
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
  }
}


template <>
inline void tti_DecompSpinSimulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}
} // namespace NEdyson
#endif // header guard
