//
// Created by tblommel on 6/8/20
//

#ifndef NE_SIMULATION_SPIN_IMPL_H
#define NE_SIMULATION_SPIN_IMPL_H

#include <fstream>
#include <boost/math/special_functions/chebyshev.hpp>

namespace NEdyson {

template <typename Repr>
SpinSimulation<Repr>::SpinSimulation(const gfmol::HartreeFock &hf,
                             const gfmol::RepresentationBase<Repr> &frepr,
                             const gfmol::RepresentationBase<Repr> &brepr,
                             int nt, int ntau, int k, double dt, int nw, double wmax,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping) : 
                                 SimulationBase(hf, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps),
                                 hmf(2, nt+1, nao_, nao_), 
                                 h0(hf.hcore()), 
                                 rho(2,nao_,nao_)
{
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::SpinSimulation<Repr> >(new gfmol::SpinSimulation<Repr>(hf, frepr, brepr, mode, 0.));
      beta_ = p_MatSim_->frepr().beta();
      dtau_ = beta_/ntau;
      p_NEgf2_ = std::unique_ptr<molGF2SolverSpin>(new molGF2SolverSpin(hf.uchem()));
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
void SpinSimulation<Repr>::free_gf() {
  G0_from_h0(Gup, p_MatSim_->mu()[0], h0, p_MatSim_->frepr().beta(), dt_);
  G0_from_h0(Gdown, p_MatSim_->mu()[1], h0, p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void SpinSimulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}


template <typename Repr>
void SpinSimulation<Repr>::do_boot() {
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
    err = dyson_start(I, Gup, Sup, hmf.data(), p_MatSim_->mu()[0], beta_, dt_);
    err = dyson_start(I, Gdown, Sdown, hmf.data() + (nt_+1)*nao2, p_MatSim_->mu()[1], beta_, dt_);

    std::cout<<"Bootstrapping iteration : "<<iter<<" | Error = "<<err<<std::endl;
    if(err<BootTol_){
      bootstrap_converged = true;
      break;
    }
  }
}


template <typename Repr>
void SpinSimulation<Repr>::do_tstp(int tstp) {
  int nao2 = nao_ * nao_;
  // Predictor
  Extrapolate(I, Gup, tstp);
  Extrapolate(I, Gdown, tstp);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    Gup.get_dm(tstp, rho.data());
    Gdown.get_dm(tstp, rho.data() + nao2);

    ZMatrixMap(hmf.data() + tstp*nao2, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
    ZMatrixMap(hmf.data() + (nt_+1)*nao2 + tstp*nao2, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);

    p_NEgf2_->solve_HF(tstp, hmf, rho);
    p_NEgf2_->solve(tstp, Sigma, G);

    dyson_step(tstp, I, Gup, Sup, hmf.data(), p_MatSim_->mu()[0], beta_, dt_);
    dyson_step(tstp, I, Gdown, Sdown, hmf.data() + (nt_+1)*nao2, p_MatSim_->mu()[1], beta_, dt_);
  }
}


template <typename Repr>
void SpinSimulation<Repr>::save(h5::File &file, const std::string &path) {
  Gup.print_to_file("/home/thomas/Libraries/NEdyson/up", dt_, dtau_);
  Gdown.print_to_file("/home/thomas/Libraries/NEdyson/down", dt_, dtau_);
  std::ofstream ofile;
  ofile.open("/home/thomas/Libraries/NEdyson/h2H.dat",std::ofstream::out);
  ofile<<nao_<<std::endl;
  ofile<<p_MatSim_->mu()[0]<<std::endl;
  ofile<<p_MatSim_->mu()[1]<<std::endl;
  for(int i = 0; i < nao_; i++){
    for(int j = 0; j < nao_; j++){
      ofile<<h0(i,j)<<std::endl;
    }
  }
  for(int i = 0; i < nao_; i++){
    for(int j = 0; j < nao_; j++){
      for(int k = 0; k < nao_; k++){
        for(int l = 0; l < nao_; l++){
          ofile<<p_NEgf2_->Uijkl()(i,j,k,l)<<std::endl;
        }
      }
    }
  }
  ofile.close();
  A.print_to_file("/home/thomas/Libraries/NEdyson/Aup");
  A.AfromG(Gdown,nw_,wmax_,dt_);
  A.print_to_file("/home/thomas/Libraries/NEdyson/Adown");
}


template <typename Repr>
void SpinSimulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}


template <typename Repr>
void SpinSimulation<Repr>::do_spectral() {
  A.AfromG(Gup,nw_,wmax_,dt_);
}


template <>
inline void SpinSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
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
inline void SpinSimulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}

//////////////////////////
//////////////////////////
//////////////////////////




template <typename Repr>
tti_SpinSimulation<Repr>::tti_SpinSimulation(const gfmol::HartreeFock &hf,
                             const gfmol::RepresentationBase<Repr> &frepr,
                             const gfmol::RepresentationBase<Repr> &brepr,
                             int nt, int ntau, int k, double dt, int nw, double wmax,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping) : 
                                 SimulationBase(hf, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps),
                                 hmf(2, nao_, nao_), 
                                 h0(hf.hcore()), 
                                 rho(2,nao_,nao_)
{
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::SpinSimulation<Repr> >(new gfmol::SpinSimulation<Repr>(hf, frepr, brepr, mode, 0.));
      beta_ = p_MatSim_->frepr().beta();
      dtau_ = beta_/ntau;
      p_NEgf2_ = std::unique_ptr<tti_molGF2SolverSpin>(new tti_molGF2SolverSpin(hf.uchem()));
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
void tti_SpinSimulation<Repr>::free_gf() {
  G0_from_h0(Gup, p_MatSim_->mu()[0], h0, p_MatSim_->frepr().beta(), dt_);
  G0_from_h0(Gdown, p_MatSim_->mu()[1], h0, p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void tti_SpinSimulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}


template <typename Repr>
void tti_SpinSimulation<Repr>::do_boot() {
  int nao2 = nao_*nao_;
  for(int iter = 0; iter <= BootMax_; iter++){
    double err = 0;

    // Update mean field & self energy
    Gup.get_dm(0, rho.data());
    Gdown.get_dm(0, rho.data() + nao2);

    ZMatrixMap(hmf.data(), nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
    ZMatrixMap(hmf.data() + nao2, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);

    p_NEgf2_->solve_HF(hmf, rho);

    for(int tstp = 0; tstp <= k_; tstp++){
      p_NEgf2_->solve(tstp, Sigma, G);
    }

    // Solve G Equation of Motion
    err = dyson_start(I, Gup, Sup, hmf.data(), p_MatSim_->mu()[0], beta_, dt_);
    err = dyson_start(I, Gdown, Sdown, hmf.data() + nao2, p_MatSim_->mu()[1], beta_, dt_);

    std::cout<<"Bootstrapping iteration : "<<iter<<" | Error = "<<err<<std::endl;
    if(err<BootTol_){
      bootstrap_converged = true;
      break;
    }
  }
}


template <typename Repr>
void tti_SpinSimulation<Repr>::do_tstp(int tstp) {
  int nao2 = nao_ * nao_;
  // Predictor
  Extrapolate(I, Gup, tstp);
  Extrapolate(I, Gdown, tstp);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    p_NEgf2_->solve(tstp, Sigma, G);

    dyson_step(tstp, I, Gup, Sup, hmf.data(), p_MatSim_->mu()[0], beta_, dt_);
    dyson_step(tstp, I, Gdown, Sdown, hmf.data() + nao2, p_MatSim_->mu()[1], beta_, dt_);
  }
}


template <typename Repr>
void tti_SpinSimulation<Repr>::save(h5::File &file, const std::string &path) {
  Gup.print_to_file("/home/thomas/Libraries/NEdyson/up", dt_, dtau_);
  Gdown.print_to_file("/home/thomas/Libraries/NEdyson/down", dt_, dtau_);
  std::ofstream ofile;
  ofile.open("/home/thomas/Libraries/NEdyson/h2H.dat",std::ofstream::out);
  ofile<<nao_<<std::endl;
  ofile<<p_MatSim_->mu()[0]<<std::endl;
  ofile<<p_MatSim_->mu()[1]<<std::endl;
  for(int i = 0; i < nao_; i++){
    for(int j = 0; j < nao_; j++){
      ofile<<h0(i,j)<<std::endl;
    }
  }
  for(int i = 0; i < nao_; i++){
    for(int j = 0; j < nao_; j++){
      for(int k = 0; k < nao_; k++){
        for(int l = 0; l < nao_; l++){
          ofile<<p_NEgf2_->Uijkl()(i,j,k,l)<<std::endl;
        }
      }
    }
  }
  ofile.close();
  A.print_to_file("/home/thomas/Libraries/NEdyson/Aup");
  A.AfromG(Gdown,nw_,wmax_,dt_);
  A.print_to_file("/home/thomas/Libraries/NEdyson/Adown");
}


template <typename Repr>
void tti_SpinSimulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}


template <typename Repr>
void tti_SpinSimulation<Repr>::do_spectral() {
  A.AfromG(Gup,nw_,wmax_,dt_);
}


template <>
inline void tti_SpinSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
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
inline void tti_SpinSimulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}

} // namespace NEdyson
#endif // header guard
