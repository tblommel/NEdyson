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
  int nl = frepr.nl();
  const size_t size_MB = 1024*1024;
  size_t mem = 0;
  // Members of gfmol::sim
  mem += 5*2*nao_*nao_*sizeof(double); 
  mem += 5*2*nao_*nao_*nl*sizeof(double);
  mem += 3*2*nao_*nao_*nl*sizeof(cplx);
  // gfmol::SESolver
  mem += 2*nao_*nao_*nao_*nao_*sizeof(double); // interaction tensor
  mem += 3*nao_*nao_*nao_*sizeof(double); // temp objects
  mem += 2*nao_*nao_*sizeof(double); // rho
  // gfmol::repn
  mem += 2*2*nl*sizeof(double);
  mem += 2*2*nl*sizeof(int);
  mem += 2*3*nl*nl*sizeof(double);
  mem += 2*2*nl*nl*sizeof(cplx);
  // Members of NEdyson::sim
  mem += 2*(nt+1)*(nao_*nao_+2)*sizeof(double); // hmf, rho
  mem += 4*(nt+1)*(nt+2)/2*nao_*nao_*sizeof(cplx); // G,S, R,<
  mem += 4*(ntau+1)*(nt+1)*nao_*nao_*sizeof(cplx); // G,S, tv
  mem += 4*(ntau+1)*nao_*nao_*sizeof(cplx);        // G,S, M
  mem += nw*(nt+1)*nao_*nao_*sizeof(cplx);         // A
  // molNEgf2
  mem += 5*nao_*nao_*(nao_+1)*sizeof(double);      // temp objects
  // dyson
  mem += k*nao_*nao_*(k+2)*sizeof(cplx);
  mem += (2*nt+ntau)*nao_*nao_*sizeof(cplx); // temporary integral storage
  
  std::cout<< " Approximate memory needed for simulation : " << std::ceil(mem / (double)size_MB) << " MB"<<std::endl;

  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::SpinSimulation<Repr> >(new gfmol::SpinSimulation<Repr>(hf, frepr, brepr, mode, 0.));
      beta_ = p_MatSim_->frepr().beta();
      dtau_ = beta_/ntau;
      p_NEgf2_ = std::unique_ptr<molGF2SolverSpin>(new molGF2SolverSpin(hf.uchem(), p_MatSim_->u_exch()));
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
  Dyson.G0_from_h0(Gup, p_MatSim_->mu()[0], p_MatSim_->fock().data(), p_MatSim_->frepr().beta(), dt_);
  Dyson.G0_from_h0(Gdown, p_MatSim_->mu()[1], p_MatSim_->fock().data() + nao_*nao_, p_MatSim_->frepr().beta(), dt_);
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
void SpinSimulation<Repr>::do_energy() {
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
void SpinSimulation<Repr>::do_tstp(int tstp) {
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
  ofile.open("/home/thomas/Libraries/NEdyson/energy.dat", std::ofstream::out);
  ofile << nt_ << std::endl;
  ofile << p_MatSim_->etot()-p_MatSim_->enuc() << std::endl;
  for(int t=0; t<=nt_; t++) ofile << ePot_(t)+eKin_(t) << std::endl;
  ofile.close();
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
                                 h0(hf.hcore())
{
  int nl = frepr.nl();
  const size_t size_MB = 1024*1024;
  size_t mem = 0;
  // Members of gfmol::sim
  mem += 5*2*nao_*nao_*sizeof(double); 
  mem += 5*2*nao_*nao_*nl*sizeof(double);
  mem += 3*2*nao_*nao_*nl*sizeof(cplx);
  // gfmol::SESolver
  mem += 2*nao_*nao_*nao_*nao_*sizeof(double); // interaction tensor
  mem += 3*nao_*nao_*nao_*sizeof(double); // temp objects
  mem += 2*nao_*nao_*sizeof(double); // rho
  // gfmol::repn
  mem += 2*2*nl*sizeof(double);
  mem += 2*2*nl*sizeof(int);
  mem += 2*3*nl*nl*sizeof(double);
  mem += 2*2*nl*nl*sizeof(cplx);
  // Members of NEdyson::sim
  mem += 2*(nao_*nao_)*sizeof(double); // rho
  mem += 4*(nt+1)*nao_*nao_*sizeof(cplx); // G,S, R,<
  mem += 4*(ntau+1)*(nt+1)*nao_*nao_*sizeof(cplx); // G,S, tv
  mem += 4*(ntau+1)*nao_*nao_*sizeof(cplx);        // G,S, M
  mem += nw*(nt+1)*nao_*nao_*sizeof(cplx);         // A
  // molNEgf2
  mem += 5*nao_*nao_*(nao_+1)*sizeof(double);      // temp objects
  // dyson
  mem += k*nao_*nao_*(k+2)*sizeof(cplx);
  mem += (2*nt+ntau)*nao_*nao_*sizeof(cplx); // temporary integral storage
  
  std::cout<< " Approximate memory needed for simulation : " << std::ceil(mem / (double)size_MB) << " MB"<<std::endl;

  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::SpinSimulation<Repr> >(new gfmol::SpinSimulation<Repr>(hf, frepr, brepr, mode, 0.));
      beta_ = p_MatSim_->frepr().beta();
      dtau_ = beta_/ntau;
      p_NEgf2_ = std::unique_ptr<tti_molGF2SolverSpin>(new tti_molGF2SolverSpin(hf.uchem(), p_MatSim_->u_exch()));
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
  Dyson.G0_from_h0(Gup, p_MatSim_->mu()[0], p_MatSim_->fock().data(), p_MatSim_->frepr().beta(), dt_);
  Dyson.G0_from_h0(Gdown, p_MatSim_->mu()[1], p_MatSim_->fock().data() + nao_*nao_, p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void tti_SpinSimulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
}

template <typename Repr>
void tti_SpinSimulation<Repr>::do_energy() {
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
void tti_SpinSimulation<Repr>::do_boot() {
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
void tti_SpinSimulation<Repr>::do_tstp(int tstp) {
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
void tti_SpinSimulation<Repr>::save(h5::File &file, const std::string &path) {
  Gup.print_to_file("/home/thomas/Libraries/NEdyson/up", dt_, dtau_);
  Gdown.print_to_file("/home/thomas/Libraries/NEdyson/down", dt_, dtau_);
  Sup.print_to_file("/home/thomas/Libraries/NEdyson/Sup", dt_, dtau_);
  Sdown.print_to_file("/home/thomas/Libraries/NEdyson/Sdown", dt_, dtau_);
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
  ofile.open("/home/thomas/Libraries/NEdyson/energy.dat", std::ofstream::out);
  ofile << nt_ << std::endl;
  ofile << p_MatSim_->etot()-p_MatSim_->enuc() << std::endl;
  for(int t=0; t<=nt_; t++) ofile << ePot_(t)+eKin_(t) << std::endl;
  ofile.close();
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
