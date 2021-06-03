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
                             int nt, int ntau, int k, double dt,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping, bool hfbool, bool boolPumpProbe, 
                             std::string PumpProbeInp, std::string MolInp,
                             double lPumpProbe, double nPumpProbe) : 
                                 SimulationBase(hf, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, hfbool, boolPumpProbe, PumpProbeInp, MolInp, lPumpProbe, nPumpProbe),
                                 hmf(nt+1, nao_, nao_), 
                                 h0(hf.hcore()), 
                                 rho(nao_, nao_)
{
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::Simulation<Repr> >(new gfmol::Simulation<Repr>(hf, frepr, brepr, mode, 0., hfbool));
      beta_ = p_MatSim_->frepr().beta();
      p_NEgf2_ = std::unique_ptr<molGF2Solver>(new molGF2Solver(hf.uchem(), p_MatSim_->u_exch()));
  }

  Sigma = GREEN(nt, ntau, nao_, -1);
  G = GREEN(nt, ntau, nao_, -1);
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
      G.get_dm(tstp, rho);
      ZMatrixMap(hmf.data() + tstp*nao_*nao_, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
      p_NEgf2_->solve_HF(tstp, hmf, rho);
      if(!hfbool_)  p_NEgf2_->solve(tstp, Sigma, G);
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
void Simulation<Repr>::do_tstp(int tstp) {
  // Predictor
  Dyson.Extrapolate(tstp, G);

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
  std::ofstream out;
  std::string data_dir = std::string(DATA_DIR);
  out.open(data_dir + "tottiming.dat" + "," + std::to_string(G.size1()) + "," + std::to_string(G.nt()) + "," + std::to_string(G.ntau()), std::ofstream::app);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    G.get_dm(tstp, rho);
    ZMatrixMap(hmf.data() + tstp*nao_*nao_, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);

    // HF Contractions
    start = std::chrono::system_clock::now();
    p_NEgf2_->solve_HF(tstp, hmf, rho);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    out << elapsed_seconds.count() << " ";

    // 2B Contractions
    start = std::chrono::system_clock::now();
    if(!hfbool_) p_NEgf2_->solve(tstp, Sigma, G);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    out << elapsed_seconds.count() << " ";

    // Efield contractions
    if(boolPumpProbe_) {
      Dyson.dipole_field(tstp, dfield_, G, G, dipole_, lPumpProbe_, nPumpProbe_, dt_);
      Ed_contractions(tstp);
    }

    start = std::chrono::system_clock::now();
    Dyson.dyson_step(tstp, G, Sigma, hmf, p_MatSim_->mu(), beta_, dt_);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    out << elapsed_seconds.count() << std::endl;
  }
}


template <typename Repr>
void Simulation<Repr>::save(h5::File &file, const std::string &path) {
  save_base(file, path);

  G.print_to_file(file, path + "/G");
  Sigma.print_to_file(file, path + "/Sigma");

  //FUCK
//  h5e::File fin("/pauli-storage/tblommel/He-VB2PP/UGU/H2-pvdz1.4out_nt100dt01_HF.h5");
//  G.read_from_file(fin, "/G");
  for(int t = 0; t<=nt_; t++) {
    p_NEgf2_->solve(t, Sigma, G);
  }
  Sigma.print_to_file(file, path + "/Sigma2BHF");

  GREEN Ints = GREEN(nt_, ntau_, nao_, -1);
  std::memset(Ints.lesptr(0,0), 0, ((nt_+1)*(nt_+2))/2*nao_*nao_*sizeof(cplx));
  for(int n = 0; n <= nt_; n++) {
    for(int m = 0; m <= n; m++) {
      Dyson.Cles2_tstp(m, m, Sigma, Sigma, G, G, n, dt_, Ints.lesptr(m,n));
    }
  }
  Ints.print_to_file(file, path + "/IL2");
  std::memset(Ints.lesptr(0,0), 0, ((nt_+1)*(nt_+2))/2*nao_*nao_*sizeof(cplx));
  for(int n = 0; n <= nt_; n++) {
    for(int m = 0; m <= n; m++) {
      Dyson.Cles3_tstp(m, m, Sigma, Sigma, G, G, n, beta_, Ints.lesptr(m,n));
    }
  }
  Ints.print_to_file(file, path + "/IL3");
  std::memset(Ints.lesptr(0,0), 0, ((nt_+1)*(nt_+2))/2*nao_*nao_*sizeof(cplx));
  for(int m = 0; m <= nt_; m++) {
    for(int n = 0; n <= m; n++) {
      if(n > k_) {
        for(int j = 0; j <= n; j++) {
          ZMatrixMap(Ints.lesptr(n,m), nao_, nao_) += dt_ * Dyson.II().gregory_weights(n,j) *              ZMatrixMap(Sigma.retptr(n,j), nao_, nao_) * ZMatrixMap(G.lesptr(j,m), nao_, nao_);
        }
      }
      else if(n <= k_) {
        for(int j = 0; j <= k_; j++) {  
          if(n >= j && j <= m){
            ZMatrixMap(Ints.lesptr(n,m), nao_, nao_) += dt_ * Dyson.II().gregory_weights(n,j) *              ZMatrixMap(Sigma.retptr(n,j), nao_, nao_) * ZMatrixMap(G.lesptr(j,m), nao_, nao_);
          }
          else if(n >= j && j > m){
            ZMatrixMap(Ints.lesptr(n,m), nao_, nao_) -= dt_ * Dyson.II().gregory_weights(n,j) *              ZMatrixMap(Sigma.retptr(n,j), nao_, nao_) * ZMatrixMap(G.lesptr(m,j), nao_, nao_).adjoint();
          }
          else if(n < j && j <= m){
            ZMatrixMap(Ints.lesptr(n,m), nao_, nao_) -= dt_ * Dyson.II().gregory_weights(n,j) *              ZMatrixMap(Sigma.retptr(j,n), nao_, nao_).adjoint() * ZMatrixMap(G.lesptr(j,m), nao_, nao_);
          }
          else if(n < j && j > m){
            ZMatrixMap(Ints.lesptr(n,m), nao_, nao_) += dt_ * Dyson.II().gregory_weights(n,j) *              ZMatrixMap(Sigma.retptr(j,n), nao_, nao_).adjoint() * ZMatrixMap(G.lesptr(m,j), nao_, nao_).adjoint();
          }
        }
      }
    }
  }
  Ints.print_to_file(file, path + "/IL1");
  std::memset(Ints.retptr(0,0), 0, ((nt_+1)*(nt_+2))/2*nao_*nao_*sizeof(cplx));
  for(int n = 0; n <= nt_; n++) {
    for(int m = 0; m <= n; m++) {
      if(n > k_ && n-m > k_){
        for(int j = m; j <= n; j++) {
          ZMatrixMap(Ints.retptr(n,m), nao_, nao_) += dt_ * Dyson.II().gregory_weights(n-m, j-m) * ZMatrixMap(G.retptr(n,j), nao_, nao_) * ZMatrixMap(Sigma.retptr(j,m), nao_, nao_);
        }
      }
      else if(n > k_ && n-m <= k_){
        for(int j = 0; j <= k_; j++) {
          if(n-j >= m) {
            ZMatrixMap(Ints.retptr(n,m), nao_, nao_) += dt_ * Dyson.II().gregory_weights(n-m, j) * ZMatrixMap(G.retptr(n,n-j), nao_, nao_) * ZMatrixMap(Sigma.retptr(n-j,m), nao_, nao_);
          }
          else if(n-j < m) {
            ZMatrixMap(Ints.retptr(n,m), nao_, nao_) -= dt_ * Dyson.II().gregory_weights(n-m, j) * ZMatrixMap(G.retptr(n,n-j), nao_, nao_) * ZMatrixMap(Sigma.retptr(m,n-j), nao_, nao_).adjoint();
          }
        }
      }
      else if(n <= k_) {
        for(int j=0; j <= k_; j++) {
          if(n >= j && j >= m){
            ZMatrixMap(Ints.retptr(n,m), nao_, nao_) += dt_ * Dyson.II().poly_integ(m,n,j) * ZMatrixMap(G.retptr(n,j), nao_, nao_) * ZMatrixMap(Sigma.retptr(j,m), nao_, nao_); 
          }
          if(n >= j && j < m){
            ZMatrixMap(Ints.retptr(n,m), nao_, nao_) -= dt_ * Dyson.II().poly_integ(m,n,j) * ZMatrixMap(G.retptr(n,j), nao_, nao_) * ZMatrixMap(Sigma.retptr(m,j), nao_, nao_).adjoint(); 
          }
          if(n < j && j >= m){
            ZMatrixMap(Ints.retptr(n,m), nao_, nao_) -= dt_ * Dyson.II().poly_integ(m,n,j) * ZMatrixMap(G.retptr(j,n), nao_, nao_).adjoint() * ZMatrixMap(Sigma.retptr(j,m), nao_, nao_); 
          }
          if(n < j && j < m){
            ZMatrixMap(Ints.retptr(n,m), nao_, nao_) += dt_ * Dyson.II().poly_integ(m,n,j) * ZMatrixMap(G.retptr(j,n), nao_, nao_).adjoint() * ZMatrixMap(Sigma.retptr(m,j), nao_, nao_).adjoint(); 
          }
        }
      }
    }
  }
  Ints.print_to_file(file, path + "/IR");
  
  //FUCK
  
  h5e::dump(file, path + "/params/beta", beta_);
  h5e::dump(file, path + "/mu", p_MatSim_->mu());
  h5e::dump(file, path + "/hmf", hmf);
  h5e::dump(file, path + "/params/h0", h0);
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
                             int nt, int ntau, int k, double dt,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping, bool hfbool) : 
                                 SimulationBase(hf, nt, ntau, k, dt, MatMax, MatTol, BootMax, BootTol, CorrSteps, hfbool, false, "", "", 0, 0),
                                 h0(hf.hcore()) 
{
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::Simulation<Repr> >(new gfmol::Simulation<Repr>(hf, frepr, brepr, mode, 0., hfbool));
      beta_ = p_MatSim_->frepr().beta();
      p_NEgf2_ = std::unique_ptr<tti_molGF2Solver>(new tti_molGF2Solver(hf.uchem(), p_MatSim_->u_exch()));
  }

  Sigma = TTI_GREEN(nt, ntau, nao_, -1);
  G = TTI_GREEN(nt, ntau, nao_, -1);
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

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
  std::ofstream out;
  std::string data_dir = std::string(DATA_DIR);
  out.open(data_dir + "tti_tottiming.dat" + "," + std::to_string(G.size1()) + "," + std::to_string(G.nt()) + "," + std::to_string(G.ntau()), std::ofstream::app);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {

    start = std::chrono::system_clock::now();
    if(!hfbool_) p_NEgf2_->solve(tstp, Sigma, G);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    out << elapsed_seconds.count() << " ";

    start = std::chrono::system_clock::now();
    Dyson.dyson_step(tstp, G, Sigma, p_MatSim_->fock(), p_MatSim_->mu(), beta_, dt_);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    out << elapsed_seconds.count() << std::endl;
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

  std::string data_dir = std::string(DATA_DIR);
  std::ofstream ofile;

  ofile.open(data_dir + "/tti_LegCoeff.dat" + "," + std::to_string(nao_) + "," + std::to_string(nt_) + "," + std::to_string(ntau_), std::ofstream::out);
  auto c = Dyson.Convolution().collocation();
  ZTensor<3> G_naa(ntau_+1, nao_, nao_);
  for(int t = 0; t <= nt_; t++) {
    ZTensorView<3> GTVTV(G.tvptr(t,0), ntau_+1, nao_, nao_);
    c.to_spectral(G_naa, GTVTV);
    for(int tau = 0; tau <= ntau_; tau++) {
      ZTensorView<1> GetMaxOf(G_naa.data() + tau*nao_*nao_, nao_*nao_);
      Sigma.tvptr(t,tau)[0] = CoeffMax<Repr>(GetMaxOf);
    }
  }

  for(int tau = 0; tau <= ntau_; tau++) {
    for(int t = 0; t <= nt_; t++) {
      ofile << Sigma.tvptr(t, tau)[0].real() << " ";
    }
    ofile << std::endl;
  }

  ofile.close();

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
inline void tti_Simulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}

} // namespace NEdyson
#endif // header guard
