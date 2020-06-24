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
                             int nt, int ntau, int k, double dt, int nw, double wmax,
                             int MatMax, double MatTol, int BootMax, double BootTol, int CorrSteps,
                             gfmol::Mode mode,
                             double damping,
                             double decomp_prec) : 
                                 SimulationBase(hf, nt, ntau, k, dt, nw, wmax, MatMax, MatTol, BootMax, BootTol, CorrSteps),
                                 hmf(nt+1, nao_, nao_), 
                                 h0(hf.hcore()), 
                                 rho(nao_,nao_)
{
  switch (mode) {
    case gfmol::Mode::GF2:
      p_MatSim_ = std::unique_ptr<gfmol::DecompSimulation<Repr> >(new gfmol::DecompSimulation<Repr>(hf, frepr, brepr, mode, 0., decomp_prec));
      beta_ = p_MatSim_->frepr().beta();
      dtau_ = beta_/ntau;
      p_NEgf2_ = std::unique_ptr<molGF2SolverDecomp>(new molGF2SolverDecomp(p_MatSim_->Vija()));
  }

  Sigma = GREEN(nt, ntau, nao_, -1);
  G = GREEN(nt, ntau, nao_, -1);
  A = SPECT();
}


template <typename Repr>
void DecompSimulation<Repr>::free_gf() {
  G0_from_h0(G, p_MatSim_->mu(), h0, p_MatSim_->frepr().beta(), dt_);
}


template <typename Repr>
void DecompSimulation<Repr>::do_mat() {
  p_MatSim_->run(MatMax_, MatTol_, nullptr);
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
      p_NEgf2_->solve(tstp, Sigma, G);
    }

    // Solve G Equation of Motion
    err = dyson_start(I, G, Sigma, hmf, p_MatSim_->mu(), beta_, dt_);

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
  Extrapolate(I,G,tstp);

  // Corrector
  for(int iter = 0; iter < CorrSteps_; iter++) {
    G.get_dm(tstp, rho);
    ZMatrixMap(hmf.data() + tstp*nao_*nao_, nao_, nao_) = DMatrixConstMap(h0.data(),nao_,nao_);
    p_NEgf2_->solve_HF(tstp, hmf, rho);
    p_NEgf2_->solve(tstp, Sigma, G);

    dyson_step(tstp, I, G, Sigma, hmf, p_MatSim_->mu(), beta_, dt_);
  }
}


template <typename Repr>
void DecompSimulation<Repr>::save(h5::File &file, const std::string &path) {
  G.print_to_file("/home/thomas/Libraries/NEdyson/", dt_, dtau_);
  std::ofstream ofile;
  ofile.open("/home/thomas/Libraries/NEdyson/h2H.dat",std::ofstream::out);
  ofile<<nao_<<std::endl;
  ofile<<p_MatSim_->mu()<<std::endl;
  ofile<<p_MatSim_->mu()<<std::endl;
  for(int i = 0; i < nao_; i++){
    for(int j = 0; j < nao_; j++){
      ofile<<h0(i,j)<<std::endl;
    }
  }
/*  for(int i = 0; i < nao_; i++){
    for(int j = 0; j < nao_; j++){
      for(int k = 0; k < nao_; k++){
        for(int l = 0; l < nao_; l++){
          ofile<<p_NEgf2_->Uijkl()(i,j,k,l)<<std::endl;
        }
      }
    }
  }
  ofile.close();*/
  A.print_to_file("/home/thomas/Libraries/NEdyson/");
}


template <typename Repr>
void DecompSimulation<Repr>::load(const h5::File &file, const std::string &path) {
  int a=0;
}


template <typename Repr>
void DecompSimulation<Repr>::do_spectral() {
  A.AfromG(G,nw_,wmax_,dt_);
}


template <>
inline void DecompSimulation<gfmol::ChebyshevRepr>::L_to_Tau(){
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
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(G.matptr(0)+i, ntau_+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(p_MatSim_->gl().data()+i, nL, Eigen::InnerStride<>(nao2));
  }
}


template <>
inline void DecompSimulation<gfmol::IntermediateRepr>::L_to_Tau(){
  std::cout<<"Not yet implemented"<<std::endl;
}

} // namespace NEdyson
#endif // header guard
