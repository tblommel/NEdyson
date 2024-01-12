//
// Created by Thomas Blommel on 5/21/20
// This defines teh functions that solve for the HF self energy using the non-decomp interaction tensor
//

#ifndef HF_SIGMASC_DEFN
#define HF_SIGMASC_DEFN

#include "greens.h"
#include "utils.h"

namespace NEdyson{

class SCHFSolver {
protected:
  ZMatrix Sigma3_;
  ZMatrix tmp_ij_;
  ZMatrix tmp2_ij_;

public:
  SCHFSolver() {
    Sigma3_ = ZMatrix::Identity(2,2);
    Sigma3_(1,1) = -1;
    tmp_ij_ = ZMatrix::Zero(2,2);
    tmp2_ij_ = ZMatrix::Zero(2,2);
  }

  void solve_Delta(int tstp, cplx *t0, GREEN &G, GREEN &Sigma) {
    int ntau = G.ntau();

    for(int t = 0; t <= tstp; t++) {
      tmp2_ij_ = Sigma3_ * ZMatrixMap(G.retptr(tstp,t), 2, 2) * Sigma3_;
      ZMatrixMap(Sigma.retptr(tstp,t), 2, 2) += 0.5 * ZMatrixMap(t0+tstp*4, 2, 2) * tmp2_ij_ * ZMatrixMap(t0+t*4, 2, 2).conjugate();
      ZMatrixMap(Sigma.retptr(tstp,t), 2, 2) += 0.5 * ZMatrixMap(t0+tstp*4, 2, 2).conjugate() * tmp2_ij_ * ZMatrixMap(t0+t*4, 2, 2);
      
      tmp2_ij_ = Sigma3_ * ZMatrixMap(G.lesptr(t,tstp), 2, 2) * Sigma3_;
      ZMatrixMap(Sigma.lesptr(t,tstp), 2, 2) += 0.5 * ZMatrixMap(t0+t*4, 2, 2) * tmp2_ij_ * ZMatrixMap(t0+tstp*4, 2, 2).conjugate();
      ZMatrixMap(Sigma.lesptr(t,tstp), 2, 2) += 0.5 * ZMatrixMap(t0+t*4, 2, 2).conjugate() * tmp2_ij_ * ZMatrixMap(t0+tstp*4, 2, 2);
    }

    for(int tau = 0; tau <= ntau; tau++) {
      tmp2_ij_ = Sigma3_ * ZMatrixMap(G.tvptr(tstp,tau), 2, 2) * Sigma3_;
      ZMatrixMap(Sigma.tvptr(tstp,tau), 2, 2) += 0.5 * ZMatrixMap(t0+tstp*4, 2, 2) * tmp2_ij_ * ZMatrixMap(t0, 2, 2).conjugate();
      ZMatrixMap(Sigma.tvptr(tstp,tau), 2, 2) += 0.5 * ZMatrixMap(t0+tstp*4, 2, 2).conjugate() * tmp2_ij_ * ZMatrixMap(t0, 2, 2);
    }
  }

  void solve_Sigma(int tstp, double U, GREEN &G, GREEN &Sigma) {
    int ntau = G.ntau();

    for(int t = 0; t <= tstp; t++) {
      // Greater component, evaluated on Retarded part of plane
      tmp_ij_ = ZMatrixMap(G.retptr(tstp, t), 2, 2) - ZMatrixMap(G.lesptr(t, tstp), 2, 2).adjoint();
      cplx *SLptr = Sigma.lesptr(t,tstp);
      cplx *SRptr = Sigma.retptr(tstp,t);
      cplx *GLptr = G.lesptr(t,tstp);
      cplx *GRptr = G.retptr(tstp,t);
      cplx *GGptr = tmp_ij_.data();

      SLptr[0] += U * U * GLptr[0] * GLptr[3] * GGptr[3];
      SLptr[1] += U * U * GLptr[1] * GLptr[2] * GGptr[1];
      SLptr[2] += U * U * GLptr[2] * GLptr[1] * GGptr[2];
      SLptr[3] += U * U * GLptr[3] * GLptr[0] * GGptr[0];

      SLptr[0] -= U * U * GLptr[1] * GGptr[3] * GLptr[2];
      SLptr[1] -= U * U * GLptr[0] * GGptr[1] * GLptr[3];
      SLptr[2] -= U * U * GLptr[3] * GGptr[2] * GLptr[0];
      SLptr[3] -= U * U * GLptr[2] * GGptr[0] * GLptr[1];

      SRptr[0] += U * U * (GGptr[0] * GGptr[3] * GLptr[3] + std::conj(GLptr[0] * GLptr[3] * GGptr[3]));
      SRptr[1] += U * U * (GGptr[1] * GGptr[2] * GLptr[1] + std::conj(GLptr[2] * GLptr[1] * GGptr[2]));
      SRptr[2] += U * U * (GGptr[2] * GGptr[1] * GLptr[2] + std::conj(GLptr[1] * GLptr[2] * GGptr[1]));
      SRptr[3] += U * U * (GGptr[3] * GGptr[0] * GLptr[0] + std::conj(GLptr[3] * GLptr[0] * GGptr[0]));

      SRptr[0] -= U * U * (GGptr[1] * GLptr[3] * GGptr[2] + std::conj(GLptr[2] * GGptr[3] * GLptr[1]));
      SRptr[1] -= U * U * (GGptr[0] * GLptr[1] * GGptr[3] + std::conj(GLptr[0] * GGptr[2] * GLptr[3]));
      SRptr[2] -= U * U * (GGptr[3] * GLptr[2] * GGptr[0] + std::conj(GLptr[3] * GGptr[1] * GLptr[0]));
      SRptr[3] -= U * U * (GGptr[2] * GLptr[0] * GGptr[1] + std::conj(GLptr[1] * GGptr[0] * GLptr[2]));
    }

    for(int tau = 0; tau <= ntau; tau++) {
      cplx *GTVptr = G.tvptr(tstp, tau);
      cplx *GVTptr = G.tvptr(tstp, ntau-tau);
      cplx *STVptr = Sigma.tvptr(tstp, tau);

      STVptr[0] += U * U * GTVptr[0] * GTVptr[3] * std::conj(GVTptr[3]);
      STVptr[1] += U * U * GTVptr[1] * GTVptr[2] * std::conj(GVTptr[2]);
      STVptr[2] += U * U * GTVptr[2] * GTVptr[1] * std::conj(GVTptr[1]);
      STVptr[3] += U * U * GTVptr[3] * GTVptr[0] * std::conj(GVTptr[0]);

      STVptr[0] -= U * U * GTVptr[1] * std::conj(GVTptr[3]) * GTVptr[2];
      STVptr[1] -= U * U * GTVptr[0] * std::conj(GVTptr[2]) * GTVptr[3];
      STVptr[2] -= U * U * GTVptr[3] * std::conj(GVTptr[1]) * GTVptr[0];
      STVptr[3] -= U * U * GTVptr[2] * std::conj(GVTptr[0]) * GTVptr[1];
    }

  }

  void solve_Sigma_Fock(int tstp, double U, GREEN &G, cplx *H) {
    cplx cplxi(0.,1.);

    cplx *les = G.lesptr(tstp, tstp);
    H[tstp*4 + 0] = 0;
    H[tstp*4 + 1] = cplxi * U * les[1];
    H[tstp*4 + 2] = cplxi * U * les[2];
    H[tstp*4 + 3] = 0;
  }

  void solve_2b(int tstp, double U, cplx *t0, GREEN &G, GREEN &Sigma) {
    Sigma.set_tstp_zero(tstp);
    solve_Delta(tstp, t0, G, Sigma);
    solve_Sigma(tstp, U, G, Sigma);
  }
};

} // namespace

#endif // HF_SIGMA_DEFN
