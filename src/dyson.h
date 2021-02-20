#ifndef DYSON_DECL
#define DYSON_DECL

#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>

#include "greens.h"
#include "integration.h"
#include "utils.h"
#include "cpsitop/nonequilibrium.h"
#include "chrono"

namespace NEdyson{


class dyson {
private:
  int nt_;
  int ntau_;
  int nao_;
  int es_;
  int k_;
  bool hfbool_;

  mutable ZTensor<1> ex_weights;
  mutable ZTensor<1> tmp;
  mutable ZTensor<1> tmp2;
  mutable ZTensor<1> iden;
  mutable ZTensor<1> M;
  mutable ZTensor<1> Q;
  mutable ZTensor<1> X;
  mutable ZTensor<1> NTauTmp;

  const INTEG I;

  mutable cpsitop::nonequilibrium::convolution Conv;

public:
  dyson(int nt, int ntau, int nao, int k, bool hfbool_);

  // return convolution class
  const cpsitop::nonequilibrium::convolution Convolution() const { return Conv; }
  
  // Non constant free GF ===============================================
  void G0_from_h0(GREEN &G, double mu, const double *hM, const double *ht, double beta, double dt) const ;

  // Constant free GF ==================================================
  void G0_from_h0(GREEN &G, double mu, const DTensor<2> &H0, double beta, double h) const ;
  void G0_from_h0(GREEN &G, double mu, const double *H0, double beta, double h) const ;
  void G0_from_h0(TTI_GREEN &G, double mu, const DTensor<2> &H0, double beta, double h) const ;
  void G0_from_h0(TTI_GREEN &G, double mu, const double *H0, double beta, double h) const ;

  // Start functions ===================================================
  double dyson_start(GREEN &G, const GREEN &Sig, const cplx* hmf, double mu, double beta, double dt) const ;
  double dyson_start(GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt) const ;

  // Start functions individual components ============================
  double dyson_start_les(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const ;
  double dyson_start_ret(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double dt) const ;
  double dyson_start_tv(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const ;

  // Start functions HF individual components =========================
  double dyson_start_les_hf(GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;
  double dyson_start_ret_hf(GREEN &G, const cplx *hmf, double mu, double dt) const ;
  double dyson_start_tv_hf(GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;

  // TTI start functions ===============================================
  double dyson_start(TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const ;
  double dyson_start(TTI_GREEN &G, const TTI_GREEN &Sig, const DTensor<2> &hmf, double mu, double beta, double dt) const ;

  // TTI start functions individual components =========================
  double dyson_start_les(TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const ;
  double dyson_start_ret(TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double dt) const ;
  double dyson_start_tv(TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const ;

  // TTI start functions HF individual components =====================
  double dyson_start_les_hf(TTI_GREEN &G, const double *hmf, double mu, double beta, double dt) const ;
  double dyson_start_ret_hf(TTI_GREEN &G, const double *hmf, double mu, double dt) const ;
  double dyson_start_tv_hf(TTI_GREEN &G, const double *hmf, double mu, double beta, double dt) const ;

  // Step functions ===================================================
  void dyson_step(int n, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const ;
  void dyson_step(int n, GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt) const ;

  // Step functions individual components ============================
  double dyson_step_les(int n, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const ;
  void dyson_step_tv(int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const ;
  void dyson_step_ret(int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double dt) const ;

  // Step functions HF individual components ==========================
  double dyson_step_les_hf(int n, GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;
  void dyson_step_tv_hf(int tstp, GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;
  void dyson_step_ret_hf(int tstp, GREEN &G, const cplx *hmf, double mu, double dt) const ;

  // TTI step functions ==============================================
  void dyson_step(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const ;
  void dyson_step(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const DTensor<2> &hmf, double mu, double beta, double dt) const ;

  // TTI step functions individual components =======================
  void dyson_step_les(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const ;
  void dyson_step_tv(int tstp, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const ;
  void dyson_step_ret(int tstp, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double dt) const ;

  // TTI step functions HF individual components ====================
  void dyson_step_les_hf(int n, TTI_GREEN &G, const double *hmf, double mu, double beta, double dt) const ;
  void dyson_step_tv_hf(int tstp, TTI_GREEN &G, const double *hmf, double mu, double beta, double dt) const ;
  void dyson_step_ret_hf(int tstp, TTI_GREEN &G, const double *hmf, double mu, double dt) const ;

  // Extrapolation ==================================================
  void Extrapolate(int n, GREEN &G) const;
  void Extrapolate(int n, TTI_GREEN &G) const;

  // TV convolutions ================================================
  void Ctv_tstp(int tstp, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, double beta, double dt) const ;
  void CTV1(cplx *ctv, const GREEN &A, const GREEN &Acc, const GREEN &B, int n, double dt) const ;
  void CTV2(const GREEN &A, const GREEN &B, int n, int m, double beta, cplx *res) const ;
  void CTV3(const GREEN &A, const GREEN &B, int n, int m, double beta, cplx *res) const ;

  // TTI TV convolutions ============================================
  void Ctv_tstp(int tstp, TTI_GREEN &C, const TTI_GREEN &A, const TTI_GREEN &Acc, const TTI_GREEN &B, const TTI_GREEN &Bcc, double beta, double dt) const ;
  void CTV1(cplx *ctv, const TTI_GREEN &A, const TTI_GREEN &Acc, const TTI_GREEN &B, int n, double dt) const ;
  void CTV2(const TTI_GREEN &A, const TTI_GREEN &B, int n, int m, double beta, cplx *res) const ;
  void CTV3(const TTI_GREEN &A, const TTI_GREEN &B, int n, int m, double beta, cplx *res) const ;

  // Les convolutions ===============================================
  void Cles2_tstp(int j1, int j2, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double dt, cplx *res) const ;
  void Cles2_tstp(const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double dt, cplx *res) const ;
  void Cles3_tstp(int j1, int j2, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double beta, cplx *res) const ;
  void Cles3_tstp(const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double beta, cplx *res) const ;

  // Energy calculation =============================================
  double energy_conv(int tstp, const GREEN &Sig, const GREEN &G, double beta, double dt) const ;
  double energy_conv(int tstp, const TTI_GREEN &Sig, const TTI_GREEN &G, double beta, double dt) const ;

  // Calculate dipole field =========================================
  void dipole_field(int tstp, ZTensor<2> &dfield, const GREEN &Gu, const GREEN &Gd, const DTensor<3> &dipole, double l, double n, double dt) const ;

}; // class dyson


} // namespace NEdyson

#endif // header guard
