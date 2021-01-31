#ifndef DYSON_HF_DECL
#define DYSON_HF_DECL

#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>

#include "greens.h"
#include "integration.h"
#include "utils.h"
#include "cpsitop/nonequilibrium.h"
#include "chrono"

namespace NEdyson{


class dyson_hf {
private:
  int nt_;
  int ntau_;
  int nao_;
  int es_;
  int k_;

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
  dyson_hf(int nt, int ntau, int nao, int k);

  // return convolution class
  const cpsitop::nonequilibrium::convolution Convolution() const { return Conv; }
  
  // Non constant free GF
  void G0_from_h0(GREEN &G, double mu, const double *hM, const double *ht, double beta, double dt) const ;

  // Start functions
  double dyson_start(GREEN &G, const cplx* hmf, double mu, double beta, double dt) const ;
  double dyson_start(GREEN &G, const ZTensor<3> &hmf, double mu, double beta, double dt) const ;
  double dyson_start_les(GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;
  double dyson_start_ret(GREEN &G, const cplx *hmf, double mu, double dt) const ;
  double dyson_start_tv(GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;

  //Step functions
  void dyson_step(int n, GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;
  void dyson_step(int n, GREEN &G, const ZTensor<3> &hmf, double mu, double beta, double dt) const ;
  double dyson_step_les(int n, GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;
  void dyson_step_tv(int tstp, GREEN &G, const cplx *hmf, double mu, double beta, double dt) const ;
  void dyson_step_ret(int tstp, GREEN &G, const cplx *hmf, double mu, double dt) const ;

  // Extrapolation
  void Extrapolate(int n, GREEN &G) const;

  // Energy calculation
  double energy_conv(int tstp, const GREEN &Sig, const GREEN &G, double beta, double dt) const ;

}; // class dyson


} // namespace NEdyson

#endif // header guard
