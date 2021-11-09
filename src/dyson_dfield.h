#ifndef DYSON_DFIELD_IMPL
#define DYSON_DFIELD_IMPL


namespace NEdyson{

void dyson::dipole_field(int tstp, ZTensor<2> &dfield, const GREEN &Gu, const GREEN &Gd, const DTensor<3> &dipole, double l, double n, double dt) const {
  assert(Gu.nao() == dipole.shape()[1]);
  assert(Gu.nao() == dipole.shape()[2]);
  assert(dfield.shape()[0] == dipole.shape()[0]);
  assert(tstp <= Gu.nt());
  assert(tstp <= dfield.shape()[1]);

  if( tstp <= k_) { // Need to use poly_diff coeff
    for(int d = 0; d < 3; d++) {
      dfield(d, tstp) = 0;

      for(int i = 0; i <= k_; i++) {
        cplx yi = DMatrixConstMap(dipole.data() + d * es_, nao_, nao_).cwiseProduct(
                      ZMatrixConstMap(Gu.lesptr(i,i), nao_, nao_).transpose() 
                    + ZMatrixConstMap(Gd.lesptr(i,i), nao_, nao_).transpose()).sum();
        dfield(d, tstp) += 1./dt * I.poly_diff(tstp, i) * yi;
      }
      dfield(d, tstp) *= cplx(0., Gu.sig() * 1.) * l * n * 2. * PI / 137.035999206;
    }
  }
  else { // Use bd_weights
    for(int d = 0; d < 3; d++) {
      dfield(d,tstp) = 0;
      for(int i = 0; i <= k_+1; i++) {
        cplx yi = DMatrixConstMap(dipole.data() + d * es_, nao_, nao_).cwiseProduct(
                      ZMatrixConstMap(Gu.lesptr(tstp-i,tstp-i), nao_, nao_).transpose() 
                    + ZMatrixConstMap(Gd.lesptr(tstp-i,tstp-i), nao_, nao_).transpose()).sum();
        dfield(d, tstp) += 1./dt * I.bd_weights(i) * yi;
      }
      dfield(d, tstp) *= cplx(0.,Gu.sig() * 1.) * l * n * 2. * PI / 137.035999206;
    }
  }
}


} // Namespace

#endif
