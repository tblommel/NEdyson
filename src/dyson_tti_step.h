#ifndef DYSON_TTI_STEP_IMPL
#define DYSON_TTI_STEP_IMPL

namespace NEdyson {

void dyson::dyson_step_ret(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double dt) const {
  ZMatrixMap(G.retptr(n), nao_, nao_).noalias() = G.sig() * ZMatrixMap(G.tvptr(n, G.ntau()), nao_, nao_) - ZMatrixMap(G.tvptr(n, 0), nao_, nao_);
}

void dyson::dyson_step_tv(int tstp, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const {
  // Counters and sizes
  int m, l, n, i;

  cplx cplxi = cplx(0.,1.);
  auto IMap = ZMatrixMap(iden.data(), nao_, nao_);
  auto QMap = ZMatrixMap(Q.data(), nao_, nao_);
  auto MMap = ZMatrixMap(M.data(), nao_, nao_);

  memset(G.tvptr(tstp,0),0,(ntau_+1)*es_*sizeof(cplx));

  // Do integrals
  Ctv_tstp(tstp, G, Sig, Sig, G, G, beta, dt);

  // Put derivatives into GRM(tstp,m)
  for(l=1; l<=k_+1; l++) {
    auto GTVMap = ZColVectorMap(G.tvptr(tstp,0), (ntau_+1)*es_);
    GTVMap.noalias() += -cplxi/dt*I.bd_weights(l) * ZColVectorMap(G.tvptr(tstp-l,0), (ntau_+1)*es_);
  }

  // Make M
  MMap.noalias() = (cplxi/dt*I.bd_weights(0) + mu) * IMap
                                             - DMatrixConstMap(hmf, nao_, nao_)
                                             - dt*I.omega(0) * ZMatrixMap(Sig.retptr(0), nao_, nao_);
  // Solve MX=Q
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  for(m=0; m<=ntau_; m++) {
    QMap.noalias() = ZMatrixMap(G.tvptr(tstp, m), nao_, nao_);
    ZMatrixMap(G.tvptr(tstp, m), nao_, nao_).noalias() = lu.solve(QMap);
  }
}

void dyson::dyson_step_les(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const {
  ZMatrixMap(G.lesptr(-n), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();
}

void dyson::dyson_step(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == nao_);
  assert(G.nt() == nt_);
  assert(G.ntau() == ntau_);
  assert(G.nt() > k_);
  assert(n > k_);
  assert(n <= G.nt());

  if(mode_ == gfmol::Mode::GF2) {
    assert(G.sig() == Sig.sig());
    assert(G.ntau() == Sig.ntau());
    assert(G.nt() == Sig.nt());
    assert(G.size1() == Sig.size1());

    dyson_step_tv(n, G, Sig, hmf, mu, beta, dt);
    dyson_step_ret(n, G, Sig, hmf, mu, dt);
    dyson_step_les(n, G, Sig, hmf, mu, beta, dt);
  }
  else {
    dyson_step_tv_hf(n, G, hmf, mu, beta, dt);
    dyson_step_ret_hf(n, G, hmf, mu, dt);
    dyson_step_les_hf(n, G, hmf, mu, beta, dt); 
  }
}

void dyson::dyson_step(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const DTensor<2> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[0]);
  assert(G.size1() == hmf.shape()[1]);

  dyson_step(n, G, Sig, hmf.data(), mu, beta, dt);
}

} // namespace NEdyson

#endif // DYSON_TTI_STEP_IMPL
