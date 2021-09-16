#ifndef DYSON_TTI_STEP_HF_IMPL
#define DYSON_TTI_STEP_HF_IMPL

namespace NEdyson {

void dyson::dyson_step_ret_hf(int tstp, TTI_GREEN &G, const double *hmf, double mu, double dt) const {
  int l;

  cplx ncplxi = cplx(0,-1);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);
  ZMatrixMap QMap = ZMatrixMap(Q.data(), nao_, nao_);
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_, nao_);

  memset(M.data(), 0, es_*sizeof(cplx));
  memset(Q.data(), 0, es_*sizeof(cplx));

  // Set up Q
  for(l=1; l<=k_+1; l++){
    QMap.noalias() += I.bd_weights(l)*ncplxi/dt * ZMatrixMap(G.retptr(tstp-l), nao_, nao_).transpose();
  }

  // Set up mm
  MMap.noalias() = -DMatrixConstMap(hmf, nao_, nao_).transpose();
  MMap.noalias() += (mu - I.bd_weights(0)*ncplxi/dt) * IMap;

  // Solve XM=Q for X
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  ZMatrixMap(G.retptr(tstp), nao_, nao_) = lu.solve(QMap).transpose();
}


void dyson::dyson_step_tv_hf(int tstp, TTI_GREEN &G, const double *hmf, double mu, double beta, double dt) const {
  // Counters and sizes
  int m, l, n, i;

  cplx cplxi = cplx(0.,1.);
  auto IMap = ZMatrixMap(iden.data(), nao_, nao_);
  auto QMap = ZMatrixMap(Q.data(), nao_, nao_);
  auto MMap = ZMatrixMap(M.data(), nao_, nao_);

  memset(G.tvptr(tstp,0),0,(ntau_+1)*es_*sizeof(cplx));

  // Put derivatives into GRM(tstp,m)
  for(l=1; l<=k_+1; l++) {
    auto GTVMap = ZColVectorMap(G.tvptr(tstp,0), (ntau_+1)*es_);
    GTVMap.noalias() += -cplxi/dt*I.bd_weights(l) * ZColVectorMap(G.tvptr(tstp-l,0), (ntau_+1)*es_);
  }

  // Make M
  MMap.noalias() = (cplxi/dt*I.bd_weights(0) + mu) * IMap
                                             - DMatrixConstMap(hmf, nao_, nao_);

  // Solve MX=Q
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  for(m=0; m<=ntau_; m++) {
    QMap.noalias() = ZMatrixMap(G.tvptr(tstp, m), nao_, nao_);
    ZMatrixMap(G.tvptr(tstp, m), nao_, nao_).noalias() = lu.solve(QMap);
  }
}

void dyson::dyson_step_les_hf(int n, TTI_GREEN &G, const double *hmf, double mu, double beta, double dt) const {
  ZMatrixMap(G.lesptr(-n), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();
}

} // namespace NEdyson

#endif // DYSON_TTI_STEP_IMPL
