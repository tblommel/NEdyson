#ifndef DYSON_TTI_STEP_IMPL
#define DYSON_TTI_STEP_IMPL

namespace NEdyson {

void dyson::dyson_step_ret(int tstp, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double dt) const {
  int l;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
  start = std::chrono::system_clock::now();

  cplx ncplxi = cplx(0,-1);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);
  ZMatrixMap QMap = ZMatrixMap(Q.data(), nao_, nao_);
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_, nao_);

  memset(M.data(), 0, es_*sizeof(cplx));
  memset(Q.data(), 0, es_*sizeof(cplx));

  start = std::chrono::system_clock::now();

  // doing the integral
  // Q = \sum_{l=0}^{n-1} w_{nl} GR(t,t-l) SR(t-l,t-n) for n=k+1...tstp
  for(l=0; l<tstp; l++) {
    QMap.noalias() += I.gregory_weights(tstp,l) * (ZMatrixMap(G.retptr(l), nao_, nao_)
                                        * ZMatrixMap(Sig.retptr(tstp-l), nao_, nao_)).transpose();
  }

  QMap *= dt;
  for(l=1; l<=k_+1; l++){
    QMap.noalias() += I.bd_weights(l)*ncplxi/dt * ZMatrixMap(G.retptr(tstp-l), nao_, nao_).transpose();
  }

  // Set up mm
  MMap.noalias() = -DMatrixConstMap(hmf, nao_, nao_).transpose();
  MMap.noalias() -= dt*I.omega(0) * ZMatrixMap(Sig.retptr(0), nao_, nao_).transpose();
  MMap.noalias() += (mu - I.bd_weights(0)*ncplxi/dt) * IMap;

  // Solve XM=Q for X
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  ZMatrixMap(G.retptr(tstp), nao_, nao_) = lu.solve(QMap).transpose();

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  double runtime = elapsed_seconds.count();

  std::ofstream out;
  std::string data_dir = std::string(DATA_DIR);
  out.open(data_dir + "tti_dystiming.dat" + "," + std::to_string(G.size1()) + "," + std::to_string(G.nt()) + "," + std::to_string(G.ntau()), std::ofstream::app);
  out<<runtime<<" ";
}


void dyson::dyson_step_tv(int tstp, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const {
  // Counters and sizes
  int m, l, n, i;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
  start = std::chrono::system_clock::now();

  cplx cplxi = cplx(0.,1.);
  auto IMap = ZMatrixMap(iden.data(), nao_, nao_);
  auto QMap = ZMatrixMap(Q.data(), nao_, nao_);
  auto MMap = ZMatrixMap(M.data(), nao_, nao_);

  memset(G.tvptr(tstp,0),0,(ntau_+1)*es_*sizeof(cplx));

  start = std::chrono::system_clock::now();

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
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  double runtime = elapsed_seconds.count();

  std::ofstream out;
  std::string data_dir = std::string(DATA_DIR);
  out.open(data_dir + "tti_dystiming.dat" + "," + std::to_string(G.size1()) + "," + std::to_string(G.nt()) + "," + std::to_string(G.ntau()), std::ofstream::app);
  out<<runtime<<std::endl;
}

void dyson::dyson_step_les(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const {
  ZMatrixMap(G.lesptr(-n), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();
}

void dyson::dyson_step(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == Sig.size1());
  assert(G.size1() == nao_);
  assert(G.nt() == Sig.nt());
  assert(G.nt() == nt_());
  assert(G.ntau() == Sig.ntau());
  assert(G.ntau() == ntau_);
  assert(G.nt() > k_);
  assert(G.sig() == Sig.sig());
  assert(n > k_);
  assert(n <= G.nt());

  dyson_step_ret(n, G, Sig, hmf, mu, dt);
  dyson_step_tv(n, G, Sig, hmf, mu, beta, dt);
  dyson_step_les(n, G, Sig, hmf, mu, beta, dt);
}

void dyson::dyson_step(int n, TTI_GREEN &G, const TTI_GREEN &Sig, const DTensor<2> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[0]);
  assert(G.size1() == hmf.shape()[1]);
  assert(G.size1() == Sig.size1());
  assert(G.size1() == nao_);
  assert(G.nt() == Sig.nt());
  assert(G.nt() == nt_());
  assert(G.ntau() == Sig.ntau());
  assert(G.ntau() == ntau_);
  assert(G.nt() > k_);
  assert(G.sig() == Sig.sig());
  assert(n > k_);
  assert(n <= G.nt());

  dyson_step_ret(n, G, Sig, hmf.data(), mu, dt);
  dyson_step_tv(n, G, Sig, hmf.data(), mu, beta, dt);
  dyson_step_les(n, G, Sig, hmf.data(), mu, beta, dt);
}

} // namespace NEdyson

#endif // DYSON_TTI_STEP_IMPL
