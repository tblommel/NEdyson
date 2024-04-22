void hubGF2Solver::solve_les(int tstp, GREEN &Sigma, GREEN &G) const {
  int nao2 = nao_ * nao_;
  int nao3 = nao2 * nao_;

  for(int t=0; t<=tstp; ++t) {
    ZMatrixMap(A1_aaa.data(), nao_, nao_) = (ZMatrixMap(G.retptr(tstp,t), nao_, nao_) - ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).adjoint());
    ZMatrixMap(Sigma.lesptr(t,tstp), nao_, nao_) = U_*U_*ZMatrixMap(G.lesptr(t,tstp), nao_, nao_)
                                     .cwiseProduct(ZMatrixMap(G.lesptr(t,tstp), nao_, nao_))
                                     .cwiseProduct(ZMatrixMap(A1_aaa.data(),    nao_, nao_).transpose());
  }
}

void hubGF2Solver::solve_ret(int tstp, GREEN &Sigma, GREEN &G) const {
  ZMatrixMap GGTt = ZMatrixMap(A1_aaa.data(), nao_, nao_);
  ZMatrixMap GGtT = ZMatrixMap(B1_aaa.data(), nao_, nao_);
  ZMatrixMap GLTt = ZMatrixMap(B2_aaa.data(), nao_, nao_);
  ZMatrixMap GLtT = ZMatrixMap(rho_T.data(), nao_, nao_);

  for(int t=0; t<=tstp; ++t) {
    GGTt = ZMatrixMap(G.retptr(tstp,t), nao_, nao_) - ZMatrixMap(G.lesptr(t,tstp), nao_, nao_).adjoint();
    GGtT = -GGTt.adjoint();
    GLtT = ZMatrixMap(G.lesptr(t,tstp), nao_, nao_);
    GLTt = -GLtT.adjoint();

    ZMatrixMap(Sigma.retptr(tstp,t), nao_, nao_) = U_*U_* GGTt
                                            .cwiseProduct(GGTt)
                                            .cwiseProduct(GLtT.transpose());
    ZMatrixMap(Sigma.retptr(tstp,t), nao_, nao_) -= U_*U_* GLTt
                                             .cwiseProduct(GLTt)
                                             .cwiseProduct(GGtT.transpose());
  }
}


void hubGF2Solver::solve(int tstp, GREEN &Sigma, GREEN &G) const {
  assert(G.sig() == Sigma.sig());

  assert(G.ntau() == Sigma.ntau());

  assert(tstp <= G.nt());
  assert(tstp <= Sigma.nt());

  assert(G.size1() == nao_);
  assert(Sigma.size1() == nao_);

  // Set self energy to be zero
  Sigma.set_tstp_zero(tstp);
  
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;

  // Perform contractions
  start = std::chrono::system_clock::now();
  solve_les(tstp, Sigma, G);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;

  start = std::chrono::system_clock::now();
  solve_ret(tstp, Sigma, G);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
}

