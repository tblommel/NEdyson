#ifndef DYSON_TTI_BOOT_IMPL
#define DYSON_TTI_BOOT_IMPL

namespace NEdyson {

double dyson::dyson_start_ret(TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double dt) const {
  double err = 0;
  ZMatrixMap QMap = ZMatrixMap(Q.data(), nao_, nao_);

  for(int n=0; n<=k_; n++) {
    QMap = G.sig() * ZMatrixMap(G.tvptr(n, G.ntau()), nao_, nao_) - ZMatrixMap(G.tvptr(n, 0), nao_, nao_);
    ZMatrixMap retMap = ZMatrixMap(G.retptr(n), nao_, nao_);

    err += (QMap - retMap).lpNorm<2>();
    retMap.noalias() = QMap;
  }

  return err;
}


double dyson::dyson_start_tv(TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const {
  // Counters and sizes
  int m, l, n, i;
  double err=0;

  cplx cplxi = cplx(0,1);
  ZMatrixMap QMap = ZMatrixMap(Q.data(), nao_*k_, nao_);
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*k_, nao_*k_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);

  // Boundary Conditions
  for(m=0; m<=ntau_; m++) {
    auto tvmap = ZColVectorMap(G.tvptr(0,m), es_);
    auto matmap = ZColVectorMap(G.matptr(ntau_-m), es_);
    err += (tvmap - (double)G.sig()*cplxi*matmap).norm();
    tvmap.noalias() = (double)G.sig()*cplxi*matmap;
  }


  // Do the integrals
  for(n=1; n<=k_; n++) {
    auto ZTVNTT = ZTensorView<3>(NTauTmp.data() + (n-1)*(ntau_+1)*nao_*nao_, ntau_+1, nao_, nao_);
    Conv.mixing(ZTVNTT,
                ZTensorView<3>(Sig.tvptr(n,0), ntau_+1, nao_, nao_),
                ZTensorView<3>(G.matptr(0), ntau_+1, nao_, nao_),
                beta, (double)G.sig());
  }

  // At each m, get n=1...k
  for(m=0; m<=ntau_; m++) {
    memset(M.data(),0,k_*k_*es_*sizeof(cplx));
    memset(Q.data(),0,k_*es_*sizeof(cplx));

    // Set up the kxk linear problem MX=Q
    for(n=1; n<=k_; n++) {
      auto QMapBlock = ZMatrixMap(Q.data() + (n-1)*es_, nao_, nao_);

      for(l=0; l<=k_; l++) {
        auto MMapBlock = MMap.block((n-1)*nao_, (l-1)*nao_, nao_, nao_);

        // Derivative term
        if(l == 0){ // Put into Q
          QMapBlock.noalias() -= cplxi*I.poly_diff(n,l)/dt * ZMatrixMap(G.tvptr(0,m), nao_, nao_);
        }
        else{ // Put into M
          MMapBlock.noalias() += cplxi*I.poly_diff(n,l)/dt * IMap;
        }

        // Delta energy term
        if(l==n){
          MMapBlock.noalias() += mu*IMap - DMatrixConstMap(hmf, nao_, nao_);
        }

        // Integral term
        if(l==0){ // Put into Q
          QMapBlock.noalias() += dt*I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(n-l), nao_, nao_) * ZMatrixMap(G.tvptr(l,m), nao_, nao_);
        }
        else{ // Put into M
          if(n>=l){ // Have Sig
            MMapBlock.noalias() -= dt*I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(n-l), nao_, nao_);
          }
          else{ // Dont have Sig
            MMapBlock.noalias() += dt*I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(l-n), nao_, nao_).adjoint();
          }
        }
      }

      // Add in the integral
      QMapBlock.noalias() += ZMatrixMap(NTauTmp.data() + (n-1)*(ntau_+1)*nao_*nao_ + m*nao_*nao_, nao_, nao_);
    }

    // Solve MX=Q
    Eigen::FullPivLU<ZMatrix> lu(MMap);
    ZMatrixMap(X.data(), k_*nao_, nao_).noalias() = lu.solve(QMap);

    for(l=0; l<k_; l++){
      err += (ZColVectorMap(G.tvptr(l+1,m), es_) - ZColVectorMap(X.data() + l*es_, es_)).norm();
      ZMatrixMap(G.tvptr(l+1,m), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_);
    }
  }

  return err;
} 

double dyson::dyson_start_les(TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double beta, double dt) const {
  double err = 0;
  for(int n=0; n<=k_; n++) {
    ZMatrixMap tvMap = ZMatrixMap(G.tvptr(n,0), nao_, nao_);
    ZMatrixMap lesMap= ZMatrixMap(G.lesptr(-n), nao_, nao_);
    err += (-tvMap.adjoint() - lesMap).lpNorm<2>();
    lesMap.noalias() = -tvMap.adjoint();
  }
  return err;
}

double dyson::dyson_start(TTI_GREEN &G, const TTI_GREEN &Sig, const double* hmf, double mu, double beta, double dt) const {
  assert(G.size1() == nao_);
  assert(G.nt() == nt_);
  assert(G.nt() >= k_);
  assert(G.ntau() == ntau_);

  double err=0;
  if(mode_ == gfmol::Mode::GF2) {
    assert(G.sig() == Sig.sig());
    assert(G.ntau() == Sig.ntau());
    assert(G.nt() == Sig.nt());
    assert(G.size1() == Sig.size1());
    err += dyson_start_tv(G, Sig, hmf, mu, beta, dt);
    err += dyson_start_ret(G, Sig, hmf, mu, dt);
    err += dyson_start_les(G, Sig, hmf, mu, beta, dt);
  }
  else {
    err += dyson_start_tv_hf(G, hmf, mu, beta, dt);
    err += dyson_start_ret_hf(G, hmf, mu, dt);
    err += dyson_start_les_hf(G, hmf, mu, beta, dt);
  }
  return err;
}

double dyson::dyson_start(TTI_GREEN &G, const TTI_GREEN &Sig, const DTensor<2> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[1]);
  assert(G.size1() == hmf.shape()[0]);

  double err=0;
  err += dyson_start(G, Sig, hmf.data(), mu, beta, dt);
  return err;
}






} // namespace NEdyson

#endif // DYSON_TTI_BOOT_IMPL
