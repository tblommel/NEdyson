#ifndef DYSON_BOOT_IMPL
#define DYSON_BOOT_IMPL

namespace NEdyson {

double dyson::dyson_start_ret(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double dt) const {

  double err = 0;
/*  
  // Counters
  int m, l, n, i;

  cplx ncplxi = cplx(0, -1);

  ZMatrixMap QMap = ZMatrixMap(Q.data(), nao_*k_, nao_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);

  // Initial condition
  for(i=0; i<=k_; i++){
    ZMatrixMap(G.retptr(i,i), nao_, nao_).noalias() = ncplxi*IMap;
  }

  // Fill GR(:k+1,0)
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*k_, nao_*k_);
  memset(M.data(), 0, k_*k_*es_*sizeof(cplx));
  memset(Q.data(), 0, k_*es_*sizeof(cplx));

  for(n=1; n<=k_; n++) {
    auto QMapBlock = QMap.block((n-1)*nao_, 0, nao_, nao_);

    for(l=0; l<=k_; l++) {
      auto MMapBlock = MMap.block((n-1)*nao_, (l-1)*nao_, nao_, nao_);

      if(l==0){ // We know these G's. Put into Q
        QMapBlock.noalias() += ncplxi/dt * I.poly_diff(n,l) * -1.*ZMatrixMap(G.retptr(0,l), nao_, nao_).adjoint()
                            + dt*I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_) * -1*ZMatrixMap(G.retptr(0,l), nao_, nao_).adjoint();
      }
      else{ // Don't have these. Put into M

        // Derivative term
        MMapBlock.noalias() = -ncplxi/dt * I.poly_diff(n,l) * IMap;

        // Delta energy term
        if(n==l){
          MMapBlock.noalias() += mu*IMap - ZMatrixConstMap(hmf + l*es_, nao_, nao_);
        }

        // Integral term
        if(n>=l){ // We have Sig
          MMapBlock.noalias() -= dt*I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_);
        }
        else{ // Don't have it
          MMapBlock.noalias() += dt*I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(l,n), nao_, nao_).adjoint();
        }
      }
    }
  }

  // Solve MX=Q for X
  Eigen::FullPivLU<ZMatrix> lu(ZMatrixMap(M.data(), k_*nao_, k_*nao_));
  ZMatrixMap(X.data(), k_*nao_, nao_) = lu.solve(ZMatrixMap(Q.data(), k_*nao_, nao_));
  
  // Put X into G
  for(l=0; l<k_; l++){
    err += (ZColVectorMap(G.retptr(l+1,0), es_) - ZColVectorMap(X.data() + l*es_, es_)).norm();
    ZMatrixMap(G.retptr(l+1,0), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_);
  }

  // Enforce TTI for bootstrapping routine
  for(l = 1; l <= k_; l++) {
    for(n = l; n <= k_; n++) {
      err += (ZColVectorMap(G.retptr(n,l), es_) - ZColVectorMap(G.retptr(n-l,0), es_)).norm();
      ZMatrixMap(G.retptr(n,l), nao_, nao_).noalias() = ZMatrixMap(G.retptr(n-l,0), nao_, nao_);
    }
  }
*/  

  for(int n = 0; n <= k_; n++) {
    ZMatrixMap GR = ZMatrixMap(Q.data(), nao_, nao_);
    GR = G.sig() * ZMatrixMap(G.tvptr(n, G.ntau()), nao_, nao_) - ZMatrixMap(G.tvptr(n,0), nao_, nao_);
    if(n==0) GR = cplx(0.,-1.) * ZMatrixMap(iden.data(), nao_, nao_);

    for(int l = 0; l <= k_-n; l++) {
      ZMatrixMap retMap = ZMatrixMap(G.retptr(n+l,l), nao_, nao_);
      err += (GR - retMap).lpNorm<2>();
      retMap = GR;
    }
  }

  return err;
}

double dyson::dyson_start_tv(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  // Counters and sizes
  int m, l, n, i;
  double err=0;

  cplx cplxi = cplx(0,1);
  ZMatrixMap QMap = ZMatrixMap(Q.data(), nao_*k_, nao_);
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*k_, nao_*k_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);

  // Boundary Conditions
  for(m=0; m<=ntau_; m++) {
    auto tvmap = ZMatrixMap(G.tvptr(0,m), nao_, nao_);
    auto matmap = ZMatrixMap(G.matptr(ntau_-m), nao_, nao_);
    err += (tvmap - (double)G.sig()*cplxi*matmap).norm();
    tvmap.noalias() = (double)G.sig()*cplxi*matmap;
  }

  memset(NTauTmp.data(),0,k_*(ntau_+1)*es_*sizeof(cplx));

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
      tv_it_conv(m, n, Sig, G, Q.data() + (n-1)*es_);

      auto QMapBlock = ZMatrixMap(Q.data() + (n-1)*es_, nao_, nao_);

      for(l=0; l<=k_; l++) {
        auto MMapBlock = MMap.block((n-1)*nao_, ((l==0?1:l)-1)*nao_, nao_, nao_);

        // Derivative term
        if(l == 0){ // Put into Q
          QMapBlock.noalias() -= cplxi*I.poly_diff(n,l)/dt * ZMatrixMap(G.tvptr(0,m), nao_, nao_);
        }
        else{ // Put into M
          MMapBlock.noalias() += cplxi*I.poly_diff(n,l)/dt * IMap;
        }

        // Delta energy term
        if(l==n){
          MMapBlock.noalias() += mu*IMap - ZMatrixConstMap(hmf + l*es_, nao_, nao_);
        }

        // Integral term
        if(l==0){ // Put into Q
          QMapBlock.noalias() += dt*I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_) * ZMatrixMap(G.tvptr(l,m), nao_, nao_);
        }
        else{ // Put into M
          if(n>=l){ // Have Sig
            MMapBlock.noalias() -= dt*I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_);
          }
          else{ // Dont have Sig
            MMapBlock.noalias() += dt*I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(l,n), nao_, nao_).adjoint();
          }
        }
      }

      // Add in the integral
//      QMapBlock.noalias() += ZMatrixMap(NTauTmp.data() + (n-1)*(ntau_+1)*nao_*nao_ + m*nao_*nao_, nao_, nao_);
    }

    // Solve MX=Q
    Eigen::FullPivLU<ZMatrix> lu(MMap);
    ZMatrixMap(X.data(), k_*nao_, nao_) = lu.solve(QMap);

    for(l=0; l<k_; l++){
      err += (ZColVectorMap(G.tvptr(l+1,m), es_) - ZColVectorMap(X.data() + l*es_, es_)).norm();
      ZMatrixMap(G.tvptr(l+1,m), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_);
    }
  }

  return err;
}

double dyson::dyson_start_les(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  double err=0;
  cplx cplxi = cplx(0,1); 

  for(int l = 0; l <= k_; l++) {
    for(int n = 0; n <= l; n++) {
      err += (ZMatrixMap(G.lesptr(n,l), nao_, nao_) + ZMatrixMap(G.tvptr(l-n,0), nao_, nao_).adjoint()).norm();
      ZMatrixMap(G.lesptr(n,l), nao_, nao_) = -ZMatrixMap(G.tvptr(l-n,0), nao_, nao_).adjoint();
    }
  }
  return err;
}


double dyson::dyson_start(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == nao_);
  assert(G.nt() == nt_);
  assert(G.nt() >= k_);
  assert(G.ntau() == ntau_);

  double err=0;
  if(mode_ == gfmol::Mode::GF2) {
    assert(G.ntau() == Sig.ntau());
    assert(G.nt() == Sig.nt());
    assert(G.size1() == Sig.size1());
    assert(G.sig() == Sig.sig());

    err += dyson_start_tv(G, Sig, hmf, mu, beta, dt);
    std::cout << err << std::endl;
    err += dyson_start_ret(G, Sig, hmf, mu, dt);
    std::cout << err << std::endl;
    err += dyson_start_les(G, Sig, hmf, mu, beta, dt);
    std::cout << err << std::endl;
  }
  else {
    err += dyson_start_tv_hf(G, hmf, mu, beta, dt);
    err += dyson_start_ret_hf(G, hmf, mu, dt);
    err += dyson_start_les_hf(G, hmf, mu, beta, dt);
  }
  return err;
}


double dyson::dyson_start(GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[2]);
  assert(G.size1() == hmf.shape()[1]);
  assert(G.nt() == hmf.shape()[0]-1);

  return dyson_start(G, Sig, hmf.data(), mu, beta, dt);
}

} // namespace NEdyson

#endif // DYSON_BOOT_IMPL
