#ifndef DYSON_BOOT_HF_IMPL
#define DYSON_BOOT_HF_IMPL

namespace NEdyson {

double dyson::dyson_start_ret_hf(GREEN &G, const cplx *hmf, double mu, double dt) const {

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

  // Fill the first k timesteps
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*(k_), nao_*(k_));
  memset(M.data(), 0, k_*k_*es_*sizeof(cplx));
  memset(Q.data(), 0, k_*es_*sizeof(cplx));

  for(n=1; n<=k_; n++) {
    auto QMapBlock = QMap.block((n-1)*nao_, 0, nao_, nao_);

    for(l=0; l<=k_; l++) {
      auto MMapBlock = MMap.block((n-1)*nao_, (l-1)*nao_, nao_, nao_);

      if(l<=0){ // We know these G's. Put into Q
        QMapBlock.noalias() += ncplxi/dt * I.poly_diff(n,l) * -1.*ZMatrixMap(G.retptr(0,l), nao_, nao_).adjoint();
      }
      else{ // Don't have these. Put into M

        // Derivative term
        MMapBlock.noalias() = -ncplxi/dt * I.poly_diff(n,l) * IMap;

        // Delta energy term
        if(n==l){
          MMapBlock.noalias() += mu*IMap - ZMatrixConstMap(hmf + l*es_, nao_, nao_);
        }
      }
    }
  }

  //solve MX=Q for X
  Eigen::FullPivLU<ZMatrix> lu(ZMatrixMap(M.data(), (k_)*nao_, (k_)*nao_));
  ZMatrixMap(X.data(), (k_)*nao_, nao_) = lu.solve(ZMatrixMap(Q.data(), (k_)*nao_, nao_));

  //put X into G
  for(l=0; l<k_; l++){
    err += (ZColVectorMap(G.retptr(l+1,0), es_) - ZColVectorMap(X.data() + l*es_, es_)).norm();
    ZMatrixMap(G.retptr(l+1,0), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_);
  }

  // Enforce TTI for bootstrapping routine
  for(l = 1; l <= k_; l++) {
    for(n = l; n <= k_; n++) {
      err += (ZColVectorMap(G.retptr(n,l), es_) - ZColVectorMap(G.retptr(n-l,0), es_)).norm();
      ZMatrixMap(G.retptr(n,l), nao_, nao_) = ZMatrixMap(G.retptr(n-l,0), nao_, nao_);
    }
  }*/

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

double dyson::dyson_start_tv_hf(GREEN &G, const cplx *hmf, double mu, double beta, double dt) const {
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
          MMapBlock.noalias() += mu*IMap - ZMatrixConstMap(hmf + l*es_, nao_, nao_);
        }
      }
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

double dyson::dyson_start_les_hf(GREEN &G, const cplx *hmf, double mu, double beta, double dt) const {
  double err=0;
  for(int l = 0; l <= k_; l++) {
    for(int n = 0; n <= l; n++) {
      err += (ZMatrixMap(G.lesptr(n,l), nao_, nao_) + ZMatrixMap(G.tvptr(l-n,0), nao_, nao_).adjoint()).norm();
      ZMatrixMap(G.lesptr(n,l), nao_, nao_) = -ZMatrixMap(G.tvptr(l-n,0), nao_, nao_).adjoint();
    }
  }
  return err;
}

} // namespace NEdyson

#endif // DYSON_BOOT_IMPL
