#ifndef DYSON_BOOT_HF_IMPL
#define DYSON_BOOT_HF_IMPL

namespace NEdyson {

double dyson::dyson_start_ret_hf(GREEN &G, const cplx *hmf, double mu, double dt) const {
  // Counters
  int m, l, n, i;

  double err = 0;
  cplx ncplxi = cplx(0, -1);

  ZMatrixMap QMap = ZMatrixMap(Q.data(), nao_*k_, nao_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);

  // Initial condition
  for(i=0; i<=k_; i++){
    ZMatrixMap(G.retptr(i,i), nao_, nao_).noalias() = ncplxi*IMap;
  }

  // Fill the first k timesteps
  for(m=0; m<k_; m++) {
    ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*(k_-m), nao_*(k_-m));
    memset(M.data(), 0, k_*k_*es_*sizeof(cplx));
    memset(Q.data(), 0, k_*es_*sizeof(cplx));

    for(n=m+1; n<=k_; n++) {
      auto QMapBlock = QMap.block((n-m-1)*nao_, 0, nao_, nao_);

      for(l=0; l<=k_; l++) {
        auto MMapBlock = MMap.block((n-m-1)*nao_, (l-m-1)*nao_, nao_, nao_);

        if(l<=m){ // We know these G's. Put into Q
          QMapBlock.noalias() += ncplxi/dt * I.poly_diff(n,l) * -1.*ZMatrixMap(G.retptr(m,l), nao_, nao_).adjoint();
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
    Eigen::FullPivLU<ZMatrix> lu(ZMatrixMap(M.data(), (k_-m)*nao_, (k_-m)*nao_));
    ZMatrixMap(X.data(), (k_-m)*nao_, nao_) = lu.solve(ZMatrixMap(Q.data(), (k_-m)*nao_, nao_));

    //put X into G
    for(l=0; l<k_-m; l++){
      err += (ZColVectorMap(G.retptr(l+m+1,m), es_) - ZColVectorMap(X.data() + l*es_, es_)).norm();
      ZMatrixMap(G.retptr(l+m+1,m), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_);
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
  for(int n=0; n<=k_; n++) err += dyson_step_les_hf(n,G,hmf,mu,beta,dt);
  return err;
}

} // namespace NEdyson

#endif // DYSON_BOOT_IMPL
