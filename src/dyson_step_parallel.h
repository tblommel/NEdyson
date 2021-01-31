#ifndef DYSON_OMP_IMPL
#define DYSON_OMP_IMPL

namespace NEdyson{

void dyson::dyson_step_omp_ret(int threads, int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double dt) const {
  int l, n, i;
  cplx ncplxi = cplx(0,-1);

  // Set timestep to zero
  memset(G.retptr(tstp,0), 0, (tstp+1)*es_*sizeof(cplx));
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);
  ZMatrixMap QMap = ZMatrixMap(Q.data(), k_*nao_, nao_);
  ZMatrixMap MMap = ZMatrixMap(M.data(), k_*nao_, k_*nao_);
  ZMatrixMap XMap = ZMatrixMap(X.data(), k_*nao_, nao_);

  // Initial condition
  ZMatrixMap(G.retptr(tstp,tstp), nao_, nao_) = ncplxi * IMap;

  // Do first k values
  // Exact same as serial implementation
  for(n=1; n<=k; n++) {
    ZMatrixMap QMapBlock = ZMatrixMap(Q.data() + (n-1)*es_, nao_, nao_);

    for(l=0; l<=k; l++) {
      auto MMapBlock = MMap.block((n-1)*nao_, (l-1)*nao_, nao_, nao_);

      // Derivative terms 
      if(l==0){ // Goes into Q
        QMapBlock.noalias() += ncplxi/dt*I.poly_diff(n,l) * ZMatrixMap(G.retptr(tstp, tstp), nao_, nao_).transpose();
      }
      else{ // Goes into M
        MMapBlock.noalias() -= ncplxi/dt*I.poly_diff(n,l) * IMap;
      }

      // Delta Energy term
      if(l==n){
        MMapBlock.noalias() += mu * IMap - ZMatrixConstMap(hmf+(tstp-n)*es_, nao_, nao_).transpose();
      }

      // Integral term
      if(l==0){ // Goes into Q
        QMapBlock.noalias() += dt*I.gregory_weights(n,l) * (ZMatrixMap(G.retptr(tstp,tstp), nao_, nao_)
                                               * ZMatrixMap(Sig.retptr(tstp,tstp-n), nao_, nao_)).transpose();
      }
      else{ // Goes into M
        if(tstp-l >= tstp-n) { // We have Sig
          MMapBlock.noalias() -= dt*I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(tstp-l,tstp-n), nao_, nao_).transpose();
        }
        else{ // Don't have Sig
          MMapBlock.noalias() += dt*I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(tstp-n,tstp-l), nao_, nao_).conjugate();
        }
      }
    }
  }

  // Solve XM=Q for X
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  XMap = lu.solve(QMap);

  // Put X into G
  for(l=0; l<k_; l++){
    ZMatrixMap(G.retptr(tstp, tstp-l-1), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_).transpose();
  }


  // Now we will use the non-conjugate equation to evaluate the remaining values of G
  // I d/dt GR(t,t') = ...
#pragma omp parallel num_threads(threads)
  {
    int j,m;
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    std::vector<bool> mask_ret(tstp+1, false);

    ZMatrix QQ(nao_, nao_);
    ZMatrix MM(nao_, nao_);

    // Each thread is assigned every nthreads^th timeslice
    for(j=0; j<tstp-k; j++) {
      if(j%nthreads == tid) {
        mask_ret[j] = true;
      }
    }
    incr_convolution_ret(tstp, mask_ret, G, Sig, Sig, G, G, dt);
    
    // Solve
    for(m=0; m<tstp-k; m++) {
      if(mask_ret[m]) {
        QQ.noalias() = ZMatrixMap(G.retptr(tstp, m));
        for(j=1; j<=k+1; j++) Q -= -ncplxi/dt * I.bd_weights(j) * ZMatrixMap(G.retptr(tstp-j, m), nao_, nao_);

        MM.noalias() = (-ncplxi/dt * I.bd_weights(0) + mu) * IMap;
        MM.noalias() -= ZMatrixConstMap(hmf + tstp * es_, nao_, nao_);
        MM.noalias() -= dt * I.omega(0) * ZMatrixMap(Sig.retptr(tstp, tstp), nao_, nao_);

        Eigen::FullPivLU<ZMatrix> lu2(MM);
        ZMatrixMap(G.retptr(tstp, m), nao_, nao_).noalias() = lu2.solve(QQ);
      }
    }
  }
  return;
}


void dyson::dyson_step_omp_tv(int threads, int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  cplx cplxi = cplx(0,1);
  auto IMap = ZMatrixMap(iden.data(), nao_, nao_);

  memset(G.tvptr(tstp, 0), 0, (ntau_+1)*es_*sizeof(cplx));

#pragma omp parallel num_threads(threads)
  {
    int j,i,l;
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    std::vector<bool> mask_tv(ntau+1, false);

    for(i = 0; i <= ntau; i++) {
      if(i%nthreads == tid) {
        mask_tv[i] = true;
      }
    }
    
    // Do integrals
    incr_convolution_tv(tstp, mask_tv, G, Sig, Sig, G, G, I, beta, dt);
    
    // Solve
    ZMatrix MM(nao_, nao_);
    ZMatrix QQ(nao_, nao_);

    MM.noalias() = (cplxi/dt*I.bd_weights(0) + mu) * IMap
                    - ZMatrixConstMap(hmf + tstp*es_, nao_, nao_)
                    - dt*I.omega(0) * ZMatrixMap(Sig.retptr(tstp, tstp), nao_, nao_);
    Eigen::FullPivLU<ZMatrix> lu(MM);

    for(j=0; j<=ntau; j++) {
      if(mask_tv[j]) {
        // integrals
        QQ.noalias() = ZMatrixMap(G.tvptr(tstp,j), nao_, nao_);

        // derivatives
        for(l=1; l<=k+1; l++) {
          QQ.noalias() += -cplxi/dt * I.bd_weights(l) * ZMatrixMap(G.tvptr(tstp-l,j), nao_, nao_);
        }
        
        // solve
        ZMatrixMap(G.tvptr(tstp,j), nao_, nao_).noalias() = lu.solve(QQ);
      }
    }
  }
  return;
}


void dyson::dyson_step_omp_les(int threads, int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  // Parallelization starts only for tstp >= 2*k+1
  if(tstp < 2*k+1){
    dyson_step_les(tstp, G, Sig, hmf, mu, beta, dt);
    return;
  }

  auto IMap = ZMatrixMap(iden.data(), nao_, nao_);

  // Set timestep to zero
  memset(G.lesptr(0,tstp), 0, (tstp+1)*es_*sizeof(cplx));

  // This solves for G^L(j,tstp) j=0...n-k-1.  Rest of the timestep done serially due to unavailablity of derivative approximation points
#pragma omp parallel num_threads(threads)
  {
    int i,j;
    cplx cplxi = cplx(0,1);
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    std::vector<bool> mask_les(tstp+1, false);

    ZMatrix QQ(nao_, nao_);
    ZMatrix MM(nao_, nao_);

    // Each thread gets assigned every nthreads^th timeslice
    for(i=0; i<tstp-k; i++) {
      if(i%nthreads == tid) { 
        mask_les[i] = true;
      }
    }
    incr_convolution_les(tstp, mask_les, G, G, G, Sig, Sig, beta, dt);

    // Solve
    for(j=0; j<tstp-k; j++) {
      if(mask_les[j]) {
        QQ.noalias() = ZMatrixMap(G.lesptr(j,tstp));
        for(i=0; i<=k+1; i++) QQ.noalias() += I.bd_weights(i) * cplxi / dt * ZMatrixMap(G.lesptr(j,tstp-i));

        MM.noalias() = (-cplxi / dt * I.bd_weights(0) + mu) * IMap;
        MM.noalias() -= dt * I.omega(0) * ZMatrixMap(Sig.retptr(tstp, tstp), nao_, nao_).adjoint();
        MM.noalias() -= ZMatrixConstMap(hmf + tstp * es_, nao_, nao_);
        
        Eigen::FullPivLU<ZMatrix> lu2(MM);
        ZMatrixMap(G.lestptr(j, tstp), nao_, nao_).noalias() = lu2.solve(QQ);
      }
    }
  }

  //Parallelize the remaining integrals
#pragma omp parallel num_threads(threads)
  {
    int c;
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    for(c = tstp-k; c<=tstp; c++) {
      if(c%nthreads == tid) {
        element_set_zero(size1,gles+(c-tstp+k)*es);
        memset
        Cles2_tstp(c,c,I,Sig,Sig,G,G,tstp,dt,gles+(c-tstp+k)*es);
        Cles3_tstp(c,c,I,Sig,Sig,G,G,tstp,beta,gles+(c-tstp+k)*es);
      }
    }
  }
  
  int j,l;
  for(j=tstp-k;j<=tstp;j++){
    element_set(size1,Q,gles+(j-tstp+k)*es);
    for(l=1;l<=k+1;l++) element_incr(size1, Q, -cplxi/dt*I.bd_weights(l),G.lesptr(j-l,tstp));
    for(l=0;l<j;l++) element_incr(size1, Q, dt*I.gregory_weights(j,l),Sig.retptr(j,l), G.lesptr(l,tstp));
    element_iden(size1, M, mu+cplxi/dt*I.bd_weights(0));
    element_incr(size1, M, -1, hmf.ptr(j));
    element_incr(size1, M, -dt*I.gregory_weights(j,j), Sig.retptr(j,j));

    element_linsolve_left(size1, size1, size1, M, G.lesptr(j,tstp), Q);
  }

  return;
}

void dyson::dyson_step_omp(int omp_threads, int n, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == Sig.size1());
  assert(G.size1() == nao_);
  assert(G.nt() == Sig.nt());
  assert(G.nt() == nt_);
  assert(G.ntau() == Sig.ntau());
  assert(G.ntau() == ntau_);
  assert(G.nt() >= k_);
  assert(G.sig() == Sig.sig());
  assert(n > k_);
  assert(n <= G.nt());

  int threads = (omp_threads == -1 ? omp_get_max_threads() : omp_threads);
  dyson_step_omp_ret(threads, n, I, G, Sig, hmf, mu, dt);
  dyson_step_omp_tv(threads, n, I, G, Sig, hmf, mu, beta, dt);
  dyson_step_omp_les(threads, n, I, G, Sig, hmf, mu, beta, dt);
}

void dyson::dyson_step_omp(int omp_threads, int n, GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[2]);
  assert(G.size1() == hmf.shape()[1]);
  assert(G.size1() == (hmf.shape()[0]-1));

  dyson_step_omp(omp_threads, n, G, Sig, hmf.data(), mu, beta, dt);
}
} // namespace


#endif // DYSON_OMP_IMPL
