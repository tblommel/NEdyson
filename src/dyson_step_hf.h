#ifndef DYSON_STEP_HF_IMPL
#define DYSON_STEP_HF_IMPL

namespace NEdyson{

// i dt' GR(t,t-t') - GR(t-t')hmf(t-t') - \int_0^t' GR(t,t-s) SR(t-s,t-t') = 0
void dyson_hf::dyson_step_ret(int tstp, GREEN &G, const cplx *hmf, double mu, double dt) const {
  int m, l, n, i;

  std::chrono::time_point<std::chrono::system_clock> intstart, intend, start, end;
  std::chrono::duration<double> elapsed_seconds, inttime;
  start = std::chrono::system_clock::now();

  cplx ncplxi = cplx(0,-1);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);
  ZMatrixMap QMap = ZMatrixMap(Q.data(), k_*nao_, nao_);
  ZMatrixMap MMap = ZMatrixMap(M.data(), k_*nao_, k_*nao_);
  ZMatrixMap XMap = ZMatrixMap(X.data(), k_*nao_, nao_);

  memset(M.data(), 0, k_*k_*es_*sizeof(cplx));
  memset(Q.data(), 0, (tstp+1)*es_*sizeof(cplx));

  start = std::chrono::system_clock::now();
            
  // Initial condition
  ZMatrixMap(G.retptr(tstp,tstp), nao_, nao_) = ncplxi * IMap;
        
  // Fill the first k timesteps 
  for(n=1; n<=k_; n++) {
    ZMatrixMap QMapBlock = ZMatrixMap(Q.data() + (n-1)*es_, nao_, nao_);
        
    for(l=0; l<=k_; l++) {
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
    }
  }

  // Solve XM=Q for X
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  XMap = lu.solve(QMap);

  // Put X into G
  for(l=0; l<k_; l++){
    ZMatrixMap(G.retptr(tstp, tstp-l-1), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_).transpose();
  }

  // Do k+1 ... tstp
  ZMatrixMap MMapSmall = ZMatrixMap(M.data(), nao_, nao_);
  for(n=k_+1; n<=tstp; n++) {
    ZMatrixMap QMapBlock = ZMatrixMap(Q.data() + n*es_, nao_, nao_);
    QMapBlock *= dt;
    // derivative info goes into Q
    for(l=1; l<=k_+1; l++){
      QMapBlock.noalias() += I.bd_weights(l)*ncplxi/dt * ZMatrixMap(G.retptr(tstp,tstp-n+l), nao_, nao_).transpose();
    }

    // Set up M
    MMapSmall.noalias() = -ZMatrixConstMap(hmf+(tstp-n)*es_, nao_, nao_).transpose();
    MMapSmall.noalias() += (mu - I.bd_weights(0)*ncplxi/dt)*IMap;

    // Solve XM=Q for X
    Eigen::FullPivLU<ZMatrix> lu2(MMapSmall);
    ZMatrixMap(G.retptr(tstp,tstp-n), nao_, nao_).noalias() = lu2.solve(QMapBlock).transpose();
  }

  // Output timing information
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  double runtime = elapsed_seconds.count();
  double intruntime = inttime.count();

  std::ofstream out;
  std::string data_dir = std::string(DATA_DIR);
  out.open(data_dir + "dystiming.dat" + "," + std::to_string(G.size1()) + "," + std::to_string(G.nt()) + "," + std::to_string(G.ntau()), std::ofstream::app);
  out << intruntime << " " << runtime << " ";
}

// idt GRM(tstp,m) - hmf(tstp) GRM(t,m) - \int_0^t dT SR(t,T) GRM(T,m) = \int_0^{beta} dT SRM(t,T) GM(T-m)
void dyson_hf::dyson_step_tv(int tstp, GREEN &G, const cplx *hmf, double mu, double beta, double dt) const {
  int m, l, n, i;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;

  cplx cplxi = cplx(0.,1.);
  auto IMap = ZMatrixMap(iden.data(), nao_, nao_);
  auto QMap = ZMatrixMap(Q.data(), nao_, nao_);
  auto MMap = ZMatrixMap(M.data(), nao_, nao_);

  memset(G.tvptr(tstp,0),0,(ntau_+1)*es_*sizeof(cplx));

  start = std::chrono::system_clock::now();

  // Put derivatives into GRM(tstp,m)
  for(l=1; l<=k_+1; l++) {
    auto GTVMap = ZColVectorMap(G.tvptr(tstp,0), (ntau_+1)*es_);
    GTVMap.noalias() += -cplxi/dt*I.bd_weights(l) * ZColVectorMap(G.tvptr(tstp-l,0), (ntau_+1)*es_);
  }

  // Make M
  MMap.noalias() = (cplxi/dt*I.bd_weights(0) + mu) * IMap
                                             - ZMatrixConstMap(hmf+tstp*es_, nao_, nao_);

  Eigen::FullPivLU<ZMatrix> lu(MMap);
  // Solve MX=Q
  for(m=0; m<=ntau_; m++) {
    QMap.noalias() = ZMatrixMap(G.tvptr(tstp, m), nao_, nao_);
    ZMatrixMap(G.tvptr(tstp,m), nao_, nao_).noalias() = lu.solve(QMap);
  }

  // Output timing information
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  double runtime = elapsed_seconds.count();

  std::ofstream out;
  std::string data_dir = std::string(DATA_DIR);
  out.open(data_dir + "dystiming.dat" + "," + std::to_string(G.size1()) + "," + std::to_string(G.nt()) + "," + std::to_string(G.ntau()), std::ofstream::app);
  out<<runtime<<" ";
}

double dyson_hf::dyson_step_les(int n, GREEN &G, const cplx *hmf, double mu, double beta, double dt) const {
  // iterators
  int m, l, i;
  int num = n>=k_ ? n : k_;
  double err=0;

  std::chrono::time_point<std::chrono::system_clock> start, end, intstart, intend;
  std::chrono::duration<double> elapsed_seconds, int1, int2, int3;

  // Matricies
  cplx cplxi = cplx(0,1);
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*k_, nao_*k_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);

  start = std::chrono::system_clock::now();

  // Initial condition
  err += (ZMatrixMap(G.lesptr(0,n), nao_, nao_) + ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint()).lpNorm<2>();
  ZMatrixMap(G.lesptr(0,n), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();

  // Integrals go into Q
  memset(Q.data(),0,sizeof(cplx)*(num+1)*es_);
  memset(M.data(),0,sizeof(cplx)*k_*k_*es_);

  // Set up the kxk linear problem MX=Q
  for(m=1; m<=k_; m++) {
    auto QMapBlock = ZMatrixMap(Q.data() + m*es_, nao_, nao_);

    for(l=0; l<=k_; l++) {
      auto MMapBlock = MMap.block((m-1)*nao_, (l-1)*nao_, nao_, nao_);

      // Derivative term
      if(l==0){ // We put this in Q
        QMapBlock.noalias() -= cplxi/dt*I.poly_diff(m,l) * ZMatrixMap(G.lesptr(l,n), nao_, nao_);
      }
      else{ // It goes into M
        MMapBlock.noalias() += cplxi/dt*I.poly_diff(m,l) * IMap;
      }

      // Delta energy term
      if(m==l){
        MMapBlock.noalias() += mu*IMap - ZMatrixConstMap(hmf+l*es_, nao_, nao_);
      }
    }
  }

  // Solve Mx=Q
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  ZMatrixMap(X.data() + es_, k_*nao_, nao_).noalias() = lu.solve(ZMatrixMap(Q.data()+es_, k_*nao_, nao_));

  ZMatrixMap(X.data(), nao_, nao_).noalias() = ZMatrixMap(G.lesptr(0,n), nao_, nao_);

  // Timestepping
  ZMatrixMap MMapSmall = ZMatrixMap(M.data(), nao_, nao_);
  for(m=k_+1; m<=n; m++) {
    auto QMapBlock = ZMatrixMap(Q.data() + m*es_, nao_, nao_);

    // Set up M
    MMapSmall.noalias() = -ZMatrixConstMap(hmf+m*es_, nao_, nao_) + (cplxi/dt*I.bd_weights(0) + mu)*IMap;

    // Derivatives into Q
    for(l=1; l<=k_+1; l++) {
      QMapBlock.noalias() -= cplxi/dt*I.bd_weights(l) * ZMatrixMap(X.data() + (m-l)*es_, nao_, nao_);
    }

    // Solve MX=Q
    Eigen::FullPivLU<ZMatrix> lu2(MMapSmall);
    ZMatrixMap(X.data()+m*es_, nao_, nao_) = lu2.solve(ZMatrixMap(Q.data() + m*es_, nao_, nao_));
  }

  // Write elements into G
  for(l=1; l<=n; l++) {
    err += (ZColVectorMap(G.lesptr(l,n), es_) - ZColVectorMap(X.data() + l*es_, es_)).norm();
    ZMatrixMap(G.lesptr(l,n), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_);
  }

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  double runtime = elapsed_seconds.count();
  if(n>k_){
    std::ofstream out;
    std::string data_dir = std::string(DATA_DIR);
    out.open(data_dir + "dystiming.dat" + "," + std::to_string(G.size1()) + "," + std::to_string(G.nt()) + "," + std::to_string(G.ntau()), std::ofstream::app);
    out << int1.count() << " " << int2.count() << " " << int3.count() << " " << runtime << " " << std::endl;
  }
  return err;
}


void dyson_hf::dyson_step(int n, GREEN &G, const cplx *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == nao_);
  assert(G.nt() == nt_);
  assert(G.ntau() == ntau_);
  assert(G.nt() > k_);
  assert(n > k_);
  assert(n <= G.nt());

  dyson_step_ret(n, G, hmf, mu, dt);
  dyson_step_tv(n, G, hmf, mu, beta, dt);
  dyson_step_les(n, G, hmf, mu, beta, dt);
}



void dyson::dyson_step(int n, GREEN &G, const ZTensor<3> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[2]);
  assert(G.size1() == hmf.shape()[1]);
  assert(G.nt() == (hmf.shape()[0]-1));

  dyson_step(n, G, hmf.data(), mu, beta, dt);
}


} // namespace
#endif // DYSON_STEP_IMPL