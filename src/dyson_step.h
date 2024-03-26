#ifndef DYSON_STEP_IMPL
#define DYSON_STEP_IMPL

namespace NEdyson{

// i dt' GR(t,t-t') - GR(t-t')hmf(t-t') - \int_0^t' GR(t,t-s) SR(t-s,t-t') = 0
double dyson::dyson_step_ret(int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double dt) const {

  cplx cplxi = cplx(0,1);
  double err = 0;
  ZMatrixMap(Q.data(), (tstp+1)*nao_*nao_, 1) = ZMatrix::Zero((tstp+1)*nao_*nao_, 1);
  ZMatrixMap IMap(iden.data(), nao_, nao_);

  // Initial condition
  ZMatrixMap(G.retptr(tstp,tstp), nao_, nao_) = -cplxi * IMap;

/*
  for(int tp = 0; tp < tstp; tp++) {
    // M
    double weight = std::min(tp, tstp-k_) == tstp-k_ ? I.poly_integ(tp-tstp+k_, k_, k_) : I.omega(0);
    ZMatrixMap(M.data(), nao_, nao_) = ZMatrixConstMap(hmf+tstp*nao_*nao_, nao_, nao_) - (mu + cplxi/dt*I.bd_weights(0)) * ZMatrixMap(iden.data(), nao_, nao_) + dt * weight * ZMatrixMap(Sig.retptr(tstp, tstp), nao_, nao_);

    ZMatrixMap XMap(X.data() + tp*nao_*nao_, nao_, nao_);
    ZMatrixMap QMap(Q.data() + tp*nao_*nao_, nao_, nao_);

    if(tp >= tstp-k_ ) {
      ZMatrix BD_M(k_+2, k_+2);
      ZMatrix BD_P(k_+2, nao_*nao_);
      ZMatrix BD_I = ZMatrix::Identity(k_+2, k_+2);
      ZMatrix BD_MI(k_+2, k_+2);

      // first rows are for f'
      for(int i = 0; i < tp-tstp+k_+1; i++) {
        for(int j = 0; j <= k_+1; j++) {
          BD_M(i,j) = j==0 ? 0 : j*std::pow(i*dt,j-1);
        }
      }
      // last rows are for f
      for(int i = tp-tstp+k_+1; i <= k_+1; i++) {
        for(int j = 0; j <= k_+1; j++) {
          BD_M(i,j) = std::pow(i*dt,j);
        }
      }
      Eigen::FullPivLU<ZMatrix> BD_lu(BD_M);
      BD_MI = BD_lu.solve(BD_I);

      // Fill in P
      for(int i = 0; i < tp-tstp+k_+1; i++) {
        ZMatrixMap(BD_P.data() + i*nao_*nao_, nao_, nao_) = -cplxi * ZMatrixConstMap(hmf+(tstp-k_-1+i)*nao_*nao_, nao_, nao_) * -ZMatrixMap(G.retptr(tp, tstp-k_-1+i), nao_, nao_).adjoint();
        for(int l = 0; l <= k_; l++) {
          if(tstp-k_-1+i >= tstp-k_-1+l && tstp-k_-1+l >= tp) {
            ZMatrixMap(BD_P.data() + i*nao_*nao_, nao_, nao_) -= cplxi * dt * I.poly_integ(i, k_+1-(tstp-tp), l) * ZMatrixMap(Sig.retptr(tstp-k_-1+i, tstp-k_-1+l), nao_, nao_)
                                                                                                                 * ZMatrixMap(G  .retptr(tstp-k_-1+l, tp)         , nao_, nao_);
          }
          else if(tstp-k_-1+i < tstp-k_-1+l && tstp-k_-1+l >= tp) {
            ZMatrixMap(BD_P.data() + i*nao_*nao_, nao_, nao_) += cplxi * dt * I.poly_integ(i, k_+1-(tstp-tp), l) * ZMatrixMap(Sig.retptr(tstp-k_-1+l, tstp-k_-1+i), nao_, nao_).adjoint()
                                                                                                                 * ZMatrixMap(G  .retptr(tstp-k_-1+l, tp)         , nao_, nao_);
          }
          else if(tstp-k_-1+i >= tstp-k_-1+l && tstp-k_-1+l < tp) {
            ZMatrixMap(BD_P.data() + i*nao_*nao_, nao_, nao_) += cplxi * dt * I.poly_integ(i, k_+1-(tstp-tp), l) * ZMatrixMap(Sig.retptr(tstp-k_-1+i, tstp-k_-1+l), nao_, nao_)
                                                                                                                 * ZMatrixMap(G  .retptr(tp, tstp-k_-1+l)         , nao_, nao_).adjoint();
          }
          else if(tstp-k_-1+i < tstp-k_-1+l && tstp-k_-1+l < tp) {
            ZMatrixMap(BD_P.data() + i*nao_*nao_, nao_, nao_) -= cplxi * dt * I.poly_integ(i, k_+1-(tstp-tp), l) * ZMatrixMap(Sig.retptr(tstp-k_-1+l, tstp-k_-1+i), nao_, nao_).adjoint()
                                                                                                                 * ZMatrixMap(G  .retptr(tp, tstp-k_-1+l)         , nao_, nao_).adjoint();
          }
        }
      }
      for(int i = tp-tstp+k_+1; i <= k_+1; i++) {
        ZMatrixMap(BD_P.data() + i*nao_*nao_, nao_, nao_) = ZMatrixMap(G.retptr(tstp-k_-1+i, tp), nao_, nao_);
      }

      // remove bd_w(0) from M
      ZMatrixMap(M.data(), nao_, nao_) += cplxi/dt*I.bd_weights(0) * ZMatrixMap(iden.data(), nao_, nao_);

      // add to M
      for(int i = 0; i <= k_+1; i++) {
        ZMatrixMap(M.data(), nao_, nao_) += -cplxi * (double)i * BD_MI(i,k_+1) * std::pow((k_+1)*dt, i-1) * ZMatrix::Identity(nao_, nao_);
      }
      
      // add to Q
      for(int i = 0; i <= k_+1; i++) {
        for(int j = 0; j < k_+1; j++) {
          QMap += cplxi * (double)i * BD_MI(i,j) * std::pow((k_+1)*dt, i-1) * ZMatrixMap(BD_P.data() + j*nao_*nao_, nao_, nao_);
        }
      }

    }
    else {
      for(int l = 1; l <= k_+1; l++) {
        if(tstp-l >= tp)  QMap +=  cplxi/dt*I.bd_weights(l) * ZMatrixMap(G.retptr(tstp-l,tp), nao_, nao_);
        else              QMap += -cplxi/dt*I.bd_weights(l) * ZMatrixMap(G.retptr(tp,tstp-l), nao_, nao_).adjoint();
      }
    }

    int min = std::min(tp, tstp-k_);
    if(min == tstp-k_) {
      for(int l = min; l<tstp; l++) {
        if(l>=tp)  QMap -=  dt * I.poly_integ( tp-min, k_, l-min) * ZMatrixMap(Sig.retptr(tstp, l), nao_, nao_) * ZMatrixMap(G.retptr(l,tp), nao_, nao_);
        else       QMap -= -dt * I.poly_integ( tp-min, k_, l-min) * ZMatrixMap(Sig.retptr(tstp, l), nao_, nao_) * ZMatrixMap(G.retptr(tp,l), nao_, nao_).adjoint();
      }
    }
    else {
      for(int l = min; l<tstp; l++) {
        if(l>=tp)  QMap -=  dt * I.gregory_weights(tstp-tp, l-min) * ZMatrixMap(Sig.retptr(tstp, l), nao_, nao_) * ZMatrixMap(G.retptr(l,tp), nao_, nao_);
        else       QMap -= -dt * I.gregory_weights(tstp-tp, l-min) * ZMatrixMap(Sig.retptr(tstp, l), nao_, nao_) * ZMatrixMap(G.retptr(tp,l), nao_, nao_).adjoint();
      }
    }
    ZMatrixMap MMap2 = ZMatrixMap(M.data(), nao_, nao_);
    Eigen::FullPivLU<ZMatrix> lu2(MMap2);
    XMap = lu2.solve(QMap);
  }

//  err = (ZMatrixMap(X.data(), tstp*nao_*nao_, 1) - ZMatrixMap(G.retptr(tstp,0), tstp*nao_*nao_, 1)).norm();
//  ZMatrixMap(G.retptr(tstp,0), tstp*nao_*nao_, 1) = ZMatrixMap(X.data(), tstp*nao_*nao_, 1);
  cplx *check_conj = new cplx[tstp*nao_*nao_];
  cplx *G_old = new cplx[tstp*nao_*nao_];
  ZMatrixMap(check_conj, tstp*nao_*nao_, 1) = ZMatrixMap(X.data(), tstp*nao_*nao_, 1);
  ZMatrixMap(G_old, tstp*nao_*nao_, 1) = ZMatrixMap(G.retptr(tstp,0), tstp*nao_*nao_, 1);

//  return err;

*/






  int m, l, n, i;

//  double err = 0;

  std::chrono::time_point<std::chrono::system_clock> intstart, intend;
  std::chrono::duration<double> inttime;

  cplx ncplxi = cplx(0,-1);
//  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);
  ZMatrixMap QMap = ZMatrixMap(Q.data(), k_*nao_, nao_);
  ZMatrixMap MMap = ZMatrixMap(M.data(), k_*nao_, k_*nao_);
  ZMatrixMap XMap = ZMatrixMap(X.data(), k_*nao_, nao_);

  memset(M.data(), 0, k_*k_*es_*sizeof(cplx));
  memset(Q.data(), 0, (tstp+1)*es_*sizeof(cplx));

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
    err += (ZMatrixMap(G.retptr(tstp, tstp-l-1), nao_, nao_) - ZMatrixMap(X.data() + l*es_, nao_, nao_).transpose()).norm();
    ZMatrixMap(G.retptr(tstp, tstp-l-1), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_).transpose();
  }

  // Start doing the integrals
  // Q_n = \sum_{l=0}^{n-1} w_{nl} GR(t,t-l) SR(t-l,t-n) for n=k+1...tstp
  intstart = std::chrono::system_clock::now();
  for(n = k_+1; n <= tstp; n++) {
    ZMatrixMap QMapBlock = ZMatrixMap(Q.data() + n*es_, nao_, nao_);
    for(l=0; l<=k_; l++) {
      QMapBlock.noalias() += I.gregory_weights(n,l) * (ZMatrixMap(G.retptr(tstp,tstp-l), nao_, nao_)
                                          * ZMatrixMap(Sig.retptr(tstp-l,tstp-n), nao_, nao_)).transpose();
    }
  }
  intend = std::chrono::system_clock::now();
  inttime += intend-intstart;

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
    MMapSmall.noalias() -= dt*I.omega(0) * ZMatrixMap(Sig.retptr(tstp-n,tstp-n), nao_, nao_).transpose();
    MMapSmall.noalias() += (mu - I.bd_weights(0)*ncplxi/dt)*IMap;

    // Solve XM=Q for X
    Eigen::FullPivLU<ZMatrix> lu2(MMapSmall);
    ZMatrixMap(X.data(), nao_, nao_) = lu2.solve(QMapBlock).transpose();   
    err += (ZMatrixMap(G.retptr(tstp, tstp-n), nao_, nao_) - ZMatrixMap(X.data(), nao_, nao_)).norm();
    ZMatrixMap(G.retptr(tstp,tstp-n), nao_, nao_).noalias() = ZMatrixMap(X.data(), nao_, nao_);

    // Add this newly computed value to the integrals which need it
    intstart = std::chrono::system_clock::now();
    for(m=n+1; m<=tstp; m++) {
      ZMatrixMap(Q.data() + m*es_, nao_, nao_).noalias() += I.gregory_weights(m,n)
                                          *(ZMatrixMap(G.retptr(tstp,tstp-n), nao_, nao_)
                                          * ZMatrixMap(Sig.retptr(tstp-n,tstp-m), nao_, nao_)).transpose();
    }
    intend = std::chrono::system_clock::now();
    inttime += intend-intstart;
  }

//  std::cout << "RET ERR CONJ & NOCONJ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
//  for(int i = 0; i < tstp; i++) {
//    std::cout << tstp << " " << i << " " << (ZMatrixMap(G.retptr(tstp,i), nao_, nao_) - ZMatrixMap(check_conj+i*nao_*nao_, nao_, nao_)).norm()<< std::endl;
//  }

//  if(tstp >= 50 && false) {
//    err = (ZMatrixMap(G_old, tstp*nao_*nao_, 1) - ZMatrixMap(check_conj, tstp*nao_*nao_, 1)).norm();
//    ZMatrixMap(G.retptr(tstp,0), tstp*nao_*nao_, 1) = ZMatrixMap(check_conj, tstp*nao_*nao_, 1);
//  }

  //TIMING
//  std::ofstream out;
//  std::string timing_data_dir = std::string(TIMING_DATA_DIR);
//  out.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "ret_int.dat", std::ofstream::app);
//  out << inttime.count() << "\n" ;
//  out.close();
  // TIMING

  return err;

}

// idt GRM(tstp,m) - hmf(tstp) GRM(t,m) - \int_0^t dT SR(t,T) GRM(T,m) = \int_0^{beta} dT SRM(t,T) GM(T-m)
double dyson::dyson_step_tv(int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  int m, l, n, i;

  double err = 0;

  cplx cplxi = cplx(0.,1.);
  auto IMap = ZMatrixMap(iden.data(), nao_, nao_);
  auto QMap = ZMatrixMap(Q.data(), nao_, nao_);
  auto MMap = ZMatrixMap(M.data(), nao_, nao_);
  auto XMap = ZMatrixMap(X.data(), nao_, nao_);

  // Copy previous G to temporary storage
  std::memcpy(NTauTmp.data() + (ntau_+1)*es_, G.tvptr(tstp,0), (ntau_+1)*es_*sizeof(cplx));
  memset(G.tvptr(tstp,0),0,(ntau_+1)*es_*sizeof(cplx));

  // Do integrals, results go into G.tv(tstp,:), not by increment
  Ctv_tstp(tstp, G, Sig, Sig, G, G, beta, dt);

  // Put derivatives into GRM(tstp,m)
  for(l=1; l<=k_+1; l++) {
    auto GTVMap = ZColVectorMap(G.tvptr(tstp,0), (ntau_+1)*es_);
    GTVMap.noalias() += -cplxi/dt*I.bd_weights(l) * ZColVectorMap(G.tvptr(tstp-l,0), (ntau_+1)*es_);
  }

  // Make M
  MMap.noalias() = (cplxi/dt*I.bd_weights(0) + mu) * IMap
                                             - ZMatrixConstMap(hmf+tstp*es_, nao_, nao_)
                                             - dt*I.omega(0) * ZMatrixMap(Sig.retptr(tstp, tstp), nao_, nao_);

  Eigen::FullPivLU<ZMatrix> lu(MMap);
  // Solve MX=Q
  for(m=0; m<=ntau_; m++) {
    QMap.noalias() = ZMatrixMap(G.tvptr(tstp, m), nao_, nao_);
    XMap = lu.solve(QMap);
    err += (ZMatrixMap(NTauTmp.data() + (ntau_+1+m)*es_, nao_, nao_) - XMap).norm();
    ZMatrixMap(G.tvptr(tstp,m), nao_, nao_).noalias() = XMap;
  }

  return err;
}

double dyson::dyson_step_les(int n, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {

  int tstp = n;
  cplx cplxi = cplx(0,1);
  double err = 0;
  ZMatrixMap(Q.data(), (tstp+1)*nao_*nao_, 1) = ZMatrix::Zero((tstp+1)*nao_*nao_, 1);

/*
  // M
  ZMatrixMap(M.data(), nao_, nao_) = (ZMatrixConstMap(hmf+tstp*nao_*nao_, nao_, nao_) - (mu - cplxi/dt*I.bd_weights(0)) * ZMatrixMap(iden.data(), nao_, nao_) + dt*I.omega(0)*ZMatrixMap(Sig.retptr(tstp,tstp), nao_, nao_).adjoint()).transpose();
  ZMatrixMap MMap2 = ZMatrixMap(M.data(), nao_, nao_);
  Eigen::FullPivLU<ZMatrix> lu3(MMap2);

  for(int t = 0; t < tstp; t++) {
    ZMatrixMap XMap(X.data() + t*nao_*nao_, nao_, nao_);
    ZMatrixMap QMap(Q.data() + t*nao_*nao_, nao_, nao_);

    for(int l = 1; l <= k_+1; l++) {
      if(tstp-l>=t) QMap += -cplxi/dt*I.bd_weights(l) * ZMatrixMap(G.lesptr(t,tstp-l), nao_, nao_);
      else          QMap +=  cplxi/dt*I.bd_weights(l) * ZMatrixMap(G.lesptr(tstp-l,t), nao_, nao_).adjoint();
    }

    for(int l = 0; l < tstp; l++) {
      if(l>=t) QMap += -dt * I.gregory_weights(tstp, l) * ZMatrixMap(G.lesptr(t,l), nao_, nao_)           * ZMatrixMap(Sig.retptr(tstp,l), nao_, nao_).adjoint();
      else     QMap +=  dt * I.gregory_weights(tstp, l) * ZMatrixMap(G.lesptr(l,t), nao_, nao_).adjoint() * ZMatrixMap(Sig.retptr(tstp,l), nao_, nao_).adjoint();
    }

    int max = std::max(k_,t);
    for(int l = 0; l <= max; l++) {
      if(t>=l) QMap += -dt * I.gregory_weights(t, l) * ZMatrixMap(G.retptr(t,l), nao_, nao_)           * ZMatrixMap(Sig.lesptr(l,tstp), nao_, nao_);
      else     QMap +=  dt * I.gregory_weights(t, l) * ZMatrixMap(G.retptr(l,t), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(l,tstp), nao_, nao_);
    }

    XMap = QMap.transpose();
    QMap = lu3.solve(XMap).transpose();
    XMap = QMap;

  }


  int t = tstp;
  ZMatrixMap XMap(X.data() + t*nao_*nao_, nao_, nao_);
  XMap = ZMatrixMap(G.lesptr(t,t), nao_, nao_);
  ZMatrixMap QMap(Q.data() + t*nao_*nao_, nao_, nao_);

  for(int l = 0; l <= tstp; l++) {
    QMap += -dt * I.gregory_weights(tstp, l) * ZMatrixMap(X.data() + l*nao_*nao_, nao_, nao_)           * ZMatrixMap(Sig.retptr(tstp,l), nao_, nao_).adjoint();
  }

  int max = std::max(k_,t);
  for(int l = 0; l <= max; l++) {
    if(t>=l) QMap += -dt * I.gregory_weights(t, l) * ZMatrixMap(G.retptr(t,l), nao_, nao_)           * ZMatrixMap(Sig.lesptr(l,tstp), nao_, nao_);
    else     QMap +=  dt * I.gregory_weights(t, l) * ZMatrixMap(G.retptr(l,t), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(l,tstp), nao_, nao_);
  }

  QMap += -ZMatrixMap(G.lesptr(tstp,tstp), nao_, nao_) * (ZMatrixConstMap(hmf+tstp*nao_*nao_, nao_, nao_) - mu * ZMatrixMap(iden.data(), nao_, nao_));
  
  QMap = cplxi*QMap;
  XMap = QMap-QMap.adjoint();
  QMap = XMap;
  for(int l = 1; l <= k_+1; l++) {
    QMap += 1./dt*I.bd_weights(l) * ZMatrixMap(G.lesptr(tstp-l,tstp-l), nao_, nao_);
  }
  ZMatrixMap(M.data(), nao_, nao_) =  -1./dt*I.bd_weights(0) * ZMatrixMap(iden.data(), nao_, nao_);
  Eigen::FullPivLU<ZMatrix> lud(MMap2);
  XMap = lud.solve(QMap);


//  err = (ZMatrixMap(X.data(), (tstp+1)*nao_*nao_, 1) - ZMatrixMap(G.lesptr(0,tstp), (tstp+1)*nao_*nao_, 1)).norm();
//  ZMatrixMap(G.lesptr(0,tstp), (tstp+1)*nao_*nao_, 1) = ZMatrixMap(X.data(), (tstp+1)*nao_*nao_, 1);
  cplx *check_conj = new cplx[(tstp+1)*nao_*nao_];
  cplx *G_old = new cplx[(tstp+1)*nao_*nao_];
  ZMatrixMap(check_conj, (tstp+1)*nao_*nao_, 1) = ZMatrixMap(X.data(), (tstp+1)*nao_*nao_, 1);
  ZMatrixMap(G_old, (tstp+1)*nao_*nao_, 1) = ZMatrixMap(G.lesptr(0,tstp), (tstp+1)*nao_*nao_, 1);

//  return err;

*/



  // iterators
  int m, l, i;
  int num = n>=k_ ? n : k_;

//  double err=0;

  std::chrono::time_point<std::chrono::system_clock> intstart, intend;
  std::chrono::duration<double> int1, int2, int3;

  // Matricies
//  cplx cplxi = cplx(0,1);
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*k_, nao_*k_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);

  // Integrals go into Q via increment.  Q must be 0.
  memset(Q.data(),0,sizeof(cplx)*(num+1)*es_);
  memset(M.data(),0,sizeof(cplx)*k_*k_*es_);

  intstart = std::chrono::system_clock::now();
//  Cles3_tstp(Sig,Sig,G,G,n,beta,Q.data());
  intend = std::chrono::system_clock::now();
  int2 = intend-intstart;
  //TIMING
//  std::ofstream out;
//  std::string timing_data_dir = std::string(TIMING_DATA_DIR);
//  out.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "les_int_tvvt.dat", std::ofstream::app);
//  out << int2.count() << "\n" ;
//  out.close();
  // TIMING

  intstart = std::chrono::system_clock::now();
  Cles2_tstp(Sig,Sig,G,G,n,dt,Q.data());
  intend = std::chrono::system_clock::now();
  int1 = intend-intstart;
  //TIMING
//  std::ofstream out2;
//  out2.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "les_int_la.dat", std::ofstream::app);
//  out2 << int1.count() << "\n" ;
//  out2.close();
  // TIMING


  // Initial condition
//  err += (ZMatrixMap(G.lesptr(0,n), nao_, nao_) + ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint()).lpNorm<2>();
//  ZMatrixMap(G.lesptr(0,n), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();
//  ZMatrixMap(X.data(), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();


// ===================== CASE FOR NO MATSUBARA BRANCH =========================
  ZMatrix QIC = ZMatrix::Zero(nao_,nao_);
  ZMatrix MIC = ZMatrix::Zero(nao_,nao_);
  MIC = (-cplxi/dt*I.bd_weights(0) * IMap - ZMatrixConstMap(hmf + n*nao_*nao_, nao_, nao_) + mu*IMap - dt*I.gregory_weights(n, n) * ZMatrixMap(Sig.retptr(n,n), nao_, nao_).adjoint()).transpose();
  for(int i = 0; i < n; i++) {
    QIC += dt * I.gregory_weights(n, i) * ZMatrixMap(Sig.retptr(n,i), nao_, nao_).conjugate() * ZMatrixMap(G.lesptr(0,i), nao_, nao_).transpose();
  }
  for(int l = 1; l <= k_+1; l++) {
    QIC += cplxi/dt * I.bd_weights(l) * ZMatrixMap(G.lesptr(0,n-l),nao_,nao_).transpose();
  }
  Eigen::FullPivLU<ZMatrix> luIC(MIC);
  ZMatrixMap(X.data(), nao_, nao_).noalias() = luIC.solve(QIC).transpose();
  ZMatrixMap(G.lesptr(0,n), nao_, nao_).noalias() = ZMatrixMap(X.data(), nao_, nao_);

// ===================== CASE FOR NO MATSUBARA BRANCH =========================


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
      // Integral term
      if(l==0){ // Goes into Q
        QMapBlock.noalias() += dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(m,l), nao_, nao_) * ZMatrixMap(G.lesptr(l,n), nao_, nao_);
      }
      else{ // Goes into M
        if(m>=l){ // Have Sig
          MMapBlock.noalias() -= dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(m,l), nao_, nao_);
        }
        else{ // Dont have Sig
          MMapBlock.noalias() += dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(l,m), nao_, nao_).adjoint();
        }
      }
    }
  }

  // Solve Mx=Q
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  ZMatrixMap(X.data() + es_, k_*nao_, nao_).noalias() = lu.solve(ZMatrixMap(Q.data()+es_, k_*nao_, nao_));

  // Timestepping
  ZMatrixMap MMapSmall = ZMatrixMap(M.data(), nao_, nao_);
//  for(m=k_+1; m<=n; m++) {
  for(m=k_+1; m<n; m++) {
    auto QMapBlock = ZMatrixMap(Q.data() + m*es_, nao_, nao_);
    // Set up M
    MMapSmall.noalias() = -ZMatrixConstMap(hmf+m*es_, nao_, nao_) + (cplxi/dt*I.bd_weights(0) + mu)*IMap - dt*I.omega(0)*ZMatrixMap(Sig.retptr(m,m), nao_, nao_);
    // Derivatives into Q
    for(l=1; l<=k_+1; l++) {
      QMapBlock.noalias() -= cplxi/dt*I.bd_weights(l) * ZMatrixMap(X.data() + (m-l)*es_, nao_, nao_);
    }
    // Rest of the retles integral
    intstart = std::chrono::system_clock::now();
    for(l=0; l<m; l++) {
      QMapBlock.noalias() += dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(m,l), nao_, nao_)
                                             * ZMatrixMap(X.data() + l*es_, nao_, nao_);
    }
    intend = std::chrono::system_clock::now();
    int3 += intend-intstart;
    // Solve MX=Q
    Eigen::FullPivLU<ZMatrix> lu2(MMapSmall);
    ZMatrixMap(X.data()+m*es_, nao_, nao_) = lu2.solve(ZMatrixMap(Q.data() + m*es_, nao_, nao_));
  }

  // Timestepping for diagonal part (density matrix)
  // Extrapolate for G(t,t)
  memset(ex_weights.data(), 0, (k_+1)*sizeof(cplx));
  for(int ll=0; ll<=k_; ll++) {
    for(int jj=0; jj<=k_; jj++) {
      ex_weights(ll)+=I.poly_interp(jj,ll)*(1-2*(jj%2));
    }
  }
  memset(X.data()+n*es_, 0, nao_*nao_*sizeof(cplx));
  for(int jj = 0; jj <= k_; jj++) {
    ZMatrixMap(X.data()+n*es_, nao_, nao_) += ex_weights(jj) * ZMatrixMap(G.lesptr(n-jj-1,n-jj-1), nao_, nao_);
  }
  auto QMapBlock = ZMatrixMap(Q.data() + n*es_, nao_, nao_);
  QMapBlock = ZMatrix::Zero(nao_, nao_);
  // do first integral
  for(l = 0; l <= n; l++) {
    QMapBlock.noalias() -= cplxi * dt * I.gregory_weights(n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_) * ZMatrixMap(X.data()+l*es_, nao_, nao_);
  }
  // do second integral
  for(l = 0; l <= n; l++) {
    QMapBlock.noalias() += cplxi * dt * I.gregory_weights(n,l) * ZMatrixMap(Sig.lesptr(l,n), nao_, nao_).adjoint() * ZMatrixMap(G.retptr(n,l), nao_, nao_).adjoint();
  }
  // add in hamiltonian term
  QMapBlock.noalias() -= cplxi * (ZMatrixConstMap(hmf+n*es_, nao_, nao_) - mu*IMap) * ZMatrixMap(X.data()+n*es_, nao_, nao_);
  // for diagonal timestepping
  ZMatrixMap(Q.data(), nao_, nao_) = QMapBlock - QMapBlock.adjoint();
  // add in derivative term
  for(l=1; l<=k_+1; l++) {
    ZMatrixMap(Q.data(), nao_, nao_) -= 1./dt * I.bd_weights(l) * ZMatrixMap(G.lesptr(n-l,n-l), nao_, nao_);
  }
  
  MMapSmall.noalias() = 1./dt*I.bd_weights(0) * IMap;

  Eigen::FullPivLU<ZMatrix> lu2(MMapSmall);
  ZMatrixMap(X.data()+n*es_, nao_, nao_) = lu2.solve(ZMatrixMap(Q.data(), nao_, nao_));

  // Write elements into G
  for(l=0; l<=n; l++) {
    err += (ZColVectorMap(G.lesptr(l,n), es_) - ZColVectorMap(X.data() + l*es_, es_)).norm();
    ZMatrixMap(G.lesptr(l,n), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_);
  }

  //TIMING
//  std::ofstream out3;
//  out3.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "les_int_rl.dat", std::ofstream::app);
//  out3 << int3.count() << "\n" ;
//  out3.close();

/*
  std::cout << std::endl << "LES ERR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  for(int i = 0; i <= tstp; i++) {
    std::cout << tstp << " " << i << " " << (ZMatrixMap(G.lesptr(i,tstp), nao_, nao_) - ZMatrixMap(check_conj+i*nao_*nao_, nao_, nao_)).norm()<< std::endl;
  }

  if(tstp >= 50) {
    int nconj = tstp-k_+5;
    err = (ZMatrixMap(G_old          , nconj*nao_*nao_,          1) - ZMatrixMap(check_conj          , nconj*nao_*nao_,          1)).norm();
    err = (ZMatrixMap(G_old+nconj*es_, (tstp+1-nconj)*nao_*nao_, 1) - ZMatrixMap(X.data() + nconj*es_, (tstp+1-nconj)*nao_*nao_, 1)).norm();

    ZMatrixMap(G.lesptr(0,tstp)      , nconj*nao_*nao_,          1) = ZMatrixMap(check_conj          , nconj*nao_*nao_,          1);
    ZMatrixMap(G.lesptr(nconj,tstp)  , (tstp+1-nconj)*nao_*nao_, 1) = ZMatrixMap(X.data() + nconj*es_, (tstp+1-nconj)*nao_*nao_, 1);
  }
*/
  // TIMING
  return err;














/*

  // iterators
  int m, l, i;
  int num = n>=k_ ? n : k_;
  double err=0;

  std::chrono::time_point<std::chrono::system_clock> intstart, intend;
  std::chrono::duration<double> int1, int2, int3;

  // Matricies
  cplx cplxi = cplx(0,1);
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*k_, nao_*k_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);

  // Initial condition
  err += (ZMatrixMap(G.lesptr(0,n), nao_, nao_) + ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint()).lpNorm<2>();
  ZMatrixMap(G.lesptr(0,n), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();

  // Integrals go into Q via increment.  Q must be 0.
  memset(Q.data(),0,sizeof(cplx)*(num+1)*es_);
  memset(M.data(),0,sizeof(cplx)*k_*k_*es_);

  intstart = std::chrono::system_clock::now();
  Cles3_tstp(Sig,Sig,G,G,n,beta,Q.data());
  intend = std::chrono::system_clock::now();
  int2 = intend-intstart;
  //TIMING
  std::ofstream out;
  std::string timing_data_dir = std::string(TIMING_DATA_DIR);
  out.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "les_int_tvvt.dat", std::ofstream::app);
  out << int2.count() << "\n" ;
  out.close();
  // TIMING

  intstart = std::chrono::system_clock::now();
  Cles2_tstp(Sig,Sig,G,G,n,dt,Q.data());
  intend = std::chrono::system_clock::now();
  int1 = intend-intstart;
  //TIMING
  std::ofstream out2;
  out2.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "les_int_la.dat", std::ofstream::app);
  out2 << int1.count() << "\n" ;
  out2.close();
  // TIMING


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

      // Integral term
      if(l==0){ // Goes into Q
        QMapBlock.noalias() += dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(m,l), nao_, nao_) * ZMatrixMap(G.lesptr(l,n), nao_, nao_);
      }
      else{ // Goes into M
        if(m>=l){ // Have Sig
          MMapBlock.noalias() -= dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(m,l), nao_, nao_);
        }
        else{ // Dont have Sig
          MMapBlock.noalias() += dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(l,m), nao_, nao_).adjoint();
        }
      }
    }
  }

  // Solve Mx=Q
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  ZMatrixMap(X.data() + es_, k_*nao_, nao_).noalias() = lu.solve(ZMatrixMap(Q.data()+es_, k_*nao_, nao_));
  ZMatrixMap(X.data(), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();


  // Timestepping
  ZMatrixMap MMapSmall = ZMatrixMap(M.data(), nao_, nao_);
  for(m=k_+1; m<=n; m++) {
    auto QMapBlock = ZMatrixMap(Q.data() + m*es_, nao_, nao_);

    // Set up M
    MMapSmall.noalias() = -ZMatrixConstMap(hmf+m*es_, nao_, nao_) + (cplxi/dt*I.bd_weights(0) + mu)*IMap - dt*I.omega(0)*ZMatrixMap(Sig.retptr(m,m), nao_, nao_);

    // Derivatives into Q
    for(l=1; l<=k_+1; l++) {
      QMapBlock.noalias() -= cplxi/dt*I.bd_weights(l) * ZMatrixMap(X.data() + (m-l)*es_, nao_, nao_);
    }

    // Rest of the retles integral
    intstart = std::chrono::system_clock::now();
    for(l=0; l<m; l++) {
      QMapBlock.noalias() += dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(m,l), nao_, nao_)
                                             * ZMatrixMap(X.data() + l*es_, nao_, nao_);
    }
    intend = std::chrono::system_clock::now();
    int3 += intend-intstart;

    // Solve MX=Q
    Eigen::FullPivLU<ZMatrix> lu2(MMapSmall);
    ZMatrixMap(X.data()+m*es_, nao_, nao_) = lu2.solve(ZMatrixMap(Q.data() + m*es_, nao_, nao_));
  }

  // Write elements into G
  for(l=0; l<=n; l++) {
    err += (ZColVectorMap(G.lesptr(l,n), es_) - ZColVectorMap(X.data() + l*es_, es_)).norm();
    ZMatrixMap(G.lesptr(l,n), nao_, nao_).noalias() = ZMatrixMap(X.data() + l*es_, nao_, nao_);
  }
  //TIMING
  std::ofstream out3;
  out3.open(timing_data_dir + "Nao" + std::to_string(G.size1()) + "Nt" + std::to_string(G.nt()) + "Ntau" + std::to_string(G.ntau()) + "les_int_rl.dat", std::ofstream::app);
  out3 << int3.count() << "\n" ;
  out3.close();
  // TIMING

  return err;
*/
}


double dyson::dyson_step(int n, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == nao_);
  assert(G.nt() == nt_);
  assert(G.ntau() == ntau_);
  assert(G.nt() > k_);
  assert(n > k_);
  assert(n <= G.nt());

  double err = 0;

  if(mode_ == gfmol::Mode::GF2) {
    assert(G.nt() == Sig.nt());
    assert(G.size1() == Sig.size1());
    assert(G.ntau() == Sig.ntau());
    assert(G.sig() == Sig.sig());
    err += dyson_step_ret(n, G, Sig, hmf, mu, dt);
    std::cout << "err = " << err;
//    err += dyson_step_tv(n, G, Sig, hmf, mu, beta, dt);
    std::cout << " " << err;
    err += dyson_step_les(n, G, Sig, hmf, mu, beta, dt);
    std::cout << " " << err << std::endl;
  }
  else {
    err += dyson_step_ret_hf(n, G, hmf, mu, dt);
    err += dyson_step_tv_hf(n, G, hmf, mu, beta, dt);
    err += dyson_step_les_hf(n, G, hmf, mu, beta, dt);
  }
  return err; 
}



double dyson::dyson_step(int n, GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[2]);
  assert(G.size1() == hmf.shape()[1]);
  assert(G.nt() == (hmf.shape()[0]-1));

  double err = dyson_step(n, G, Sig, hmf.data(), mu, beta, dt);
  return err;
}


} // namespace
#endif // DYSON_STEP_IMPL
