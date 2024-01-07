#ifndef DYSON_BOOT_IMPL
#define DYSON_BOOT_IMPL

namespace NEdyson {

double dyson::dyson_start_ret(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double dt) const {
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

  // Fill GR(:k+1,0)
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*k_, nao_*k_);
  memset(M.data(), 0, k_*k_*es_*sizeof(cplx));
  memset(Q.data(), 0, k_*es_*sizeof(cplx));

  for(m=0; m<k_; m++) {
    ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*(k_-m), nao_*(k_-m));
    memset(M.data(), 0, k_*k_*es_*sizeof(cplx));
    memset(Q.data(), 0, k_*es_*sizeof(cplx));
      
    for(n=m+1; n<=k_; n++) {
      auto QMapBlock = QMap.block((n-m-1)*nao_, 0, nao_, nao_);
  
      for(l=0; l<=m; l++) {
//        QMapBlock.noalias() -= -ncplxi/dt * I.poly_diff(n,l) * -1.*ZMatrixMap(G.retptr(m,l), nao_, nao_).adjoint() - dt * I.poly_integ(m,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_) * -1. * ZMatrixMap(G.retptr(m,l), nao_, nao_).adjoint();
        QMapBlock.noalias() -= -ncplxi/dt * I.poly_diff(n,l) * -1.*ZMatrixMap(G.retptr(m,l), nao_, nao_).adjoint();// - dt * I.poly_integ(m,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_) * -1. * ZMatrixMap(G.retptr(m,l), nao_, nao_).adjoint();
      }

      for(l = m+1; l <= k_; l++) {
        MMap.block((n-m-1)*nao_, (l-m-1)*nao_, nao_, nao_) += -ncplxi/dt * I.poly_diff(n,l) * IMap;
        if(n==l) MMap.block((n-m-1)*nao_, (l-m-1)*nao_, nao_, nao_) += mu*IMap - ZMatrixConstMap(hmf + l*es_, nao_, nao_);
//        if(n>=l) MMap.block((n-m-1)*nao_, (l-m-1)*nao_, nao_, nao_) += -dt*I.poly_integ(m,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_);
//        else     MMap.block((n-m-1)*nao_, (l-m-1)*nao_, nao_, nao_) -= -dt*I.poly_integ(m,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_).adjoint();
      }

    }

    // Solve MX=Q for X
    Eigen::FullPivLU<ZMatrix> lu(ZMatrixMap(M.data(), (k_-m)*nao_, (k_-m)*nao_));
    ZMatrixMap(X.data(), (k_-m)*nao_, nao_) = lu.solve(ZMatrixMap(Q.data(), (k_-m)*nao_, nao_));
    
    // Put X into G
    for(n=m+1; n<=k_; n++){
      err += (ZColVectorMap(G.retptr(n,m), es_) - ZColVectorMap(X.data() + (n-m-1)*es_, es_)).norm();
      ZMatrixMap(G.retptr(n,m), nao_, nao_).noalias() = ZMatrixMap(X.data() + (n-m-1)*es_, nao_, nao_);
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
    auto tvmap = ZColVectorMap(G.tvptr(0,m), es_);
    auto matmap = ZColVectorMap(G.matptr(ntau_-m), es_);
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
      QMapBlock.noalias() += ZMatrixMap(NTauTmp.data() + (n-1)*(ntau_+1)*nao_*nao_ + m*nao_*nao_, nao_, nao_);
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
  cplx cplxi = cplx(0,1);
  ZMatrixMap(G.lesptr(0,0), nao_, nao_) = -(-cplxi*ZMatrixMap(G.matptr(ntau_), nao_, nao_)).adjoint();

/*
  G.lesptr(0,0)[0]  = cplxi * 0.5;
  G.lesptr(0,0)[1]  = cplxi * 0.447214;
  G.lesptr(0,0)[2]  = cplxi * 4.7358e-16;
  G.lesptr(0,0)[3]  = cplxi * -0.223607;
  G.lesptr(0,0)[4]  = cplxi * 0.447214;
  G.lesptr(0,0)[5]  = cplxi * 0.5;
  G.lesptr(0,0)[6]  = cplxi * 0.223607;
  G.lesptr(0,0)[7]  = cplxi * -4.7358e-16;
  G.lesptr(0,0)[8]  = cplxi * 4.7358e-16;
  G.lesptr(0,0)[9]  = cplxi * 0.223607;
  G.lesptr(0,0)[10] = cplxi * 0.5;
  G.lesptr(0,0)[11] = cplxi * 0.447214;
  G.lesptr(0,0)[12] = cplxi * -0.223607;
  G.lesptr(0,0)[13] = cplxi * -4.7358e-16;
  G.lesptr(0,0)[14] = cplxi * 0.447214;
  G.lesptr(0,0)[15] = cplxi * 0.5;


  G.lesptr(0,0)[0]  = cplxi * 0.5;
  G.lesptr(0,0)[1]  = cplxi * 0.444505;
  G.lesptr(0,0)[2]  = cplxi * 3.72966e-17;
  G.lesptr(0,0)[3]  = cplxi * -0.219471;
  G.lesptr(0,0)[4]  = cplxi * 0.444505;
  G.lesptr(0,0)[5]  = cplxi * 0.5;
  G.lesptr(0,0)[6]  = cplxi * 0.22256;
  G.lesptr(0,0)[7]  = cplxi * -1.30104e-17;
  G.lesptr(0,0)[8]  = cplxi * 3.72966e-17;
  G.lesptr(0,0)[9]  = cplxi * 0.22256;
  G.lesptr(0,0)[10] = cplxi * 0.5;
  G.lesptr(0,0)[11] = cplxi * 0.444505;
  G.lesptr(0,0)[12] = cplxi * -0.219471;
  G.lesptr(0,0)[13] = cplxi * -1.30104e-17;
  G.lesptr(0,0)[14] = cplxi * 0.444505;
  G.lesptr(0,0)[15] = cplxi * 0.5;

*/




  ZMatrix QIC = ZMatrix::Zero(k_*nao_, nao_);
  ZMatrix XIC = ZMatrix::Zero(k_*nao_, nao_);
  ZMatrix MIC = ZMatrix::Zero(k_*nao_, k_*nao_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);
  double err = 0;
  int m = 0;

  // FIRST COLUMN, T=0
  for(int n = 1; n <= k_; n++) {
    auto QBlock = QIC.block((n-1)*nao_, 0, nao_, nao_);
    for(int l = 0; l <= k_; l++) {
      if(m>=l && n>=l) QBlock += (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(m,l), nao_, nao_) * ZMatrixMap(Sig.lesptr(l,n), nao_, nao_)).transpose();
      else if(m<l && n>=l) QBlock -= (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(l,m), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(l,n), nao_, nao_)).transpose();
      else if(m>=l && n<l) QBlock -= (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(m,l), nao_, nao_) * ZMatrixMap(Sig.lesptr(n,l), nao_, nao_).adjoint()).transpose();
      else if(m<l && n<l) QBlock += (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(l,m), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(n,l), nao_, nao_).adjoint()).transpose();
    }
    int l = 0;
    QBlock += (dt * I.poly_integ(0,n,l) * ZMatrixMap(G.lesptr(m,l), nao_, nao_) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_).adjoint()).transpose();
    QBlock += (cplxi/dt * I.poly_diff(n,l) * ZMatrixMap(G.lesptr(m,l), nao_, nao_)).transpose();
  }

  for(int n = 1; n <= k_; n++) {
    for(int l = 1; l <= k_; l++) {
      auto MBlock = MIC.block((n-1)*nao_, (l-1)*nao_, nao_, nao_);
      MBlock += -cplxi/dt * I.poly_diff(n,l) * IMap;
      if(n==l) MBlock -= ZMatrixConstMap(hmf+l*es_, nao_, nao_).transpose() - mu*IMap;
      if(n>=l) MBlock -= (dt * I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_).adjoint()).transpose();
      else if(n<l) MBlock += (dt * I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(l,n), nao_, nao_)).transpose();
    }
  }
  Eigen::FullPivLU<ZMatrix> lu(MIC);
  XIC = lu.solve(QIC);
  for(int i = 1; i <= k_; i++) {
    err += (ZMatrixMap(G.lesptr(0,i), nao_, nao_) - ZMatrixMap(XIC.data() + (i-1)*es_, nao_, nao_).transpose()).norm();
    ZMatrixMap(G.lesptr(0,i), nao_, nao_) = ZMatrixMap(XIC.data() + (i-1)*es_, nao_, nao_).transpose();
  }


  ZMatrix DIC = ZMatrix::Zero(k_*nao_, nao_);
  for(int i = 1; i <= k_; i++) {
    ZMatrixMap(DIC.data() + (i-1)*es_, nao_, nao_) = ZMatrixMap(G.lesptr(i,i), nao_, nao_);
  }

  // REST OF THE COLUMNS
  for(m = 1; m <= k_; m++) {
    MIC = ZMatrix::Zero(k_*nao_, k_*nao_);
    XIC = ZMatrix::Zero(k_*nao_, nao_);
    QIC = ZMatrix::Zero(k_*nao_, nao_);

    for(int n = 1; n <= k_; n++) {
      auto QBlock = QIC.block((n-1)*nao_, 0, nao_, nao_);
      for(int l = 0; l <= k_; l++) {
        if(m>=l && n>=l) QBlock += (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(m,l), nao_, nao_) * ZMatrixMap(Sig.lesptr(l,n), nao_, nao_)).transpose();
        else if(m<l && n>=l) QBlock -= (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(l,m), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(l,n), nao_, nao_)).transpose();
        else if(m>=l && n<l) QBlock -= (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(m,l), nao_, nao_) * ZMatrixMap(Sig.lesptr(n,l), nao_, nao_).adjoint()).transpose();
        else if(m<l && n<l) QBlock += (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(l,m), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(n,l), nao_, nao_).adjoint()).transpose();
      }
      QBlock -= (dt * I.poly_integ(0,n,0) * ZMatrixMap(G.lesptr(0,m), nao_, nao_).adjoint() * ZMatrixMap(Sig.retptr(n,0), nao_, nao_).adjoint()).transpose();
      QBlock -= (cplxi/dt * I.poly_diff(n,0) * ZMatrixMap(G.lesptr(0,m), nao_, nao_).adjoint()).transpose();
    }
    for(int n = 1; n <= k_; n++) {
      for(int l = 1; l <= k_; l++) {
        auto MBlock = MIC.block((n-1)*nao_, (l-1)*nao_, nao_, nao_);
        MBlock += -cplxi/dt * I.poly_diff(n,l) * IMap;
        if(n==l) MBlock -= ZMatrixConstMap(hmf+l*es_, nao_, nao_).transpose() - mu*IMap;
        if(n>=l) MBlock -= (dt * I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_).adjoint()).transpose();
        else if(n<l) MBlock += (dt * I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(l,n), nao_, nao_)).transpose();
      }
    }
    Eigen::FullPivLU<ZMatrix> lu2(MIC);
    XIC = lu2.solve(QIC);
    for(int i = m; i <= k_; i++) {
      err += i==m ? 0 : (ZMatrixMap(G.lesptr(m,i), nao_, nao_) - ZMatrixMap(XIC.data() + (i-1)*es_, nao_, nao_).transpose()).norm();
      ZMatrixMap(G.lesptr(m,i), nao_, nao_) = ZMatrixMap(XIC.data() + (i-1)*es_, nao_, nao_).transpose();
    }
  }

  // redo diagonal
  MIC = ZMatrix::Zero(k_*nao_, k_*nao_);
  XIC = ZMatrix::Zero(k_*nao_, nao_);
  QIC = ZMatrix::Zero(k_*nao_, nao_);

  for(int i = 1; i <= k_; i++) {
    for(int j = 1; j <= k_; j++) {
      auto MBlock = MIC.block((i-1)*nao_, (j-1)*nao_, nao_, nao_);
      MBlock = 1./dt * I.poly_diff(i,j) * IMap;
    }
  }

  for(int i = 1; i <= k_; i++) {
    auto QBlock = QIC.block((i-1)*nao_, 0, nao_, nao_);
    QBlock += cplxi * ZMatrixMap(G.lesptr(i,i), nao_, nao_) * (ZMatrixConstMap(hmf+i*es_, nao_, nao_) - mu*IMap);

    for(int l = 0; l <= i; l++) {
      QBlock += cplxi * I.poly_integ(0,i,l) * dt * ZMatrixMap(G.retptr(i,l), nao_, nao_) * ZMatrixMap(Sig.lesptr(l,i), nao_, nao_);
    }
    for(int l = i+1; l <= k_; l++) {
      QBlock += cplxi * I.poly_integ(0,i,l) * dt * ZMatrixMap(G.retptr(l,i), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(i,l), nao_, nao_).adjoint();
    }
    for(int l = 0; l <= i; l++) {
      QBlock -= cplxi * I.poly_integ(0,i,l) * dt * ZMatrixMap(G.lesptr(l,i), nao_, nao_).adjoint() * ZMatrixMap(Sig.retptr(i,l), nao_, nao_).adjoint();
    }
    for(int l = i+1; l <= k_; l++) {
      QBlock -= cplxi * I.poly_integ(0,i,l) * dt * ZMatrixMap(G.lesptr(i,l), nao_, nao_) * ZMatrixMap(Sig.retptr(l,i), nao_, nao_);
    }

    ZMatrixMap(XIC.data(), nao_, nao_) = QBlock-QBlock.adjoint();
    QBlock = ZMatrixMap(XIC.data(), nao_, nao_);

    QBlock -= 1./dt * I.poly_diff(i,0) * ZMatrixMap(G.lesptr(0,0), nao_, nao_);
  }

  Eigen::FullPivLU<ZMatrix> lu3(MIC);
  XIC = lu3.solve(QIC);
  for(int i = 1; i <= k_; i++) {
    err += (ZMatrixMap(XIC.data()+(i-1)*es_, nao_, nao_) - ZMatrixMap(DIC.data() + (i-1)*es_, nao_, nao_)).norm();
    ZMatrixMap(G.lesptr(i,i), nao_, nao_) = ZMatrixMap(XIC.data() + (i-1)*es_, nao_, nao_);
  }

  
  
/*
  for(int m = 1; m <= k_; m++) {
    ZMatrixMap MMap = ZMatrixMap(MIC.data(), (k_-m+1)*nao_, (k_-m+1)*nao_);
    ZMatrixMap QMap = ZMatrixMap(QIC.data(), (k_-m+1)*nao_, nao_);
    ZMatrixMap XMap = ZMatrixMap(XIC.data(), (k_-m+1)*nao_, nao_);
   sdfasdfdsf MMap = ZMatrix::Zero((k_-m+1)*nao_, (k_-m+1)*nao_);
    QMap = ZMatrix::Zero((k_-m+1)*nao_, nao_);
    XMap = ZMatrix::Zero((k_-m+1)*nao_, nao_);
    for(int n = m; n <= k_; n++) {
      auto QBlock = QMap.block((n-m)*nao_, 0, nao_, nao_);
      for(int l = 0; l <= k_; l++) {
        if(m>=l && n>=l) QBlock += (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(m,l), nao_, nao_) * ZMatrixMap(Sig.lesptr(l,n), nao_, nao_)).transpose();
        else if(m<l && n>=l) QBlock -= (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(l,m), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(l,n), nao_, nao_)).transpose();
        else if(m>=l && n<l) QBlock -= (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(m,l), nao_, nao_) * ZMatrixMap(Sig.lesptr(n,l), nao_, nao_).adjoint()).transpose();
        else if(m<l && n<l) QBlock += (dt * I.poly_integ(0,m,l) * ZMatrixMap(G.retptr(l,m), nao_, nao_).adjoint() * ZMatrixMap(Sig.lesptr(n,l), nao_, nao_).adjoint()).transpose();
      }
      for(int l = 0; l < m; l++) {
        QBlock -= (cplxi/dt * I.poly_diff(n,l) * ZMatrixMap(G.lesptr(l,m), nao_, nao_).adjoint()).transpose();
      }
      for(int l = 0; l < m; l++) {
        QBlock -= (dt * I.poly_integ(0,n,l) * ZMatrixMap(G.lesptr(l,m), nao_, nao_).adjoint()*ZMatrixMap(Sig.retptr(n,l), nao_, nao_).adjoint()).transpose();
      }
    }
    for(int n = m; n <= k_; n++) {
      for(int l = m; l <= k_; l++) {
        auto MBlock = MMap.block((n-m)*nao_, (l-m)*nao_, nao_, nao_);
        MBlock += -cplxi/dt * I.poly_diff(n,l) * IMap;
        if(n==l) MBlock -= ZMatrixConstMap(hmf+l*es_, nao_, nao_).transpose() - mu*IMap;
        if(n>=l) MBlock -= (dt * I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_).adjoint()).transpose();
        else if(n<l) MBlock += (dt * I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(l,n), nao_, nao_)).transpose();
      }
    }

    Eigen::FullPivLU<ZMatrix> lu(MMap);
    XIC = lu.solve(QMap);
    if(m==2) std::cout<< "OLD M VS NEW M: " << (MMap-MIC2).norm() << std::endl;
    if(m==2) std::cout<< "OLD Q VS NEW Q: " << (QMap-QIC2).norm() << std::endl;
    if(m==2) std::cout<< "OLD X VS NEW X: " << (XMap-XIC2).norm() << std::endl;
    if(m==2) std::cout<< MMap-MIC2 << std::endl;
    for(int i = m; i <= k_; i++) {
      err += (ZMatrixMap(G.lesptr(m,i), nao_, nao_) - ZMatrixMap(XIC.data() + (i-m)*es_, nao_, nao_).transpose()).norm();
      ZMatrixMap(G.lesptr(m,i), nao_, nao_) = ZMatrixMap(XIC.data() + (i-m)*es_, nao_, nao_).transpose();
    }
  }
*/
  return err;
}


double dyson::dyson_start(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == nao_);
  assert(G.nt() == nt_);
  assert(G.nt() >= k_);
  assert(G.ntau() == ntau_);

  double err=0;
  if(mode_ == gfmol::Mode::GF2) {
    assert(G.sig() == Sig.sig());
    assert(G.nt() == Sig.nt());
    assert(G.ntau() == Sig.ntau());
    assert(G.size1() == Sig.size1());
    err += dyson_start_ret(G, Sig, hmf, mu, dt);
    err += dyson_start_tv(G, Sig, hmf, mu, beta, dt);
    err += dyson_start_les(G, Sig, hmf, mu, beta, dt);
  }
  else {
    err += dyson_start_ret_hf(G, hmf, mu, dt);
    err += dyson_start_tv_hf(G, hmf, mu, beta, dt);
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
