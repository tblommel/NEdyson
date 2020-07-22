#ifndef DYSON_IMPL
#define DYSON_IMPL

#include "dyson.h"
#include "chrono"
#include "dyson_integrals.h"
#include "dyson_free.h"
#include "dyson_tti_integrals.h"

namespace NEdyson{

dyson::dyson(int nt, int ntau, int nao, int k) : nt_(nt),
                                                 ntau_(ntau),
                                                 nao_(nao),
                                                 es_(nao*nao),
                                                 k_(k),
                                                 ex_weights(k+1),
                                                 tmp(es_),
                                                 tmp2(es_),
                                                 iden(es_),
                                                 M(k_*k_*es_),
                                                 Q((nt_+1)*es_),
                                                 X((nt_+1)*es_),
                                                 NTauTmp((ntau_+1)*es_),
                                                 I(k_)
{
  ZMatrixMap(iden.data(), nao_, nao_).noalias() = ZMatrix::Identity(nao_, nao_);
}

// Extrapolates a Green's function object from [n-k-1,n-1] to n
void dyson::Extrapolate(int n, GREEN &G) const {
  assert(n>k_);
  assert(n<=G.nt());
  assert(G.nt() == nt_);
  assert(G.ntau() == ntau_);
  assert(G.size1() == nao_);
  

  int l, j, jcut;
  memset(ex_weights.data(), 0, (k_+1)*sizeof(cplx));

  for(l=0; l<=k_; l++) {
    for(j=0; j<=k_; j++) {
      ex_weights(l)+=I.poly_interp(j,l)*(1-2*(j%2));
    }
  }

  //right mixing
  memset(G.tvptr(n,0), 0, es_*(ntau_+1)*sizeof(cplx));
  for(l=0; l<=ntau_; l++) {
    ZMatrixMap resMap = ZMatrixMap(G.tvptr(n,l), nao_, nao_);
    for(j=0; j<=k_; j++) { 
      resMap.noalias() += ex_weights(j) * ZMatrixMap(G.tvptr(n-j-1,l), nao_, nao_);
    }
  }

  //retarded
  memset(G.retptr(n,0), 0, (n+1)*es_*sizeof(cplx));
  for(l=0; l<n-k_; l++) {
    ZMatrixMap resMap = ZMatrixMap(G.retptr(n,n-l), nao_, nao_);
    for(j=0; j<=k_; j++) {
      resMap.noalias() += ex_weights(j) * ZMatrixMap(G.retptr(n-j-1,n-l-j-1), nao_, nao_);
    }
  }
  for(l=0; l<=k_; l++) {
    jcut = (l<=n-k_-1)?k_:(n-l-1);
    ZMatrixMap resMap = ZMatrixMap(G.retptr(n,l), nao_, nao_);

    for(j=0; j<=jcut; j++) {
      resMap.noalias() += ex_weights(j) * ZMatrixMap(G.retptr(n-j-1,l), nao_, nao_);
    }

    for(j=jcut+1; j<=k_; j++) {
      resMap.noalias() -= ex_weights(j) * ZMatrixMap(G.retptr(l,n-j-1), nao_, nao_).adjoint();
    }
  }

  //less
  memset(G.lesptr(0,n), 0, (n+1)*es_*sizeof(cplx));
  ZMatrixMap(G.lesptr(0,n), nao_, nao_).noalias() = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();
  for(l=1; l<=k_; l++) {
    jcut=(k_>n-l-1)?(n-l-1):k_;

    ZMatrixMap resMap = ZMatrixMap(G.lesptr(l,n), nao_, nao_);

    for(j=0; j<=jcut; j++) {
      resMap.noalias() += ex_weights(j) * ZMatrixMap(G.lesptr(l,n-1-j), nao_, nao_);
    }
    for(j=jcut+1; j<=k_; j++) {
      resMap.noalias() -= ex_weights(j) * ZMatrixMap(G.lesptr(n-1-j,l), nao_, nao_).adjoint();
    }
  }
  for(l=k_+1; l<=n; l++) {
    ZMatrixMap resMap = ZMatrixMap(G.lesptr(l,n), nao_, nao_);

    for(j=0; j<=k_; j++) {
      resMap.noalias() += ex_weights(j) * ZMatrixMap(G.lesptr(l-j-1, n-j-1), nao_, nao_);
    }
  }
}





// Start Functions =======================================================================================
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
          QMapBlock.noalias() += ncplxi/dt * I.poly_diff(n,l) * -1.*ZMatrixMap(G.retptr(m,l), nao_, nao_).adjoint() 
                              + dt*I.poly_integ(m,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_) * -1*ZMatrixMap(G.retptr(m,l), nao_, nao_).adjoint();
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
            MMapBlock.noalias() -= dt*I.poly_integ(m,n,l) * ZMatrixMap(Sig.retptr(n,l), nao_, nao_);
          }
          else{ // Don't have it
            MMapBlock.noalias() += dt*I.poly_integ(m,n,l) * ZMatrixMap(Sig.retptr(l,n), nao_, nao_).adjoint();
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

      // Add in the integrals
      CTV2(Sig, G, n, m, beta, tmp.data());
      CTV3(Sig, G, n, m, beta, tmp2.data());
      QMapBlock.noalias() += ZMatrixMap(tmp.data(), nao_, nao_) + ZMatrixMap(tmp2.data(), nao_, nao_);
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



// Step Functions ==========================================================================================

// i dt' GR(t,t-t') - GR(t-t')hmf(t-t') - \int_0^t' GR(t,t-s) SR(t-s,t-t') = 0
void dyson::dyson_step_ret(int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double dt) const {
  int m, l, n, i;
 
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
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

  // Start doing the integrals
  // Q_n = \sum_{l=0}^{n-1} w_{nl} GR(t,t-l) SR(t-l,t-n) for n=k+1...tstp
  for(n=k_+1; n<=tstp; n++) {
    ZMatrixMap QMapBlock = ZMatrixMap(Q.data() + n*es_, nao_, nao_);
    for(l=0; l<=k_; l++) {
      QMapBlock.noalias() += I.gregory_weights(n,l) * (ZMatrixMap(G.retptr(tstp,tstp-l), nao_, nao_)
                                          * ZMatrixMap(Sig.retptr(tstp-l,tstp-n), nao_, nao_)).transpose();
    }
  }

  ZMatrixMap MMapSmall = ZMatrixMap(M.data(), nao_, nao_);
  for(n=k_+1; n<=tstp; n++) {
    ZMatrixMap QMapBlock = ZMatrixMap(Q.data() + n*es_, nao_, nao_);
    QMapBlock *= dt;
    for(l=1; l<=k_+1; l++){
      QMapBlock.noalias() += I.bd_weights(l)*ncplxi/dt * ZMatrixMap(G.retptr(tstp,tstp-n+l), nao_, nao_).transpose();
    }

    // Set up mm
    MMapSmall.noalias() = -ZMatrixConstMap(hmf+(tstp-n)*es_, nao_, nao_).transpose();
    MMapSmall.noalias() -= dt*I.omega(0) * ZMatrixMap(Sig.retptr(tstp-n,tstp-n), nao_, nao_).transpose();
    MMapSmall.noalias() += (mu - I.bd_weights(0)*ncplxi/dt)*IMap;

    // Solve XM=Q for X
    Eigen::FullPivLU<ZMatrix> lu2(MMapSmall);
    ZMatrixMap(G.retptr(tstp,tstp-n), nao_, nao_).noalias() = lu2.solve(QMapBlock).transpose();
    
    // Add this newly computed value to the integrals which need it
    for(m=n+1; m<=tstp; m++) {
      ZMatrixMap(Q.data() + m*es_, nao_, nao_).noalias() += I.gregory_weights(m,n)
                                          *( ZMatrixMap(G.retptr(tstp,tstp-n), nao_, nao_)
                                          * ZMatrixMap(Sig.retptr(tstp-n,tstp-m), nao_, nao_)).transpose();
    }
  }

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  double runtime = elapsed_seconds.count();

  std::ofstream out;
  out.open("dystiming.dat", std::ofstream::app);
  out<<runtime<<" ";
}


// idt GRM(tstp,m) - hmf(tstp) GRM(t,m) - \int_0^t dT SR(t,T) GRM(T,m) = \int_0^{beta} dT SRM(t,T) GM(T-m)
void dyson::dyson_step_tv(int tstp, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
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
                                             - ZMatrixConstMap(hmf+tstp*es_, nao_, nao_) 
                                             - dt*I.omega(0) * ZMatrixMap(Sig.retptr(tstp, tstp), nao_, nao_);

  Eigen::FullPivLU<ZMatrix> lu(MMap);
  // Solve MX=Q
  for(m=0; m<=ntau_; m++) {
    QMap.noalias() = ZMatrixMap(G.tvptr(tstp, m), nao_, nao_);
    ZMatrixMap(G.tvptr(tstp,m), nao_, nao_).noalias() = lu.solve(QMap);
  }

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  double runtime = elapsed_seconds.count();

  std::ofstream out;
  out.open("dystiming.dat", std::ofstream::app);
  out<<runtime<<" ";
}

double dyson::dyson_step_les(int n, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  // iterators
  int m, l, i;
  int num = n>=k_ ? n : k_;
  double err=0;

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;
  
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
  Cles2_tstp(Sig,Sig,G,G,n,dt,Q.data());
  Cles3_tstp(Sig,Sig,G,G,n,beta,Q.data());

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
  
  ZMatrixMap(X.data(), nao_, nao_).noalias() = ZMatrixMap(G.lesptr(0,n), nao_, nao_);

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
    for(l=0; l<m; l++) {
      QMapBlock.noalias() += dt*I.gregory_weights(m,l) * ZMatrixMap(Sig.retptr(m,l), nao_, nao_)
                                             * ZMatrixMap(X.data() + l*es_, nao_, nao_);
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
    out.open("dystiming.dat", std::ofstream::app);
    out<<runtime<<std::endl;
  }
  return err;
}


double dyson::dyson_start_les(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  double err=0;
  for(int n=0; n<=k_; n++) err += dyson_step_les(n,G,Sig,hmf,mu,beta,dt);
  return err;
}


double dyson::dyson_start(GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == nao_);
  assert(G.size1() == Sig.size1());
  assert(G.nt() == nt_);
  assert(G.nt() == Sig.nt());
  assert(G.nt() >= k_);
  assert(G.ntau() == Sig.ntau());
  assert(G.ntau() == ntau_);
  assert(G.sig() == Sig.sig());

  double err=0;
  err += dyson_start_ret(G, Sig, hmf, mu, dt);
  err += dyson_start_tv(G, Sig, hmf, mu, beta, dt);
  err += dyson_start_les(G, Sig, hmf, mu, beta, dt);
  return err;
}


double dyson::dyson_start(GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[2]);
  assert(G.size1() == hmf.shape()[1]);
  assert(G.nt() == hmf.shape()[0]-1);

  return dyson_start(G, Sig, hmf.data(), mu, beta, dt);
}



void dyson::dyson_step(int n, GREEN &G, const GREEN &Sig, const cplx *hmf, double mu, double beta, double dt) const {
  assert(G.size1() == Sig.size1());
  assert(G.size1() == nao_);
  assert(G.nt() == Sig.nt());
  assert(G.nt() == nt_());
  assert(G.ntau() == Sig.ntau());
  assert(G.ntau() == ntau());
  assert(G.nt() > k_);
  assert(G.sig() == Sig.sig());
  assert(n > k_);
  assert(n <= G.nt());

  dyson_step_ret(n, G, Sig, hmf, mu, dt);
  dyson_step_tv(n, G, Sig, hmf, mu, beta, dt);
  dyson_step_les(n, G, Sig, hmf, mu, beta, dt);
}



void dyson::dyson_step(int n, GREEN &G, const GREEN &Sig, const ZTensor<3> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[2]);
  assert(G.size1() == hmf.shape()[1]);
  assert(G.nt() == (hmf.shape()[0]-1));

  dyson_step(n, G, Sig, hmf.data(), mu, beta, dt);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Extrapolates a Green's function object from [n-k-1,n-1] to n
void dyson::Extrapolate(int n, TTI_GREEN &G) const {
  assert(n > k_);
  assert(n <= G.nt());
  assert(G.nt() == nt_);
  assert(G.ntau() == ntau_);
  assert(G.size1() == nao_);

  int l, j, jcut;
  memset(ex_weights.data(), 0, (k_+1)*sizeof(cplx));

  for(l=0; l<=k_; l++) {
    for(j=0; j<=k_; j++) {
      ex_weights(l) += I.poly_interp(j,l)*(1-2*(j%2));
    }
  }

  //right mixing
  memset(G.tvptr(n,0), 0, es_*(ntau_+1)*sizeof(cplx));
  for(l=0; l<=ntau_; l++) {
    ZMatrixMap resMap = ZMatrixMap(G.tvptr(n,l), nao_, nao_);
    for(j=0; j<=k_; j++) { 
      resMap.noalias() += ex_weights(j) * ZMatrixMap(G.tvptr(n-j-1,l), nao_, nao_);
    }
  }
  //retarded
  memset(G.retptr(n), 0, es_*sizeof(cplx));
  ZMatrixMap resMap = ZMatrixMap(G.retptr(n), nao_, nao_);
  for(j=0; j<=k_; j++) {
    resMap.noalias() += ex_weights(j) * ZMatrixMap(G.retptr(n-j-1), nao_, nao_);
  }
  //less
  ZMatrixMap(G.lesptr(-n), nao_, nao_) = -ZMatrixMap(G.tvptr(n,0), nao_, nao_).adjoint();
}


double dyson::dyson_start_ret(TTI_GREEN &G, const TTI_GREEN &Sig, const double *hmf, double mu, double dt) const {
  // Counters
  int m, l, n, i;
  
  double err=0;
  cplx ncplxi = cplx(0, -1);
  
  ZMatrixMap QMap = ZMatrixMap(Q.data(), nao_*k_, nao_);
  ZMatrixMap MMap = ZMatrixMap(M.data(), nao_*k_, nao_*k_);
  ZMatrixMap IMap = ZMatrixMap(iden.data(), nao_, nao_);

  // Boundary Condition
  ZMatrixMap(G.retptr(0), nao_, nao_) = ncplxi*IMap;
  
  memset(M.data(), 0, k_*k_*es_*sizeof(cplx));
  memset(Q.data(), 0, k_*es_*sizeof(cplx));

  for(n=1; n<=k_; n++) {
    auto QMapBlock = ZMatrixMap(Q.data() + (n-1)*es_, nao_, nao_);

    for(l=0; l<=k_; l++) {
      auto MMapBlock = MMap.block((n-1)*nao_, (l-1)*nao_, nao_, nao_);

      if(l == 0){ // We know these G's. Put into Q
        QMapBlock.noalias() += ncplxi/dt * I.poly_diff(n,l) * -1.*ZMatrixMap(G.retptr(0), nao_, nao_).adjoint() 
                            + dt*I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(n), nao_, nao_) * -1*ZMatrixMap(G.retptr(0), nao_, nao_).adjoint();
      }
      else{ // Don't have these. Put into M
        // Derivative term
        MMapBlock.noalias() = -ncplxi/dt * I.poly_diff(n,l) * IMap;

        // Delta energy term
        if(n==l){
          MMapBlock.noalias() += mu*IMap - DMatrixConstMap(hmf, nao_, nao_);
        }

        // Integral term
        if(n>=l){ // We have Sig
          MMapBlock.noalias() -= dt*I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(n-l), nao_, nao_);
        }
        else{ // Don't have it
          MMapBlock.noalias() += dt*I.poly_integ(0,n,l) * ZMatrixMap(Sig.retptr(l-n), nao_, nao_).adjoint();
        }
      }
    }
  }

  //solve MX=Q for X
  Eigen::FullPivLU<ZMatrix> lu(MMap);
  ZMatrixMap(X.data(), k_*nao_, nao_) = lu.solve(QMap);

  //put X into G
  for(l=0; l<k_; l++) {
    ZColVectorMap retMap = ZColVectorMap(G.retptr(l+1), es_);
    ZColVectorMap XMap = ZColVectorMap(X.data() + l*es_, es_);
    err += (retMap - XMap).norm();
    retMap = XMap;
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

      // Add in the integrals
      CTV2(Sig, G, n, m, beta, tmp.data());
      CTV3(Sig, G, n, m, beta, tmp2.data());
      QMapBlock.noalias() += ZMatrixMap(tmp.data(), nao_, nao_) + ZMatrixMap(tmp2.data(), nao_, nao_);
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
  assert(G.size1() == Sig.size1());
  assert(G.nt() == nt_);
  assert(G.nt() == Sig.nt());
  assert(G.nt() >= k_);
  assert(G.ntau() == Sig.ntau());
  assert(G.ntau() == ntau_);
  assert(G.sig() == Sig.sig());

  double err=0;
  err += dyson_start_ret(G, Sig, hmf, mu, dt);
  err += dyson_start_tv(G, Sig, hmf, mu, beta, dt);
  err += dyson_start_les(G, Sig, hmf, mu, beta, dt);
  return err;
}

double dyson::dyson_start(TTI_GREEN &G, const TTI_GREEN &Sig, const DTensor<2> &hmf, double mu, double beta, double dt) const {
  assert(G.size1() == hmf.shape()[1]);
  assert(G.size1() == hmf.shape()[0]);

  double err=0;
  err += dyson_start_ret(G, Sig, hmf.data(), mu, dt);
  err += dyson_start_tv(G, Sig, hmf.data(), mu, beta, dt);
  err += dyson_start_les(G, Sig, hmf.data(), mu, beta, dt);
  return err;
}

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
  out.open("dystiming.dat", std::ofstream::app);
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
  out.open("tti_dystiming.dat", std::ofstream::app);
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
  assert(G.ntau() == ntau());
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

  dyson_step_ret(n, G, Sig, hmf.data(), mu, dt);
  dyson_step_tv(n, G, Sig, hmf.data(), mu, beta, dt);
  dyson_step_les(n, G, Sig, hmf.data(), mu, beta, dt);
}

} // Namespace

#endif
