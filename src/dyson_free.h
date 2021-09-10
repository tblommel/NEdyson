#ifndef FREE_DYSON_IMPL
#define FREE_DYSON_IMPL

namespace NEdyson{

#define EXPMAX 100

// 1/(1+exp(ometa*beta))
double fermi(double beta, double omega) {
    double arg = omega * beta;
    if (fabs(arg) > EXPMAX) {
        return (arg > 0.0 ? 0.0 : 1.0);
    } else {
        return 1.0 / (1.0 + exp(arg));
    }
}


// 1/(1+exp(ometa*beta))
DColVector fermi(double beta, DColVector &omega){
   int size=omega.size();
   DColVector tmp(size);
   for(int i=0;i<size;i++){
      tmp(i)=fermi(beta,omega(i));
   }
   return tmp;
}


// exp(w*t)/(1+exp(w*b))
double fermi_exp(double beta, double tau, double omega) {
    if (omega < 0) {
        return exp(omega * tau) *
               fermi(beta, omega); // exp(w*t)/(1+exp(b*w)) always OK for w<0
    } else {
        return exp((tau - beta) * omega) *
               fermi(beta, -omega); // exp((t-b)*w)/(1+exp(-w*b))
    }
}


// exp(w*t)/(1+exp(w*b))
DColVector fermi_exp(double beta,double tau,DColVector &omega){
   int size=omega.size();
   DColVector tmp(size);

   for(int i=0;i<size;i++){
      tmp(i)=fermi_exp(beta,tau,omega(i));
   }
   return tmp;
}


// 1/(exp(w*b)-1)
double bose(double beta, double omega) {
    double arg = omega * beta;
    if (arg < 0)
        return (-1.0 - bose(beta, -omega));
    if (fabs(arg) > EXPMAX) {
        return 0.0;
    } else if (arg < 1e-10) {
        return 1.0 / arg;
    } else {
        return 1.0 / (exp(arg) - 1.0);
    }
}


// 1/(exp(w*b)-1)
DColVector bose(double beta,DColVector &omega){
   int size=omega.size();
   DColVector tmp(size);
   for(int i=0;i<size;i++){
      tmp(i)=bose(beta,omega(i));
   }
   return tmp;
}


// exp(w*t)/(exp(w*b)-1)
double bose_exp(double beta, double tau, double omega) {
    if (omega < 0)
        return exp(tau * omega) * bose(beta, omega);
    else
        return -exp((tau - beta) * omega) * bose(beta, -omega);
}


// exp(w*t)/(exp(w*b)-1)
DColVector bose_exp(double beta,double tau,DColVector &omega){
   int size=omega.size();
   DColVector tmp(size);

   for(int i=0;i<size;i++){
      tmp(i)=bose_exp(beta,tau,omega(i));
   }
   return tmp;
}



// gives an interpolation of a function on the domain [tstp-k,tstp]
// eps(tinterp) = sum_{l=0}^k eps(tstp-k+l) * sum_{p=0}^k (tinterp+k-tstp)^p P_{pl}
void interpolate(int tstp, double tinterp, const INTEG &I, const double *eps, int size, cplx *res){
  int p,l,k=I.k(),es=size*size;
  double timk=-tstp+k+tinterp;
  
  memset(res, 0, es*sizeof(cplx));  
  ZMatrixMap resMap = ZMatrixMap(res, size, size);
  
  double weight,pref;

  for(l=0; l<=k; l++) {
    pref = 1.;
    weight = I.poly_interp(0,l);
    for(p=1; p<=k; p++) {
      pref *= timk;
      weight += pref*I.poly_interp(p,l);
    }
    resMap.noalias() += weight * DMatrixConstMap(eps+(tstp-k+l)*es, size, size);
  }
}



// evaluates Ut(tstp) by the approximation exp(-I*(a1*eps(tstp-1+c1)+a2*eps(tstp-1+c2)))
// *exp(-I*(a1*eps(tstp-1+c2)+a2*eps(tstp-1+c1)))*Ut(tstp-1)
// a1=(3-2sqrt(3))/12
// a2=(3+2sqrt(3))/12
// c1=(1-sqrt(3)/3)/2
// c2=(1+sqrt(3)/3)/2
void propagator_exp(int tstp, const INTEG &I, cplx *Ut, const double *ht, double dt, int size){
  assert(tstp > 0);
  assert(tstp <= I.k());
  
  ZMatrix eps1(size, size), eps2(size, size);
  ZMatrix arg1(size, size), arg2(size, size);
  
  int n = (tstp<I.k())?(I.k()):(tstp);

  interpolate(n, tstp-1./2-sqrt(3.)/6., I, ht, size, eps1.data());
  interpolate(n, tstp-1./2+sqrt(3.)/6., I, ht, size, eps2.data());
  double a1 = (3.-2*sqrt(3.))/12.;
  double a2 = (3.+2*sqrt(3.))/12.;
  arg1=-cplx(0.,dt)*(a1*eps1+a2*eps2);
  arg2=-cplx(0.,dt)*(a2*eps1+a1*eps2);

  ZMatrixMap(Ut + tstp*size*size, size, size).noalias() = arg1.exp() * arg2.exp() * ZMatrixMap(Ut + (tstp-1)*size*size, size, size);
}

// Gives the free greens function from a non-constant hamiltonian
void dyson::G0_from_h0(GREEN &G, double mu, const double *hM, const double *ht, double beta, double dt) const {
  assert(G.nt() == nt_);
  assert(G.ntau() == ntau_);
  assert(G.size1() == nao_);
  assert(G.nt() >= k_);

  int sig=G.sig(), m, n;
  double tau, t;
  
  ZMatrix mHmu = mu*ZMatrixMap(iden.data(), nao_, nao_) - DMatrixConstMap(hM, nao_, nao_);
  Eigen::SelfAdjointEigenSolver<ZMatrix> eigensolver(mHmu);
  
  ZMatrix evec0(nao_, nao_);
  DColVector eval0(nao_), eval0m(nao_);

  evec0 = eigensolver.eigenvectors();
  eval0 = eigensolver.eigenvalues();
  eval0m= -eval0;

  // Matsubara
  for(m=0; m<=ntau_; m++){
    tau = beta/2. * (Convolution().collocation().x_i()(m) + 1.);
    if(sig==-1)     ZMatrixMap(G.matptr(m), nao_, nao_) = -evec0*fermi_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    else if(sig==1) ZMatrixMap(G.matptr(m), nao_, nao_) =  evec0* bose_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
  }

  // Fill the Propagator Ut
  ZMatrixMap(M.data(), nao_, nao_) = ZMatrixMap(iden.data(), nao_, nao_);
  for(int tstp=1; tstp<=k_; tstp++){
    propagator_exp(tstp, I, M.data(), ht, dt, nao_);
  }
  // Multiply by exp(i*mu*t)
  for(int tstp=1; tstp<=k_; tstp++) {
    ZMatrixMap(M.data() + tstp*es_, nao_, nao_) *= cplx(cos(mu*dt*tstp), sin(mu*dt*tstp));
  }
  
  
  // TV
  ZMatrixMap UtMap = ZMatrixMap(tmp.data(), nao_, nao_);
  ZMatrixMap tmp2Map = ZMatrixMap(tmp2.data(), nao_, nao_);
  ZMatrixMap UtMap2 = ZMatrixMap(NTauTmp.data(), nao_, nao_);
  
  for(n=0; n<=k_; n++) {
    UtMap = ZMatrixMap(M.data() + n*es_, nao_, nao_);

    for(m=0; m<=ntau_; m++) {
      tau = beta/2. * (Convolution().collocation().x_i()(m) + 1.);

      if(sig==-1) ZMatrixMap(G.tvptr(n,m), nao_, nao_) = cplx(0.,1.) * UtMap * evec0 * fermi_exp(beta,tau,eval0m).asDiagonal() * evec0.adjoint();
      else if(sig==1) ZMatrixMap(G.tvptr(n,m), nao_, nao_) = cplx(0.,1.) * UtMap * evec0 * bose_exp(beta,tau,eval0m).asDiagonal() * evec0.adjoint();
    }
  }
  
  // Ret and Les
  if(sig==-1)     tmp2Map = evec0 * fermi(beta,eval0m).asDiagonal() * evec0.adjoint();
  else if(sig==1) tmp2Map =-evec0 * bose(beta,eval0m).asDiagonal() * evec0.adjoint();

  for(m=0; m<=k_; m++) {
    UtMap = ZMatrixMap(M.data()+ m*es_, nao_, nao_);
    
    for(n=0; n<=m; n++) {
      UtMap2 = ZMatrixMap(M.data() + n*es_, nao_, nao_);
      ZMatrixMap(G.lesptr(n,m), nao_, nao_) = cplx(0.,1.) * UtMap2 * tmp2Map * UtMap.adjoint();
      ZMatrixMap(G.retptr(m,n), nao_, nao_) = cplx(0.,-1.) * UtMap * UtMap2.adjoint();
    }
  }
}


// Gives free green's function from constant hamiltonian
// G^M(\tau) = s f_s(mu-h) exp((mu-h)\tau)
// G^{TV}(n,\tau) = -is U_{n,0} f_s(h-mu) exp((h-mu)\tau)
// G^R(n,j) = -iU_{n,j} = U_{n,0} (U_{j,0})^\dagger
// G^L(j,n) = -siU_{j,0} f_s(h-mu) (U_{n,0})^\dagger
void dyson::G0_from_h0(GREEN &G, double mu, const double *H0, double beta, double h) const {
  assert(G.nt() == nt_);
  assert(G.nt() >= k_);
  assert(G.ntau() == ntau_);
  assert(G.size1() == nao_);

  int sign=G.sig();
  double tau,t;
  
  // Make Hamiltonian and solve eigen problem
  ZMatrix Hmu = mu*ZMatrixMap(iden.data(), nao_, nao_) - DMatrixConstMap(H0, nao_, nao_);
  Eigen::SelfAdjointEigenSolver<ZMatrix> eigensolver(Hmu);
  
  ZMatrix evec0(nao_, nao_);
  DColVector eval0(nao_), eval0m(nao_);
  
  evec0 = eigensolver.eigenvectors();
  eval0 = eigensolver.eigenvalues();
  eval0m = -eval0;

  // Matsubara
  for(int m=0; m<=ntau_; m++) {
    tau = beta/2. * (Convolution().collocation().x_i()(m) + 1.);
    if(sign==-1){
      ZMatrixMap(G.matptr(m), nao_, nao_) = -evec0 * fermi_exp(beta,tau,eval0).asDiagonal() * evec0.adjoint();
    } else if(sign==1){
      ZMatrixMap(G.matptr(m), nao_, nao_) =  evec0 * bose_exp(beta,tau,eval0).asDiagonal() * evec0.adjoint();
    }
  }

  // Ut
  ZMatrixMap IHdt = ZMatrixMap(tmp.data(), nao_, nao_);
  ZMatrixMap Udt = ZMatrixMap(tmp2.data(), nao_, nao_);
  IHdt = std::complex<double>(0,1.0) * h * Hmu;
  Udt = IHdt.exp();

  ZMatrixMap(M.data(), nao_, nao_) = ZMatrixMap(iden.data(), nao_, nao_);
  for(int n=1; n<=k_; n++) {
    ZMatrixMap(M.data() + n*es_, nao_, nao_) = ZMatrixMap(M.data() + (n-1)*es_, nao_, nao_) * Udt;
  }

  // TV
  for(int m=0; m <= ntau_; m++) {
    tau = beta/2. * (Convolution().collocation().x_i()(m) + 1.);

    for(int n=0; n<=k_; n++) {
      ZMatrixMap UtMap = ZMatrixMap(M.data() + n*es_, nao_, nao_);
      if(sign==-1){
        ZMatrixMap(G.tvptr(n,m), nao_, nao_) = std::complex<double>(0,1.0) *  UtMap * evec0 * fermi_exp(beta,tau,eval0m).asDiagonal() * evec0.adjoint();
      } else if(sign==1){
        ZMatrixMap(G.tvptr(n,m), nao_, nao_) = std::complex<double>(0,-1.0) * UtMap * evec0 * bose_exp(beta,tau,eval0m).asDiagonal() * evec0.adjoint();
      } 
    }
  } 
  
  // Ret and Less
  ZMatrixMap value = ZMatrixMap(tmp.data(), nao_, nao_);
  if(sign==-1){
    value = evec0*fermi(beta,eval0m).asDiagonal()*evec0.adjoint();
  }else if(sign==1){
    value = -1.0*evec0*bose(beta,eval0m).asDiagonal()*evec0.adjoint();
  }

  for(int m=0; m<=k_; m++) {
    for(int n=0; n<=m; n++) {
      ZMatrixMap Ut1 = ZMatrixMap(M.data() + m*es_, nao_, nao_);
      ZMatrixMap Ut2 = ZMatrixMap(M.data() + n*es_, nao_, nao_);
      ZMatrixMap(G.retptr(m,n), nao_, nao_) = std::complex<double>(0,-1.0) * Ut1 * Ut2.adjoint();
      ZMatrixMap(G.lesptr(n,m), nao_, nao_) = std::complex<double>(0,1.0) * Ut2 * value * Ut1.adjoint();
    }
  }
}


void dyson::G0_from_h0(GREEN &G, double mu, const DTensor<2> &H0, double beta, double h) const {
  assert(H0.shape()[0] == nao_);
  assert(H0.shape()[1] == nao_);
  assert(G.size1() == nao_);

  G0_from_h0(G, mu, H0.data(), beta, h);
}


// Gives free green's function from constant hamiltonian
// G^M(\tau) = s f_s(mu-h) exp((mu-h)\tau)
// G^{TV}(n,\tau) = -is U_{n,0} f_s(h-mu) exp((h-mu)\tau)
// G^R(n,j) = -iU_{n,j} = U_{n,0} (U_{j,0})^\dagger
// G^L(j,n) = -siU_{j,0} f_s(h-mu) (U_{n,0})^\dagger
void dyson::G0_from_h0(TTI_GREEN &G, double mu, const double *H0, double beta, double h) const {
  assert(G.nt() == nt_);
  assert(G.nt() >= k_);
  assert(G.ntau() == ntau_);
  assert(G.size1() == nao_);

  int sign=G.sig();
  double tau, t;
  
  ZMatrix Hmu = mu*ZMatrixMap(iden.data(), nao_, nao_) - DMatrixConstMap(H0, nao_, nao_);
  Eigen::SelfAdjointEigenSolver<ZMatrix> eigensolver(Hmu);  
  ZMatrix evec0(nao_, nao_);
  DColVector eval0(nao_), eval0m(nao_);
  
  evec0 = eigensolver.eigenvectors();
  eval0 = eigensolver.eigenvalues();
  eval0m = -eval0;

  for(int m=0; m<=ntau_; m++) {
    tau = beta/2. * (Convolution().collocation().x_i()(m) + 1.);
    if(sign==-1){
      ZMatrixMap(G.matptr(m), nao_, nao_) = -evec0 * fermi_exp(beta,tau,eval0).asDiagonal() * evec0.adjoint();
    } else if(sign==1){
      ZMatrixMap(G.matptr(m), nao_, nao_) =  evec0 * bose_exp(beta,tau,eval0).asDiagonal() * evec0.adjoint();
    }
  }

  ZMatrixMap IHdt = ZMatrixMap(tmp.data(), nao_, nao_);
  ZMatrixMap Udt = ZMatrixMap(tmp2.data(), nao_, nao_);
  IHdt = std::complex<double>(0,1.0) * h * Hmu;
  Udt = IHdt.exp();

  ZMatrixMap(M.data(), nao_, nao_) = ZMatrixMap(iden.data(), nao_, nao_);
  for(int n=1; n<=k_; n++) {
    ZMatrixMap(M.data() + n*es_, nao_, nao_) = ZMatrixMap(M.data() + (n-1)*es_, nao_, nao_) * Udt;
  }

  for(int m=0; m<=ntau_; m++) {
    tau = beta/2. * (Convolution().collocation().x_i()(m) + 1.);

    for(int n=0; n<=k_; n++) {
      ZMatrixMap UtMap = ZMatrixMap(M.data() + n*es_, nao_, nao_);
      if(sign==-1){
        ZMatrixMap(G.tvptr(n,m), nao_, nao_) = std::complex<double>(0,1.0) *  UtMap * evec0 * fermi_exp(beta,tau,eval0m).asDiagonal() * evec0.adjoint();
      } else if(sign==1){
        ZMatrixMap(G.tvptr(n,m), nao_, nao_) = std::complex<double>(0,-1.0) * UtMap * evec0 * bose_exp(beta,tau,eval0m).asDiagonal() * evec0.adjoint();
      } 
    }
  } 
  
  ZMatrixMap value = ZMatrixMap(tmp.data(), nao_, nao_);
  if(sign == -1) {
    value = evec0*fermi(beta,eval0m).asDiagonal()*evec0.adjoint();
  }else if(sign == 1) {
    value = -1.0*evec0*bose(beta,eval0m).asDiagonal()*evec0.adjoint();
  }

  for(int m=0; m<=k_; m++) {
      ZMatrixMap Ut1 = ZMatrixMap(M.data() + m*es_, nao_, nao_);
      ZMatrixMap(G.retptr(m),  nao_, nao_) = std::complex<double>(0,-1.0) * Ut1;
      ZMatrixMap(G.lesptr(-m), nao_, nao_) = std::complex<double>(0,1.0) * value * Ut1.adjoint();
  }
}

void dyson::G0_from_h0(TTI_GREEN &G, double mu, const DTensor<2> &H0, double beta, double h) const {
  G0_from_h0(G, mu, H0.data(), beta, h);
}

} // Namespace

#endif
