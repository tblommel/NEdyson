#ifndef DYSON_IMPL
#define DYSON_IMPL

#include <iostream>

#include "dyson.h"

namespace NEdyson{

#define EXPMAX 100

// Free functions=================================================================================

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
dvector fermi(double beta, dvector &omega){
   int size=omega.size();
   dvector tmp(size);
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
dvector fermi_exp(double beta,double tau,dvector &omega){
   int size=omega.size();
   dvector tmp(size);

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
dvector bose(double beta,dvector &omega){
   int size=omega.size();
   dvector tmp(size);
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
dvector bose_exp(double beta,double tau,dvector &omega){
   int size=omega.size();
   dvector tmp(size);

   for(int i=0;i<size;i++){
      tmp(i)=bose_exp(beta,tau,omega(i));
   }
   return tmp;
}



// gives an interpolation of a function on the domain [tstp-k,tstp]
// eps(tinterp) = sum_{l=0}^k eps(tstp-k+l) * sum_{p=0}^k (tinterp+k-tstp)^p P_{pl}
cdmatrix interpolate(int tstp, double tinterp, const INTEG &I, const function &eps){
  int size = eps.size1(), p,l,k=I.k();
  double timk=-tstp+k+tinterp;
  cdmatrix res(size,size), tmp(size,size);
  res.setZero();
  double weight,pref;

  for(l=0;l<=k;l++){
    pref = 1.;
    weight = I.poly_interp(0,l);
    for(p=1;p<=k;p++){
      pref*=timk;
      weight+=pref*I.poly_interp(p,l);
    }
    eps.get_value(tstp-k+l,tmp);
    res += weight*tmp;
  }
  return res;
}


// evaluates Ut(tstp) by the approximation exp(-I*(a1*eps(tstp-1+c1)+a2*eps(tstp-1+c2)))
// *exp(-I*(a1*eps(tstp-1+c2)+a2*eps(tstp-1+c1)))*Ut(tstp-1)
// a1=(3-2sqrt(3))/12
// a2=(3+2sqrt(3))/12
// c1=(1-sqrt(3)/3)/2
// c2=(1+sqrt(3)/3)/2
void propagator_exp(int tstp, const INTEG &I, function &Ut, const function &eps, double dt){
  int size = Ut.size1();
  cdmatrix eps1(size,size), eps2(size,size);
  cdmatrix arg1(size,size), arg2(size,size);
  
  int n = (tstp<I.k())?(I.k()):(tstp);

  eps1 = interpolate(n,tstp-1./2-sqrt(3.)/6.,I,eps);
  eps2 = interpolate(n,tstp-1./2+sqrt(3.)/6.,I,eps);
  double a1 = (3.-2*sqrt(3.))/12.;
  double a2 = (3.+2*sqrt(3.))/12.;
  arg1=-cplx(0.,dt)*(a1*eps1+a2*eps2);
  arg2=-cplx(0.,dt)*(a2*eps1+a1*eps2);

  cdmatrix res(size,size);
  Ut.get_value(tstp-1,res);
  res=arg1.exp()*arg2.exp()*res;
  Ut.set_value(tstp,res);
}


// Gives the free greens function from a non constant hamiltonian
void G0_from_h0(GREEN &G, const INTEG &I, double mu, const function &eps, double beta, double dt){
  int nt=G.nt(), ntau=G.ntau(), size=G.size1(), sig=G.sig(),m,n;
  double tau, t, dtau=beta/ntau;
  
  cdmatrix iden = cdmatrix::Identity(size,size);
  cdmatrix tmp(size,size), mHmu(size,size), tmp1(size,size), tmp2(size,size);
  eps.get_value(-1,tmp);
  mHmu=mu*iden-tmp;
  function Ut(nt,size);

  Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(mHmu);
  cdmatrix evec0(size,size),value(size,size);
  dvector eval0(size), eval0m(size);
  evec0=eigensolver.eigenvectors();
  eval0=eigensolver.eigenvalues();
  eval0m=-eval0;

  // Matsubara
  for(m=0;m<=ntau;m++){
    tau = m*dtau;
    if(sig==-1) value = -evec0*fermi_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    else if(sig==1) value = evec0*bose_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    G.set_mat(m,value);
  }

  // Fill the Propagator Ut
  Ut.set_value(-1,iden);
  Ut.set_value(0,iden);
  for(int tstp=1;tstp<=nt;tstp++){
    propagator_exp(tstp,I,Ut,eps,dt);
  }
  // Multiply by exp(i*mu*t)
  for(int tstp=1;tstp<=nt;tstp++){
    Ut.get_value(tstp,value);
    value=value*cplx(cos(mu*dt*tstp),sin(mu*dt*tstp));
    Ut.set_value(tstp,value);
  }
  
  // TV
  for(n=0;n<=nt;n++){
    Ut.get_value(n,tmp);
    for(m=0;m<=ntau;m++){
      tau=m*dtau;
      if(sig==-1) value = cplx(0.,1.)*tmp*evec0*fermi_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
      else if(sig==1) value = cplx(0.,1.)*tmp*evec0*bose_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
      G.set_tv(n,m,value);
    }
  }
  
  // Ret and Les
  if(sig==-1) value = evec0*fermi(beta,eval0m).asDiagonal()*evec0.adjoint();
  else if(sig==1) value = -evec0*bose(beta,eval0m).asDiagonal()*evec0.adjoint();
  for(m=0;m<=nt;m++){
    Ut.get_value(m,tmp1);
    for(n=0;n<=m;n++){
      Ut.get_value(n,tmp2);
      tmp = cplx(0.,1.)*tmp2*value*tmp1.adjoint();
      G.set_les(n,m,tmp);
      tmp = cplx(0.,-1.)*tmp1*tmp2.adjoint();
      G.set_ret(m,n,tmp);
    }
  }
}


// Gives free green's function from constant hamiltonian
// G^M(\tau) = s f_s(mu-h) exp((mu-h)\tau)
// G^{TV}(n,\tau) = -is U_{n,0} f_s(h-mu) exp((h-mu)\tau)
// G^R(n,j) = -iU_{n,j} = U_{n,0} (U_{j,0})^\dagger
// G^L(j,n) = -siU_{j,0} f_s(h-mu) (U_{n,0})^\dagger
void G0_from_h0(GREEN &G, double mu, const cdmatrix &H0, double beta, double h){
  int nt=G.nt(),ntau=G.ntau(), size = G.size1();
  int sign=G.sig();
  double tau,t,dtau=beta/ntau;
  cdmatrix idm(size,size);
  cdmatrix Udt(size,size);
  cdmatrix IHdt(size,size);
  cdmatrix Hmu(size,size);
  cdmatrix evec0(size,size),value(size,size);
  dvector eval0(size),eval0m(size);

  idm = Eigen::MatrixXcd::Identity(size,size);
  Hmu = -H0 + mu * idm;

  Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(Hmu);
  evec0=eigensolver.eigenvectors();
  eval0=eigensolver.eigenvalues();
  eval0m=(-1.0)*eval0;

  for(int m=0;m<=ntau;m++){
    tau=m*dtau;
    if(sign==-1){
      value=(-1.0)*evec0*fermi_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    }else if(sign==1){
      value=(1.0)*evec0*bose_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    }
    G.set_mat(m,value);
  }

  if(nt >=0 ){
    IHdt = std::complex<double>(0,-1.0) * h * H0;
    Udt = IHdt.exp();

    NEdyson::function Ut(nt,size);
    cdmatrix Un(size,size);
    Ut.set_value(-1,idm);
    Ut.set_value(0,idm);
    for(int n=1;n<=nt;n++){
      Ut.get_value(n-1,Un);
      Un = Un * Udt;
      Ut.set_value(n,Un);
    }

    cdmatrix expp(size,size);
    for(int m=0;m<=ntau;m++){
      tau=m*dtau;
      for(int n=0;n<=nt;n++){
        Ut.get_value(n,expp);
        if(sign==-1){
          value=std::complex<double>(0,1.0)*expp*evec0*fermi_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }else if(sign==1){
          value=std::complex<double>(0,-1.0)*expp*evec0*bose_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        } 
        G.set_tv(n,m,value);
      } 
    } 
    
    if(sign==-1){
      value=evec0*fermi(beta,eval0m).asDiagonal()*evec0.adjoint();
    }else if(sign==1){
      value=-1.0*evec0*bose(beta,eval0m).asDiagonal()*evec0.adjoint();
    } 
    cdmatrix exppt1(size,size);
    cdmatrix exppt2(size,size);
    for(int m=0;m<=nt;m++){
      for(int n=0;n<=m;n++){
        cdmatrix tmp(size,size);
        Ut.get_value(m,exppt1);
        Ut.get_value(n,exppt2);
        tmp = std::complex<double>(0,-1.0)*exppt1*exppt2.adjoint();
        G.set_ret(m,n,tmp);
        tmp=std::complex<double>(0,1.0)*exppt2*value*exppt1.adjoint();
        G.set_les(n,m,tmp);
      }
    }
  }
}



// Matsubara solver =====================================================================================
// Solves Dysons equation in frequency space
double mat_fourier(GREEN &G, const GREEN &Sigma, double mu, const cplx *h0, double beta){
  cplx *sigmadft, *sigmaiomn, *z1, *z2, *one;
  cplx *gmat, *hj, iomn, *zinv;
  int ntau, m, r, pcf, p, m2, sg, ss, l, sig, size1=G.size1();
  double dtau;

  assert(G.ntau()==Sigma.ntau());
  sig=G.sig();
  assert(G.sig()==Sigma.sig());
  sg=G.element_size();
  ss=Sigma.element_size();
  ntau=G.ntau();
  dtau=beta/ntau;
  
  if(ntau%2==1){
    std::cerr << "must have ntau even" <<std::endl;
    abort();
  }
  sigmadft=new cplx[(ntau+1)*ss];
  sigmaiomn = new cplx[ss];
  gmat = new cplx[(ntau + 1) * sg];
  z1 = new cplx[sg];
  z2 = new cplx[sg];
  hj = new cplx[sg];
  one = new cplx[sg];
  zinv = new cplx[sg];
  element_iden(size1, one);
  element_set(size1,hj,h0);
  for(l=0;l<sg;l++) hj[l]-=mu*one[l];

  pcf=10;
  m2=ntau/2;
  matsubara_dft(sigmadft,Sigma,sig);
  set_first_order_tail(gmat,one,beta,sg,ntau,sig,size1);
  for(p=-pcf;p<=pcf;p++){
    int l=0,h=0;
    if(sig==1&&p==0) h=1;
    else if(sig==1&&p>0){
      h=1;
      l=-1;
    }
    for(m=-m2-l;m<=m2-1+h;m++){
      iomn=cplx(0,get_omega(m+p*ntau,beta,sig));
      matsubara_ft(sigmaiomn,m+p*ntau,Sigma,sigmadft,sig,beta);
      element_set(size1,z1,sigmaiomn);
      element_incr(size1,z1,hj);
      if(sig==1&&m+p*ntau==0){//zero frequency for bosons
        element_inverse(size1,zinv,z1);
        element_set(size1,zinv,z1);
        element_smul(size1,z2,-1.0);
      }
      else{
        for(l=0;l<sg;l++){
          z2[l]=iomn*one[l]-z1[l];
        }
        element_inverse(size1,zinv,z2);
        for(l=0;l<sg;l++){
          zinv[l]-=1./iomn*one[l];
        }
      }
      element_smul(size1,zinv,1./beta);
      for(r=0;r<=ntau;r++){
        cplx expfac=std::exp(-get_tau(r,beta,ntau)*iomn);
        for(l=0;l<sg;l++){
          gmat[r*sg+l]+=zinv[l]*expfac;
        }
      }
    }
  }
  double err=0;
  for(r=0;r<=ntau;r++){
    err += element_diff(size1,G.matptr(r),gmat+r*sg);
    element_set(size1,G.matptr(r),gmat+r*sg);
  }

  force_mat_herm(G);

  delete[] sigmadft;
  delete[] sigmaiomn;
  delete[] gmat;
  delete[] z1;
  delete[] z2;
  delete[] hj;
  delete[] one;
  delete[] zinv;
  return err;
}



// Extrapolation ==========================================================================================
// Extrapolates a Green's function object from [n-k-1,n-1] to n
void Extrapolate(const INTEG &I, GREEN &G, int n){
  int k=I.k(), size1=G.size1(), ntau=G.ntau(),l,j,es=size1*size1,jcut;
  double *pref= new double[k+1];
  for(l=0;l<=k;l++) pref[l]=0;
  for(l=0;l<=k;l++){
    for(j=0;j<=k;j++){
      pref[l]+=I.poly_interp(j,l)*(1-2*(j%2));
    }
  }
  cplx *sav;
  cplx *tmp=new cplx[es]; 
  //right mixing
  for(l=0;l<=ntau;l++){
    sav=G.tvptr(n,l);
    element_set_zero(size1,sav);
    for(j=0;j<=k;j++) element_incr(size1,sav,pref[j],G.tvptr(n-j-1,l));
  }
  //retarded
  for(l=0;l<n-k;l++){
    sav=G.retptr(n,n-l);
    element_set_zero(size1,sav);
    for(j=0;j<=k;j++){
      element_incr(size1,sav,pref[j],G.retptr(n-j-1,n-l-j-1));
    }
  }
  for(l=0;l<=k;l++){
    jcut=(l<=n-k-1)?k:(n-l-1);
    sav=G.retptr(n,l);
    element_set_zero(size1,sav);
    for(j=0;j<=jcut;j++){
      element_incr(size1,sav,pref[j],G.retptr(n-j-1,l));
    }
    for(j=jcut+1;j<=k;j++){
      element_conj(size1,tmp,G.retptr(l,n-j-1));
      element_incr(size1,sav,-1.*pref[j],tmp);
    }
  }
  //less
  element_set(size1,G.lesptr(0,n),G.tvptr(n,0));
  for(l=0;l<=k;l++){
    jcut=(k>n-l-1)?(n-l-1):(k);
    sav=G.lesptr(l,n);
    element_set_zero(size1,sav);
    for(j=0;j<=jcut;j++){
      element_incr(size1,sav,pref[j],G.lesptr(l,n-1-j));
    }
    for(j=jcut+1;j<=k;j++){
      element_conj(size1,tmp,G.lesptr(n-1-j,l));
      element_incr(size1,sav,-1.*pref[j],tmp);
    }
  }
  for(l=k+1;l<=n;l++){
    sav=G.lesptr(l,n);
    element_set_zero(size1,sav);
    for(j=0;j<=k;j++){
      element_incr(size1,sav,pref[j],G.lesptr(l-j-1,n-j-1));
    }
  }
  delete[] pref;
  delete[] tmp;
}





// Start Functions =======================================================================================
double GRstart(const INTEG &I, GREEN &G, const GREEN &Sig, const CFUNC &hmf, double mu, double dt){
  // Counters and sizes
  int k=I.k(), size1=G.size1(), es=G.element_size(),m,l,n,i;
  
  // Matricies
  cplx *M = new cplx[k*k*es];
  cplx *Q = new cplx[k*es];
  cplx *X = new cplx[k*es];
  cplx *tmp = new cplx[es];
  cplx *stmp = new cplx[es];
  cplx *iden = new cplx[es];
  cplx weight;
  cplx ncplxi = cplx(0,-1);
  element_iden(size1,iden);

  // Initial condition
  for(i=0;i<=k;i++){
    element_iden(size1,ncplxi,G.retptr(i,i));
  }
  double err=0;
  // Fill the first k timesteps
  for(m=0;m<k;m++){
    for(l=0;l<k*k*es;l++){
      M[l]=0.;
    }
    for(l=0;l<k*es;l++){
      Q[l]=0.;
    }
    for(n=m+1;n<=k;n++){
      for(l=0;l<=k;l++){
        if(l<=m){ // We know these G's. Put into Q
          element_conj(size1,tmp,G.retptr(m,l));
          element_smul(size1,tmp,-1);
          element_mult(size1,stmp,Sig.retptr(n,l),tmp);
          for(i=0;i<es;i++){
            Q[(n-m-1)*es+i]+=ncplxi/dt*I.poly_diff(n,l)*tmp[i]+dt*I.poly_integ(m,n,l)*stmp[i];
          }
        }
        else{ // Don't have these. Put into M
          // Derivative term
          for(i=0;i<es;i++) M[es*((n-m-1)*(k-m))+(i/size1)*(k-m)*(size1)+(l-m-1)*size1+i%size1] = -ncplxi/dt*I.poly_diff(n,l)*iden[i];

          // Delta energy term
          if(n==l){
            element_set(size1,tmp,hmf.ptr(l));
            for(i=0;i<es;i++) M[es*((n-m-1)*(k-m))+(i/size1)*(k-m)*(size1)+(l-m-1)*size1+i%size1] += mu*iden[i]-tmp[i];
          }

          // Integral term
          weight=dt*I.poly_integ(m,n,l);
          if(n>=l){ // We have Sig
            element_set(size1,stmp,Sig.retptr(n,l));
          }
          else{ // Don't have it
            element_set(size1,stmp,Sig.retptr(l,n));
            element_conj(size1,stmp);
            weight *= -1;
          }
          for(i=0;i<es;i++){
            M[es*((n-m-1)*(k-m))+(i/size1)*(k-m)*(size1)+(l-m-1)*size1+i%size1] -= weight*stmp[i];
          }
        }
      }
    }
    //solve MX=Q for X
    element_linsolve_left((k-m)*size1,(k-m)*size1,size1,M,X,Q);
    //put X into G
    for(l=0;l<k-m;l++){
      err += element_diff(size1,G.retptr(l+m+1,m),X+l*es);
      element_set(size1,G.retptr(l+m+1,m),X+l*es);
    }
  }

  delete[] M;
  delete[] X;
  delete[] Q;
  delete[] tmp;
  delete[] stmp;
  delete[] iden;
  return err;
}


double GTVstart(const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  // Counters and sizes
  int k=I.k(), size1=G.size1(), es=G.element_size(),ntau=G.ntau(),m,l,n,i;
  cplx weight;

  // Matricies
  cplx *M = new cplx[k*k*es];
  cplx *Q = new cplx[k*es];
  cplx *X = new cplx[k*es];
  cplx *tmp = new cplx[es];
  cplx *stmp = new cplx[es];
  cplx *iden = new cplx[es];

  cplx cplxi = cplx(0,1);
  element_iden(size1,iden);

  // Boundary Conditions
  for(m=0;m<=ntau;m++){
    for(i=0;i<es;i++){
      G.tvptr(0,m)[i] = (double)G.sig()*cplxi*G.matptr(ntau-m)[i];
    }
  }
  // Do the Q integrals and store in Gtv(n,m)
  for(n=1;n<=k;n++){
    for(m=0;m<=ntau;m++){
      CTV2(I, Sig, G, n, m, beta, tmp);
      CTV3(I, Sig, G, n, m, beta, stmp);
      for(i=0;i<es;i++){
        G.tvptr(n,m)[i] = tmp[i] + stmp[i];
      }
    }
  }

  // At each m, get n=1...k
  for(m=0;m<=ntau;m++){
    memset(M,0,k*k*es*sizeof(cplx));
    memset(Q,0,k*es*sizeof(cplx));

    // Set up the kxk linear problem MX=Q
    for(n=1;n<=k;n++){
      for(l=0;l<=k;l++){

        // Derivative term
        weight = cplxi*I.poly_diff(n,l)/dt;
        if(l==0){ // Put into Q
          for(i=0;i<es;i++){
            Q[(n-1)*es+i] -= weight*G.tvptr(0,m)[i];
          }
        }
        else{ // Put into M
          for(i=0;i<es;i++){
            M[(n-1)*es*k+(i/size1)*k*size1+(l-1)*size1+i%size1] += weight*iden[i];
          }
        }

        // Delta energy term
        if(l==n){
          element_set(size1,tmp,hmf.ptr(l));
          for(i=0;i<es;i++) M[es*(n-1)*k+(i/size1)*k*size1+(l-1)*size1+i%size1] += mu*iden[i]-tmp[i];
        }

        // Integral term
        weight = -dt*I.gregory_weights(n,l);
        if(l==0){ // Put into Q
          element_incr(size1,Q+(n-1)*es,-weight,Sig.retptr(n,l),G.tvptr(l,m));
        }
        else{ // Put into M
          if(n>=l){ // Have Sig
            element_set(size1,stmp,Sig.retptr(n,l));
          }
          else{ // Dont have Sig
            element_set(size1,stmp,Sig.retptr(l,n));
            element_conj(size1,stmp);
            element_smul(size1,stmp,-1);
          }
          for(i=0;i<es;i++){
            M[es*(n-1)*k+(i/size1)*k*size1+(l-1)*size1+i%size1] += weight*stmp[i];
          }
        }
      }
      // Add in the integrals
      CTV2(I, Sig, G, n, m, beta, tmp);
      CTV3(I, Sig, G, n, m, beta, stmp);
      element_incr(size1,Q+(n-1)*es,tmp);
      element_incr(size1,Q+(n-1)*es,stmp);
    }
    // Solve MX=Q
    element_linsolve_left(k*size1,k*size1,size1,M,X,Q);
    for(l=0;l<k;l++){
      element_set(size1,G.tvptr(l+1,m),X+l*es);
    }
  }

  delete[] M;
  delete[] Q;
  delete[] X;
  delete[] tmp;
  delete[] stmp;
  delete[] iden;
}






// Step Functions ==========================================================================================

// i dt' GR(t,t-t') - GR(t-t')hmf(t-t') - \int_0^t' GR(t,t-s) SR(t-s,t-t') = 0
void GRstep(int tstp, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double dt){
  
  // Counters and sizes
  int k=I.k(), size1=G.size1(), es=G.element_size(),m,l,n,i;
  
  // Matricies
  cplx *M = new cplx[k*k*es];
  cplx *Q = new cplx[k*es];
  cplx *X = new cplx[k*es];
  cplx *tmp = new cplx[es];
  cplx *stmp = new cplx[es];
  cplx *iden = new cplx[es];
  cplx *qqint = new cplx[es*(tstp+1)];
  cplx weight;
  cplx ncplxi = cplx(0,-1);
  element_iden(size1,iden);

  // Initial condition
  element_iden(size1,ncplxi,G.retptr(tstp,tstp));

  // Fill the first k timesteps
  for(l=0;l<k*k*es;l++){
    M[l]=0.;
  }
  for(l=0;l<k*es;l++){
    Q[l]=0.;
  }
  for(n=1;n<=k;n++){
    for(l=0;l<=k;l++){
      // Derivative terms
      weight = -ncplxi/dt*I.poly_diff(n,l);
      if(l==0){ // Goes into Q
        for(i=0;i<es;i++) Q[(i/size1)*k*size1+(n-1)*size1+i%size1] += -weight*G.retptr(tstp,tstp)[i];
      }
      else{ // Goes into M
        for(i=0;i<es;i++) M[es*k*(l-1)+(i/size1)*k*size1+(n-1)*size1+i%size1] += weight*iden[i];
      }

      // Delta Energy term
      if(l==n){
        element_set(size1,tmp,hmf.ptr(tstp-n));
        for(i=0;i<es;i++) M[es*k*(l-1)+(i/size1)*k*size1+(n-1)*size1+i%size1] += mu*iden[i] - tmp[i];
      }

      // Integral term
      weight = dt*I.gregory_weights(n,l);
      if(tstp-l>=tstp-n){ // We have Sig
        element_set(size1,stmp,Sig.retptr(tstp-l,tstp-n));
      }
      else{ // Don't have it
        element_set(size1,stmp,Sig.retptr(tstp-n,tstp-l));
        element_conj(size1,stmp);
        weight*=-1;
      }
      if(l==0){ // Goes into Q
        element_mult(size1,tmp,G.retptr(tstp,tstp),stmp);
        for(i=0;i<es;i++) Q[(i/size1)*k*size1+(n-1)*size1+i%size1] += weight*tmp[i];
      }
      else{ // Into M
        for(i=0;i<es;i++){
          M[es*k*(l-1)+(i/size1)*k*size1+(n-1)*size1+i%size1] -= weight*stmp[i];
        }
      }
    }
  }
  
  // Solve XM=Q for X
  element_linsolve_right(size1,k*size1,k*size1,X,M,Q);
  
      
  // Put X into G
  cplx *Gp;
  for(l=0;l<k;l++){
    Gp=G.retptr(tstp,tstp-l-1);
    for(i=0;i<es;i++){
      Gp[i] = X[(i/size1)*k*size1+l*size1+i%size1];
    }
  }

  // Start doing the integrals
  // qqint_n = \sum_{l=0}^{n-1} w_{nl} GR(t,t-l) SR(t-l,t-n) for n=k+1...tstp
  memset(qqint,0,sizeof(cplx)*es*(tstp+1));
  for(n=k+1;n<=tstp;n++){
    for(l=0;l<=k;l++){
      element_incr(size1,qqint+n*es,I.gregory_weights(n,l),G.retptr(tstp,tstp-l),Sig.retptr(tstp-l,tstp-n));
    }
  }

  cplx *bdweight = new cplx[k+2];
  for(l=0;l<k+2;l++) bdweight[l] = I.bd_weights(l)*-ncplxi/dt;
  double w0 = dt*I.omega(0);

  for(n=k+1;n<=tstp;n++){
    // Set up qq
    for(i=0;i<es;i++){
      qqint[n*es+i]*=dt;
      for(l=1;l<=k+1;l++){
        qqint[n*es+i]-=bdweight[l]*G.retptr(tstp,tstp-n+l)[i];
      }
    }

    // Set up mm
    element_set(size1,M,hmf.ptr(tstp-n));
    element_smul(size1,M,-1);
    element_incr(size1,M,-w0,Sig.retptr(tstp-n,tstp-n));
    for(i=0;i<es;i++) M[i] += (mu+bdweight[0])*iden[i];

    // Solve XM=Q for X
    element_linsolve_right(size1,size1,size1,X,M,qqint+n*es);
    
    // Put X into G
    element_set(size1,G.retptr(tstp,tstp-n),X);

    // Add this newly computed value to the integrals which need it
    for(m=n+1;m<=tstp;m++) element_incr(size1,qqint+m*es,I.gregory_weights(m,n),G.retptr(tstp,tstp-n),Sig.retptr(tstp-n,tstp-m));
  }

  delete[] bdweight;
  delete[] qqint;
  delete[] M;
  delete[] X;
  delete[] Q;
  delete[] tmp;
  delete[] stmp;
  delete[] iden;
  return;
}


// idt GRM(tstp,m) - hmf(tstp) GRM(t,m) - \int_0^t dT SR(t,T) GRM(T,m) = \int_0^{beta} dT SRM(t,T) GM(T-m)
void GTVstep(int tstp, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  
  // Counters and sizes
  int k=I.k(), size1=G.size1(), es=G.element_size(), ntau=G.ntau(), m, l, n, i;
  cplx weight;
  cplx weight2;

  // Matricies
  cplx *iden = new cplx[es];
  cplx *M = new cplx[es];
  cplx *Q = new cplx[es];
  cplx *tmp = new cplx[es];
  cplx *stmp = new cplx[es];
  cplx *resptr;
  cplx *Gptr;

  element_iden(size1,iden);
  cplx cplxi = cplx(0.,1.);
  
  // First do the SRM GM integral.  Put it into GRM(tstp,m)
  for(m=0;m<=ntau;m++){
    CTV2(I, Sig, G, tstp, m, beta, tmp);
    CTV3(I, Sig, G, tstp, m, beta, stmp);
    for(i=0;i<es;i++){
      G.tvptr(tstp,m)[i] = tmp[i] + stmp[i];
    }
  }

  // Next do the SR GRM integral.  Adds to GRM(tstp,m)
  CTV1(I, Sig, Sig, G, tstp, dt);

  // Put derivatives into GRM(tstp,m)
  for(l=1;l<=k+1;l++){
    weight = -cplxi/dt*I.bd_weights(l);
    resptr=G.tvptr(tstp,0);
    Gptr=G.tvptr(tstp-l,0);
    for(m=0;m<(ntau+1)*es;m++){
      resptr[m] += weight*Gptr[m];
    }
  }
  
  // Make M
  weight = cplxi/dt*I.bd_weights(0);
  resptr = hmf.ptr(tstp);
  weight2= -dt*I.omega(0);
  Gptr = Sig.retptr(tstp,tstp);
  for(i=0;i<es;i++) M[i] = (weight+mu)*iden[i] - resptr[i] + weight2*Gptr[i];

  // Solve MX=Q
  for(m=0;m<=ntau;m++){
    element_set(size1,Q,G.tvptr(tstp,m));
    element_linsolve_left(size1,size1,size1,M,G.tvptr(tstp,m),Q);
  }

  delete[] stmp;
  delete[] tmp;
  delete[] iden;
  delete[] M;
  delete[] Q;
}


void GLstep(int n, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  // Sizes and iterators
  int k=I.k(), size1=G.size1(), es=G.element_size(), m,l,i;
  cplx weight;
  int num = (n>=k)?(n):(k);

  // Matricies
  cplx *M = new cplx[k*k*es];
  cplx *X = new cplx[(num+1)*es];
  cplx *Q = new cplx[(num+1)*es];
  cplx *tmp = new cplx[es];
  cplx *stmp = new cplx[es];
  cplx *iden = new cplx[es];
  cplx cplxi = cplx(0,1);
  element_iden(size1,iden);

  // Initial condition
  element_set(size1,G.lesptr(0,n),G.tvptr(n,0));
  element_conj(size1,G.lesptr(0,n));
  element_smul(size1,G.lesptr(0,n),-1);
  
  // Integrals go into Q
  memset(Q,0,sizeof(cplx)*num*es);
  Cles2_tstp(I,Sig,Sig,G,n,dt,Q);
  Cles3_tstp(I,Sig,G,n,beta,Q);

  // Set up the kxk linear problem MX=Q
  for(m=1;m<=k;m++){
    for(l=0;l<=k;l++){
      
      // Derivative term
      weight = cplxi/dt*I.poly_diff(m,l);
      if(l==0){ // We put this in Q
        for(i=0;i<es;i++){
          Q[(m)*es+i] -= weight*G.lesptr(l,n)[i];
        }
      }
      else{ // It goes into M
        for(i=0;i<es;i++){
          M[(m-1)*es*k+(i/size1)*k*size1+(l-1)*size1+i%size1] += weight*iden[i];
        }
      }

      // Delta energy term
      if(m==l){
        element_set(size1,tmp,hmf.ptr(l));
        for(i=0;i<es;i++){
          M[es*(m-1)*k+(i/size1)*k*size1+(l-1)*size1+i%size1] += mu*iden[i]-tmp[i];
        }
      }

      // Integral term
      weight = -dt*I.gregory_weights(m,l);
      if(l==0){ // Goes into Q
        element_incr(size1,Q+(m)*es,-weight,Sig.retptr(m,l),G.lesptr(l,n));
      }
      else{ // Goes into M
        if(m>=l){ // Have Sig
          element_set(size1,stmp,Sig.retptr(m,l));
        }
        else{ // Dont have Sig
          element_set(size1,stmp,Sig.retptr(l,m));
          element_conj(size1,stmp);
          element_smul(size1,stmp,-1);
        }
        for(i=0;i<es;i++){
          M[es*(m-1)*k+(i/size1)*k*size1+(l-1)*size1+i%size1] += weight*stmp[i];
        }
      }
    }
  }

  // Solve Mx=Q
  element_linsolve_left(k*size1,k*size1,size1,M,X+es,Q+es);
  element_set(size1,X,G.lesptr(0,n));

  cplx bd0 = cplxi/dt*I.bd_weights(0) + mu;
  cplx bdl;
  double nhw = -dt*I.omega(0);
  cplx *sigmm;

  // Timestepping
  for(m=k+1;m<=n;m++){
    // Set up M
    sigmm=Sig.retptr(m,m);
    element_set(size1,M,hmf.ptr(m));
    element_smul(size1,M,-1);
    for(i=0;i<es;i++) M[i] += bd0*iden[i] + nhw*sigmm[i];

    // Derivatives into Q
    for(l=1;l<=k+1;l++){
      bdl = -cplxi/dt*I.bd_weights(l);
      element_incr(size1,Q+(m)*es,bdl,X+(m-l)*es);
    }

    // Rest of the retles integral
    for(l=0;l<m;l++){
      element_incr(size1,Q+(m)*es,dt*I.gregory_weights(m,l),Sig.retptr(m,l),X+l*es);
    }

    // Solve MX=Q
    element_linsolve_left(size1,size1,size1,M,X+m*es,Q+(m)*es);
  }

  // Write elements into G
  for(l=1;l<=n;l++) element_set(size1,G.lesptr(l,n),X+l*es);

  
  delete[] M;
  delete[] Q;
  delete[] X;
  delete[] tmp;
  delete[] stmp;
  delete[] iden;
}


void GLstart(const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  int k=I.k();
  for(int n=0;n<=k;n++) GLstep(n,I,G,Sig,hmf,mu,beta,dt);
}

double dyson_start(const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  double err=0;
  err += GRstart(I, G, Sig, hmf, mu, dt);
  err += GTVstart(I, G, Sig, hmf, mu, beta, dt);
  err += GLstart(I, G, Sig, hmf, mu, beta, dt);
  return err;
}



void dyson_step(int n, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  GRstep(n, I, G, Sig, hmf, mu, dt);
  GTVstep(n, I, G, Sig, hmf, mu, beta, dt);
  GLstep(n, I, G, Sig, hmf, mu, beta, dt);
}




}//namespace
#endif
