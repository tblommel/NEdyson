#ifndef DYSON_OMP_IMPL
#define DYSON_OMP_IMPL

#include "dyson_omp.h"

namespace NEdyson{

#if USE_OMP == 1

void dyson_step_omp_ret(int threads, int tstp, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double dt){
  assert(G.size1()==Sig.size1());
  assert(G.size1()==hmf.size1());
  assert(G.nt()==Sig.nt());
  assert(G.nt()==hmf.nt());
  assert(G.nt()>=I.k());
  assert(G.sig()==Sig.sig());
  assert(tstp>I.k());
  
  int k=I.k(), size1=G.size1(), es=G.element_size(), l,n,i;
  cplx ncplxi = cplx(0,-1), weight;

  // Set timestep to zero
  memset(G.retptr(tstp,0),0,(tstp+1)*es*sizeof(cplx));

  // Initial condition
  element_iden(size1, G.retptr(tstp,tstp), ncplxi);

  // Do first k values
  // Exact same as serial implementation
  {
    cplx *M = new cplx[k*k*es];
    cplx *Q = new cplx[k*es];
    cplx *X = new cplx[k*es];
    cplx *tmp = new cplx[es];
    cplx *stmp = new cplx[es];
    cplx *iden = new cplx[es];
    element_iden(size1,iden);

    memset(M,0,k*k*es*sizeof(cplx));
    memset(Q,0,k*es*sizeof(cplx));
    
    for(n=1;n<=k;n++){
      for(l=0;l<=k;l++){
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


    delete[] M;
    delete[] Q;
    delete[] X;
    delete[] tmp;
    delete[] stmp;
    delete[] iden;
  }

  // Now we will use the non-conjugate equation to evaluate the remaining values of G
  // I d/dt GR(t,t') = ...

#pragma omp parallel num_threads(threads)
  {
    int j,m;
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    std::vector<bool> mask_ret(tstp+1, false);
    cplx *bd = new cplx[k+2];
    cplx *Q = new cplx[es];
    cplx *M = new cplx[es];
    for(j=0;j<=k+1;j++) bd[j] = I.bd_weights(j)*ncplxi/dt;

    // Each thread is assigned every nthreads^th timeslice
    for(j=0;j<tstp-k;j++){
      if(j%nthreads == tid){
        mask_ret[j] = true;
      }
    }
    incr_convolution_ret(tstp, mask_ret, G, Sig, Sig, G, G, I, dt);
    
    // Solve
    for(m=0;m<tstp-k;m++){
      if(mask_ret[m]){
        element_set(size1, Q, G.retptr(tstp, m));
        for(j=1;j<=k+1;j++) element_incr(size1, Q, bd[j], G.retptr(tstp-j,m));

        element_iden(size1, M, -bd[0]+mu);
        element_incr(size1, M, -1, hmf.ptr(tstp));
        element_incr(size1, M, -dt*I.omega(0), Sig.retptr(tstp, tstp));

        element_linsolve_left(size1, size1, size1, M, G.retptr(tstp, m), Q);
      }
    }
    delete[] M;
    delete[] Q;
    delete[] bd;
  }
  return;
}


void dyson_step_omp_tv(int threads, int tstp, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  assert(G.size1()==Sig.size1());
  assert(G.size1()==hmf.size1());
  assert(G.nt()==Sig.nt());
  assert(G.nt()==hmf.nt());
  assert(G.nt()>=I.k());
  assert(G.sig()==Sig.sig());
  assert(tstp>I.k());
  assert(G.ntau()==Sig.ntau());
  
  int size1 = G.size1(), es = size1*size1, k = I.k(), ntau = G.ntau();
  cplx cplxi = cplx(0,1);
  for(int j=0;j<=ntau;j++) element_set_zero(size1, G.tvptr(tstp,j));

#pragma omp parallel num_threads(threads)
  {
    int j,i,l;
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    std::vector<bool> mask_tv(ntau+1, false);
    cplx *bd = new cplx[k+2];
    cplx *Q = new cplx[es];
    cplx *M = new cplx[es];
    for(j=0;j<=k+1;j++) bd[j] = -I.bd_weights(j)*cplxi/dt;
    for(i=0;i<=ntau;i++)
      if(i%nthreads==tid)
        mask_tv[i] = true;
    incr_convolution_tv(tstp, mask_tv, G, Sig, Sig, G, G, I, beta, dt);
    
    // Solve
    element_iden(size1, M, -bd[0]+mu);
    element_incr(size1, M, -1, hmf.ptr(tstp));
    element_incr(size1, M, -dt*I.omega(0), Sig.retptr(tstp,tstp));
    for(j=0;j<=ntau;j++){
      if(mask_tv[j]){
        element_set(size1, Q, G.tvptr(tstp,j));
        for(l=1;l<=k+1;l++){
          element_incr(size1, Q, bd[l], G.tvptr(tstp-l,j));
        }
        element_linsolve_left(size1,size1,size1, M,G.tvptr(tstp, j),  Q);
      }
    }

    delete[] M;
    delete[] Q;
    delete[] bd;
  }
  return;
}


void dyson_step_omp_les(int threads, int tstp, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  assert(G.size1()==Sig.size1());
  assert(G.size1()==hmf.size1());
  assert(G.nt()==Sig.nt());
  assert(G.nt()==hmf.nt());
  assert(G.nt()>=I.k());
  assert(G.sig()==Sig.sig());
  assert(tstp>I.k());
  
  cplx cplxi = cplx(0,1);
  int k = I.k(), size1 = G.size1(), es = size1*size1, ntop = (tstp>k?tstp:k);

  // Parallelization starts only for tstp>=2*k+1
  if(tstp < 2*k+1){
    dyson_step_les(tstp, I, G, Sig, hmf.ptr(0), mu, beta, dt);
    return;
  }
  for(int j=0;j<=tstp;j++) element_set_zero(size1, G.lesptr(j,tstp));

  // This solves for G^L(j,tstp) j=0...n-k-1.  Rest of the timestep done serially
#pragma omp parallel num_threads(threads)
  {
    int i,j,l;
    cplx cplxi = cplx(0,1);
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    std::vector<bool> mask_les(tstp+1, false);
    for(i=0;i<tstp-k;i++){
      if(i%nthreads==tid) mask_les[i]=true;
    }
    cplx *bd = new cplx[k+2];
    cplx *Q = new cplx[es];
    cplx *M = new cplx[es];
    for(i=0;i<=k+1;i++) bd[i] = I.bd_weights(i) * cplxi /dt;
    incr_convolution_les(tstp, mask_les, G, G, G, Sig, Sig, I, beta, dt);

    // Solve
    for(j=0;j<tstp-k;j++){
      if(mask_les[j]){
        element_iden(size1, M, -bd[0]+mu);
        element_conj(size1, Q, Sig.retptr(j,j));
        element_incr(size1, M, -dt*I.omega(0), Q);
        element_incr(size1, M, -1, hmf.ptr(tstp));
        
        element_set(size1, Q, G.lesptr(j,tstp));
        for(l=1;l<=k+1;l++) element_incr(size1, Q, bd[l], G.lesptr(j,tstp-l));

        element_linsolve_right(size1, size1, size1, G.lesptr(j, tstp), M, Q);
      }
    }

    delete[] Q;
    delete[] M;
    delete[] bd;
  }

  //Parallelize the remaining integrals
  cplx *gles = new cplx[(k+1)*es];
  cplx *Q = new cplx[es];
  cplx *M = new cplx[es];
#pragma omp parallel num_threads(threads)
  {
    int c;
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    for(c=tstp-k;c<=tstp;c++){
      if(c%nthreads==tid){
        element_set_zero(size1,gles+(c-tstp+k)*es);
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

  delete[] gles;
  delete[] M;
  delete[] Q;
  return;
}



void dyson_step_omp(int omp_threads, int n, const INTEG &I, GREEN &G, const GREEN &Sig, const function &hmf, double mu, double beta, double dt){
  assert(G.size1()==Sig.size1());
  assert(G.size1()==hmf.size1());
  assert(G.nt()==Sig.nt());
  assert(G.nt()==hmf.nt());
  assert(G.nt()>=I.k());
  assert(G.sig()==Sig.sig());
  assert(n>I.k());
  assert(G.ntau()==Sig.ntau());

  int threads = (omp_threads == -1 ? omp_get_max_threads() : omp_threads);
  dyson_step_omp_ret(threads, n, I, G, Sig, hmf, mu, dt);
  dyson_step_omp_tv(threads, n, I, G, Sig, hmf, mu, beta, dt);
  dyson_step_omp_les(threads, n, I, G, Sig, hmf, mu, beta, dt);
}
  #endif // USE_OMP

} // namespace

  #endif // DYSON_OMP_IMPL
