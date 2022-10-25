#ifndef DYSON_IMPL
#define DYSON_IMPL

#include "dyson.h"

#include "dyson_free.h"

#include "dyson_integrals.h"
#include "dyson_tti_integrals.h"

#include "dyson_step.h"
#include "dyson_tti_step.h"
#include "dyson_step_hf.h"
#include "dyson_tti_step_hf.h"

#include "dyson_boot.h"
#include "dyson_tti_boot.h"
#include "dyson_boot_hf.h"
#include "dyson_tti_boot_hf.h"

#include "dyson_dfield.h"

namespace NEdyson{

dyson::dyson(int nt, int ntau, int nao, int k, bool hfbool) : nt_(nt),
                                                 ntau_(ntau),
                                                 nao_(nao),
                                                 es_(nao*nao),
                                                 k_(k),
                                                 hfbool_(hfbool),
                                                 ex_weights(k+1),
                                                 tmp(es_),
                                                 tmp2(es_),
                                                 iden(es_),
                                                 M(k_*k_*es_),
                                                 Q((nt_+1)*es_),
                                                 X((nt_+1)*es_),
                                                 NTauTmp((k+1)*(ntau_+1)*es_),
                                                 I(k_),
                                                 Conv(nao, ntau+1)
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
//  for(j = 0; j <= k_; j++) {
//    ZMatrixMap(G.lesptr(0,n), nao_, nao_).noalias() += ex_weights(j) * ZMatrixMap(G.lesptr(0,n-1-j), nao_, nao_);
//  }

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


} // Namespace

#endif
