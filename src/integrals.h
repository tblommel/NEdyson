#ifndef INTEGRALS_DECL
#define INTEGRALS_DECL

#include "elementops.h"
#include "greens.h"
#include "integration.h"
#include "utils.h"
#include <complex>

namespace NEdyson{

void CTV1(const INTEG &I, const GREEN &A, const GREEN &Acc, GREEN &B, int n, double dt);

void CTV2(const INTEG &I, const GREEN &A, const GREEN &B, int n, int m, double beta, cplx *res);

void CTV3(const INTEG &I, const GREEN &A, const GREEN &B, int n, int m, double beta, cplx *res);

void Cles2_tstp(const INTEG &I, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double dt, cplx *Q);

void Cles3_tstp(const INTEG &I,const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double beta, cplx *Q);

void Cles2_tstp(int j1, int j2, const INTEG &I, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double dt, cplx *Q);

void Cles3_tstp(int j1, int j2, const INTEG &I,const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, int m, double beta, cplx *Q);

void incr_convolution_ret(int tstp, const std::vector<bool> &mask_ret, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I , double dt);

void incr_convolution_tv(int tstp, const std::vector<bool> &mask_tv, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I , double dt);

void incr_convolution_les(int tstp, const std::vector<bool> &mask_les, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I , double dt);

}//namespace

#endif
