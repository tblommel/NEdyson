#ifndef MAT_UTILS_DECL
#define MAT_UTILS_DECL

#include <complex>
#include "greens.h"
#include "elementops.h"
#include "mat_utils.h"
#include "utils.h"

namespace NEdyson{

void get_dftcorr_cubic(double th, double *corfac, cplx *endcor);
void matsubara_ft(cplx *res, int m, const green_func &Sig, const cplx *sigdft, int sig, double beta);
double get_tau(int r, double beta, int ntau);
double get_omega(int m, double beta, int sig);
void matsubara_dft(cplx *mdft,const green_func &G, int sig);
void set_first_order_tail(cplx *gmat, cplx *one, double beta, int sg, int ntau, int sig, int size1);
void force_mat_herm(green_func &G);













}//namespace
#endif
