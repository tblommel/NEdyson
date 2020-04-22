#ifndef MAT_UTILS_DECL
#define MAT_UTILS_DECL

namespace NEdyson{

void get_dftcorr_cubic(double th, double *corfac, cplx *endcor);
void matsubara_ft(cplx *res, int m, green_func &Sig, cplx *sigdft, int sig, double beta);
double get_tau(int r, double beta, int ntau);
double get_omega(int m, double beta, int sig);
void matsubara_dft(cplx *mdft, green_func &G, int sig);
void set_first_order_tail(cplx *gmat, cplx *one, double beta, int sg, int ntau, int sig, int size1);
void force_mat_herm(green_func &G);













}//namespace
#endif
