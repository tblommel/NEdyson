#ifndef ELEMENTS_DECL
#define ELEMENTS_DECL

#include <Eigen/Eigen>
#include <string.h>
#include <iostream>

#include "utils.h"

namespace NEdyson{
void element_set_zero(int size1,cplx *z);
void element_set(int size1, cplx *z, cplx *z1);
void element_set(int size1, cplx *z, const cplx *z1);
void element_iden(int size1, cplx *z);
void element_iden(int size1, cplx *z, cplx a);
void element_conj(int size1, cplx *z, const cplx *z1);
void element_conj(int size1, cplx *z);
void element_smul(int size1, cplx *z, double z0);
void element_smul(int size1, double *z, double z0);
void element_smul(int size1, cplx *z, cplx z0);
void element_smul(int size1, cplx *z, double z0, const cplx *z1);
void element_smul(int size1, cplx *z, cplx z0, const cplx *z1);
void element_incr(int size1, cplx *z, cplx *z1);
void element_incr(int size1, cplx *z, cplx z0, cplx *z1);
void element_incr(int size1, cplx *z, cplx *z0, cplx *z1);
double element_diff(int size1, const cplx *A, const cplx *B);
void element_incr(int size1, cplx *z, cplx z0, cplx *z1, cplx *z2);
void element_mult_small(cplx *z, cplx *z1, cplx *z2);
void element_mult_large(int size1, cplx *z,  cplx *z1,  cplx *z2);
void element_mult(int size1, cplx *z,  cplx *z1,  cplx *z2);
void element_mult(int a, int b, int c, cplx *A,  cplx *B,  cplx *C);
void element_inverse_large(int size1, cplx *z,  cplx *z1);
void element_inverse_small(cplx *z,  cplx *z1);
void element_inverse(int size1, cplx *z,  cplx *z1);
void element_linsolve_left_large(int a, int b, int c,  cplx *M, cplx *X,  cplx *Q);
void element_linsolve_left_small( cplx *M, cplx *X,  cplx *Q);
void element_linsolve_left(int a, int b, int c,  cplx *M, cplx *X,  cplx *Q);
void element_linsolve_right_large(int a, int b, int c, cplx *X,  cplx *M,  cplx*Q);
void element_linsolve_right(int a, int b, int c, cplx *X,  cplx *M,  cplx *Q);
double element_norm2(int size1, cplx *A);
void element_imag_incr(int size1, double *z0, cplx weight, const cplx *z1);
}

#endif
