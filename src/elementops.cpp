#ifndef ELEMENTS_IMPL
#define ELEMENTS_IMPL

#include "elementops.h"

namespace NEdyson{

// Sets Z to be zero
void element_set_zero(int size1,cplx *z){
	std::memset(z,0,sizeof(cplx)*size1*size1);
}

// Copies Z1 to Z
void element_set(int size1, cplx *z, cplx *z1){
	std::memcpy(z,z1,sizeof(cplx)*size1*size1);
}

// Copies Z1 to Z
void element_set(int size1, cplx *z, const cplx *z1){
	std::memcpy(z,z1,sizeof(cplx)*size1*size1);
}

// Sets Z to be the identity matrix
void element_iden(int size1, cplx *z){
	element_set_zero(size1, z);
	for(int i=0;i<size1;i++) z[i*size1+i]=1.;
}

// Sets Z to be a*I
void element_iden(int size1, cplx *z, cplx a){
	element_set_zero(size1,z);
	for(int i=0;i<size1;i++) z[i*size1+i]=a;
}

// puts Hermitian Conj of Z1 into Z
void element_conj(int size1, cplx *z, const cplx *z1){
	int l,m;
	cplx tmp;
	for(int l=0;l<size1;l++){
		for(int m=0;m<size1;m++){
			tmp=z1[m*size1+l];
			z[l*size1+m]=cplx(tmp.real(),-tmp.imag());
		}
	}
}

// Conjugates Z
void element_conj(int size1, cplx *z){
	cplx *z1 = new cplx[size1*size1];
	element_set(size1,z1,z);
	element_conj(size1,z,z1);
	delete[] z1;
}

// Multiplies Z by z0
void element_smul(int size1, cplx *z, double z0){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]*=z0;
}

// Multiplies Z by z0
void element_smul(int size1, double *z, double z0){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]*=z0;
}

// Multiplies Z by z0
void element_smul(int size1, cplx *z, cplx z0){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]*=z0;
}

// Z=z0*Z1
void element_smul(int size1, cplx *z, double z0, const cplx *z1){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]=z0*z1[i];
}

// Z=z0*Z1
void element_smul(int size1, cplx *z, cplx z0, const cplx *z1){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]=z0*z1[i];
}

// Z+=Z1
void element_incr(int size1, cplx *z, cplx *z1){
	int es=size1*size1,i;
	for(int i=0;i<es;i++) z[i]+=z1[i];
}

// Z+=z0*Z1
void element_incr(int size1, cplx *z, cplx z0, cplx *z1){
	int es=size1*size1,i;
	for(int i=0;i<es;i++) z[i]+=z0*z1[i];
}

// Z0+=Imag(weight*Z1)
void element_imag_incr(int size1, double *z0, cplx weight, const cplx *z1){
	int es = size1*size1,i;
	for(int i=0;i<es;i++) z0[i] += (weight*z1[i]).imag();	
}

// z += Z0*Z1
void element_incr(int size1, cplx *z, cplx *z0, cplx *z1){
	switch(size1){
		case 2:
			if(z!=z0 && z!=z1)
			ZMatrixMap2(z,2,2).noalias()+=ZMatrixMap2(z0,2,2)*ZMatrixMap2(z1,2,2);
			else
			ZMatrixMap2(z,2,2)+=ZMatrixMap2(z0,2,2)*ZMatrixMap2(z1,2,2);
			break;
		case 3:
			if(z!=z0&&z!=z1)
			ZMatrixMap3(z,3,3).noalias()+=ZMatrixMap3(z0,3,3)*ZMatrixMap3(z1,3,3);
			else
			ZMatrixMap3(z,3,3)+=ZMatrixMap3(z0,3,3)*ZMatrixMap3(z1,3,3);
			break;
		case 4:
			if(z!=z0&&z!=z1)
			ZMatrixMap4(z,4,4).noalias()+=ZMatrixMap4(z0,4,4)*ZMatrixMap4(z1,4,4);
			else
			ZMatrixMap4(z,4,4)+=ZMatrixMap4(z0,4,4)*ZMatrixMap4(z1,4,4);
			break;
		default:
			if(z!=z0&&z!=z1)
			ZMatrixMap(z,size1,size1).noalias()+=ZMatrixMap(z0,size1,size1)*ZMatrixMap(z1,size1,size1);
			else
			ZMatrixMap(z,size1,size1)+=ZMatrixMap(z0,size1,size1)*ZMatrixMap(z1,size1,size1);
			break;
	}
}

// Z+=z0*Z1*Z2
void element_incr(int size1, cplx *z, cplx z0, cplx *z1, cplx *z2){
	switch(size1){
		case 2:
			if(z!=z1&&z!=z2)
			ZMatrixMap2(z,2,2).noalias()+=z0*ZMatrixMap2(z1,2,2)*ZMatrixMap2(z2,2,2);
			else
			ZMatrixMap2(z,2,2)+=z0*ZMatrixMap2(z1,2,2)*ZMatrixMap2(z2,2,2);
			break;
		case 3:
			if(z!=z1&&z!=z2)
			ZMatrixMap3(z,3,3).noalias()+=z0*ZMatrixMap3(z1,3,3)*ZMatrixMap3(z2,3,3);
			else
			ZMatrixMap3(z,3,3)+=z0*ZMatrixMap3(z1,3,3)*ZMatrixMap3(z2,3,3);
			break;
		case 4:
			if(z!=z1&&z!=z2)
			ZMatrixMap4(z,4,4).noalias()+=z0*ZMatrixMap4(z1,4,4)*ZMatrixMap4(z2,4,4);
			else
			ZMatrixMap4(z,4,4)+=z0*ZMatrixMap4(z1,4,4)*ZMatrixMap4(z2,4,4);
			break;
		default:
			if(z!=z1&&z!=z2)
			ZMatrixMap(z,size1,size1).noalias()+=z0*ZMatrixMap(z1,size1,size1)*ZMatrixMap(z2,size1,size1);
			else
			ZMatrixMap(z,size1,size1)+=z0*ZMatrixMap(z1,size1,size1)*ZMatrixMap(z2,size1,size1);
			break;
	}
}

// z+=z1*z2
void element_mult_small(cplx *z, cplx *z1, cplx *z2){
	*z=*z1*(*z2);
}

// Z=Z1*Z2
void element_mult_large(int size1, cplx *z, cplx *z1, cplx *z2){
	switch(size1){
		case 2:
			if(z!=z1&&z!=z2)
			ZMatrixMap2(z,2,2).noalias()=ZMatrixMap2(z1,2,2)*ZMatrixMap2(z2,2,2);
			else
			ZMatrixMap2(z,2,2)=ZMatrixMap2(z1,2,2)*ZMatrixMap2(z2,2,2);
			break;
		case 3:
			if(z!=z1&&z!=z2)
			ZMatrixMap3(z,3,3).noalias()=ZMatrixMap3(z1,3,3)*ZMatrixMap3(z2,3,3);
			else
			ZMatrixMap3(z,3,3)=ZMatrixMap3(z1,3,3)*ZMatrixMap3(z2,3,3);
			break;
		case 4:
			if(z!=z1&&z!=z2)
			ZMatrixMap4(z,4,4).noalias()=ZMatrixMap4(z1,4,4)*ZMatrixMap4(z2,4,4);
			else
			ZMatrixMap4(z,4,4)=ZMatrixMap4(z1,4,4)*ZMatrixMap4(z2,4,4);
			break;
		default:
			if(z!=z1&&z!=z2)
			ZMatrixMap(z,size1,size1).noalias()=ZMatrixMap(z1,size1,size1)*ZMatrixMap(z2,size1,size1);
			else
			ZMatrixMap(z,size1,size1)=ZMatrixMap(z1,size1,size1)*ZMatrixMap(z2,size1,size1);
			break;
	}
}

// Z=Z1*Z2
void element_mult(int size1, cplx *z,cplx *z1, cplx *z2){
	if(size1==1) element_mult_small(z, z1, z2);
	else element_mult_large(size1, z, z1, z2);
}

// B*C=A. A is axc.  B is axb. C is bxc
void element_mult(int a, int b, int c, cplx *A, cplx *B, cplx *C){
	ZMatrixMap(A,a,c)=ZMatrixMap(B,a,b)*ZMatrixMap(C,b,c);
}


// Z = Z1^{-1}
void element_inverse_large(int size1, cplx *z, cplx *z1){
	int i,es=size1*size1;
	ZMatrix A(ZMatrixMap(z1,size1,size1));
	ZMatrix X(size1,size1);
	Eigen::FullPivLU<ZMatrix> lu(A);
	X=lu.inverse();
	ZMatrixMap(z,size1,size1)=X;
}

// z = 1/z1
void element_inverse_small(cplx *z, cplx *z1){
	*z=1./(*z1);
}

// Z = Z1^{-1}
void element_inverse(int size1, cplx *z, cplx *z1){
	if(size1==1) element_inverse_small(z,z1);
	else element_inverse_large(size1, z, z1);
}

// solve MX=Q for X. M is axb. X is bxc. Q is axc.
void element_linsolve_left_large(int a, int b, int c, cplx *M, cplx *X, cplx *Q){
	Eigen::FullPivLU<ZMatrix> lu(ZMatrixMap(M,a,b));
	ZMatrixMap(X,b,c)=lu.solve(ZMatrixMap(Q,a,c));
}

// solve mx=q for x
void element_linsolve_left_small(cplx *M, cplx *X, cplx *Q){
	*X=(*Q)/(*M);
}

// solve MX=Q for X. M is axb. X is bxc. Q is axc.
void element_linsolve_left(int a, int b, int c, cplx *M, cplx *X, cplx *Q){
	if(a==1&&b==1&&c==1) element_linsolve_left_small(M,X,Q);
	else element_linsolve_left_large(a, b, c, M, X, Q);
}

//Solve XM=Q for X. M is bxc. X is axb. Q is axc.
void element_linsolve_right_large(int a, int b, int c, cplx *X, cplx *M, cplx*Q){
	ZMatrix MT(c,b);
	MT=ZMatrixMap(M,b,c).transpose();
	ZMatrix QT(c,a);
	QT=ZMatrixMap(Q,a,c).transpose();
	ZMatrix XT(b,a);
	Eigen::FullPivLU<ZMatrix> lu(MT);
	XT=lu.solve(QT);
	ZMatrix XTT(a,b);	
	XTT=XT.transpose();
	ZMatrixMap(X,a,b)=XTT;
}


//Solve XM=Q for X. M is bxc. X is axb. Q is axc.
void element_linsolve_right(int a, int b, int c, cplx *X, cplx *M, cplx *Q){
	if(a==1&&b==1&&c==1) element_linsolve_left_small(M,X,Q);
	else element_linsolve_right_large(a,b,c,X,M,Q);
}


double element_norm2(int size1, cplx *A){
	double ret = 0;
	int es=size1*size1;
	for(int i=0;i<es;i++) ret+=(A[i].real()*A[i].real()+A[i].imag()*A[i].imag());
	return sqrt(ret);
}

double element_diff(int size1, const cplx *A, const cplx *B){
  double ret = 0, imag, real;
  int es=size1*size1;
  for(int i=0;i<es;i++){
    real = A[i].real()-B[i].real();
    imag = A[i].imag()-B[i].imag();
    ret += real*real+imag*imag;
  }
  return sqrt(ret);
}


}//nameespace
#endif
