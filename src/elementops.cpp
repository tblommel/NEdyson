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
			cdmmap2(z,2,2).noalias()+=cdmmap2(z0,2,2)*cdmmap2(z1,2,2);
			else
			cdmmap2(z,2,2)+=cdmmap2(z0,2,2)*cdmmap2(z1,2,2);
			break;
		case 3:
			if(z!=z0&&z!=z1)
			cdmmap3(z,3,3).noalias()+=cdmmap3(z0,3,3)*cdmmap3(z1,3,3);
			else
			cdmmap3(z,3,3)+=cdmmap3(z0,3,3)*cdmmap3(z1,3,3);
			break;
		case 4:
			if(z!=z0&&z!=z1)
			cdmmap4(z,4,4).noalias()+=cdmmap4(z0,4,4)*cdmmap4(z1,4,4);
			else
			cdmmap4(z,4,4)+=cdmmap4(z0,4,4)*cdmmap4(z1,4,4);
			break;
		default:
			if(z!=z0&&z!=z1)
			cdmmap(z,size1,size1).noalias()+=cdmmap(z0,size1,size1)*cdmmap(z1,size1,size1);
			else
			cdmmap(z,size1,size1)+=cdmmap(z0,size1,size1)*cdmmap(z1,size1,size1);
			break;
	}
}

// Z+=z0*Z1*Z2
void element_incr(int size1, cplx *z, cplx z0, cplx *z1, cplx *z2){
	switch(size1){
		case 2:
			if(z!=z1&&z!=z2)
			cdmmap2(z,2,2).noalias()+=z0*cdmmap2(z1,2,2)*cdmmap2(z2,2,2);
			else
			cdmmap2(z,2,2)+=z0*cdmmap2(z1,2,2)*cdmmap2(z2,2,2);
			break;
		case 3:
			if(z!=z1&&z!=z2)
			cdmmap3(z,3,3).noalias()+=z0*cdmmap3(z1,3,3)*cdmmap3(z2,3,3);
			else
			cdmmap3(z,3,3)+=z0*cdmmap3(z1,3,3)*cdmmap3(z2,3,3);
			break;
		case 4:
			if(z!=z1&&z!=z2)
			cdmmap4(z,4,4).noalias()+=z0*cdmmap4(z1,4,4)*cdmmap4(z2,4,4);
			else
			cdmmap4(z,4,4)+=z0*cdmmap4(z1,4,4)*cdmmap4(z2,4,4);
			break;
		default:
			if(z!=z1&&z!=z2)
			cdmmap(z,size1,size1).noalias()+=z0*cdmmap(z1,size1,size1)*cdmmap(z2,size1,size1);
			else
			cdmmap(z,size1,size1)+=z0*cdmmap(z1,size1,size1)*cdmmap(z2,size1,size1);
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
			cdmmap2(z,2,2)=cdmmap2(z1,2,2)*cdmmap2(z2,2,2);
			break;
		case 3:
			cdmmap3(z,3,3)=cdmmap3(z1,3,3)*cdmmap3(z2,3,3);
			break;
		case 4:
			cdmmap4(z,4,4)=cdmmap4(z1,4,4)*cdmmap4(z2,4,4);
			break;
		default:
			cdmmap(z,size1,size1)=cdmmap(z1,size1,size1)*cdmmap(z2,size1,size1);
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
	cdmmap(A,a,c)=cdmmap(B,a,b)*cdmmap(C,b,c);
}


// Z = Z1^{-1}
void element_inverse_large(int size1, cplx *z, cplx *z1){
	int i,es=size1*size1;
	cdmatrix A(cdmmap(z1,size1,size1));
	cdmatrix X(size1,size1);
	Eigen::FullPivLU<cdmatrix> lu(A);
	X=lu.inverse();
	cdmmap(z,size1,size1)=X;
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
	Eigen::FullPivLU<cdmatrix> lu(cdmmap(M,a,b));
	cdmmap(X,b,c)=lu.solve(cdmmap(Q,a,c));
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
	cdmatrix MT(c,b);
	MT=cdmmap(M,b,c).transpose();
	cdmatrix QT(c,a);
	QT=cdmmap(Q,a,c).transpose();
	cdmatrix XT(b,a);
	Eigen::FullPivLU<cdmatrix> lu(MT);
	XT=lu.solve(QT);
	cdmatrix XTT(a,b);	
	XTT=XT.transpose();
	cdmmap(X,a,b)=XTT;
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
