#ifndef ELEMENTS_IMPL
#define ELEMENTS_IMPL

#include "elementops.hpp"
#define map Eigen::Map<Eigen::Matrix<cplx,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >

namespace NEdyson{


void element_set_zero(int size1,cplx *z){
	std::memset(z,0,sizeof(cplx)*size1*size1);
}


void element_set(int size1, cplx *z, cplx *z1){
	std::memcpy(z,z1,sizeof(cplx)*size1*size1);
}


void element_iden(int size1, cplx *z){
	element_set_zero(size1, z);
	for(int i=0;i<size1;i++) z[i*size1+i]=1.;
}


void element_iden(int size1, cplx &a, cplx *z){
	element_set_zero(size1,z);
	for(int i=0;i<size1;i++) z[i*size1+i]=a;
}


void element_conj(int size1, cplx *z, cplx *z1){
	int l,m;
	cplx tmp;
	for(int l=0;l<size1;l++){
		for(int m=0;m<size1;m++){
			tmp=z1[m*size1+l];
			z[l*size1+m]=cplx(tmp.real(),-tmp.imag());
		}
	}
}


void element_conj(int size1, cplx *z){
	cplx *z1 = new cplx[size1*size1];
	element_set(size1,z1,z);
	element_conj(size1,z,z1);
	delete[] z1;
}


void element_smul(int size1, cplx *z, double z0){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]*=z0;
}

void element_smul(int size1, double *z, double z0){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]*=z0;
}


void element_smul(int size1, cplx *z, cplx z0){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]*=z0;
}



void element_smul(int size1, cplx *z, double z0, cplx *z1){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]=z0*z1[i];
}



void element_smul(int size1, cplx *z, cplx z0, cplx *z1){
	int es=size1*size1,i;
	for(i=0;i<es;i++) z[i]=z0*z1[i];
}



void element_incr(int size1, cplx *z, cplx *z1){
	int es=size1*size1,i;
	for(int i=0;i<es;i++) z[i]+=z1[i];
}



void element_incr(int size1, cplx *z, cplx z0, cplx *z1){
	int es=size1*size1,i;
	for(int i=0;i<es;i++) z[i]+=z0*z1[i];
}

void element_imag_incr(int size1, double *z0, cplx weight, cplx *z1){
	int es = size1*size1,i;
	for(int i=0;i<es;i++) z0[i] += (weight*z1[i]).imag();	
}



void element_incr(int size1, cplx *z, cplx *z0, cplx *z1){
	switch(size1){
		case 2:
			if(z!=z0&&z!=z1)
			Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z,2,2).noalias()+=Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z0,2,2)*Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z1,2,2);
			else
			Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z,2,2)+=Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z0,2,2)*Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z1,2,2);
			break;
		case 3:
			if(z!=z0&&z!=z1)
			Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z,3,3).noalias()+=Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z0,3,3)*Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z1,3,3);
			else
			Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z,3,3)+=Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z0,3,3)*Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z1,3,3);
			break;
		case 4:
			if(z!=z0&&z!=z1)
			Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z,4,4).noalias()+=Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z0,4,4)*Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z1,4,4);
			else
			Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z,4,4)+=Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z0,4,4)*Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z1,4,4);
			break;
		default:
			if(z!=z0&&z!=z1)
			Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z,size1,size1).noalias()+=Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z0,size1,size1)*Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z1,size1,size1);
			else
			Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z,size1,size1)+=Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z0,size1,size1)*Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z1,size1,size1);
			break;
	}
}



void element_incr(int size1, cplx *z, cplx z0, cplx *z1, cplx *z2){
	switch(size1){
		case 2:
			if(z!=z1&&z!=z2)
			Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z,2,2).noalias()+=z0*Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z1,2,2)*Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z2,2,2);
			else
			Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z,2,2)+=z0*Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z1,2,2)*Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z2,2,2);
			break;
		case 3:
			if(z!=z1&&z!=z2)
			Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z,3,3).noalias()+=z0*Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z1,3,3)*Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z2,3,3);
			else
			Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z,3,3)+=z0*Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z1,3,3)*Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z2,3,3);
			break;
		case 4:
			if(z!=z1&&z!=z2)
			Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z,4,4).noalias()+=z0*Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z1,4,4)*Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z2,4,4);
			else
			Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z,4,4)+=z0*Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z1,4,4)*Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z2,4,4);
			break;
		default:
			if(z!=z1&&z!=z2)
			Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z,size1,size1).noalias()+=z0*Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z1,size1,size1)*Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z2,size1,size1);
			else
			Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z,size1,size1)+=z0*Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z1,size1,size1)*Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z2,size1,size1);
			break;
	}
}



void element_mult_small(cplx *z, cplx *z1, cplx *z2){
	*z=*z1*(*z2);
}


void element_mult_large(int size1, cplx *z, cplx *z1, cplx *z2){
	switch(size1){
		case 2:
			Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z,2,2)=Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z1,2,2)*Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >(z2,2,2);
			break;
		case 3:
			Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z,3,3)=Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z1,3,3)*Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >(z2,3,3);
			break;
		case 4:
			Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z,4,4)=Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z1,4,4)*Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >(z2,4,4);
			break;
		default:
			Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z,size1,size1)=Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z1,size1,size1)*Eigen::Map<Eigen::Matrix<cplx,-1,-1,Eigen::RowMajor> >(z2,size1,size1);
			break;
	}
}


void element_mult(int size1, cplx *z, cplx *z1, cplx *z2){
	if(size1==1) element_mult_small(z, z1, z2);
	else element_mult_large(size1, z, z1, z2);
}

//A=B*C. A is axc.  B is axb. C is bxc

void element_mult(int a, int b, int c, cplx *A, cplx *B, cplx *C){
	map(A,a,c)=map(B,a,b)*map(C,b,c);
}



void element_inverse_large(int size1, cplx *z, cplx *z1){
	int i,es=size1*size1;
	cdmatrix A(size1,size1);
	cdmatrix X(size1,size1);
	A=map(z1,size1,size1);
	Eigen::FullPivLU<cdmatrix> lu(A);
	X=lu.inverse();
	map(z,size1,size1)=X;
}


void element_inverse_small(cplx *z, cplx *z1){
	*z=1./(*z1);
}


void element_inverse(int size1, cplx *z, cplx *z1){
	if(size1==1) element_inverse_small(z,z1);
	else element_inverse_large(size1, z, z1);
}

//solve MX=Q for X. M is axb. X is bxc. Q is axc.

void element_linsolve_left_large(int a, int b, int c, cplx *M, cplx *X, cplx *Q){
	Eigen::FullPivLU<cdmatrix> lu(map(M,a,b));
	map(X,b,c)=lu.solve(map(Q,a,c));
}


void element_linsolve_left_small(cplx *M, cplx *X, cplx *Q){
	*X=(*Q)/(*M);
}


void element_linsolve_left(int a, int b, int c, cplx *M, cplx *X, cplx *Q){
	if(a==1&&b==1&&c==1) element_linsolve_left_small(M,X,Q);
	else element_linsolve_left_large(a, b, c, M, X, Q);
}

//Solve XM=Q for X. M is bxc. X is axb. Q is axc.

void element_linsolve_right_large(int a, int b, int c, cplx *X, cplx *M, cplx*Q){
	cdmatrix MT(c,b);
	MT=map(M,b,c).transpose();
	cdmatrix QT(c,a);
	QT=map(Q,a,c).transpose();
	cdmatrix XT(b,a);
	Eigen::FullPivLU<cdmatrix> lu(MT);
	XT=lu.solve(QT);
//	XT=MT.colPivHouseholderQr().solve(QT);
	cdmatrix XTT(a,b);	
	XTT=XT.transpose();
	map(X,a,b)=XTT;
}


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


}//nameespace
#endif
