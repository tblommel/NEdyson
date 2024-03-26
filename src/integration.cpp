#ifndef INTEGRATION_IMPL
#define INTEGRATION_IMPL

#include "integration.h"

namespace Integration{

void make_poly_interp(int k, long double *P){
	Eigen::Matrix<long double, -1, -1> M(k+1,k+1);
	for(int i=0;i<=k;i++){
		M(i,0)=1.;
		for(int n=1;n<=k;n++){
			M(i,n)=M(i,n-1)*i;
		}
	}
	Eigen::FullPivLU<Eigen::Matrix<long double, -1, -1> > lu(M);
	M=lu.inverse();
	for(int i=0;i<=k;i++){
		for(int j=0;j<=k;j++){
			P[i*(k+1)+j]=M(i,j);
		}
	}
}



void make_poly_interp(int k, double *P){
	Eigen::MatrixXd M(k+1,k+1);
	for(int i=0;i<=k;i++){
		M(i,0)=1.;
		for(int n=1;n<=k;n++){
			M(i,n)=M(i,n-1)*i;
		}
	}
	Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
	M=lu.inverse();
	for(int i=0;i<=k;i++){
		for(int j=0;j<=k;j++){
			P[i*(k+1)+j]=M(i,j);
		}
	}
}


void make_poly_diff(int k, long double *P, long double *D){
	long double npow,x;
	int k1=k+1;
	for(int n=0;n<=k;n++){
		for(int l=0;l<=k;l++){
			npow=1.;
			x=0.;
			for(int a=1;a<=k;a++){
				if(a>1) npow*=n;
				x+=(long double)P[a*k1+l]*a*npow;
			}
			D[n*k1+l]=x;
		}
	}
}


void make_poly_diff(int k, double *P, double *D){
	long double npow,x;
	int k1=k+1;
	for(int n=0;n<=k;n++){
		for(int l=0;l<=k;l++){
			npow=1.;
			x=0.;
			for(int a=1;a<=k;a++){
				if(a>1) npow*=n;
				x+=(long double)P[a*k1+l]*a*npow;
			}
			D[n*k1+l]=x;
		}
	}
}

void make_poly_integ(int k, double *P, double *I){
	long double x, npow, mpow;
  int k2 = k+2;
	long double *PP = new long double [k2*k2];
	make_poly_interp(k,PP);
	int k1=k+1;
	for(int m=0;m<=k;m++){
		for(int n=0;n<=k;n++){
			for(int l=0;l<=k;l++){
				mpow=m;
				npow=n;
				x=0.;
				for(int a=0;a<=k;a++){
					x+=(long double)PP[a*k1+l]*(npow-mpow)/(a+1.);
					mpow*=m;
					npow*=n;
				}
				I[m*k1*k1+n*k1+l]=x;
			}
		}
	}
}

void make_bd_weights(int k, double *BD){
	int k2=k+2;
	long double *P = new long double [k2*k2];
  long double *D = new long double [k2*k2];
	make_poly_interp(k+1, P);
	make_poly_diff(k+1,P,D);
	for(int i=0;i<k2;i++){
		BD[i]=-D[i];
	}
	delete [] P;
	delete [] D;
}

void make_start(int k, double *P, double *s){
	int k1=k+1,k2=k+2;
	long double x, npow;
	long double *PP = new long double [k2*k2];
	make_poly_interp(k,PP);
	for(int n=0;n<=k;n++){
		for(int l=0;l<=k;l++){
			x=0.;
			npow=n;
			for(int a=0;a<=k;a++){
				x+=(long double)PP[a*k1+l]*npow/(a+1.);
				npow*=n;
			}
			s[n*k1+l]=x;
		}
	}
  delete [] PP;
}

int fact(int n){
	if(n>1) return n*fact(n-1);
	else return 1;
}

int pascal(int row, int col){
	return fact(row)/(fact(row-col)*fact(col));
}

void make_Omega(int k, double *O){
	int k1=k+1;
	Eigen::Matrix<long double,-1,-1> A(k1,k1);
	A.setZero();
	Eigen::Matrix<long double,-1,-1> y(k1,1);
	Eigen::Matrix<long double,-1,-1> res(k1,1);
	for(int i=0;i<=k;i++){
		y(i,0)=(1-2*((i+1)%2))/(i+2.);
		for(int j=0;j<=i;j++){
			A(i,j)=(1-2*((i+j)%2))/(i-j+1.);
		}
	}
	Eigen::FullPivLU<Eigen::Matrix<long double,-1,-1> > lu(A);
	res=lu.solve(y);
	long double x;
	for(int n=0;n<=k;n++){
		Eigen::Matrix<long double,-1,-1> I=Eigen::Matrix<long double,-1,-1>::Ones(k+n+2,1);
		for(int i=0;i<=k;i++){
			for(int j=0;j<=i;j++){
				x=(1-2*((i+j)%2))*res(i,0)*pascal(i,j);
				I(j,0)+=x;
				I(k+n-j+1,0)+=x;
			}
		}
		for(int i=0;i<=k;i++){
			O[n*k1+i]=I(i,0);
		}
	}
}

void make_rcorr(int k, double *P, double *R){
	int k1=k+1;
	double x;
	int mpowa,mpowb;
	for(int m=0;m<=k;m++){
		for(int r=0;r<=k;r++){
			for(int s=0;s<=k;s++){
				x=0.;
				mpowa=m;
				mpowb=1;
				for(int a=0;a<=k;a++){
					for(int b=0;b<=k;b++){
						x+=P[a*k1+r]*P[b*k1+s]*fact(a)*fact(b)*mpowa*mpowb/(double)fact(a+b+1);
						mpowb*=m;
					}
					mpowa*=m;
					mpowb=1;
				}
				R[m*k1*k1+r*k1+s]=x;
			}
		}
	}
}

}//namespace
#endif
