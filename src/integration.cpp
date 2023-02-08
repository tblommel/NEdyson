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
  long double *PL = new long double [(k+1)*(k+1)];
  make_poly_interp(k, PL);
  for(int i = 0; i < (k+1)*(k+1); i++) P[i] = PL[i];
  delete [] PL;
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


void make_poly_diff(int k, double *D){
  long double *PL = new long double [(k+1)*(k+1)];
  long double *DL = new long double [(k+1)*(k+1)];
  make_poly_interp(k, PL);
  make_poly_diff(k, PL, DL);
  for(int i = 0; i < (k+1)*(k+1); i++) D[i] = DL[i];
  delete [] PL;
  delete [] DL;
}

void make_poly_integ(int k, long double *P, long double *I){
  long double npow, mpow,x;
  int k1=k+1;
  for(int m=0;m<=k;m++){
    for(int n=0;n<=k;n++){
      for(int l=0;l<=k;l++){
        mpow=m;
        npow=n;
        x=0.;
        for(int a=0;a<=k;a++){
          x += P[a*k1+l]*(npow-mpow)/(a+1);
          mpow*=m;
          npow*=n;
        }
        I[m*k1*k1+n*k1+l]=x;
      }
    }
  }
}

void make_poly_integ(int k, double *I){
  long double *PL = new long double [(k+1)*(k+1)];
  long double *IL = new long double [(k+1)*(k+1)*(k+1)];
  make_poly_interp(k, PL);
  make_poly_integ(k, PL, IL);
  for(int i = 0; i < (k+1)*(k+1)*(k+1); i++) I[i] = IL[i];
  delete [] PL;
  delete [] IL;
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

void make_start(int k, double *s){
  int k1=k+1,k2=k+2;
  long double x, npow;
  long double *PP = new long double [k2*k2];
  make_poly_interp(k,PP);
  for(int n=0;n<=k;n++){
    for(int l=0;l<=k;l++){
      x=0.;
      npow=n;
      for(int a=0;a<=k;a++){
        x+=PP[a*k1+l]*npow/(a+1.);
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

void make_ex_weights(int k, double *P, double *E) {
  memset(E, 0, (k+1)*sizeof(double));

  for(int l=0; l<=k; l++) {
    for(int j=0; j<=k; j++) {
      E[l] += P[j*(k+1)+l] * (1 - 2*(j%2));
    }
  }
}

}//namespace
#endif
