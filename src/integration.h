#ifndef INTEGRATION_DECL
#define INTEGRATION_DECL

#include <cmath>
#include <iostream>
#include <complex>
#include <stdlib.h>
#include <Eigen/Eigen>

#include "utils.h"

namespace Integration{
	
	void make_poly_interp(int k, double *P);
	void make_poly_diff(int k,double *P, double *D);
	void make_poly_integ(int k, double *P, double *I);
	void make_bd_weights(int k, double *BD);
	void make_start(int k, double *P, double *s);
	void make_Omega(int k, double *O);
	void make_rcorr(int k, double *P, double *R); 


	class Integrator{
		public:

		Integrator(int k){
			int k1=k+1;
			if(k<0){std::cout << "Integrator: k out of range"<<std::endl; abort();}
			poly_interp_= new double [k1*k1];
			poly_interp1= new double [2*2];
			poly_interp2= new double [3*3];
			poly_interp3= new double [4*4];
			poly_interp4= new double [5*5];
			poly_interp5= new double [6*6];
			poly_diff_= new double [k1*k1];
			poly_integ_= new double [k1*k1*k1];
			bd_weights_= new double [k1+1];
			gregory_start_= new double [k1*k1];
			gregory_start1= new double [2*2];
			gregory_start2= new double [3*3];
			gregory_start3= new double [4*4];
			gregory_start4= new double [5*5];
			gregory_start5= new double [6*6];
			gregory_Omega_= new double [k1*k1];
			gregory_Omega1= new double [2*2];
			gregory_Omega2= new double [3*3];
			gregory_Omega3= new double [4*4];
			gregory_Omega4= new double [5*5];
			gregory_Omega5= new double [6*6];
			gregory_omega_= gregory_Omega_+(k*k1);
			gregory_omega1= gregory_Omega1+(1*2);
			gregory_omega2= gregory_Omega2+(2*3);
			gregory_omega3= gregory_Omega3+(3*4);
			gregory_omega4= gregory_Omega4+(4*5);
			gregory_omega5= gregory_Omega5+(5*6);
			rcorr_= new double [k1*k1*k1];
			k_=k;
			Integration::make_poly_interp(k_,poly_interp_);
			Integration::make_poly_interp(1,poly_interp1);
			Integration::make_poly_interp(2,poly_interp2);
			Integration::make_poly_interp(3,poly_interp3);
			Integration::make_poly_interp(4,poly_interp4);
			Integration::make_poly_interp(5,poly_interp5);
			Integration::make_poly_diff(k_,poly_interp_,poly_diff_);
			Integration::make_poly_integ(k_,poly_interp_,poly_integ_);
			Integration::make_bd_weights(k_,bd_weights_);
			Integration::make_start(k_,poly_interp_,gregory_start_);
			Integration::make_start(1,poly_interp1,gregory_start1);
			Integration::make_start(2,poly_interp2,gregory_start2);
			Integration::make_start(3,poly_interp3,gregory_start3);
			Integration::make_start(4,poly_interp4,gregory_start4);
			Integration::make_start(5,poly_interp5,gregory_start5);
			Integration::make_Omega(k_,gregory_Omega_);
			Integration::make_Omega(1,gregory_Omega1);
			Integration::make_Omega(2,gregory_Omega2);
			Integration::make_Omega(3,gregory_Omega3);
			Integration::make_Omega(4,gregory_Omega4);
			Integration::make_Omega(5,gregory_Omega5);
			Integration::make_rcorr(k_,poly_interp_,rcorr_);
		}

		Integrator& operator=(const Integrator& I){
			int k1;
			if(this->k_==I.k_) return *this;
			delete [] poly_interp_;
			delete [] poly_diff_;
			delete [] poly_integ_;
			delete [] bd_weights_;
			delete [] gregory_start_;
			delete [] gregory_Omega_;
			delete [] rcorr_;
			k_=I.k_;
			k1=k_+1;
			poly_interp_= new double [k1*k1];
			poly_diff_= new double [k1*k1];
			poly_integ_= new double [k1*k1*k1];
			bd_weights_= new double [k1+1];
			gregory_start_= new double [k1*k1];
			gregory_Omega_= new double [k1*k1];
			gregory_omega_= gregory_Omega_+(k_*k1);
			rcorr_= new double [k1*k1*k1];
			Integration::make_poly_interp(k_,poly_interp_);
			Integration::make_poly_diff(k_,poly_interp_,poly_diff_);
			Integration::make_poly_integ(k_,poly_interp_,poly_integ_);
			Integration::make_bd_weights(k_,bd_weights_);
			Integration::make_start(k_,poly_interp_,gregory_start_);
			Integration::make_Omega(k_,gregory_Omega_);
			Integration::make_rcorr(k_,poly_interp_,rcorr_);
			return *this;
		}

		~Integrator(){
			delete [] poly_interp_;
			delete [] poly_diff_;
			delete [] poly_integ_;
			delete [] bd_weights_;
			delete [] gregory_start_;
			delete [] gregory_Omega_;
			delete [] rcorr_;
		}
		
		Integrator(const Integrator &I){
			int k1;
			k_=I.k_;
			k1=k_+1;
			poly_interp_= new double [k1*k1];
			poly_diff_= new double [k1*k1];
			poly_integ_= new double [k1*k1*k1];
			bd_weights_= new double [k1+1];
			gregory_start_= new double [k1*k1];
			gregory_Omega_= new double [k1*k1];
			gregory_omega_= gregory_Omega_+(k_*k1);
			rcorr_= new double [k1*k1*k1];
			Integration::make_poly_interp(k_,poly_interp_);
			Integration::make_poly_diff(k_,poly_interp_,poly_diff_);
			Integration::make_poly_integ(k_,poly_interp_,poly_integ_);
			Integration::make_bd_weights(k_,bd_weights_);
			Integration::make_start(k_,poly_interp_,gregory_start_);
			Integration::make_Omega(k_,gregory_Omega_);
			Integration::make_rcorr(k_,poly_interp_,rcorr_);
		}
		
		int k(void) const {return k_;}
		double poly_interp(int a, int l) const {
      assert(a<=k_);
      assert(l<=k_);
      return poly_interp_[a*(k_+1)+l];
    }
		double poly_diff(int m, int l) const {
      assert(m<=k_);
      assert(l<=k_);
      return poly_diff_[m*(k_+1)+l];
    }
		double poly_integ(int m, int n, int l) const {
      assert(m<=k_);
      assert(n<=k_);
      assert(l<=k_);
      return poly_integ_[m*(k_+1)*(k_+1)+n*(k_+1)+l];
    }
		double bd_weights(int l) const {
      assert(l<=k_+1);
      return bd_weights_[l];
    }
		double gregory_weights(int n, int j) const {
			if(n<=k_&&j<=k_){
				return gregory_start_[n*(k_+1)+j];
			}
			else if(n<=2*k_+1){
				if(j<=k_) return gregory_Omega_[(n-k_-1)*(k_+1)+j];
				else return gregory_omega_[n-j];
			}
			else{
				if(j<=k_) return gregory_omega_[j];
				else if(j<n-k_) return 1;
				else return gregory_omega_[n-j];
			}
		}
		double omega(int j) const {
      assert(j<=k_);
      return gregory_omega_[j];
    }
		double start(int i, int j) const {
      assert(i<=k_);
      assert(j<=k_);
      return gregory_start_[i*(k_+1)+j];
    }
		double Omega(int i, int j) const {
      assert(i<=k_);
      assert(j<=k_);
      return gregory_Omega_[i*(k_+1)+j];
    }
		double rcorr(int m, int j, int l) const {
      assert(m<=k_);
      assert(j<=k_);
      assert(l<=k_);
      return rcorr_[m*(k_+1)*(k_+1)+j*(k_+1)+l];
    }


		double gregory_weights(int k, int n, int j) const {
			if(n <= k && j <= k ) {
        if(k == 1) {
				  return gregory_start1[n*(k+1)+j];
        }
        if(k == 2) {
				  return gregory_start2[n*(k+1)+j];
        }
        if(k == 3) {
				  return gregory_start3[n*(k+1)+j];
        }
        if(k == 4) {
				  return gregory_start4[n*(k+1)+j];
        }
        if(k == 5) {
				  return gregory_start5[n*(k+1)+j];
        }
			}
			else if(n <= 2*k + 1){
				if(j<=k) {
          if( k == 1) {
            return gregory_Omega1[(n-k-1)*(k+1)+j];
          }
          if( k == 2) {
            return gregory_Omega2[(n-k-1)*(k+1)+j];
          }
          if( k == 3) {
            return gregory_Omega3[(n-k-1)*(k+1)+j];
          }
          if( k == 4) {
            return gregory_Omega4[(n-k-1)*(k+1)+j];
          }
          if( k == 5) {
            return gregory_Omega5[(n-k-1)*(k+1)+j];
          }
        }
				else {
          if( k == 1) {
            return gregory_omega1[n-j];
          }
          if( k == 2) {
            return gregory_omega2[n-j];
          }
          if( k == 3) {
            return gregory_omega3[n-j];
          }
          if( k == 4) {
            return gregory_omega4[n-j];
          }
          if( k == 5) {
            return gregory_omega5[n-j];
          }
        }
			}
			else{
				if(j<=k) {
          if(k==1){
            return gregory_omega1[j];
          }
          if(k==2){
            return gregory_omega2[j];
          }
          if(k==3){
            return gregory_omega3[j];
          }
          if(k==4){
            return gregory_omega4[j];
          }
          if(k==5){
            return gregory_omega5[j];
          }
        }
				else if(j<n-k) return 1;
				else {  
          if(k==1){
            return gregory_omega1[n-j];
          }
          if(k==2){
            return gregory_omega2[n-j];
          }
          if(k==3){
            return gregory_omega3[n-j];
          }
          if(k==4){
            return gregory_omega4[n-j];
          }
          if(k==5){
            return gregory_omega5[n-j];
          }
        }
			}
		}
		

		private:
			int k_;
			double *poly_interp_;
			double *poly_interp1;
			double *poly_interp2;
			double *poly_interp3;
			double *poly_interp4;
			double *poly_interp5;
			double *poly_diff_;
			double *poly_integ_;
			double *bd_weights_;
			double *gregory_start_;
			double *gregory_start1;
			double *gregory_start2;
			double *gregory_start3;
			double *gregory_start4;
			double *gregory_start5;
			double *gregory_Omega_;
			double *gregory_Omega1;
			double *gregory_Omega2;
			double *gregory_Omega3;
			double *gregory_Omega4;
			double *gregory_Omega5;
			double *gregory_omega_;
			double *gregory_omega1;
			double *gregory_omega2;
			double *gregory_omega3;
			double *gregory_omega4;
			double *gregory_omega5;
			double *rcorr_;
	};
}//Namespace

#endif


