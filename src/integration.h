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
	void make_poly_diff(int k, double *D);
	void make_poly_integ(int k, double *I);
	void make_bd_weights(int k, double *BD);
	void make_start(int k, double *s);
	void make_Omega(int k, double *O);
	void make_ex_weights(int k, double *P, double *E);


	class Integrator{
		public:

		Integrator(int k){
			int k1=k+1;
			if(k<0){std::cout << "Integrator: k out of range"<<std::endl; abort();}
			poly_interp_= new double [k1*k1];
			poly_diff_= new double [k1*k1];
			poly_integ_= new double [k1*k1*k1];
			bd_weights_= new double [k1+1];
			ex_weights_= new double [k1];
			gregory_start_= new double [k1*k1];
			gregory_Omega_= new double [k1*k1];
			gregory_omega_= gregory_Omega_+(k*k1);
			k_=k;
			Integration::make_poly_interp(k_, poly_interp_);
			Integration::make_poly_diff(k_, poly_diff_);
			Integration::make_poly_integ(k_, poly_integ_);
			Integration::make_bd_weights(k_,bd_weights_);
			Integration::make_start(k_, gregory_start_);
			Integration::make_Omega(k_, gregory_Omega_);
			Integration::make_ex_weights(k_, poly_interp_, ex_weights_);
		}

		~Integrator(){
			delete [] poly_interp_;
			delete [] poly_diff_;
			delete [] poly_integ_;
			delete [] bd_weights_;
			delete [] gregory_start_;
			delete [] gregory_Omega_;
			delete [] ex_weights_;
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
		double ex_weights(int j) const {
      assert(j<=k_);
      return ex_weights_[j];
    }
    double *ex_weights_ptr() const {
      return ex_weights_;
    }

		private:
			int k_;
			double *poly_interp_;
			double *poly_diff_;
			double *poly_integ_;
			double *bd_weights_;
			double *ex_weights_;
			double *gregory_start_;
			double *gregory_Omega_;
			double *gregory_omega_;
	};
}//Namespace

#endif


