#ifndef INTEGRATION_DECL
#define INTEGRATION_DECL

#include <cmath>
#include <iostream>
#include <complex>
#include <stdlib.h>

namespace integration{
	
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
			poly_diff_= new double [k1*k1];
			poly_integ_= new double [k1*k1*k1];
			bd_weights_= new double [k1+1];
			gregory_start_= new double [k1*k1];
			gregory_Omega_= new double [k1*k1];
			gregory_omega_= gregory_Omega_+(k*k1);
			rcorr_= new double [k1*k1*k1];
			k_=k;
			integration::make_poly_interp(k_,poly_interp_);
			integration::make_poly_diff(k_,poly_interp_,poly_diff_);
			integration::make_poly_integ(k_,poly_interp_,poly_integ_);
			integration::make_bd_weights(k_,bd_weights_);
			integration::make_start(k_,poly_interp_,gregory_start_);
			integration::make_Omega(k_,gregory_Omega_);
			integration::make_rcorr(k_,poly_interp_,rcorr_);
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
			integration::make_poly_interp(k_,poly_interp_);
			integration::make_poly_diff(k_,poly_interp_,poly_diff_);
			integration::make_poly_integ(k_,poly_interp_,poly_integ_);
			integration::make_bd_weights(k_,bd_weights_);
			integration::make_start(k_,poly_interp_,gregory_start_);
			integration::make_Omega(k_,gregory_Omega_);
			integration::make_rcorr(k_,poly_interp_,rcorr_);
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
			integration::make_poly_interp(k_,poly_interp_);
			integration::make_poly_diff(k_,poly_interp_,poly_diff_);
			integration::make_poly_integ(k_,poly_interp_,poly_integ_);
			integration::make_bd_weights(k_,bd_weights_);
			integration::make_start(k_,poly_interp_,gregory_start_);
			integration::make_Omega(k_,gregory_Omega_);
			integration::make_rcorr(k_,poly_interp_,rcorr_);
		}
		
		int k(void) {return k_;}
		double poly_interp(int a, int l) {return poly_interp_[a*(k_+1)+l];}
		double poly_diff(int m, int l) {return poly_diff_[m*(k_+1)+l];}
		double poly_integ(int m, int n, int l) {return poly_integ_[m*(k_+1)*(k_+1)+n*(k_+1)+l];}
		double bd_weights(int l) {return bd_weights_[l];}
		double gregory_weights(int n, int j){
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
			
//			if(n>2*k_){
//				if(j>k_&&j<n-k_) return 1;
//				else if(j<=k_) return gregory_omega_[j];
//				else return gregory_omega_[n-j];
//			}
//			else if(n>k_) return gregory_Omega_[(n-k_-1)*(k_+1)+j];
//			else return gregory_start_[n*(k_+1)+j];
		}
		double omega(int j) {return gregory_omega_[j];}
		double start(int i, int j)  {return gregory_start_[i*(k_+1)+j];}
		double Omega(int i, int j)  {return gregory_Omega_[i*(k_+1)+j];}
		double rcorr(int m, int j, int l) {return rcorr_[m*(k_+1)*(k_+1)+j*(k_+1)+l];}
		

		private:
			int k_;
			double *poly_interp_;
			double *poly_diff_;
			double *poly_integ_;
			double *bd_weights_;
			double *gregory_start_;
			double *gregory_Omega_;
			double *gregory_omega_;
			double *rcorr_;
	};
}//Namespace

#endif


