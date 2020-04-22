#ifndef SPECTRAL_DECL
#define SPECTRAL_DECL

#include <complex>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include "greens.hpp"
#include "elementops.hpp"

#define cplx std::complex<double>

namespace NEdyson{

class spectral{
	private:
		int nt_;
		int nw_;
		int size1_;
		int es_;
		double wmax_;
		double dw_;
		double dt_;
		double *A_;

	public:
		// Construct and Destruct
		spectral();
		~spectral();
		spectral(int nt, int nw, int size1, double wmax, double dt);
		spectral(const spectral &A);
		spectral &operator=(const spectral &A);
		
		// Get sizes
		int nt(void) const {return nt_;}
		int nw(void) const {return nw_;}
		int es(void) const {return es_;}
		int size1(void) const {return size1_;}
		double wmax(void) const {return wmax_;}
		double dw(void) const {return dw_;}
		double dt(void) const {return dt_;}
		
		// Get pointer
		double *ptr(int t, int w);

		// Input/Output
		void print_to_file(std::string file, int precision=16);
		void read_from_file(const char *file);
		
		// Actually Calculate
		void AfromG(green_func &G, int nw, double wmax, double dt);
		void AfromG(green_func &G, int nw, double wmax, double dt, cplx *extdata, int nfit, int ntp);
};













} // namespace

#endif
