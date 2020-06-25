#ifndef SPECTRAL_DECL
#define SPECTRAL_DECL

#include <complex>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include "greens.h"
#include "elementops.h"
#include "utils.h"

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
		double *ptr(int t, int w) const ;

		// Input/Output
		void print_to_file(std::string file, int precision=16) const;
		void read_from_file(const char *file);
    void read_from_file(h5e::File File, std::string path);
    void print_to_file(h5e::File File, std::string path) const;
		
		// Actually Calculate
		void AfromG(const green_func &G, int nw, double wmax, double dt);
		void AfromG(const green_func &G, int nw, double wmax, double dt, const cplx *extdata, int nfit, int ntp);
};


class tti_spectral{
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
		tti_spectral();
		~tti_spectral();
		tti_spectral(int nt, int nw, int size1, double wmax, double dt);
		tti_spectral(const tti_spectral &A);
		tti_spectral &operator=(const tti_spectral &A);
		
		// Get sizes
		int nt(void) const {return nt_;}
		int nw(void) const {return nw_;}
		int es(void) const {return es_;}
		int size1(void) const {return size1_;}
		double wmax(void) const {return wmax_;}
		double dw(void) const {return dw_;}
		double dt(void) const {return dt_;}
		
		// Get pointer
		double *ptr(int w) const ;

		// Input/Output
		void print_to_file(std::string file, int precision=16) const;
		void read_from_file(const char *file);
    void read_from_file(h5e::File File, std::string path);
    void print_to_file(h5e::File File, std::string path) const;
		
		// Actually Calculate
		void AfromG(const green_func &G, int nw, double wmax, double dt);
		void AfromG(const green_func &G, int nw, double wmax, double dt, const cplx *extdata, int nfit, int ntp);
};

} // namespace

#endif
