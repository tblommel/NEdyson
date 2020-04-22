#ifndef HUBB_OPS_DECL
#define HUBB_OPS_DECL

#include <Eigen/Eigen>
#include <complex>
#include "greens.hpp"
#include "function.hpp"
#include "greens_tstp.hpp"
#include "integration.hpp"
#include "vie2.hpp"

#define CFUNC NEdyson::function
#define GREEN NEdyson::green_func
#define GREEN_TSTP NEdyson::green_func_tstp
#define cdmatrix Eigen::MatrixXcd
#define INTEG integration::Integrator


namespace Hubb{
	void Ham_MF(int tstp, GREEN &G, CFUNC &U, cdmatrix &h0, CFUNC &hmf);

	void Sigma_2B(int tstp, GREEN &G, CFUNC &Ut, GREEN &Sigma);

	void Polarization(int tstp, GREEN &G, GREEN_TSTP &Pol);

	void Bubble2(int tstp, GREEN &Sigma, int s1, int s2, GREEN &G, int g1, int g2, GREEN_TSTP &Pol, int p1, int p2);
	void Bubble1(int tstp, GREEN_TSTP &Pol, int p1, int p2, GREEN &A, int a1, int a2, GREEN &B, int b1, int b2);

	void GenTPP(int tstp, double dt, double beta, GREEN &G, GREEN &Phi, CFUNC &Ut, GREEN &UxPhi, GREEN &PhixU, GREEN &TPP, INTEG &I);
	
	void GenTPP(double dt, double beta, GREEN &G, GREEN &Phi, CFUNC &Ut, GREEN &PhixU, GREEN &UxPhi, GREEN &TPP, INTEG &I);

	void Sigma_TPP(int tstp, GREEN &G, CFUNC &Ut, GREEN &TPP, GREEN &Sigma);
}

#endif
