#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <chrono>
#include <boost/math/special_functions/chebyshev.hpp>

#include "utils.h"
#include "SCHFSolver.h"
#include "dyson.h"
#include "gfmol/utils.h"



int main(const int argc, char *argv[]) {
  using namespace NEdyson;

  std::string mat_file_s(argv[1]);
  std::string phi_file_s(argv[2]);
  int nt = std::stoi(argv[3]);
  int ntau = std::stoi(argv[4]);
  double beta = std::stod(argv[5]);
  double phi0 = std::stod(argv[6]);

  // objects
  int k = 5;
  GREEN G(nt, ntau, 2, -1);
  GREEN Sigma(nt, ntau, 2, -1);
  DYSON Dyson(nt, ntau, 2, k, gfmol::Mode::GF2);
	SCHFSolver SCsolver;
  cplx *H = new cplx[(nt+1)*2*2];
  double *Hr = new double[2*2];
  cplx *t0 = new cplx[(nt+1)*2*2];
  double *phi = new double[nt+1];

  std::vector<double> timings(nt+1);
  std::chrono::time_point<std::chrono::system_clock> start, end;

/*
	DMatrix x(ntau+1, 1);
	for(int i = 0; i <= ntau; i++) {
		x(i,0) = Dyson.Convolution().collocation().x_i()(i);
	}
	h5e::File out_file("/pauli-storage/tblommel/Legendre_points.h5", h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);
	h5e::dump<DMatrix>(out_file, "x", x);
*/

  // read in G, Sigma matsubara
  h5e::File mat_file(mat_file_s);
  ZMatrixMap(G.matptr(0), (ntau+1)*2*2, 1) = h5e::load<ZMatrix>(mat_file, "G");
  ZMatrixMap(Sigma.matptr(0), (ntau+1)*2*2, 1) = h5e::load<ZMatrix>(mat_file, "S");
  ZMatrixMap(H, 2, 2) = h5e::load<ZMatrix>(mat_file, "H");
  DMatrixMap(Hr, 2, 2) = ZMatrixMap(H, 2, 2).real();

  // read in phi
  h5e::File phi_file(phi_file_s);
  DMatrix phi_vec = h5e::load<DMatrix>(phi_file, "phi");
  DMatrixMap(phi, 1, nt+1) = phi0 * DMatrixMap(phi_vec.data(), 1, nt+1);
  ZMatrixMap(t0, (nt+1)*4, 1) = ZMatrix::Zero((nt+1)*4, 1);
  for(int tstp = 0; tstp <= nt; tstp++) {
    t0[tstp*4 + 0] = std::exp( cplx(0.,1.) * phi[tstp]);
    t0[tstp*4 + 3] = std::exp(-cplx(0.,1.) * phi[tstp]);
  }

  // calc GmatConvTens
  G.set_conv_tensor(Dyson.Convolution(), beta);

	// Gfree
	Dyson.G0_from_h0(G, 0, Hr, beta, 0.02);

  // bootstrap
  for(int iter = 0; iter <= 1000; iter++){
    double err = 0;

    // Update mean field & self energy
    for(int tstp = 0; tstp <= k; tstp++) { 
      ZMatrixMap(H + tstp*4, 2, 2) = ZMatrixMap(H,2,2);
      SCsolver.solve_2b(tstp, 2, t0, G, Sigma);
    }

    // Solve G Equation of Motion
    err = Dyson.dyson_start(G, Sigma, H, 0, beta, 0.02);

    std::cout<<"Bootstrapping iteration : "<<iter<<" | Error = "<<err<<std::endl;
    if(err<1e-10){
      break;
    }
  }

  // timestep
	for(int tstp = k+1; tstp <= nt; tstp++) {
    std::cout << "tstp " << tstp << std::endl;
		Dyson.Extrapolate(tstp, G);
    start = std::chrono::system_clock::now();

	  // Corrector
	  for(int iter = 0; iter < 5; iter++) {
	    // HF Contractions
			SCsolver.solve_Sigma_Fock(tstp, 2., G, H); 
	
	    // 2B Contractions
			SCsolver.solve_2b(tstp, 2., t0, G, Sigma);
	
	    Dyson.dyson_step(tstp, G, Sigma, H, 0, beta, 0.02);
	  }
    end = std::chrono::system_clock::now();
    timings[tstp] = (end-start).count();
	}

  // output
  cplx *rho_t = new cplx[(nt+1) * 4];
  for(int i = 0; i <= nt; i++) {
    G.get_dm(i, rho_t + i*4);
  }
  h5e::File out_file_res("/pauli-storage/tblommel/hodlr_SC/cubic_timing.h5", h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);
  h5e::dump<ZMatrix>(out_file_res, "rho", ZMatrixMap(rho_t, (nt+1)*4, 1));
  h5e::dump<DMatrix>(out_file_res, "time", DMatrixMap(timings.data(), nt+1, 1));

  
}
