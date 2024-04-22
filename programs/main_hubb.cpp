#include "params.h"
#include "molNEsim.h"
#include "hub2bSolver.h"
#include "gfmol/repn.h"
#include <iostream>
#include "gfmol/utils.h"

int main(const int argc, char *const *const argv)
{
  using namespace NEdyson;
  cplx cplxi(0.,1.);

  int nao = 4;
  int k = std::stoi(argv[3]);
  int nt = std::stoi(argv[4]) + 1;
  double U = std::stod(argv[1]);
  double h = std::stod(argv[2]);

  std::cout << k << std::endl;
  std::cout << nt << std::endl;
  std::cout << U << std::endl;
  std::cout << h << std::endl;


  GREEN G(nt, 8, nao, -1);
  GREEN Sigma(nt, 8, nao, -1);
  DYSON dyson_solver(nt, 8, nao, k, gfmol::Mode::GF2);
  hubGF2Solver Hub(U, nao);

  ZMatrix Hmf(nt, nao*nao);
  ZMatrix rhot(nt,nao*nao);
  ZMatrix rho(nao,nao);
  ZMatrix H0(nao,nao);

  H0(1,0) = -1;
  H0(2,0) = -1;
  H0(0,1) = -1;
  H0(0,2) = -1;
  H0(3,1) = -1;
  H0(3,2) = -1;
  H0(1,3) = -1;
  H0(2,3) = -1;
  
  H0(0,0) = -1;
  H0(1,1) = -1;
  H0(2,2) = -1;
  H0(3,3) = -1;

  G.lesptr(0,0)[0*4 + 0] = cplxi;
  G.lesptr(0,0)[1*4 + 1] = 0.6*cplxi;
  G.lesptr(0,0)[2*4 + 2] = 0.6*cplxi;
  G.lesptr(0,0)[3*4 + 3] = 0.6*cplxi;

  ZMatrixMap(G.retptr(0,0), 4, 4) = cplxi * ZMatrix::Identity(4,4);

  // Set up HF initially
  G.get_dm(0, rho);
  std::cout << rho << std::endl;
  ZMatrixMap(Hmf.data(), nao, nao) = H0;
  Hub.solve_HF(0, Hmf.data(), rho);
  for(int i = 1; i <= k; i++) {
    ZMatrixMap(Hmf.data() + i*nao*nao, nao, nao) = ZMatrixMap(Hmf.data(), nao, nao);
  }

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();


  double booterr = 1;
  while(booterr > 2e-10) {
    booterr = dyson_solver.dyson_start(G, Sigma, Hmf.data(), 0., 0., h);
    G.get_dm(0, rho);
    std::cout << rho << std::endl;
    std::cout << "boot err " << booterr << std::endl;
    
    for(int i = 0; i <= k; i++) {
      ZMatrixMap(Hmf.data() + i*nao*nao, nao, nao) = H0;
      G.get_dm(i, rho);
      Hub.solve_HF(i, Hmf.data(), rho);

      Hub.solve(i, Sigma, G);
    }
  }


  for(int t = k+1; t < nt; t++) {
    dyson_solver.Extrapolate(t, G);

    ZMatrixMap(Hmf.data() + t*nao*nao, nao, nao) = H0;
    G.get_dm(t, rho);
    Hub.solve_HF(t, Hmf.data(), rho);

    Hub.solve(t, Sigma, G);

    double steperr = 1;
    int count = 0;
    while(steperr > 1e-10 && count < 3) {
      steperr = dyson_solver.dyson_step(t, G, Sigma, Hmf.data(), 0., 0., h);

      std::cout << t << " " << steperr << std::endl;
      count++;

      ZMatrixMap(Hmf.data() + t*nao*nao, nao, nao) = H0;
      G.get_dm(t, rho);
      Hub.solve_HF(t, Hmf.data(), rho);

      Hub.solve(t, Sigma, G);
    }
  }

  end = std::chrono::system_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

  for(int i = 0; i <= nt; i++) {
    G.get_dm(i, rhot.data()+i*nao*nao);  
  }
  
  h5::File outfile(argv[5], h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);
  h5e::dump(outfile, "rho", ZColVectorMap(rhot.data(), nt*nao*nao));
  h5e::dump(outfile, "h", h);
  h5e::dump(outfile, "w", (int)(duration.count()));

  return 0;   
}
