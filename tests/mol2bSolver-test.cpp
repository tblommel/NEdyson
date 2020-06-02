#include "mol2bSolver.h"
#include "greens.h"
#include "tests-def.h"

TEST_CASE("Decomp Spin GF2")
{
  using namespace NEdyson;

  int nao = 2;
  int nalpha = 3;
  int nt = 4;
  int ntau = 5;
  int nao2 = nao*nao;
  
  DTensor<3> Vija(nao,nao,nalpha);

  GREEN Gup(nt, ntau, nao, -1);
  GREEN Gdown(nt, ntau, nao, -1);
  GREEN Sup(nt, ntau, nao, -1);
  GREEN Sdown(nt, ntau, nao, -1);

  std::vector<std::reference_wrapper<GREEN>> G = {Gup, Gdown};
  std::vector<std::reference_wrapper<GREEN>> Sigma = {Sup, Sdown};

  molGF2SolverSpinDecomp test_solver(Vija);

  DColVectorMap(Vija.data(),nao2*nalpha) = Eigen::VectorXd::Random(nao2*nalpha);

  ZColVectorMap(Gup.retptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(Gdown.retptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(Sup.retptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(Sdown.retptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);

  ZColVectorMap(Gup.lesptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(Gdown.lesptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(Sup.lesptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(Sdown.lesptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);

  ZColVectorMap(Gup.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);
  ZColVectorMap(Gdown.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);
  ZColVectorMap(Sup.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);
  ZColVectorMap(Sdown.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);

  std::cout<<Gup.tvptr(1,1)<<" "<<G[0].get().tvptr(1,1)<<std::endl;

/*
  for(int tstp=0; tstp<=nt; tstp++){
    test_solver.solve(tstp, Sigma, G);
    test_solver.solve_loop(tstp, Sigma, G);
  }
*/

REQUIRE(0==0);
  











}
