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
  GREEN SLup(nt, ntau, nao, -1);
  GREEN SLdown(nt, ntau, nao, -1);

  std::vector<std::reference_wrapper<GREEN>> G = {Gup, Gdown};
  std::vector<std::reference_wrapper<GREEN>> Sigma = {Sup, Sdown};
  std::vector<std::reference_wrapper<GREEN>> SigmaLoop = {SLup, SLdown};

  DColVectorMap(Vija.data(),nao2*nalpha) = Eigen::VectorXd::Random(nao2*nalpha);

  molGF2SolverSpinDecomp test_solver(Vija);

  ZColVectorMap(Gup.retptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(Gdown.retptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);

  ZColVectorMap(Gup.lesptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(Gdown.lesptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);

  ZColVectorMap(Gup.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);
  ZColVectorMap(Gdown.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);

  ZColVectorMap SLUtest(Sigma[0].get().lesptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SLUloop(SigmaLoop[0].get().lesptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SLDtest(Sigma[1].get().lesptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SLDloop(SigmaLoop[1].get().lesptr(0,0), (nt+1)*(nt+2)/2*nao2);

  ZColVectorMap SRUtest(Sigma[0].get().retptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SRUloop(SigmaLoop[0].get().retptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SRDtest(Sigma[1].get().retptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SRDloop(SigmaLoop[1].get().retptr(0,0), (nt+1)*(nt+2)/2*nao2);

  ZColVectorMap STVUtest(Sigma[0].get().tvptr(0,0), (nt+1)*(ntau+1)*nao2);
  ZColVectorMap STVUloop(SigmaLoop[0].get().tvptr(0,0), (nt+1)*(ntau+1)*nao2);
  ZColVectorMap STVDtest(Sigma[1].get().tvptr(0,0), (nt+1)*(ntau+1)*nao2);
  ZColVectorMap STVDloop(SigmaLoop[1].get().tvptr(0,0), (nt+1)*(ntau+1)*nao2);

  for(int tstp=0; tstp<=nt; tstp++){
    test_solver.solve(tstp, Sigma, G);
    test_solver.solve_loop(tstp, SigmaLoop, G);
  }


  REQUIRE((SLUtest-SLUloop).norm() < 1e-10);
  REQUIRE((SLDtest-SLDloop).norm() < 1e-10);
  REQUIRE((SRUtest-SRUloop).norm() < 1e-10);
  REQUIRE((SRDtest-SRDloop).norm() < 1e-10);
  REQUIRE((STVUtest-STVUloop).norm() < 1e-10);
  REQUIRE((STVDtest-STVDloop).norm() < 1e-10);

}
