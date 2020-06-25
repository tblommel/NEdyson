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


TEST_CASE("Decomp GF2")
{
  using namespace NEdyson;

  int nao = 2;
  int nalpha = 3;
  int nt = 4;
  int ntau = 5;
  int nao2 = nao*nao;
  
  DTensor<3> Vija(nao,nao,nalpha);

  GREEN G(nt, ntau, nao, -1);
  GREEN S(nt, ntau, nao, -1);
  GREEN SL(nt, ntau, nao, -1);

  DColVectorMap(Vija.data(),nao2*nalpha) = Eigen::VectorXd::Random(nao2*nalpha);

  molGF2SolverDecomp test_solver(Vija);

  ZColVectorMap(G.retptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(G.lesptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(G.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);

  ZColVectorMap SLtest(S.lesptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SLloop(SL.lesptr(0,0), (nt+1)*(nt+2)/2*nao2);

  ZColVectorMap SRtest(S.retptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SRloop(SL.retptr(0,0), (nt+1)*(nt+2)/2*nao2);

  ZColVectorMap STVtest(S.tvptr(0,0), (nt+1)*(ntau+1)*nao2);
  ZColVectorMap STVloop(SL.tvptr(0,0), (nt+1)*(ntau+1)*nao2);

  for(int tstp=0; tstp<=nt; tstp++){
    test_solver.solve(tstp, S, G);
    test_solver.solve_loop(tstp, SL, G);
  }


  REQUIRE((SLtest-SLloop).norm() < 1e-10);
  REQUIRE((SRtest-SRloop).norm() < 1e-10);
  REQUIRE((STVtest-STVloop).norm() < 1e-10);

}


TEST_CASE("Spin GF2")
{
  using namespace NEdyson;

  int nao = 2;
  int nt = 4;
  int ntau = 5;
  int nao2 = nao*nao;
  
  DTensor<4> Uint(nao,nao,nao,nao);

  GREEN Gup(nt, ntau, nao, -1);
  GREEN Gdown(nt, ntau, nao, -1);
  GREEN Sup(nt, ntau, nao, -1);
  GREEN Sdown(nt, ntau, nao, -1);
  GREEN SLup(nt, ntau, nao, -1);
  GREEN SLdown(nt, ntau, nao, -1);

  std::vector<std::reference_wrapper<GREEN>> G = {Gup, Gdown};
  std::vector<std::reference_wrapper<GREEN>> Sigma = {Sup, Sdown};
  std::vector<std::reference_wrapper<GREEN>> SigmaLoop = {SLup, SLdown};

  DColVectorMap(Uint.data(),nao2*nao2) = Eigen::VectorXd::Random(nao2*nao2);

  molGF2SolverSpin test_solver(Uint);

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


TEST_CASE("Full GF2")
{
  using namespace NEdyson;

  int nao = 2;
  int nalpha = 3;
  int nt = 4;
  int ntau = 5;
  int nao2 = nao*nao;
  
  DTensor<4> Uijkl(nao,nao,nao,nao);

  GREEN G(nt, ntau, nao, -1);
  GREEN S(nt, ntau, nao, -1);
  GREEN SL(nt, ntau, nao, -1);

  DColVectorMap(Uijkl.data(),nao2*nao2) = Eigen::VectorXd::Random(nao2*nao2);

  molGF2Solver test_solver(Uijkl);

  ZColVectorMap(G.retptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(G.lesptr(0,0),(nt+1)*(nt+2)/2*nao2) = Eigen::VectorXcd::Random((nt+1)*(nt+2)/2*nao2);
  ZColVectorMap(G.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);

  ZColVectorMap SLtest(S.lesptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SLloop(SL.lesptr(0,0), (nt+1)*(nt+2)/2*nao2);

  ZColVectorMap SRtest(S.retptr(0,0), (nt+1)*(nt+2)/2*nao2);
  ZColVectorMap SRloop(SL.retptr(0,0), (nt+1)*(nt+2)/2*nao2);

  ZColVectorMap STVtest(S.tvptr(0,0), (nt+1)*(ntau+1)*nao2);
  ZColVectorMap STVloop(SL.tvptr(0,0), (nt+1)*(ntau+1)*nao2);

  for(int tstp=0; tstp<=nt; tstp++){
    test_solver.solve(tstp, S, G);
    test_solver.solve_loop(tstp, SL, G);
  }


  REQUIRE((SLtest-SLloop).norm() < 1e-10);
  REQUIRE((SRtest-SRloop).norm() < 1e-10);
  REQUIRE((STVtest-STVloop).norm() < 1e-10);

}


TEST_CASE("TTI Decomp Spin GF2")
{
  using namespace NEdyson;

  int nao = 2;
  int nalpha = 3;
  int nt = 4;
  int ntau = 5;
  int nao2 = nao*nao;
  
  DTensor<3> Vija(nao,nao,nalpha);

  TTI_GREEN Gup(nt, ntau, nao, -1);
  TTI_GREEN Gdown(nt, ntau, nao, -1);
  TTI_GREEN Sup(nt, ntau, nao, -1);
  TTI_GREEN Sdown(nt, ntau, nao, -1);
  TTI_GREEN SLup(nt, ntau, nao, -1);
  TTI_GREEN SLdown(nt, ntau, nao, -1);

  std::vector<std::reference_wrapper<TTI_GREEN>> G = {Gup, Gdown};
  std::vector<std::reference_wrapper<TTI_GREEN>> Sigma = {Sup, Sdown};
  std::vector<std::reference_wrapper<TTI_GREEN>> SigmaLoop = {SLup, SLdown};

  DColVectorMap(Vija.data(),nao2*nalpha) = Eigen::VectorXd::Random(nao2*nalpha);

  tti_molGF2SolverSpinDecomp test_solver(Vija);

  ZColVectorMap(Gup.retptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  ZColVectorMap(Gdown.retptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);

  ZColVectorMap(Gup.lesptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  ZColVectorMap(Gdown.lesptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);

  ZColVectorMap(Gup.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);
  ZColVectorMap(Gdown.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);

  ZColVectorMap SLUtest(Sigma[0].get().lesptr(0), (nt+1)*nao2);
  ZColVectorMap SLUloop(SigmaLoop[0].get().lesptr(0), (nt+1)*nao2);
  ZColVectorMap SLDtest(Sigma[1].get().lesptr(0), (nt+1)*nao2);
  ZColVectorMap SLDloop(SigmaLoop[1].get().lesptr(0), (nt+1)*nao2);

  ZColVectorMap SRUtest(Sigma[0].get().retptr(0), (nt+1)*nao2);
  ZColVectorMap SRUloop(SigmaLoop[0].get().retptr(0), (nt+1)*nao2);
  ZColVectorMap SRDtest(Sigma[1].get().retptr(0), (nt+1)*nao2);
  ZColVectorMap SRDloop(SigmaLoop[1].get().retptr(0), (nt+1)*nao2);

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


TEST_CASE("TTI Decomp GF2")
{
  using namespace NEdyson;

  int nao = 2;
  int nalpha = 3;
  int nt = 4;
  int ntau = 5;
  int nao2 = nao*nao;
  
  DTensor<3> Vija(nao,nao,nalpha);

  TTI_GREEN G(nt, ntau, nao, -1);
  TTI_GREEN S(nt, ntau, nao, -1);
  TTI_GREEN SL(nt, ntau, nao, -1);

  DColVectorMap(Vija.data(),nao2*nalpha) = Eigen::VectorXd::Random(nao2*nalpha);

  tti_molGF2SolverDecomp test_solver(Vija);

  ZColVectorMap(G.retptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  ZColVectorMap(G.lesptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  ZColVectorMap(G.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);

  ZColVectorMap SLtest(S.lesptr(0), (nt+1)*nao2);
  ZColVectorMap SLloop(SL.lesptr(0), (nt+1)*nao2);

  ZColVectorMap SRtest(S.retptr(0), (nt+1)*nao2);
  ZColVectorMap SRloop(SL.retptr(0), (nt+1)*nao2);

  ZColVectorMap STVtest(S.tvptr(0,0), (nt+1)*(ntau+1)*nao2);
  ZColVectorMap STVloop(SL.tvptr(0,0), (nt+1)*(ntau+1)*nao2);

  for(int tstp=0; tstp<=nt; tstp++){
    test_solver.solve(tstp, S, G);
    test_solver.solve_loop(tstp, SL, G);
  }


  REQUIRE((SLtest-SLloop).norm() < 1e-10);
  REQUIRE((SRtest-SRloop).norm() < 1e-10);
  REQUIRE((STVtest-STVloop).norm() < 1e-10);

}


TEST_CASE("TTI Spin GF2")
{
  using namespace NEdyson;

  int nao = 2;
  int nt = 4;
  int ntau = 5;
  int nao2 = nao*nao;
  
  DTensor<4> Uint(nao,nao,nao,nao);

  TTI_GREEN Gup(nt, ntau, nao, -1);
  TTI_GREEN Gdown(nt, ntau, nao, -1);
  TTI_GREEN Sup(nt, ntau, nao, -1);
  TTI_GREEN Sdown(nt, ntau, nao, -1);
  TTI_GREEN SLup(nt, ntau, nao, -1);
  TTI_GREEN SLdown(nt, ntau, nao, -1);

  std::vector<std::reference_wrapper<TTI_GREEN>> G = {Gup, Gdown};
  std::vector<std::reference_wrapper<TTI_GREEN>> Sigma = {Sup, Sdown};
  std::vector<std::reference_wrapper<TTI_GREEN>> SigmaLoop = {SLup, SLdown};

  DColVectorMap(Uint.data(),nao2*nao2) = Eigen::VectorXd::Random(nao2*nao2);

  tti_molGF2SolverSpin test_solver(Uint);

  ZColVectorMap(Gup.retptr(0),  (nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  ZColVectorMap(Gdown.retptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);

  ZColVectorMap(Gup.lesptr(0),  (nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  ZColVectorMap(Gdown.lesptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);

  ZColVectorMap(Gup.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);
  ZColVectorMap(Gdown.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);

  ZColVectorMap SLUtest(Sigma[0].get().lesptr(0), (nt+1)*nao2);
  ZColVectorMap SLUloop(SigmaLoop[0].get().lesptr(0), (nt+1)*nao2);
  ZColVectorMap SLDtest(Sigma[1].get().lesptr(0), (nt+1)*nao2);
  ZColVectorMap SLDloop(SigmaLoop[1].get().lesptr(0), (nt+1)*nao2);

  ZColVectorMap SRUtest(Sigma[0].get().retptr(0), (nt+1)*nao2);
  ZColVectorMap SRUloop(SigmaLoop[0].get().retptr(0), (nt+1)*nao2);
  ZColVectorMap SRDtest(Sigma[1].get().retptr(0), (nt+1)*nao2);
  ZColVectorMap SRDloop(SigmaLoop[1].get().retptr(0), (nt+1)*nao2);

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


TEST_CASE("TTI Full GF2")
{
  using namespace NEdyson;

  int nao = 2;
  int nalpha = 3;
  int nt = 4;
  int ntau = 5;
  int nao2 = nao*nao;
  
  DTensor<4> Uijkl(nao,nao,nao,nao);

  TTI_GREEN G(nt, ntau, nao, -1);
  TTI_GREEN S(nt, ntau, nao, -1);
  TTI_GREEN SL(nt, ntau, nao, -1);

  DColVectorMap(Uijkl.data(),nao2*nao2) = Eigen::VectorXd::Random(nao2*nao2);

  tti_molGF2Solver test_solver(Uijkl);

  ZColVectorMap(G.retptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  ZColVectorMap(G.lesptr(0),(nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  ZColVectorMap(G.tvptr(0,0),(nt+1)*(ntau+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*(ntau+1)*nao2);

  ZColVectorMap SLtest(S.lesptr(0),  (nt+1)*nao2);
  ZColVectorMap SLloop(SL.lesptr(0), (nt+1)*nao2);

  ZColVectorMap SRtest(S.retptr(0),  (nt+1)*nao2);
  ZColVectorMap SRloop(SL.retptr(0), (nt+1)*nao2);

  ZColVectorMap STVtest(S.tvptr(0,0), (nt+1)*(ntau+1)*nao2);
  ZColVectorMap STVloop(SL.tvptr(0,0), (nt+1)*(ntau+1)*nao2);

  for(int tstp=0; tstp<=nt; tstp++){
    test_solver.solve(tstp, S, G);
    test_solver.solve_loop(tstp, SL, G);
  }


  REQUIRE((SLtest-SLloop).norm() < 1e-10);
  REQUIRE((SRtest-SRloop).norm() < 1e-10);
  REQUIRE((STVtest-STVloop).norm() < 1e-10);

}
