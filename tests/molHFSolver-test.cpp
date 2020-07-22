#include "molHFSolver.h"
#include "greens.h"
#include "tests-def.h"

TEST_CASE("Sigma HF")
{
  using namespace NEdyson;
  
  int nt = 3;
  int na = 2;
  int na2 = na * na;
  int na4 = na2 * na2;

  DTensor<4> U_int(na,na,na,na);
  ZTensor<2> rho(na,na);
  ZTensor<3> hmf(nt+1, na, na);
  ZTensor<3> hmf_loop(nt+1, na, na);

  DColVectorMap(U_int.data(),na4) = Eigen::VectorXd::Random(na4);
  ZColVectorMap(rho.data(), na2) = Eigen::VectorXcd::Random(na2);
  ZColVectorMap(hmf.data(), (nt+1)*na2) = Eigen::VectorXcd::Random((nt+1)*na2);
  hmf_loop = hmf;

  molHFSolver test_solver(U_int);
  
  for(int t=0;t<=nt;t++){
    test_solver.solve_HF(t,hmf,rho);
    test_solver.solve_HF_loop(t,hmf_loop,rho);
  }


  double err = (ZColVectorMap(hmf.data(), (nt+1)*na2) - ZColVectorMap(hmf_loop.data(), (nt+1)*na2)).norm();
  REQUIRE(err < 1e-10);
}


TEST_CASE("Decomp Sigma HF")
{
  using namespace NEdyson;
  
  int nt = 3;
  int nao = 2;
  int nalpha = 5;
  int nao2 = nao * nao;

  DTensor<3> Vija(nao,nao,nalpha);
  ZTensor<2> rho(nao,nao);
  ZTensor<3> hmf(nt+1, nao, nao);
  ZTensor<3> hmf_loop(nt+1, nao, nao);

  DColVectorMap(Vija.data(),nao2*nalpha) = Eigen::VectorXd::Random(nao2*nalpha);
  ZColVectorMap(rho.data(), nao2) = Eigen::VectorXcd::Random(nao2);
  ZColVectorMap(hmf.data(), (nt+1)*nao2) = Eigen::VectorXcd::Random((nt+1)*nao2);
  hmf_loop = hmf;

  molHFSolverDecomp test_solver(Vija);
  
  for(int t=0;t<=nt;t++){
    test_solver.solve_HF(t,hmf,rho);
    test_solver.solve_HF_loop(t,hmf_loop,rho);
  }

  double err = (ZColVectorMap(hmf.data(), (nt+1)*nao2) - ZColVectorMap(hmf_loop.data(), (nt+1)*nao2)).norm();
  REQUIRE(err < 1e-10);
}


TEST_CASE("Sigma HF Spin")
{
  using namespace NEdyson;
  
  int nt = 3;
  int na = 2;
  int na2 = na * na;
  int na4 = na2 * na2;

  DTensor<4> U_int(na,na,na,na);
  ZTensor<3> rho(2,na,na);
  ZTensor<4> hmf(2,nt+1, na, na);
  ZTensor<4> hmf_loop(2,nt+1, na, na);

  DColVectorMap(U_int.data(),na4) = Eigen::VectorXd::Random(na4);
  ZColVectorMap(rho.data(), 2*na2) = Eigen::VectorXcd::Random(2*na2);
  ZColVectorMap(hmf.data(), 2*(nt+1)*na2) = Eigen::VectorXcd::Random(2*(nt+1)*na2);
  hmf_loop = hmf;

  molHFSolverSpin test_solver(U_int);
  
  for(int t=0;t<=nt;t++){
    test_solver.solve_HF(t,hmf,rho);
    test_solver.solve_HF_loop(t,hmf_loop,rho);
  }

  double err = (ZColVectorMap(hmf.data(), 2*(nt+1)*na2) - ZColVectorMap(hmf_loop.data(), 2*(nt+1)*na2)).norm();
  REQUIRE(err < 1e-10);
}



TEST_CASE("Decomp Sigma HF Spin")
{
  using namespace NEdyson;
  
  int nt = 3;
  int nao = 2;
  int nalpha = 4;
  int nao2 = nao * nao;
  int nao4 = nao2 * nao2;

  DTensor<3> Vija(nao,nao,nalpha);
  ZTensor<3> rho(2,nao,nao);
  ZTensor<4> hmf(2, nt+1, nao, nao);
  ZTensor<4> hmf_loop(2, nt+1, nao, nao);

  DColVectorMap(Vija.data(),nao2*nalpha) = Eigen::VectorXd::Random(nao2*nalpha);
  ZColVectorMap(rho.data(), 2*nao2) = Eigen::VectorXcd::Random(2*nao2);
  ZColVectorMap(hmf.data(), 2*(nt+1)*nao2) = Eigen::VectorXcd::Random(2*(nt+1)*nao2);
  hmf_loop = hmf;

  molHFSolverSpinDecomp test_solver(Vija);
  
  for(int t=0;t<=nt;t++){
    test_solver.solve_HF(t,hmf,rho);
    test_solver.solve_HF_loop(t,hmf_loop,rho);
  }

  double err = (ZColVectorMap(hmf.data(), 2*(nt+1)*nao2) - ZColVectorMap(hmf_loop.data(), 2*(nt+1)*nao2)).norm();
  REQUIRE(err < 1e-10);
}
