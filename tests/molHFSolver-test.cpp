#include "molHFSolver.h"
#include "greens.h"
#include "function.h"
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
    test_solver.solve_HF(t,rho,hmf);
    test_solver.solve_HF_loop(t,rho,hmf_loop);
  }

  double err = (ZColVectorMap(hmf.data(), (nt+1)*na2) - ZColVectorMap(hmf_loop.data(), (nt+1)*na2)).norm();
  REQUIRE(err < 1e-10);
}
