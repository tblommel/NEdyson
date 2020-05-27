#include "molHFSolver.h"
#include "greens.h"
#include "function.h"
#include "tests-def.h"

TEST_CASE("Sigma HF")
{
  using namespace NEdyson;
  
  int nt = 3;
  int na = 3;
  int na2 = na * na;
  int na4 = na2 * na2;

  DTensor<4> U_int(na,na,na,na);
  ZMatrix rho(na,na);
  CFUNC hmf(nt, na);
  CFUNC hmf_loop(nt, na);

  DColVectorMap(U_int.data(),na4) = Eigen::VectorXd::Random(na4);
  ZColVectorMap(rho.data(), na2) = Eigen::VectorXcd::Random(na2);
  ZColVectorMap(hmf.ptr(-1), (nt+2)*na2) = Eigen::VectorXcd::Random((nt+2)*na2);
  hmf_loop = hmf;

  molHFSolver test_solver(U_int);
  
  for(int t=0;t<=nt;t++){
    test_solver.solve_HF(t,rho,hmf);
    test_solver.solve_HF_loop(t,rho,hmf_loop);
  }

  double err = (ZColVectorMap(hmf.ptr(-1), (nt+2)*na2) - ZColVectorMap(hmf_loop.ptr(-1), (nt+2)*na2)).norm();
  REQUIRE(err < 1e-10);
}
