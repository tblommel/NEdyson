#include "dyson.h"
//#include "omp_integrals.h"
#include "utils.h"
#include "tests-def.h"

using namespace NEdyson;

TEST_CASE("Integrals"){
  int Nt=20, Ntau=30, Nsites=2;
  double beta=9, dt=0.1;

  GREEN A = GREEN(Nt, Ntau, Nsites, -1);
  GREEN B = GREEN(Nt, Ntau, Nsites, -1);
  GREEN Acc = GREEN(Nt, Ntau, Nsites, -1);
  GREEN Bcc = GREEN(Nt, Ntau, Nsites, -1);
  GREEN Cact = GREEN(Nt, Ntau, Nsites, -1);
  GREEN Ctest= GREEN(Nt, Ntau, Nsites, -1);
  DYSON Dyson = DYSON(Nt, Ntau, Nsites, 5);

  h5e::File File(std::string(TEST_DATA_DIR)+"/int_test_data.h5");
  A.read_from_file(File, "A");
  Acc.read_from_file(File, "Acc");
  B.read_from_file(File, "B");
  Bcc.read_from_file(File, "Bcc");
  Cact.read_from_file(File, "C");

  for(int tstp=0;tstp<=Nt;tstp++)
    Dyson.Ctv_tstp(tstp, Ctest, A, Acc, B, Bcc, beta, dt);

  ZColVectorMap CTVtstpVec = ZColVectorMap(Ctest.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
  ZColVectorMap CTVtstpCVec = ZColVectorMap(Cact.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
  REQUIRE((CTVtstpVec-CTVtstpCVec).norm() < 1e-11);

  cplx *tvvt = new cplx[984];
  int count=0;
  for(int tstp=0; tstp<=Nt; tstp++) {
    Dyson.Cles3_tstp(A, Acc, B, Bcc, tstp, beta, tvvt+count);
    count+=((tstp>5)?(tstp+1)*4:6*4);
  }
  
  ZColVectorMap CTVVTtstpVec = ZColVectorMap(tvvt, 984);
  ZColVector CTVVTtstpCVec = h5e::load<ZColVector>(File, "tvvt");
  REQUIRE((CTVVTtstpVec-CTVVTtstpCVec).norm() < 1e-11);
  delete[] tvvt;

  cplx *lesadv = new cplx[984];
  memset(lesadv,0,sizeof(cplx)*984);
  count=0;
  for(int tstp=0; tstp<=Nt; tstp++) {
    Dyson.Cles2_tstp(A, Acc, B, Bcc, tstp, dt, lesadv+count);
    count+=((tstp>5)?(tstp+1)*4:6*4);
  }

  ZColVectorMap CLAtstpVec = ZColVectorMap(lesadv, 984);
  ZColVector CLAtstpCVec = h5e::load<ZColVector>(File, "lesadv");
  REQUIRE((CLAtstpVec-CLAtstpCVec).norm() < 1e-11);
  delete[] lesadv;
  
/*
  Cact.read_from_file(File, "Comp");
  for(int t=-1; t<=Nt; t++) Ctest.set_tstp_zero(t);
  for(int t=0;t<=Nt;t++){
    std::vector<bool> mask_rl(t+1,true);
    std::vector<bool> mask_tv(Ntau+1,true);
    incr_convolution_ret(t, mask_rl, Ctest, A, Acc, B, Bcc, I, dt);
    incr_convolution_tv(t, mask_tv, Ctest, A, Acc, B, Bcc, I, beta, dt);
    incr_convolution_les(t, mask_rl, Ctest, A, Acc, B, Bcc, I, beta, dt);    
  }

  ZColVectorMap CRomptstpVec = ZColVectorMap(Ctest.retptr(0,0), (Nt+1)*(Nt+2)/2*4);
  ZColVectorMap CRomptstpCVec = ZColVectorMap(Cact.retptr(0,0), (Nt+1)*(Nt+2)/2*4);
  REQUIRE((CRomptstpVec-CRomptstpCVec).norm() < 1e-11);

  ZColVectorMap CLomptstpVec = ZColVectorMap(Ctest.lesptr(0,0), (Nt+1)*(Nt+2)/2*4);
  ZColVectorMap CLomptstpCVec = ZColVectorMap(Cact.lesptr(0,0), (Nt+1)*(Nt+2)/2*4);
  REQUIRE((CLomptstpVec-CLomptstpCVec).norm() < 1e-11);

  ZColVectorMap CTVomptstpVec = ZColVectorMap(Ctest.tvptr(0,0), (Ntau+1)*(Nt+1)*4);
  ZColVectorMap CTVomptstpCVec = ZColVectorMap(Cact.tvptr(0,0), (Ntau+1)*(Nt+1)*4);
  REQUIRE((CTVomptstpVec-CTVomptstpCVec).norm() < 1e-11);
*/
}
