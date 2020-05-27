#include "tests-def.h"
#include "integration.h"
#include "hubbchainops.h"
#include "utils.h"
#include "greens.h"

using namespace NEdyson;
TEST_CASE("Hubbard Chain Self-Energy Computations")
{
  int Nt=20, Ntau=30, Nsites=2, sig=-1;
  double Beta=9, dt=0.1;
  double MuChem = 1;
  ZMatrix h0(Nsites,Nsites);
  h0(0,0)=1;
  h0(0,1)=-3;
  h0(1,0)=-3;
  h0(1,1)=2;

  GREEN G = GREEN(Nt, Ntau, Nsites, sig);
  GREEN Sig = GREEN(Nt, Ntau, Nsites, sig);
  GREEN Sigcheck = GREEN(Nt, Ntau, Nsites, sig);
  CFUNC hmf(Nt, Nsites);
  CFUNC Ut(Nt, Nsites);
  CFUNC hmfcheck(Nt, Nsites);

  h5e::File file(std::string(TEST_DATA_DIR)+"/Gfree_const.h5");
  G.read_from_file(file,"");

  h5e::File file2(std::string(TEST_DATA_DIR)+"/hubb_SE_test_data.h5");
  Sigcheck.read_from_file(file2,"Sigma");
  hmfcheck.read_from_file(file2,"epsmf");
  Ut.read_from_file(file2,"Ut");

  for(int t=-1; t<=Nt; t++){
    Hubb::Ham_MF(t,G,Ut,h0,hmf);
    Hubb::Sigma_2B(t, G, Ut, Sig);
  }

  ZColVectorMap hmfVec = ZColVectorMap(hmf.ptr(-1), (Nt+2)*Nsites*Nsites);
  ZColVectorMap hmfcheckVec = ZColVectorMap(hmfcheck.ptr(-1), (Nt+2)*Nsites*Nsites);
  REQUIRE((hmfVec-hmfcheckVec).norm() < 1e-12);

  ZColVectorMap SMVec = ZColVectorMap(Sig.matptr(0), (Ntau+1)*Nsites*Nsites);
  ZColVectorMap SMcVec = ZColVectorMap(Sigcheck.matptr(0), (Ntau+1)*Nsites*Nsites);
  REQUIRE((SMVec-SMcVec).norm() < 1e-12);
  
  ZColVectorMap SRVec = ZColVectorMap(Sig.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
  ZColVectorMap SRcVec = ZColVectorMap(Sigcheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
  REQUIRE((SRVec-SRcVec).norm() < 1e-12);
  
  ZColVectorMap SLVec = ZColVectorMap(Sig.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
  ZColVectorMap SLcVec = ZColVectorMap(Sigcheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
  REQUIRE((SLVec-SLcVec).norm() < 1e-12);

  ZColVectorMap STVVec = ZColVectorMap(Sig.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
  ZColVectorMap STVcVec = ZColVectorMap(Sigcheck.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
  REQUIRE((STVVec-STVcVec).norm() < 1e-12);

}
