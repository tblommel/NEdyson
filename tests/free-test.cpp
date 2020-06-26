#include "NEdyson.h"
#include "tests-def.h"
#include "integration.h"
#include "utils.h"

using namespace NEdyson;
TEST_CASE("Const Free Function Calculation")
{
  int Nt=20, Ntau=30, Nsites=2, sig=-1;
  double Beta=9, dt=0.1;
  GREEN G = GREEN(Nt, Ntau, Nsites, sig);
  GREEN Gcheck = GREEN(Nt, Ntau, Nsites, sig);
  double MuChem = 1;
  ZMatrix h0(Nsites,Nsites);
  h0(0,0)=1;
  h0(0,1)=-3;
  h0(1,0)=-3;
  h0(1,1)=2;

  NEdyson::G0_from_h0(G,MuChem,h0,Beta,dt);

  h5e::File file(std::string(TEST_DATA_DIR)+"/Gfree_const.h5", h5e::File::ReadWrite);
  Gcheck.read_from_file(file,"");

  ZColVector GR = ZColVectorMap(G.retptr(0,0),(Nt+1)*(Nt+2)/2*Nsites*Nsites);
  ZColVector GTV = ZColVectorMap(G.tvptr(0,0),(Nt+1)*(Ntau+1)*Nsites*Nsites);
  ZColVector GM = ZColVectorMap(G.matptr(0),(Ntau+1)*Nsites*Nsites);
  ZColVector GL = ZColVectorMap(G.lesptr(0,0),(Nt+1)*(Nt+2)/2*Nsites*Nsites);

  ZColVector GRc = ZColVectorMap(Gcheck.retptr(0,0),(Nt+1)*(Nt+2)/2*Nsites*Nsites);
  ZColVector GTVc = ZColVectorMap(Gcheck.tvptr(0,0),(Nt+1)*(Ntau+1)*Nsites*Nsites);
  ZColVector GMc = ZColVectorMap(Gcheck.matptr(0),(Ntau+1)*Nsites*Nsites);
  ZColVector GLc = ZColVectorMap(Gcheck.lesptr(0,0),(Nt+1)*(Nt+2)/2*Nsites*Nsites);
  
  REQUIRE((GR-GRc).norm() < 1e-10);
  REQUIRE((GL-GLc).norm() < 1e-10);
  REQUIRE((GTV-GTVc).norm() < 1e-10);
  REQUIRE((GM-GMc).norm() < 1e-10);
}

TEST_CASE("Non-Const Free Function Calculation")
{
  int Nt=20, Ntau=30, Nsites=2, sig=-1;
  double Beta=9, dt=0.1;
  GREEN G = GREEN(Nt, Ntau, Nsites, sig);
  GREEN Gcheck = GREEN(Nt, Ntau, Nsites, sig);
  CFUNC H = CFUNC(Nt,Nsites);
  double MuChem = 1;

  h5e::File file(std::string(TEST_DATA_DIR)+"/Gfree_nconst.h5");
  Gcheck.read_from_file(file,"");
  h5::DataSet dataset = file.getDataSet("/hmft");
  dataset.read(H.ptr(0));
  h5::DataSet dataset2 = file.getDataSet("/h0");
  dataset2.read(H.ptr(-1));
  
  INTEG I = INTEG(5);
  G0_from_h0(G,I,MuChem,H,Beta,dt);

  ZColVector GR = ZColVectorMap(G.retptr(0,0),(Nt+1)*(Nt+2)/2*Nsites*Nsites);
  ZColVector GTV = ZColVectorMap(G.tvptr(0,0),(Nt+1)*(Ntau+1)*Nsites*Nsites);
  ZColVector GM = ZColVectorMap(G.matptr(0),(Ntau+1)*Nsites*Nsites);
  ZColVector GL = ZColVectorMap(G.lesptr(0,0),(Nt+1)*(Nt+2)/2*Nsites*Nsites);

  ZColVector GRc = ZColVectorMap(Gcheck.retptr(0,0),(Nt+1)*(Nt+2)/2*Nsites*Nsites);
  ZColVector GTVc = ZColVectorMap(Gcheck.tvptr(0,0),(Nt+1)*(Ntau+1)*Nsites*Nsites);
  ZColVector GMc = ZColVectorMap(Gcheck.matptr(0),(Ntau+1)*Nsites*Nsites);
  ZColVector GLc = ZColVectorMap(Gcheck.lesptr(0,0),(Nt+1)*(Nt+2)/2*Nsites*Nsites);
  
  REQUIRE((GR-GRc).norm() < 1e-10);
  REQUIRE((GL-GLc).norm() < 1e-10);
  REQUIRE((GTV-GTVc).norm() < 1e-10);
  REQUIRE((GM-GMc).norm() < 1e-10);
}
