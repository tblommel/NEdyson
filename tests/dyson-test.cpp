#include "dyson.h"
#include "utils.h"
#include "tests-def.h"

using namespace NEdyson;

TEST_CASE("Dyson Tests"){
  int Nt=20, Ntau=30, Nsites=2;
  double mu=1, beta=9, dt=0.1;

  GREEN Sigma = GREEN(Nt, Ntau, Nsites, -1);
  DYSON Dyson = DYSON(Nt, Ntau, Nsites, 5);
  ZTensor<1> hmf((Nt+2)*Nsites*Nsites);

  
  h5e::File File(std::string(TEST_DATA_DIR)+"/dyson_test_data.h5");
  Sigma.read_from_file(File, "Sigma");
  h5::DataSet dataset = File.getDataSet("epsmf/ft");
  size_t len = dataset.getElementCount();
  dataset.read(hmf.data());
  File.flush();



  SECTION("Bootstrapping"){

    GREEN Gtest(Nt,Ntau, Nsites, -1);
    Gtest.read_from_file(File, "Free");
    GREEN Gcheck(Nt, Ntau, Nsites, -1);
    Gcheck.read_from_file(File, "Start");
    Dyson.dyson_start(Gtest,Sigma,hmf.data()+4,mu,beta,dt);

    ZColVectorMap GMVec = ZColVectorMap(Gtest.matptr(0), (Ntau+1)*Nsites*Nsites);
    ZColVectorMap GMcVec = ZColVectorMap(Gcheck.matptr(0), (Ntau+1)*Nsites*Nsites);
    REQUIRE((GMVec-GMcVec).norm() < 1e-12);
  
    ZColVectorMap GRVec = ZColVectorMap(Gtest.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap GRcVec = ZColVectorMap(Gcheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((GRVec-GRcVec).norm() < 1e-12);

    ZColVectorMap GLVec = ZColVectorMap(Gtest.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap GLcVec = ZColVectorMap(Gcheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((GLVec-GLcVec).norm() < 1e-12);
  
    ZColVectorMap GTVVec = ZColVectorMap(Gtest.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap GTVcVec = ZColVectorMap(Gcheck.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
    REQUIRE((GTVVec-GTVcVec).norm() < 1e-12);
  }

  SECTION("Timestepping"){

    GREEN Gtest(Nt,Ntau, Nsites, -1);
    Gtest.read_from_file(File, "Free");
    GREEN Gcheck(Nt, Ntau, Nsites, -1);
    Gcheck.read_from_file(File, "Timesteps");
    for(int tstp=6;tstp<=Nt;tstp++){
      Dyson.dyson_step(tstp,Gtest,Sigma,hmf.data()+4,mu,beta,dt);
    }

    ZColVectorMap GMVec = ZColVectorMap(Gtest.matptr(0), (Ntau+1)*Nsites*Nsites);
    ZColVectorMap GMcVec = ZColVectorMap(Gcheck.matptr(0), (Ntau+1)*Nsites*Nsites);
    REQUIRE((GMVec-GMcVec).norm() < 1e-12);
  
    ZColVectorMap GRVec = ZColVectorMap(Gtest.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap GRcVec = ZColVectorMap(Gcheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((GRVec-GRcVec).norm() < 1e-12);

    ZColVectorMap GLVec = ZColVectorMap(Gtest.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap GLcVec = ZColVectorMap(Gcheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((GLVec-GLcVec).norm() < 1e-12);
  
    ZColVectorMap GTVVec = ZColVectorMap(Gtest.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap GTVcVec = ZColVectorMap(Gcheck.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
    REQUIRE((GTVVec-GTVcVec).norm() < 1e-12);
  }

  SECTION("Extrapolation"){
    
    GREEN Gtest(Nt,Ntau, Nsites, -1);
    GREEN Gcheck(Nt, Ntau, Nsites, -1);
    Gcheck.read_from_file(File, "Extrap");
    Gtest.read_from_file(File, "Start");
    for(int tstp=6;tstp<=Nt;tstp++){
      Dyson.Extrapolate(tstp,Gtest);
    }

    ZColVectorMap GMVec = ZColVectorMap(Gtest.matptr(0), (Ntau+1)*Nsites*Nsites);
    ZColVectorMap GMcVec = ZColVectorMap(Gcheck.matptr(0), (Ntau+1)*Nsites*Nsites);
    REQUIRE((GMVec-GMcVec).norm() < 1e-6);
  
    ZColVectorMap GRVec = ZColVectorMap(Gtest.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap GRcVec = ZColVectorMap(Gcheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((GRVec-GRcVec).norm() < 1e-6);

    ZColVectorMap GLVec = ZColVectorMap(Gtest.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap GLcVec = ZColVectorMap(Gcheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((GLVec-GLcVec).norm() < 1e-6);
  
    ZColVectorMap GTVVec = ZColVectorMap(Gtest.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap GTVcVec = ZColVectorMap(Gcheck.tvptr(0,0), (Nt+1)*(Ntau+1)*Nsites*Nsites);
    REQUIRE((GTVVec-GTVcVec).norm() < 1e-6);
  }
}
