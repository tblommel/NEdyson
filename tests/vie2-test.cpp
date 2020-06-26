#include "vie2.h"
#include "hubbchainops.h"
#include "utils.h"
#include "tests-def.h"

using namespace NEdyson;

TEST_CASE("VIE2 Matsubara Fourier Solution"){
  int Nt=20, Ntau=30, Nsites=2;
  double mu=1, beta=9, dt=0.1;

  GREEN G = GREEN(Nt, Ntau, Nsites, -1);
  CFUNC Ut = CFUNC(Nt, Nsites);
  INTEG I(5);

  GREEN SigmaCheck = GREEN(Nt, Ntau, Nsites, -1);
  GREEN PhiCheck = GREEN(Nt, Ntau, Nsites, 1);
  GREEN UxPhiCheck = GREEN(Nt, Ntau, Nsites, 1);
  GREEN PhixUCheck = GREEN(Nt, Ntau, Nsites, 1);
  GREEN TPPCheck = GREEN(Nt, Ntau, Nsites, 1);
  
  GREEN Sigma = GREEN(Nt, Ntau, Nsites, -1);
  GREEN Phi = GREEN(Nt, Ntau, Nsites, 1);
  GREEN UxPhi = GREEN(Nt, Ntau, Nsites, 1);
  GREEN PhixU = GREEN(Nt, Ntau, Nsites, 1);
  GREEN TPP = GREEN(Nt, Ntau, Nsites, 1);
  
  h5e::File File(std::string(TEST_DATA_DIR)+"/vie2_test_data.h5");
  G.read_from_file(File, "G");
  Ut.read_from_file(File, "Ut");
  SigmaCheck.read_from_file(File, "Sigma");
  PhiCheck.read_from_file(File, "Phi");
  UxPhiCheck.read_from_file(File, "UxPHI");
  PhixUCheck.read_from_file(File, "PHIxU");
  TPPCheck.read_from_file(File, "TPP");

  SECTION("Matsubara Fourier"){

    Hubb::GenTPP(-1, dt, beta, G, Phi, Ut, PhixU, UxPhi, TPP, I);
    Hubb::Sigma_TPP(-1, G, Ut, TPP, Sigma);

    ZColVectorMap SigMVec = ZColVectorMap(Sigma.matptr(0), (Ntau+1)*Nsites*Nsites);
    ZColVectorMap SigMcVec = ZColVectorMap(SigmaCheck.matptr(0), (Ntau+1)*Nsites*Nsites);
    REQUIRE((SigMVec-SigMcVec).norm() < 1e-12);
  
    ZColVectorMap TPPMVec = ZColVectorMap(TPP.matptr(0), (Ntau+1)*Nsites*Nsites);
    ZColVectorMap TPPMcVec = ZColVectorMap(TPPCheck.matptr(0), (Ntau+1)*Nsites*Nsites);
    REQUIRE((TPPMVec-TPPMcVec).norm() < 1e-12);
  
    ZColVectorMap PhiMVec = ZColVectorMap(Phi.matptr(0), (Ntau+1)*Nsites*Nsites);
    ZColVectorMap PhiMcVec = ZColVectorMap(PhiCheck.matptr(0), (Ntau+1)*Nsites*Nsites);
    REQUIRE((PhiMVec-PhiMcVec).norm() < 1e-12);
  
    ZColVectorMap pxuMVec = ZColVectorMap(PhixU.matptr(0), (Ntau+1)*Nsites*Nsites);
    ZColVectorMap pxuMcVec = ZColVectorMap(PhixUCheck.matptr(0), (Ntau+1)*Nsites*Nsites);
    REQUIRE((pxuMVec-pxuMcVec).norm() < 1e-12);
  
    ZColVectorMap uxpMVec = ZColVectorMap(UxPhi.matptr(0), (Ntau+1)*Nsites*Nsites);
    ZColVectorMap uxpMcVec = ZColVectorMap(UxPhiCheck.matptr(0), (Ntau+1)*Nsites*Nsites);
    REQUIRE((uxpMVec-uxpMcVec).norm() < 1e-12);
  
  }

  SECTION("Bootstrapping"){

    Hubb::GenTPP(-1, dt, beta, G, Phi, Ut, PhixU, UxPhi, TPP, I);
    Hubb::GenTPP(dt, beta, G, Phi, Ut, PhixU, UxPhi, TPP, I);

    for(int tstp=0; tstp<=5; tstp++){
      Hubb::Sigma_TPP(tstp, G, Ut, TPP, Sigma);
    }

    ZColVectorMap SigRVec = ZColVectorMap(Sigma.retptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap SigRcVec = ZColVectorMap(SigmaCheck.retptr(0,0), 21*Nsites*Nsites);
    REQUIRE((SigRVec-SigRcVec).norm() < 1e-12);
  
    ZColVectorMap TPPRVec = ZColVectorMap(TPP.retptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap TPPRcVec = ZColVectorMap(TPPCheck.retptr(0,0), 21*Nsites*Nsites);
    REQUIRE((TPPRVec-TPPRcVec).norm() < 1e-12);
  
    ZColVectorMap PhiRVec = ZColVectorMap(Phi.retptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap PhiRcVec = ZColVectorMap(PhiCheck.retptr(0,0), 21*Nsites*Nsites);
    REQUIRE((PhiRVec-PhiRcVec).norm() < 1e-12);
  
    ZColVectorMap pxuRVec = ZColVectorMap(PhixU.retptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap pxuRcVec = ZColVectorMap(PhixUCheck.retptr(0,0), 21*Nsites*Nsites);
    REQUIRE((pxuRVec-pxuRcVec).norm() < 1e-12);
  
    ZColVectorMap uxpRVec = ZColVectorMap(UxPhi.retptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap uxpRcVec = ZColVectorMap(UxPhiCheck.retptr(0,0), 21*Nsites*Nsites);
    REQUIRE((uxpRVec-uxpRcVec).norm() < 1e-12);



    ZColVectorMap SigLVec = ZColVectorMap(Sigma.lesptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap SigLcVec = ZColVectorMap(SigmaCheck.lesptr(0,0), 21*Nsites*Nsites);
    REQUIRE((SigLVec-SigLcVec).norm() < 1e-12);
  
    ZColVectorMap TPPLVec = ZColVectorMap(TPP.lesptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap TPPLcVec = ZColVectorMap(TPPCheck.lesptr(0,0), 21*Nsites*Nsites);
    REQUIRE((TPPLVec-TPPLcVec).norm() < 1e-12);
  
    ZColVectorMap PhiLVec = ZColVectorMap(Phi.lesptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap PhiLcVec = ZColVectorMap(PhiCheck.lesptr(0,0), 21*Nsites*Nsites);
    REQUIRE((PhiLVec-PhiLcVec).norm() < 1e-12);
  
    ZColVectorMap pxuLVec = ZColVectorMap(PhixU.lesptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap pxuLcVec = ZColVectorMap(PhixUCheck.lesptr(0,0), 21*Nsites*Nsites);
    REQUIRE((pxuLVec-pxuLcVec).norm() < 1e-12);
  
    ZColVectorMap uxpLVec = ZColVectorMap(UxPhi.lesptr(0,0), 21*Nsites*Nsites);
    ZColVectorMap uxpLcVec = ZColVectorMap(UxPhiCheck.lesptr(0,0), 21*Nsites*Nsites);
    REQUIRE((uxpLVec-uxpLcVec).norm() < 1e-12);



    ZColVectorMap SigTVVec = ZColVectorMap(Sigma.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap SigTVcVec = ZColVectorMap(SigmaCheck.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    REQUIRE((SigTVVec-SigTVcVec).norm() < 1e-12);
  
    ZColVectorMap TPPTVVec = ZColVectorMap(TPP.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap TPPTVcVec = ZColVectorMap(TPPCheck.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    REQUIRE((TPPTVVec-TPPTVcVec).norm() < 1e-12);
  
    ZColVectorMap PhiTVVec = ZColVectorMap(Phi.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap PhiTVcVec = ZColVectorMap(PhiCheck.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    REQUIRE((PhiTVVec-PhiTVcVec).norm() < 1e-12);
  
    ZColVectorMap pxuTVVec = ZColVectorMap(PhixU.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap pxuTVcVec = ZColVectorMap(PhixUCheck.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    REQUIRE((pxuTVVec-pxuTVcVec).norm() < 1e-12);
  
    ZColVectorMap uxpTVVec = ZColVectorMap(UxPhi.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap uxpTVcVec = ZColVectorMap(UxPhiCheck.tvptr(0,0), 6*(Ntau+1)*Nsites*Nsites);
    REQUIRE((uxpTVVec-uxpTVcVec).norm() < 1e-12);

  }
  
  SECTION("TimeStepping"){

    Hubb::GenTPP(-1, dt, beta, G, Phi, Ut, PhixU, UxPhi, TPP, I);
    Hubb::GenTPP(dt, beta, G, Phi, Ut, PhixU, UxPhi, TPP, I);
    for(int t=6;t<=Nt;t++){
      Hubb::GenTPP(t, dt, beta, G, Phi, Ut, PhixU, UxPhi, TPP, I);
    }

    for(int tstp=-1; tstp<=Nt; tstp++){
      Hubb::Sigma_TPP(tstp, G, Ut, TPP, Sigma);
    }

    ZColVectorMap SigRVec = ZColVectorMap(Sigma.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap SigRcVec = ZColVectorMap(SigmaCheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((SigRVec-SigRcVec).norm() < 1e-12);
  
    ZColVectorMap TPPRVec = ZColVectorMap(TPP.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap TPPRcVec = ZColVectorMap(TPPCheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((TPPRVec-TPPRcVec).norm() < 1e-12);
  
    ZColVectorMap PhiRVec = ZColVectorMap(Phi.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap PhiRcVec = ZColVectorMap(PhiCheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((PhiRVec-PhiRcVec).norm() < 1e-12);
  
    ZColVectorMap pxuRVec = ZColVectorMap(PhixU.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap pxuRcVec = ZColVectorMap(PhixUCheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((pxuRVec-pxuRcVec).norm() < 1e-12);
  
    ZColVectorMap uxpRVec = ZColVectorMap(UxPhi.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap uxpRcVec = ZColVectorMap(UxPhiCheck.retptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((uxpRVec-uxpRcVec).norm() < 1e-12);



    ZColVectorMap SigLVec = ZColVectorMap(Sigma.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap SigLcVec = ZColVectorMap(SigmaCheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((SigLVec-SigLcVec).norm() < 1e-12);
  
    ZColVectorMap TPPLVec = ZColVectorMap(TPP.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap TPPLcVec = ZColVectorMap(TPPCheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((TPPLVec-TPPLcVec).norm() < 1e-12);
  
    ZColVectorMap PhiLVec = ZColVectorMap(Phi.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap PhiLcVec = ZColVectorMap(PhiCheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((PhiLVec-PhiLcVec).norm() < 1e-12);
  
    ZColVectorMap pxuLVec = ZColVectorMap(PhixU.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap pxuLcVec = ZColVectorMap(PhixUCheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((pxuLVec-pxuLcVec).norm() < 1e-12);
  
    ZColVectorMap uxpLVec = ZColVectorMap(UxPhi.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    ZColVectorMap uxpLcVec = ZColVectorMap(UxPhiCheck.lesptr(0,0), (Nt+1)*(Nt+2)/2*Nsites*Nsites);
    REQUIRE((uxpLVec-uxpLcVec).norm() < 1e-12);



    ZColVectorMap SigTVVec = ZColVectorMap(Sigma.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap SigTVcVec = ZColVectorMap(SigmaCheck.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    REQUIRE((SigTVVec-SigTVcVec).norm() < 1e-12);
  
    ZColVectorMap TPPTVVec = ZColVectorMap(TPP.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap TPPTVcVec = ZColVectorMap(TPPCheck.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    REQUIRE((TPPTVVec-TPPTVcVec).norm() < 1e-12);
  
    ZColVectorMap PhiTVVec = ZColVectorMap(Phi.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap PhiTVcVec = ZColVectorMap(PhiCheck.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    REQUIRE((PhiTVVec-PhiTVcVec).norm() < 1e-12);
  
    ZColVectorMap pxuTVVec = ZColVectorMap(PhixU.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap pxuTVcVec = ZColVectorMap(PhixUCheck.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    REQUIRE((pxuTVVec-pxuTVcVec).norm() < 1e-12);
  
    ZColVectorMap uxpTVVec = ZColVectorMap(UxPhi.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    ZColVectorMap uxpTVcVec = ZColVectorMap(UxPhiCheck.tvptr(0,0), Nt*(Ntau+1)*Nsites*Nsites);
    REQUIRE((uxpTVVec-uxpTVcVec).norm() < 1e-12);

  }
}
