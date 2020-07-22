#include "spectral.h"
#include "greens.h"
#include "utils.h"
#include "tests-def.h"

using namespace NEdyson;

TEST_CASE("Spectral function calculations"){
  int Nt=20, Ntau=30, Nsites=2;
  h5e::File File(std::string(TEST_DATA_DIR) + "/Gfree_const.h5",h5e::File::ReadWrite);

  GREEN G(Nt, Ntau, Nsites,-1);
  G.read_from_file(File,"");

  SPECT ATest;
  SPECT AAct;

  h5e::File File2(std::string(TEST_DATA_DIR) + "/spectral.h5",h5e::File::ReadWrite);
  ATest.AfromG(G,100,10,0.1);
  AAct.read_from_file(File2,"Anext");

  DColVectorMap ATestVM = DColVectorMap(ATest.ptr(0,0),(ATest.nt()+1)*ATest.nw()*4);
  DColVectorMap AActVM = DColVectorMap(AAct.ptr(0,0),(AAct.nt()+1)*AAct.nw()*4);

  REQUIRE((ATestVM-AActVM).norm() < 10e-12);
}
/*
TEST_CASE("Expanded Spectral function calculations"){
  int Nt=20, Ntau=30, Nsites=2;
  h5e::File File(std::string(TEST_DATA_DIR) + "/Gfree_const.h5",h5e::File::ReadWrite);

  GREEN G(Nt, Ntau, Nsites,-1);
  G.read_from_file(File,"");

  SPECT ATest;
  SPECT AAct;
  h5e::File File2(std::string(TEST_DATA_DIR)+ "/spectral.h5",h5e::File::ReadWrite);
  AAct.read_from_file(File2,"Aext");

  std::ifstream in;
  in.open(std::string(TEST_DATA_DIR) + "/Gfree_GRExp.data");
  int nfit, ntp;
  in >> nfit >> ntp;
  cplx *extdata = new cplx[G.nt()*(G.nt()-nfit-ntp)*G.size1()*G.size1()];
  memset(extdata,0,G.nt()*(G.nt()-nfit-ntp)*G.size1()*G.size1());

  int tp, t, i, j, tpmax=G.nt()-nfit-ntp, tmax, size1=G.size1(), nt=G.nt(), es = size1*size1;
  double real, imag;
  for(tp=0;tp<tpmax;tp++){
    for(t=0;t<nt-tp;t++){
      for(i=0;i<size1;i++){
        for(j=0;j<size1;j++){
          if(!(in>>real>>imag)){
            std::cout<<t<<" "<<tp<<" "<<i<<" "<<j<<" "<<real<<" "<<imag<<std::endl;
            std::cerr<<"Error in reading extension data"<<std::endl;
          }
          extdata[tp*nt*es+t*es+i*size1+j] = cplx(real,imag);
        }
      }
    }
  }
  ATest.AfromG(G, AAct.nw(), AAct.wmax(), AAct.dt(), extdata, nfit, ntp);
  delete[] extdata;


  DColVectorMap ATestVM = DColVectorMap(ATest.ptr(0,0),(ATest.nt()+1)*ATest.nw()*4);
  DColVectorMap AActVM = DColVectorMap(AAct.ptr(0,0),(AAct.nt()+1)*AAct.nw()*4);

  REQUIRE((ATestVM-AActVM).norm() < 1e-12);
}*/
