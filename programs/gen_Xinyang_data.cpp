#include "dyson.h"
#include "greens.h"
#include "utils.h"

int main(const int argc, char *const *const argv)
{
  using namespace NEdyson;
  int Nt = 1, Ntau = strtol(argv[1], NULL ,10), Nao = 2, sig = -1;
  double beta = 2, dtau = beta/Ntau;
  GREEN G(Nt, Ntau, Nao, sig);
  GREEN Sigma(Nt, Ntau, Nao, sig);
  DYSON Dyson(Nt, Ntau, Nao, 5);

  ZMatrix SigBase(2,2);
  SigBase(0,0) = 1;
  SigBase(0,1) = 2;
  SigBase(1,0) = 3;
  SigBase(1,1) = 4;
  ZTensor<3> Conv(Ntau+1, Nao, Nao);
  ZTensor<2> tmp(Nao, Nao);

  for(int tau = 0; tau <= Ntau; tau++) {
    double t = tau*dtau;
    ZMatrixMap(Sigma.tvptr(0,tau), Nao, Nao) = t * SigBase;
    ZMatrixMap(Sigma.tvptr(1,tau), Nao, Nao) = t * SigBase;

    t = (Ntau-tau)*dtau;
    ZMatrix GMat(2,2);
    GMat(0,0) = t;
    GMat(0,1) = t*t;
    GMat(1,0) = t*t*t;
    GMat(1,1) = t*t*t*t;
    GMat *= cplx(0,1);
    ZMatrixMap(G.tvptr(0,tau), Nao, Nao) = GMat.adjoint();
    ZMatrixMap(G.tvptr(1,tau), Nao, Nao) = GMat.adjoint();

    t = tau*dtau;
    GMat(0,0) = t;
    GMat(0,1) = t*t;
    GMat(1,0) = t*t*t;
    GMat(1,1) = t*t*t*t;
    ZMatrixMap(G.matptr(tau), Nao, Nao) = GMat;
  }

  ZMatrix res(2,2);
  
  Dyson.Cles3_tstp(0, 0, Sigma, Sigma, G, G, 0, beta, res.data());
  
  for(int m=0; m<=Ntau; m++){
    Dyson.CTV2(Sigma, G, 0, m, beta, Conv.data() + m*Nao*Nao);
    Dyson.CTV3(Sigma, G, 0, m, beta, tmp.data());
    ZMatrixMap(Conv.data() + m*Nao*Nao, Nao, Nao) += ZMatrixMap(tmp.data(), Nao, Nao);
  }


  for(int i=0; i<Nao; i++)
    for(int j=0; j<Nao; j++)
      std::cout<<res(i,j)<<std::endl;

  for(int t=0; t<=Ntau; t++){
    std::cout<<t*dtau<<std::endl;
    for(int i=0; i<Nao; i++)
      for(int j=0; j<Nao; j++)
        std::cout<<Conv(t,i,j)<<std::endl;
    std::cout<<std::endl;
  }
}
