#include "dyson.h"
#include "utils.h"
#include "tests-def.h"
#include "set_up_g.h"

using namespace NEdyson;

TEST_CASE("Energy Integral"){
  // Parameters
  int Nt=50, Ntau=30, Nao=2;
  double beta=1.5, dt=0.02;

  // Objects
  GREEN S = GREEN(Nt, Ntau, Nao, -1);
  GREEN G = GREEN(Nt, Ntau, Nao, -1);
  DYSON Dyson(Nt, Ntau, Nao, 5, gfmol::Mode::GF2);

  set_up_G_herm(Dyson, G, S, dt, beta);

  // Values we are comparing to are exact solutions from functional forms
  // given in the set_up_G function
  double res = Dyson.energy_conv(0, S, G, beta, dt); 
  REQUIRE(std::abs(res + 0.28125) < 1e-11);

  res = Dyson.energy_conv(3, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.2927488689451264) < 1e-10);

  res = Dyson.energy_conv(5, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.30327566284415913) < 1e-10);

  res = Dyson.energy_conv(7, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.315908825747133) < 1e-10);

  res = Dyson.energy_conv(8, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.3229469844098787) < 1e-10);

  res = Dyson.energy_conv(10, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.33829858325203477) < 1e-10);  

  res = Dyson.energy_conv(20, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.42913160446733745) < 1e-10);
}


TEST_CASE("CTV1"){
  // Parameters
  int Nt=50, Ntau=30, Nao=2;
  double beta=1.5, dt=0.02;

  // Objects
  GREEN S = GREEN(Nt, Ntau, Nao, -1);
  GREEN G = GREEN(Nt, Ntau, Nao, -1);
  GREEN Scc = GREEN(Nt, Ntau, Nao, -1);
  GREEN Gcc = GREEN(Nt, Ntau, Nao, -1);
  DYSON Dyson(Nt, Ntau, Nao, 5, gfmol::Mode::GF2);

  ZRowVector res = ZRowVector((Ntau+1) * Nao * Nao);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  
  // Values we are comparing to are exact solutions from functional forms
  // given in the set_up_G function
  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 0, dt); 
  REQUIRE(std::abs(res.sum() - cplx(0.,0.)) < 1e-11);

  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 2, dt); 
  REQUIRE(std::abs(res.sum() - cplx(1.942981791711302 ,0.1785128961274855)) < 1e-11);

  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 4, dt); 
  REQUIRE(std::abs(res.sum() - cplx(3.8094748281637,0.7143691514656275)) < 1e-11);

  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 7, dt); 
  REQUIRE(std::abs(res.sum() - cplx(6.472103613463175,2.190432440083725)) < 1e-11);

  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 9, dt); 
  REQUIRE(std::abs(res.sum() - cplx(8.16023388825266,3.62521550316101)) < 1e-11);

  res.setZero();  
  Dyson.CTV1(res.data(), S, Scc, G, 20, dt);
  REQUIRE(std::abs(res.sum() - cplx(16.37804189025634,18.11528562371781)) < 1e-11);
}

TEST_CASE("Ctv_tstp"){
  // Parameters
  int Nt=50, Ntau=30, Nao=2;
  double beta=1.5, dt=0.02;

  // Objects
  GREEN S = GREEN(Nt, Ntau, Nao, -1);
  GREEN G = GREEN(Nt, Ntau, Nao, -1);
  GREEN Scc = GREEN(Nt, Ntau, Nao, -1);
  GREEN Gcc = GREEN(Nt, Ntau, Nao, -1);
  DYSON Dyson(Nt, Ntau, Nao, 5, gfmol::Mode::GF2);
 
  // Values we are comparing to are exact solutions from functional forms
  // given in the set_up_G function
  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(0, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(0,0), (Ntau+1) * Nao * Nao).sum() - cplx(59.99682203389833,0)) < 1e-11);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(2, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(2,0), (Ntau+1) * Nao * Nao).sum() - cplx(59.11788717375197, -0.06479860320481823)) < 1e-11);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(4, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(4,0), (Ntau+1) * Nao * Nao).sum() - cplx(58.64718604482833, 0.03667453794393438)) < 1e-10);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(7, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(7,0), (Ntau+1) * Nao * Nao).sum() - cplx(58.78029162902846, 0.5003153197597525)) < 1e-10);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(9, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(9,0), (Ntau+1) * Nao * Nao).sum() - cplx(59.47542129291611, 1.016421704133382)) < 1e-10);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(20, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(20,0), (Ntau+1) * Nao * Nao).sum() - cplx(73.36535591913849,6.770356696009172)) < 1e-10);
}

TEST_CASE("Cles2_tstp"){
  // Parameters
  int Nt=50, Ntau=30, Nao=2;
  double beta=1.5, dt=0.02;

  // Objects
  GREEN S = GREEN(Nt, Ntau, Nao, -1);
  GREEN G = GREEN(Nt, Ntau, Nao, -1);
  GREEN Scc = GREEN(Nt, Ntau, Nao, -1);
  GREEN Gcc = GREEN(Nt, Ntau, Nao, -1);
  DYSON Dyson(Nt, Ntau, Nao, 5, gfmol::Mode::GF2);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);

  ZRowVector res = ZRowVector((Nt+1) * Nao * Nao);
 
  // Values we are comparing to are exact solutions from functional forms
  // given in the set_up_G function
  res.setZero();
  Dyson.Cles2_tstp(G, Gcc, S, Scc, 0, dt, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (5+1) * Nao * Nao).sum() - cplx(0.,0.)) < 1e-11);
  
  res.setZero();
  Dyson.Cles2_tstp(G, Gcc, S, Scc, 2, dt, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (5+1) * Nao * Nao).sum() - cplx(-0.9779017217310684,-0.004556524033142945)) < 1e-11);

  res.setZero();
  Dyson.Cles2_tstp(G, Gcc, S, Scc, 5, dt, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (5+1) * Nao * Nao).sum() - cplx(-2.506405029994643,-0.028363125672494226)) < 1e-11);

  res.setZero();
  Dyson.Cles2_tstp(G, Gcc, S, Scc, 6, dt, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (6+1) * Nao * Nao).sum() - cplx(-3.534156697360735, -0.04705129775222273)) < 1e-11);

  res.setZero();
  Dyson.Cles2_tstp(G, Gcc, S, Scc, 9, dt, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (9+1) * Nao * Nao).sum() - cplx(-7.715189412549618, -0.14521215747129892)) < 1e-11);

  res.setZero();
  Dyson.Cles2_tstp(G, Gcc, S, Scc, 20, dt, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (20+1) * Nao * Nao).sum() - cplx(-37.224605902254595, -1.258015308782273)) < 1e-11);
}


TEST_CASE("Cles3_tstp"){
  // Parameters
  int Nt=50, Ntau=30, Nao=2;
  double beta=1.5, dt=0.02;

  // Objects
  GREEN S = GREEN(Nt, Ntau, Nao, -1);
  GREEN G = GREEN(Nt, Ntau, Nao, -1);
  GREEN Scc = GREEN(Nt, Ntau, Nao, -1);
  GREEN Gcc = GREEN(Nt, Ntau, Nao, -1);
  DYSON Dyson(Nt, Ntau, Nao, 5, gfmol::Mode::GF2);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);

  ZRowVector res = ZRowVector((Nt+1) * Nao * Nao);
 
  // Values we are comparing to are exact solutions from functional forms
  // given in the set_up_G function
  res.setZero();
  Dyson.Cles3_tstp(S, Scc, G, Gcc, 0, beta, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (5+1) * Nao * Nao).sum() - cplx(6.427012997367941, 5.348335104615002)) < 1e-11);
  
  res.setZero();
  Dyson.Cles3_tstp(S, Scc, G, Gcc, 2, beta, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (5+1) * Nao * Nao).sum() - cplx(6.64094640155254, 5.091254584720283)) < 1e-11);
  
  res.setZero();
  Dyson.Cles3_tstp(S, Scc, G, Gcc, 5, beta, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (5+1) * Nao * Nao).sum() - cplx(6.961846507829439, 4.705633804878206)) < 1e-11);
  
  res.setZero();
  Dyson.Cles3_tstp(S, Scc, G, Gcc, 6, beta, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (6+1) * Nao * Nao).sum() - cplx(8.42467270546969, 5.318615593818669)) < 1e-11);
  
  res.setZero();
  Dyson.Cles3_tstp(S, Scc, G, Gcc, 9, beta, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (9+1) * Nao * Nao).sum() - cplx(13.32282738165199, 6.804593007976855)) < 1e-11);
  
  res.setZero();
  Dyson.Cles3_tstp(S, Scc, G, Gcc, 20, beta, res.data());
  REQUIRE(std::abs(ZRowVectorMap(res.data(), (20+1) * Nao * Nao).sum() - cplx(37.65576268035543, 6.65193545259473)) < 1e-11);  
}


TEST_CASE("TTI Energy Integral"){
  // Parameters
  int Nt=50, Ntau=30, Nao=2;
  double beta=1.5, dt=0.02;

  // Objects
  TTI_GREEN S = TTI_GREEN(Nt, Ntau, Nao, -1);
  TTI_GREEN G = TTI_GREEN(Nt, Ntau, Nao, -1);
  DYSON Dyson(Nt, Ntau, Nao, 5, gfmol::Mode::GF2);

  set_up_G_herm(Dyson, G, S, dt, beta);

  // Values we are comparing to are exact solutions from functional forms
  // given in the set_up_G function
  double res = Dyson.energy_conv(0, S, G, beta, dt); 
  REQUIRE(std::abs(res + 0.2536938991661486) < 1e-11);

  res = Dyson.energy_conv(3, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.2912600061588824) < 1e-10);

  res = Dyson.energy_conv(5, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.3205005127226488) < 1e-10);

  res = Dyson.energy_conv(7, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.3530075191812132) < 1e-10);

  res = Dyson.energy_conv(8, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.3704494118784188) < 1e-10);

  res = Dyson.energy_conv(10, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.4076169710156956) < 1e-10);  

  res = Dyson.energy_conv(20, S, G, beta, dt);
  REQUIRE(std::abs(res + 0.6325216859891173) < 1e-10);
}

TEST_CASE("tti_CTV1"){
  // Parameters
  int Nt=50, Ntau=30, Nao=2;
  double beta=1.5, dt=0.02;

  // Objects
  TTI_GREEN S = TTI_GREEN(Nt, Ntau, Nao, -1);
  TTI_GREEN G = TTI_GREEN(Nt, Ntau, Nao, -1);
  TTI_GREEN Scc = TTI_GREEN(Nt, Ntau, Nao, -1);
  TTI_GREEN Gcc = TTI_GREEN(Nt, Ntau, Nao, -1);
  DYSON Dyson(Nt, Ntau, Nao, 5, gfmol::Mode::GF2);

  ZRowVector res = ZRowVector((Ntau+1) * Nao * Nao);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  
  // Values we are comparing to are exact solutions from functional forms
  // given in the set_up_G function
  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 0, dt); 
  REQUIRE(std::abs(res.sum() - cplx(0.,0.)) < 1e-11);

  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 2, dt); 
  REQUIRE(std::abs(res.sum() - cplx(2.021612956849114 ,0.0905539885442064)) < 1e-11);

  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 4, dt); 
  REQUIRE(std::abs(res.sum() - cplx(4.121417958322803,0.3673638697634925)) < 1e-11);

  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 7, dt); 
  REQUIRE(std::abs(res.sum() - cplx(7.415708359301342,1.148154340574662)) < 1e-11);

  res.setZero();
  Dyson.CTV1(res.data(), S, Scc, G, 9, dt); 
  REQUIRE(std::abs(res.sum() - cplx(9.70720634880647,1.922829105467475)) < 1e-11);

  res.setZero();  
  Dyson.CTV1(res.data(), S, Scc, G, 20, dt);
  REQUIRE(std::abs(res.sum() - cplx(23.65573319089267,10.1286305485927)) < 1e-11);
}

TEST_CASE("tti_Ctv_tstp"){
  // Parameters
  int Nt=50, Ntau=30, Nao=2;
  double beta=1.5, dt=0.02;

  // Objects
  TTI_GREEN S = TTI_GREEN(Nt, Ntau, Nao, -1);
  TTI_GREEN G = TTI_GREEN(Nt, Ntau, Nao, -1);
  TTI_GREEN Scc = TTI_GREEN(Nt, Ntau, Nao, -1);
  TTI_GREEN Gcc = TTI_GREEN(Nt, Ntau, Nao, -1);
  DYSON Dyson(Nt, Ntau, Nao, 5, gfmol::Mode::GF2);
 
  // Values we are comparing to are exact solutions from functional forms
  // given in the set_up_G function
  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(0, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(0,0), (Ntau+1) * Nao * Nao).sum() - cplx(59.99682203389833,0)) < 1e-11);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(2, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(2,0), (Ntau+1) * Nao * Nao).sum() - cplx(58.2314602691063, 0.968418291294562)) < 1e-9);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(4, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(4,0), (Ntau+1) * Nao * Nao).sum() - cplx(56.5390842722373, 1.984545915633391)) < 1e-9);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(7, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(7,0), (Ntau+1) * Nao * Nao).sum() - cplx(54.14804891993939, 3.596471489626535)) < 1e-8);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(9, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(9,0), (Ntau+1) * Nao * Nao).sum() - cplx(52.65820798605889, 4.729121596181367)) < 1e-8);

  set_up_G(Dyson, G, Gcc, S, Scc, dt, beta);
  Dyson.Ctv_tstp(20, G, G, Gcc, S, Scc, beta, dt);
  REQUIRE(std::abs(ZRowVectorMap(G.tvptr(20,0), (Ntau+1) * Nao * Nao).sum() - cplx(46.05010602956037,11.80946496691539)) < 1e-8);
}
