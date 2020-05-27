#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <vector>
#include <Eigen/Eigen>
#include <chrono>

#include "NEdyson.h"

using namespace NEdyson;


int main(int argc, char *argv[]){

  if(argc!=2) throw std::invalid_argument("please provide input file");
  
  // Parameters ========================================================================
  int Nsites, Ntau, MatsMaxIter, oversamp, Nt, RampSite, k, SolverOrder, BootMaxIter, StepIter,nw;
  double HoppingT, HubbardU, MuChem, Beta, dt, dtau, W0, MatsMaxErr, BootMaxErr, wmax;

  std::chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot;
  std::chrono::duration<double> runtime;
  int tstp,i,j;

  start_tot = std::chrono::system_clock::now();
  // Read in ===========================================================================
  find_param(argv[1],"__Nsites=",Nsites);
  find_param(argv[1],"__HoppingT=",HoppingT);
  find_param(argv[1],"__HubbardU=",HubbardU);
  find_param(argv[1],"__MuChem=",MuChem);
  find_param(argv[1],"__Beta=",Beta);
  find_param(argv[1],"__Nt=",Nt);
  find_param(argv[1],"__dt=",dt);
  find_param(argv[1],"__Ntau=",Ntau);
  find_param(argv[1],"__MatsMaxIter=",MatsMaxIter);
  find_param(argv[1],"__oversamp=",oversamp);
  find_param(argv[1],"__RampSite=",RampSite);
  find_param(argv[1],"__RampW0=",W0);
  find_param(argv[1],"__k",k);
  dtau=Beta/Ntau;
  find_param(argv[1],"__MatsMaxErr",MatsMaxErr);
  find_param(argv[1],"__BootMaxIter",BootMaxIter);
  find_param(argv[1],"__BootMaxErr",BootMaxErr);
  SolverOrder=k;
  find_param(argv[1],"__StepIter",StepIter);
  Eigen::VectorXd ed(Nsites);
  find_param(argv[1],"__ed=",ed);
  find_param(argv[1],"__nw",nw);
  find_param(argv[1],"__wmax",wmax);

  // Define Necessary Structures =======================================================
  INTEG I = INTEG(k);
  GREEN G = GREEN(Nt, Ntau, Nsites, -1);
  GREEN Sigma = GREEN(Nt, Ntau, Nsites, -1);

  CFUNC hmf = CFUNC(Nt, Nsites);
  CFUNC h0f = CFUNC(Nt,Nsites);
  ZMatrix h0(Nsites,Nsites), DensM(Nsites,Nsites);
  CFUNC Ut = CFUNC(Nt, Nsites);
  
  SPECT A = SPECT();

  // Get free green's func ==============================================================
  h0.setZero();
  for(i=0;i<Nsites-1;i++){
    h0(i,i+1)=-HoppingT;
    h0(i+1,i)=-HoppingT;
  }
  for(i=0;i<Nsites;i++){
    h0(i,i) = ed(i);
  }
  h0f.set_constant(h0);
  Ut.set_constant(HubbardU*ZMatrix::Identity(Nsites,Nsites));
  
  NEdyson::G0_from_h0(G,MuChem,h0,Beta,dt);

  std::string str;
  str = argv[1];
  str+="/Gfree_const.h5";
  h5e::File h5Gfree(str, h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);
  G.print_to_file(h5Gfree,"");
  
  str = argv[1];
  str+="/G";  
  G.print_to_file_mat(str,dt,dtau,16);

  G.get_dm(-1,DensM);
  std::cout.precision(17);
  double npart = DensM.trace().real();
  std::cout << "number of particles = " << npart << std::endl;
  MuChem += HubbardU*npart/Nsites;

  // Matsubara Self Consistency =============================================================
  {
    start = std::chrono::system_clock::now();
    bool converged = false;
    tstp = -1;
    for(int iter=0;iter<MatsMaxIter;iter++){
      double err;
      // Update Mean field and Self Energy
      Hubb::Ham_MF(tstp,G,Ut,h0,hmf);
      Hubb::Sigma_2B(tstp,G,Ut,Sigma);
  
      // Solve Dyson Equation
      err = NEdyson::mat_fourier(G,Sigma,MuChem,hmf.ptr(-1),Beta);

      G.get_dm(-1,DensM);
      npart = DensM.trace().real();
      std::cout<<"iteration: "<<iter<<" | N =  "<<npart<<" |  Error = "<<err<<std::endl;
      if(err<MatsMaxErr){
        converged=true;
        break;
      }
    }

    if(!converged){
      std::cout<<std::endl;
      std::cout<<" Matsubara Iteration not converging! Quitting program. " <<std::endl;
      return 0;
    }
  
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Time [equilibrium calculation] = " << elapsed_seconds.count() << "s\n\n";  
  } // End Matsubara Self Consistency Loop




  // Bootstrapping ==========================================================================
  {
    start = std::chrono::system_clock::now();
    bool bootstrap_converged=false;
    
    // to represent the quench, the free Hamiltonian is updated
    h0(RampSite-1,RampSite-1) += W0;

    for(int iter=0; iter<=BootMaxIter; iter++){
      double err=0;
      // Update mean field
      for(tstp=0; tstp<=SolverOrder; tstp++){
        Hubb::Ham_MF(tstp, G, Ut, h0, hmf);
      }

      // Update self-energy
      for(tstp=0; tstp<=SolverOrder; tstp++){
        Hubb::Sigma_2B(tstp, G, Ut, Sigma);
      }

      // Solve Dyson Eqn
      err = dyson_start(I,G,Sigma,hmf,MuChem,Beta,dt);

      std::cout<<"bootstrapping iteration : "<<iter<<" |  Error = "<<err<<std::endl;
      if(err<BootMaxErr){
        bootstrap_converged=true;
        break;
      }
    }
    if(!bootstrap_converged){
      std::cout<<std::endl;
      std::cout<<" Bootstrapping did not converge. Program Ending. "<<std::endl;
      return 0;
    }
  } // End Bootstrapping

  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "Time [bootstrapping] = " << elapsed_seconds.count() << "s\n\n";


  // Timestepping ===========================================================================
  start = std::chrono::system_clock::now();
  GREEN Gtemp(Nt,Ntau,Nsites,-1);

  for(int tstp=k+1;tstp<=Nt;tstp++){
    // Extrapolate
    NEdyson::Extrapolate(I,G,tstp);
    
    for(int iter=0;iter<StepIter;iter++){
      Gtemp.set_tstp(tstp,G);
      // Update hmf
      Hubb::Ham_MF(tstp,G,Ut,h0,hmf);
      
      // Update Self Energy
      Hubb::Sigma_2B(tstp,G,Ut,Sigma);

      // Solve Dyson Eqn
      dyson_step(tstp,I,G,Sigma,hmf,MuChem,Beta,dt);
      double err = distance_norm2(tstp,G,Gtemp);
      std::cout<<tstp<<" "<<iter<<err<<" "<<std::endl;
    }
  }

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << "Time [Propagation] = " << elapsed_seconds.count() << "s\n\n";

  str = argv[1];
  str+="/G";
  
  G.print_to_file(str,dt,dtau,16);

  A.AfromG(G,nw,wmax,dt);
  str = argv[1];
  str+="/A";
  A.print_to_file(str);
}
