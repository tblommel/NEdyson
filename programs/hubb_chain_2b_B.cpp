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
  double HoppingT, HubbardU, MuChem, Beta, dt, dtau, W0, MatsMaxErr, BootMaxErr, wmax, B;

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
  find_param(argv[1],"__k=",k);
  find_param(argv[1],"__B=",B);
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
  GREEN GU = GREEN(Nt, Ntau, Nsites, -1);
  GREEN SigmaU = GREEN(Nt, Ntau, Nsites, -1);
  GREEN GD = GREEN(Nt, Ntau, Nsites, -1);
  GREEN SigmaD = GREEN(Nt, Ntau, Nsites, -1);

  CFUNC hmfU = CFUNC(Nt, Nsites);
  CFUNC hmfD = CFUNC(Nt, Nsites);
  cdmatrix h0U(Nsites,Nsites), DensMU(Nsites,Nsites), h0D(Nsites,Nsites), DensMD(Nsites,Nsites);
  CFUNC Ut = CFUNC(Nt, Nsites);
  
  SPECT A = SPECT();

  // Get free green's func ==============================================================
  h0U.setZero();
  h0D.setZero();
  for(i=0;i<Nsites-1;i++){
    h0U(i,i+1)=-HoppingT;
    h0U(i+1,i)=-HoppingT;
    h0D(i,i+1)=-HoppingT;
    h0D(i+1,i)=-HoppingT;
  }
  for(i=0;i<Nsites;i++){
    h0U(i,i) = ed(i)-B;
    h0D(i,i) = ed(i)+B;
  }
  Ut.set_constant(HubbardU*cdmatrix::Identity(Nsites,Nsites));
 
  NEdyson::G0_from_h0(GU,MuChem,h0U,Beta,dt);
  NEdyson::G0_from_h0(GD,MuChem,h0D,Beta,dt);
  
  std::string str;
  str = argv[1];
  str+="/GU_free";
  
  GU.print_to_file_mat(str,dt,dtau,16);
  str = argv[1];
  str+="/GD_free";
  GD.print_to_file_mat(str,dt,dtau,16);

  GU.get_dm(-1,DensMU);
  GD.get_dm(-1,DensMD);

  std::cout.precision(17);
  double npart = DensMU.trace().real()+DensMD.trace().real();
  std::cout << "number of particles = " << npart << std::endl;
  MuChem += HubbardU*npart/Nsites/2;

  // Matsubara Self Consistency =============================================================
  {
    start = std::chrono::system_clock::now();
    bool converged = false;
    tstp = -1;
    for(int iter=0;iter<MatsMaxIter;iter++){
      double err;
      // Update Mean field and Self Energy
      Hubb::Ham_MF(tstp,GU,GD,Ut,h0U,h0D,hmfU,hmfD);
      Hubb::Sigma_2B(tstp,GU,GD,Ut,SigmaU,SigmaD);
  
      // Solve Dyson Equation
      err = NEdyson::mat_fourier(GU,SigmaU,MuChem,hmfU.ptr(-1),Beta);
      err += NEdyson::mat_fourier(GD,SigmaD,MuChem,hmfD.ptr(-1),Beta);

      GU.get_dm(-1,DensMU);
      GD.get_dm(-1,DensMD);
      npart = DensMU.trace().real() + DensMD.trace().real();
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


  str = argv[1];
  str+="/GU";
  
  GU.print_to_file_mat(str,dt,dtau,16);
  str = argv[1];
  str+="/GD";
  
  GD.print_to_file_mat(str,dt,dtau,16);


  // Bootstrapping ==========================================================================
  {
    start = std::chrono::system_clock::now();
    bool bootstrap_converged=false;
    
    // to represent the quench, the free Hamiltonian is updated
    h0D(RampSite-1,RampSite-1) += W0;
    h0U(RampSite-1,RampSite-1) += W0;

    for(int iter=0; iter<=BootMaxIter; iter++){
      double err=0;
      // Update mean field
      for(tstp=0; tstp<=SolverOrder; tstp++){
        Hubb::Ham_MF(tstp,GU,GD,Ut,h0U,h0D,hmfU,hmfD);
      }

      // Update self-energy
      for(tstp=0; tstp<=SolverOrder; tstp++){
        Hubb::Sigma_2B(tstp,GU,GD,Ut,SigmaU,SigmaD);
      }

      // Solve Dyson Eqn
      err = dyson_start(I,GU,SigmaU,hmfU,MuChem,Beta,dt);
      err += dyson_start(I,GD,SigmaD,hmfD,MuChem,Beta,dt);

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

  for(int tstp=k+1;tstp<=Nt;tstp++){
    // Extrapolate
    NEdyson::Extrapolate(I,GD,tstp);
    NEdyson::Extrapolate(I,GU,tstp);
    
    for(int iter=0;iter<StepIter;iter++){
      // Update hmf
      Hubb::Ham_MF(tstp,GU,GD,Ut,h0U,h0D,hmfU,hmfD);
      
      // Update Self Energy
      Hubb::Sigma_2B(tstp,GU,GD,Ut,SigmaU,SigmaD);

      // Solve Dyson Eqn
      dyson_step(tstp,I,GU,SigmaU,hmfU,MuChem,Beta,dt);
      dyson_step(tstp,I,GD,SigmaD,hmfD,MuChem,Beta,dt);
    }
    std::cout<<tstp<<std::endl;
  }

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << "Time [Propagation] = " << elapsed_seconds.count() << "s\n\n";

  str = argv[1];
  str+="/GU"; 
  GU.print_to_file_ret(str,dt,dtau,16);
  str = argv[1];
  str+="/GD"; 
  GD.print_to_file_ret(str,dt,dtau,16);

  A.AfromG(GU,nw,wmax,dt);
  str = argv[1];
  str+="/AU";
  A.print_to_file(str);
  A.AfromG(GD,nw,wmax,dt);
  str = argv[1];
  str+="/AD";
  A.print_to_file(str);
}
