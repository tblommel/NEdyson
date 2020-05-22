#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <vector>
#include <Eigen/Eigen>
#include <chrono>

#include <NEdyson.h>

using namespace NEdyson;

int main(int argc, char *argv[]){
  if(argc!=2) throw std::invalid_argument("please provide input file");
  
  // Parameters ========================================================================
  int Nsites, Ntau, MatsMaxIter, oversamp, Nt, RampSite, k, SolverOrder, BootMaxIter, StepIter;
  double HoppingT, HubbardU, MuChem, Beta, dt, dtau, W0, MatsMaxErr, BootMaxErr;

  std::chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot, start1, end1;
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
  find_param(argv[1],"__ed",ed);


  // Define Necessary Structures
  INTEG I = INTEG(k);
  GREEN G = GREEN(Nt, Ntau, Nsites, -1);
  GREEN Phi = GREEN(Nt,Ntau,Nsites,1);
  GREEN UxPhi = GREEN(Nt,Ntau,Nsites,1);
  GREEN PhixU = GREEN(Nt,Ntau,Nsites,1);
  GREEN TPP = GREEN(Nt,Ntau,Nsites,1);
  GREEN Sigma = GREEN(Nt, Ntau, Nsites, -1);

  CFUNC hmf = CFUNC(Nt, Nsites);
  CFUNC h0f = CFUNC(Nt,Nsites);
  CFUNC Ut = CFUNC(Nt, Nsites);
  ZMatrix h0(Nsites,Nsites), DensM(Nsites,Nsites);
  
  SPECT A = SPECT();

  // Get free green's func ==============================================================
  h0.setZero();
  for(i=0;i<Nsites-1;i++){
    h0(i,i+1)=-HoppingT;
    h0(i+1,i)=-HoppingT;
    h0(i,i)=ed(i);
  }
  h0f.set_constant(h0);
  Ut.set_constant(HubbardU*ZMatrix::Identity(Nsites,Nsites));
  
  NEdyson::G0_from_h0(G,MuChem,h0,Beta,dt);

  G.get_dm(-1,DensM);
  std::cout.precision(17);
  double npart = DensM.trace().real();
  std::cout << "number of particles = " << npart << std::endl;
  MuChem += HubbardU*npart/Nsites;

  std::string str = argv[1];
  str+="/Gfree";
  G.print_to_file(str,dt,dtau,16);

  A.AfromG(G,101,10,dt);
  str = argv[1];
  str+="/Afree";
  A.print_to_file(str);

  // Matsubara Self Consistency =============================================================
  {
    start = std::chrono::system_clock::now();
    bool converged = false;
    tstp = -1;

    GREEN gtemp = GREEN(k,Ntau,Nsites,-1);
    gtemp.set_tstp(tstp,G);
  
    for(int iter=0;iter<MatsMaxIter;iter++){
      // Update Mean field
      Hubb::Ham_MF(tstp,G,Ut,h0,hmf);

      // Update T-matrix
      Hubb::GenTPP(tstp,dt,Beta,G,Phi,Ut,PhixU,UxPhi,TPP,I);

      // Update Self-Energy
      Hubb::Sigma_TPP(tstp,G,Ut,TPP,Sigma);

      // Solve Dyson Equation
      NEdyson::mat_fourier(G,Sigma,MuChem,hmf.ptr(-1),Beta);
      
      
      double err = distance_norm2(tstp,G,gtemp);
      G.get_dm(-1,DensM);
      npart = DensM.trace().real();
      std::cout<<"iteration: "<<iter<<" | N =  "<<npart<<" |  Error = "<<err<<std::endl;
      if(err<MatsMaxErr){
        converged=true;
        break;
      }
      gtemp.set_tstp(tstp,G);
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

    GREEN gtemp = GREEN(k,Ntau,Nsites,-1);
    for(tstp=0; tstp<=SolverOrder; tstp++) gtemp.set_tstp(tstp,G);

    for(int iter=0; iter<=BootMaxIter; iter++){
      // Update mean field
      for(tstp=0; tstp<=SolverOrder; tstp++){
        Hubb::Ham_MF(tstp, G, Ut, h0, hmf);
      }

      // Update T-matrix
      Hubb::GenTPP(dt,Beta,G,Phi,Ut,PhixU,UxPhi,TPP,I);
      
      // Update self-energy
      for(tstp=0; tstp<=SolverOrder; tstp++){
        Hubb::Sigma_TPP(tstp,G,Ut,TPP,Sigma);
      }

      // Solve Dyson Eqn
      dyson_start(I,G,Sigma,hmf,MuChem,Beta,dt);

      double err=0.0;
      for(tstp=0;tstp<=SolverOrder;tstp++){
        err += distance_norm2(tstp,G,gtemp);
      }
      for(tstp=0; tstp<=SolverOrder; tstp++) gtemp.set_tstp(tstp,G);
      std::cout<<"bootstrapping iteration : "<<iter<<" |  Error = "<<err<<std::endl;
      if(err<BootMaxErr){
        bootstrap_converged=true;
        break;
      }
    }
    if(!bootstrap_converged){
      std::cout<<std::endl;
      std::cout<<" Bootstrapping did not converge. Program Ending. "<<std::endl;
    }
  } // End Bootstrapping

  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "Time [bootstrapping] = " << elapsed_seconds.count() << "s\n\n";  



  // Timestepping ===========================================================================
  start = std::chrono::system_clock::now();

  for(int tstp=k+1;tstp<=Nt;tstp++){

    // Extrapolate
    NEdyson::Extrapolate(I,G,tstp);

    for(int iter=0;iter<StepIter;iter++){
      // Update hmf
      Hubb::Ham_MF(tstp,G,Ut,h0,hmf);

      // Update T-matrix
      Hubb::GenTPP(tstp,dt,Beta,G,Phi,Ut,PhixU,UxPhi,TPP,I);

      // Update self-energy
      Hubb::Sigma_TPP(tstp,G,Ut,TPP,Sigma);

      // Solve Dyson Eqn
      dyson_step(tstp,I,G,Sigma,hmf,MuChem,Beta,dt);

    }

    G.get_dm(tstp,DensM);
    std::cout.precision(17);
    double npart = DensM.trace().real();
    std::cout << tstp << " number of particles = " << npart << std::endl;
  } // End Timestepping

  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << "Time [Propagation] = " << elapsed_seconds.count() << "s\n\n";  

  std::cout << "Total Time = " << (end-start_tot).count() << "s\n\n";
  
  str = argv[1];
  str+="/G";
  G.print_to_file(str,dt,dtau,16);

  A.AfromG(G,101,10,dt);
  str = argv[1];
  A.print_to_file(str);
}
