#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <vector>
#include <Eigen/Eigen>
#include <chrono>

#include <headers2/NEdyson.hpp>

#define cdmatrix Eigen::MatrixXcd
#define GREEN NEdyson::green_func
#define INTEG integration::Integrator
#define FUNC NEdyson::function
#define SPECT NEdyson::spectral

using namespace std;


int main(int argc, char *argv[]){
	if(argc!=2) throw std::invalid_argument("please provide input file");
	
	// Parameters ========================================================================
	int Nsites, Ntau, MatsMaxIter, oversamp, Nt, RampSite, k, SolverOrder, BootMaxIter, StepIter;
	double HoppingT, HubbardU, MuChem, Beta, dt, dtau, W0, MatsMaxErr, BootMaxErr;

	chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot;
	int tstp,i,j;
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
	// Define Necessary Structures
	INTEG I = INTEG(k);
	GREEN G = GREEN(Nt, Ntau, Nsites, -1);
	GREEN Sigma = GREEN(Nt, Ntau, Nsites, -1);
	FUNC hmf = FUNC(Nt, Nsites);
	cdmatrix h0(Nsites,Nsites), DensM(Nsites,Nsites);
	FUNC Ut = FUNC(Nt, Nsites);
	SPECT A = SPECT();

	// Get free green's func ==============================================================
	h0.setZero();
	for(i=0;i<Nsites-1;i++){
		h0(i,i+1)=-HoppingT;
		h0(i+1,i)=-HoppingT;
	}
	h0(0,0) = -7;
	
	FUNC h0f = FUNC(Nt,Nsites);
	h0f.set_constant(h0);

	Ut.set_constant(HubbardU*cdmatrix::Identity(Nsites,Nsites));
	
	NEdyson::G0_from_h0(G,MuChem,h0,Beta,dt);

	std::string str = argv[1];
	str+="/Gfree";
	G.print_to_file(str,dt,dtau,16);

	G.get_dm(-1,DensM);
	cout.precision(17);
	double npart = DensM.trace().real();
	cout << "number of particles = " << npart << endl;
	MuChem += HubbardU*npart/Nsites;


	// Matsubara Self Consistency =============================================================
	{
		start = chrono::system_clock::now();
		bool converged = false;
		tstp = -1;
		
		
		GREEN gtemp = GREEN(k,Ntau,Nsites,-1);
		gtemp.set_tstp(tstp,G);
	
		for(int iter=0;iter<MatsMaxIter;iter++){
			// Update Mean field and Self Energy
			Hubb::Ham_MF(tstp,G,Ut,h0,hmf);
			Hubb::Sigma_2B(tstp,G,Ut,Sigma);
	
			// Solve Dyson Equation
			NEdyson::mat_fourier(G,Sigma,MuChem,hmf.ptr(-1),Beta);
			
			
			double err = distance_norm2(tstp,G,gtemp);
			G.get_dm(-1,DensM);
			npart = DensM.trace().real();
			cout<<"iteration: "<<iter<<" | N =  "<<npart<<" |  Error = "<<err<<endl;
			if(err<MatsMaxErr){
				converged=true;
				break;
			}
			gtemp.set_tstp(tstp,G);
			if(iter==40){
				str = argv[1];
				str+="/Ginter0";
				G.print_to_file(str,dt,dtau,16);
			}
			if(iter==41){
				str = argv[1];
				str+="/Ginter1";
				G.print_to_file(str,dt,dtau,16);
			}
			if(iter==42){
				str = argv[1];
				str+="/Ginter2";
				G.print_to_file(str,dt,dtau,16);
			}
			if(iter==43){
				str = argv[1];
				str+="/Ginter3";
				G.print_to_file(str,dt,dtau,16);
				return 1;
			}
		}

		if(!converged){
			cout<<endl;
			cout<<" Matsubara Iteration not converging! Quitting program. " <<endl;
			return 0;
		}
	
		end = chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Time [equilibrium calculation] = " << elapsed_seconds.count() << "s\n\n";	
	} // End Matsubara Self Consistency Loop




	// Bootstrapping ==========================================================================
	{
		start = std::chrono::system_clock::now();
		bool bootstrap_converged=false;
		
		GREEN gtemp = GREEN(k,Ntau,Nsites,-1);
		for(tstp=0; tstp<=SolverOrder; tstp++) gtemp.set_tstp(tstp,G);

		for(int iter=0; iter<=BootMaxIter; iter++){
			// Update mean field
			for(tstp=0; tstp<=SolverOrder; tstp++){
				Hubb::Ham_MF(tstp, G, Ut, h0, hmf);
			}

			// Update self-energy
			for(tstp=0; tstp<=SolverOrder; tstp++){
				Hubb::Sigma_2B(tstp, G, Ut, Sigma);
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
    cout << "Time [bootstrapping] = " << elapsed_seconds.count() << "s\n\n";	






	// Timestepping ===========================================================================
		start = std::chrono::system_clock::now();
	for(int tstp=k+1;tstp<=Nt;tstp++){
		// Extrapolate
		NEdyson::Extrapolate(I,G,tstp);
		for(int iter=0;iter<StepIter;iter++){
			// Update hmf
			Hubb::Ham_MF(tstp,G,Ut,h0,hmf);
			
			// Update Self Energy
			Hubb::Sigma_2B(tstp,G,Ut,Sigma);

			// Solve Dyson Eqn
			dyson_step(tstp,I,G,Sigma,hmf,MuChem,Beta,dt);
		}
	}

		end = std::chrono::system_clock::now();
		elapsed_seconds = end-start;
    cout << "Time [Propagation] = " << elapsed_seconds.count() << "s\n\n";
	
		str = argv[1];
		str+="/G";
		
		G.print_to_file(str,dt,dtau,16);

		A.AfromG(G,101,10,dt);
		str = argv[1];
		A.print_to_file(str);
}
