#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <vector>
#include <chrono>

#include "NEdyson.h"

using namespace NEdyson;

bool is_file_exist(std::string fileName){
	std::ifstream infile(fileName);
	return infile.good();
}

int main(int argc, char *argv[]){
  if(argc!=5) throw std::invalid_argument("please provide G input file, A input file, nw, and wmax");
	int nw = strtol(argv[3],NULL,10);
	double wmax = strtod(argv[4],NULL);

	// Define the structures
	GREEN G = GREEN();
	SPECT A = SPECT();
	double dt, dtau;
	
	// Read in the Greens function
	G.read_from_file(argv[1],dt,dtau);

	// Check to see if there exists an extension file
	std::string actfile = "";
	actfile += argv[1];
	actfile += "_GRExp.dat";
	bool extexist = is_file_exist(actfile);

	// If it does read it in
	if(extexist){
		std::ifstream in;
		in.open(actfile);
		int nfit,ntp;
		if(!(in >> nfit >> ntp)){
			std::cerr<<"expanded file broken"<<actfile<<std::endl;
			abort();
		}
		std::cout<<nfit<<ntp<<std::endl;
		
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
		A.AfromG(G, nw, wmax, dt, extdata, nfit, ntp);
		delete[] extdata;
	}
	else{
		A.AfromG(G, nw, wmax, dt);
	}



	// Print A to file
	actfile = "";
	actfile += argv[2];
	A.print_to_file(actfile);
}	
