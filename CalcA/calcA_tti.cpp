#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <vector>
#include <chrono>

#include "greens.h"
#include "spectral.h"
#include "utils.h"

using namespace NEdyson;

int main(int argc, char *argv[]){
  if(argc!=7) throw std::invalid_argument("please provide G input file, GRfilepath,  A output file, nw, wmin, wmax");
	int nw = strtol(argv[4],NULL,10);
	double wmin = strtod(argv[5],NULL);
	double wmax = strtod(argv[6],NULL);

	// Define the structures
	TTI_GREEN G = TTI_GREEN();
	TTI_SPECT A = TTI_SPECT();
	double dt;

	// Read in the Greens function
	h5e::File File(argv[1], h5e::File::ReadWrite);
	G.read_from_file_ret(File, argv[2]);
  dt = h5e::load<double>(File, "/solve/params/dt");

  // Calculate A
	A.AfromG(G, nw, wmin, wmax, dt);

	// Print A to file
	std::string data_dir = std::string(DATA_DIR);
  std::ofstream ofile;
	ofile.open(data_dir + argv[3], std::ofstream::out);
  for(int w=0; w<nw; w++) {
    ofile << w*A.dw() + wmin << " ";
    for(int i=0; i<A.size1(); i++) {
      ofile << A.ptr(w)[i*A.size1()+i] << " ";
    }
    ofile << std::endl;
  }
  ofile.close();

	
	/*
	h5e::File File2(argv[3],h5e::File::ReadWrite);
	A.print_to_file(File2,"/A");
  */
}		
