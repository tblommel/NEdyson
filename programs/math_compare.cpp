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
  GREEN G = GREEN(1, 600, 2, -1);
  SPECT A = SPECT();
  
  int nt = 1;
  int ntau = 600;
  int nsites = 2;
  int sig = -1;
  double dt = 0.015;
  double beta=20;
  double MUChem = 0;

  double t=1, U=0, mu1 = -1, mu2 = 1;
  cdmatrix h0(2,2);
  h0.setZero();
  h0(0,0)=mu1;
  h0(0,1)=-t;
  h0(1,0)=-t;
  h0(1,1)=mu2;

  NEdyson::G0_from_h0(G,MUChem,h0,beta,dt);
  double dtau=beta/600;
  cdmatrix DensM;
  G.get_dm(-1,DensM);

  std::cout.precision(17);
  double npart = DensM.trace().real();
  std::cout << "number of particles = " << npart << std::endl;
  std::string str = ".";
  str+="/Gfree";
  G.print_to_file(str,dt,dtau,16);
}
