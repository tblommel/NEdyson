#include "NEdyson.h"
#include "utils.h"
using namespace NEdyson;
int main(int argc, char *argv[]){
  GREEN G;
  double dt,dtau;
  G.read_from_file("/home/thomas/Libraries/NEdyson/dyson_test_Gextr",dt,dtau);

  h5e::File file("/home/thomas/Libraries/NEdyson/tests/data/dyson_test_data.h5", h5e::File::Overwrite);
  G.print_to_file(file,"Extrap");

  G.read_from_file("/home/thomas/Libraries/NEdyson/dyson_test_Gfree",dt,dtau);
  G.print_to_file(file,"Free");

  G.read_from_file("/home/thomas/Libraries/NEdyson/dyson_test_GMat",dt,dtau);
  G.print_to_file(file,"MatFourier");

  G.read_from_file("/home/thomas/Libraries/NEdyson/dyson_test_Gstart",dt,dtau);
  G.print_to_file(file,"Start");

  G.read_from_file("/home/thomas/Libraries/NEdyson/dyson_test_Sig",dt,dtau);
  G.print_to_file(file,"Sigma");

  G.read_from_file("/home/thomas/Libraries/NEdyson/dyson_test_GStep",dt,dtau);
  G.print_to_file(file,"Timesteps");




  CFUNC Ut(G.nt(), G.size1());
  CFUNC eps_mf(G.nt(), G.size1());

  std::string str("/home/thomas/Libraries/NEdyson/dyson_test_eps_mf.dat");
  std::ifstream inpfile;
  inpfile.open(str);
  double real, imag;
  for(int t=-1;t<=G.nt();t++){
    for(int i=0;i<4;i++){
      inpfile>>real>>imag;
      eps_mf.ptr(t)[i] = cplx(real,imag);
    }
  }
  inpfile.close();

  eps_mf.print_to_file(file,"epsmf");


  return 1;
  
}
