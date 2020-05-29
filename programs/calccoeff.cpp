#include "NEdyson.h"
#include "utils.h"
using namespace NEdyson;
int main(int argc, char *argv[]){
  GREEN G;
  double dt,dtau;
  G.read_from_file("/home/thomas/Libraries/NEdyson/TPP_test_TPP",dt,dtau);

  h5e::File file("/home/thomas/Libraries/NEdyson/tests/data/vie2_test_data.h5", h5e::File::Overwrite);
  G.print_to_file(file,"TPP");

  G.read_from_file("/home/thomas/Libraries/NEdyson/TPP_test_PHIxU",dt,dtau);
  G.print_to_file(file,"PHIxU");

  G.read_from_file("/home/thomas/Libraries/NEdyson/TPP_test_UxPHI",dt,dtau);
  G.print_to_file(file,"UxPHI");

  G.read_from_file("/home/thomas/Libraries/NEdyson/TPP_test_Phi",dt,dtau);
  G.print_to_file(file,"Phi");

  G.read_from_file("/home/thomas/Libraries/NEdyson/TPP_test_Sigma",dt,dtau);
  G.print_to_file(file,"Sigma");

  G.read_from_file("/home/thomas/Libraries/NEdyson/TPP_test_G",dt,dtau);
  G.print_to_file(file,"G");




  CFUNC Ut(G.nt(), G.size1());
  CFUNC eps_mf(G.nt(), G.size1());

  std::string str("/home/thomas/Libraries/NEdyson/TPP_test_Ut.dat");
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

  eps_mf.print_to_file(file,"Ut");


  return 1;
  
}
