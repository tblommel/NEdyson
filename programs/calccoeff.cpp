#include "NEdyson.h"
#include "utils.h"

using namespace NEdyson;
int main(int argc, char *argv[]){
  GREEN G;
  double dt,dtau;
  h5e::File readfile("/home/thomas/Libraries/NEdyson/tests/data/Gfree_const.h5",h5e::File::ReadWrite);
  G.read_from_file(readfile,"");
  G.print_to_file("/home/thomas/Libraries/NEdyson/Gfree",0.1,9./30.);
  
  return 1;
/*
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
  */
}
