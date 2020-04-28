#include<iostream>
#include"hybfermion.hpp"
int main(int argc, char** argv){
  if(argc!=2) throw std::invalid_argument("call this code with ./fermion number_of_orbitals");
  fermion_op f(atoi(argv[1]));
}
