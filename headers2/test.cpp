#include <iostream>
#include "elementops.hpp"

int main(int argc, char *argv[]){
	std::complex<double> *M = new std::complex<double> [24];
	Eigen::MatrixXd ME(4,5);
	for(int i=0;i<4;i++){
		for(int j=0;j<5;j++){
			M[i*5+j]=i+j;
			ME(i,j)=i+j;
		}
	}

	std::cout<<ME<<std::endl<<std::endl;
	std::cout<<map(M,4,5)<<std::endl;

	

	delete[] M;
}
