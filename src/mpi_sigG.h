#ifndef GREENS_FUNCTIONS_DECL
#define GREENS_FUNCTIONS_DECL

#include <complex>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <cstring>

#include "utils.h"

namespace NEdyson{

class mpi_sigG{
	public:	
		// Construct Destruct
		~green_func();
		green_func(int nt, int ntau, int size1, int sig);
		
		// Get sizes
		int es() const {return es_;}
		int nao() const {return nao_;}
		int ntau() const {return ntau_;}
		int nt() const {return nt_;}
		int tid() const {return tid_;}
		int nthreads() const {return nthreads_;}

		// Pointers
		cplx *lesptr(int i, int j) const ;
		cplx *retptr(int i, int j) const ;
		cplx *tvptr(int i, int j) const ;
		cplx *matptr(int i) const ;

    // Density Matrix
		void get_dm(int i, ZMatrix &M) const {
 		  assert(M.rows()==size1_ && M.cols()==size1_ && i>=-1 && i<=nt_);
		  if(i==-1){
        M = -ZMatrixMap(matptr(ntau_), size1_, size1_);
		  } else {
        M = cplx(0.0, 1.0*sig_)*ZMatrixMap(lesptr(i,i), size1_, size1_);
		  }
		}

    void get_dm(int i, ZTensor<2> &rho) const {
      assert(rho.shape()[0] == size1_);
      assert(rho.shape()[1] == size1_);
      assert(i>=-1 && i<=nt_);
      if(i==-1){
        ZMatrixMap(rho.data(), size1_, size1_) = -ZMatrixMap(matptr(ntau_), size1_, size1_);
      } else {
        ZMatrixMap(rho.data(), size1_, size1_) = cplx(0.0, 1.0*sig_)*ZMatrixMap(lesptr(i,i), size1_, size1_);
      }
    }

    void get_dm(int i, cplx *rho) const {
      assert(i>=-1 && i<=nt_);
      if(i==-1){
        ZMatrixMap(rho, size1_, size1_) = -ZMatrixMap(matptr(ntau_), size1_, size1_);
      } else {
        ZMatrixMap(rho, size1_, size1_) = cplx(0.0, 1.0*sig_)*ZMatrixMap(lesptr(i,i), size1_, size1_);
      }
    }

    void set_tstp_zero(int tstp){
      assert(tstp >= -1 && tstp <= nt_);
      if(tstp == -1){
        memset(matptr(0), 0, (ntau_+1)*element_size_*sizeof(cplx));
      }
      else{
        memset(retptr(tstp,0), 0, (tstp+1)*element_size_*sizeof(cplx));
        memset(lesptr(0,tstp), 0, (tstp+1)*element_size_*sizeof(cplx));
        memset(tvptr(tstp,0), 0, (ntau_+1)*element_size_*sizeof(cplx));
      }
    }

		/*Input/Output*/
		void print_to_file(std::string file, double dt, double dtau, int precision=16) const ;
		void print_to_file_mat(std::string file, double dt, double dtau, int precision=16) const ;
		void print_to_file_ret(std::string file, double dt, double dtau, int precision=16) const ;
		void print_to_file_tv(std::string file, double dt, double dtau, int precision=16) const ;
		void print_to_file_les(std::string file, double dt, double dtau, int precision=16) const ;
		void read_from_file(const char *file, double &dt, double &dtau);
		void read_from_file_ret(const char *file, double &dt, double &dtau);

    
		void print_to_file(H5Easy::File &File, std::string path) const ;
    void print_to_file_ret(H5Easy::File &File, std::string path) const ;
		void print_to_file_mat(H5Easy::File &File, std::string path) const ;
		void print_to_file_les(H5Easy::File &File, std::string path) const ;
		void print_to_file_tv(H5Easy::File &File, std::string path) const ;
    void read_from_file(h5e::File &File, std::string path) ;
    void read_from_file_mat(h5e::File &File, std::string path) ;
    void read_from_file_ret(h5e::File &File, std::string path) ;
    void read_from_file_les(h5e::File &File, std::string path) ;
    void read_from_file_tv(h5e::File &File, std::string path) ;


	private:
    // Values of Greens Func
    ZTensor<1> Gles_;
    ZTensor<1> Gret_;
    ZTensor<1> Gtv_;
    ZTensor<1> Gmat_;

    // Values of Sigma
    ZTensor<1> Sles_;
    ZTensor<1> Sret_;
    ZTensor<1> Stv_;
    ZTensor<1> Smat_;

    // Density Matrix
    ZTensor<1> Dm_;

    // Temporary timestep storage variables
    ZTensor<1> tstp_Sles_;
    ZTensor<1> tstp_Sret_;
    ZTensor<1> tstp_Stv_;

    int tid_;
    int nthreads_;

		int nt_;
		int ntau_;
		int nao_;
		int es_;
		int sig_;
};

}

#endif
