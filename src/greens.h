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

class green_func{
	public:	
		/*construct, destruct*/
		green_func();
		~green_func();
		green_func(int nt, int ntau, int size1, int sig);
		green_func(const green_func &g);
		green_func &operator=(const green_func &g);
		
		/*get sizes*/
		int element_size(void) const {return element_size_;}
		int size1(void) const {return size1_;}
		int ntau(void) const {return ntau_;}
		int nt(void) const {return nt_;}
		int sig(void) const {return sig_;}

    void resize(int nt, int ntau, int size1);

		/*get pointers*/
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
		cplx *les_;
		cplx *ret_;
		cplx *tv_;
		cplx *mat_;
		int nt_;
		int ntau_;
		int size1_;
		int element_size_;
		int sig_;
};

class tti_green_func{
	public:	
		/*construct, destruct*/
		tti_green_func();
		~tti_green_func();
		tti_green_func(int nt, int ntau, int size1, int sig);
		tti_green_func(const tti_green_func &g);
		tti_green_func &operator=(const tti_green_func &g);
		
		/*get sizes*/
		int element_size(void) const {return element_size_;}
		int size1(void) const {return size1_;}
		int ntau(void) const {return ntau_;}
		int nt(void) const {return nt_;}
		int sig(void) const {return sig_;}

    void resize(int nt, int ntau, int size1);


		/*get pointers*/
		cplx *lesptr(int i) const ;
		cplx *retptr(int i) const ;
		cplx *tvptr(int i, int j) const ;
		cplx *matptr(int i) const ;

    // Density Matrix
		void get_dm(int i, ZMatrix &M) const {
 		  assert(M.rows()==size1_ && M.cols()==size1_ && i>=-1 && i<=nt_);
		  if(i==-1){
        M = -ZMatrixMap(matptr(ntau_), size1_, size1_);
		  } else {
        M = cplx(0.0, 1.0*sig_)*ZMatrixMap(lesptr(0), size1_, size1_);
		  }
		}

    void get_dm(int i, ZTensor<2> &rho) const {
      assert(rho.shape()[0] == size1_);
      assert(rho.shape()[1] == size1_);
      assert(i>=-1 && i<=nt_);
      if(i==-1){
        ZMatrixMap(rho.data(), size1_, size1_) = -ZMatrixMap(matptr(ntau_), size1_, size1_);
      } else {
        ZMatrixMap(rho.data(), size1_, size1_) = cplx(0.0, 1.0*sig_)*ZMatrixMap(lesptr(0), size1_, size1_);
      }
    }

    void get_dm(int i, cplx *rho) const {
      assert(i>=-1 && i<=nt_);
      if(i==-1){
        ZMatrixMap(rho, size1_, size1_) = -ZMatrixMap(matptr(ntau_), size1_, size1_);
      } else {
        ZMatrixMap(rho, size1_, size1_) = cplx(0.0, 1.0*sig_)*ZMatrixMap(lesptr(0), size1_, size1_);
      }
    }

    void set_tstp_zero(int tstp){
      assert(tstp >= -1 && tstp <= nt_);
      if(tstp == -1){
        memset(matptr(0), 0, (ntau_+1)*element_size_*sizeof(cplx));
      }
      else{
        memset(retptr(tstp), 0, element_size_*sizeof(cplx));
        memset(lesptr(-tstp), 0, element_size_*sizeof(cplx));
        memset(tvptr(tstp,0), 0, (ntau_+1)*element_size_*sizeof(cplx));
      }
    }

		/*Input/Output*/
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
		cplx *les_;
		cplx *ret_;
		cplx *tv_;
		cplx *mat_;
		int nt_;
		int ntau_;
		int size1_;
		int element_size_;
		int sig_;
};

}

#endif
