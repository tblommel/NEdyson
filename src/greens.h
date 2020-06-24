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
#include "function.h"
#include "elementops.h"
#include "greens_tstp.h"

namespace NEdyson{

class green_func{
	public:	
		/*construct, destruct*/
		green_func();
		~green_func();
		green_func(int nt, int ntau, int size1, int sig);
		green_func(const green_func &g);
		green_func &operator=(const green_func &g);
		void resize(int nt, int ntau, int size1);
		
		/*get sizes*/
		int element_size(void) const {return element_size_;}
		int size1(void) const {return size1_;}
		int ntau(void) const {return ntau_;}
		int nt(void) const {return nt_;}
		int sig(void) const {return sig_;}

		/*get pointers*/
		cplx *lesptr(int i, int j) const ;
		cplx *retptr(int i, int j) const ;
		cplx *tvptr(int i, int j) const ;
		cplx *matptr(int i) const ;

		/*fill matricies*/
		#define read_element                 \
		        {                            \
		          int r,s;                   \
		          M.resize(size1_,size1_);   \
		          for(r=0;r<size1_;r++)      \
		            for(s=0;s<size1_;s++)    \
		              M(r,s)=x[r*size1_+s];  \
		        }
		
		#define read_element_minus_conj      \
		        {                            \
		          int r,s;                   \
		          cplx w;                    \
		          M.resize(size1_,size1_);   \
		          for(r=0;r<size1_;r++)      \
		            for(s=0;s<size1_;s++){   \
		              w=x[r*size1_+s];       \
		              M(s,r)=std::complex<double>(-w.real(),w.imag());  \
		            }                                              \
		        } 
		
		template <class Matrix>
		void get_les(int i, int j, Matrix &M) const {
		  assert(i<=nt_ && j<=nt_);
		  cplx *x;
		  if(i<=j){
		    x=lesptr(i,j);
		    read_element
		  }else{
		    x=lesptr(j,i);
		    read_element_minus_conj
		  }
		}
		
		
		template <class Matrix>
		void get_ret(int i, int j, Matrix &M) const {
		  assert(i<=nt_ && j<=nt_);
		  cplx *x;
		  if(i>=j){
		    x=retptr(i,j);
		    read_element
		  }else{
		    x=retptr(j,i);
		    read_element_minus_conj
		  }
		}
		
		
		template <class Matrix>
		void get_adv(int i, int j, Matrix &M) const {
		  assert(i<=nt_ && j<=nt_);
		  get_ret(j,i,M);
		  M.adjointInPlace();
		}
		
		
		template <class Matrix>
		void get_tv(int i, int j, Matrix &M) const {
		  assert(i<=nt_ && j<=ntau_);
		  cplx *x=tvptr(i,j);
		  read_element
		}
		
		
		template <class Matrix>
		void get_vt(int i, int j, Matrix &M) const {
		  assert(i<=ntau_ && j<=nt_);
		  cplx *x= tvptr(j,ntau_-i);
		  read_element_minus_conj if (sig_==-1) M=-M;
		}
		
		
		template <class Matrix>
		void get_mat(int j, Matrix &M) const {
		  assert(j<=ntau_);
		  cplx *x=matptr(j);
		  read_element
		}
		
		
		template <class Matrix>
		void get_grt(int i, int j, Matrix &M) const {
		  assert(i<=nt_ && j<=nt_);
		  Matrix M1;
		  get_ret(i,j,M);
		  get_les(i,j,M1);
		  M+=M1;
		}

		template <class Matrix>
		void get_dm(int i, Matrix &M) const {
 		  assert(M.rows()==size1_ && M.cols()==size1_ && i>=-1 && i<=nt_);
		  if(i==-1){
		    get_mat(ntau_,M);
		    M*=-1.;
		  } else {
		    get_les(i,i,M);
		    M*=cplx(0.0,1.0*sig_);
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

		#undef read_element
		#undef read_element_minus_conj


    #define set_element          \
      {                          \
        int r,s;                 \
        for(r=0;r<size1_;r++)    \
          for(s=0;s<size1_;s++)  \
            x[r*size1_+s]=M(r,s);\
      }
      template <class Matrix>
      void set_ret(int i, int j, Matrix &M) {
        assert(i <= nt_ && j <= i);
        cplx *x = retptr(i, j);
        set_element
      }

      template <class Matrix>
      void set_les(int i, int j, Matrix &M) {
        assert(j <= nt_ && i <= j);
        cplx *x = lesptr(i, j);
        set_element
      }

      template <class Matrix>
      void set_tv(int i, int j, Matrix &M) {
        assert(i <= nt_ && j <= ntau_);
        cplx *x = tvptr(i, j);
        set_element
      }

      template <class Matrix>
      void set_mat(int j, Matrix &M) {
        assert(j <= ntau_);
        cplx *x = matptr(j);
        set_element
      }
    #undef set_element
      

    // set local timestep from G
    void set_tstp(int tstp, const green_func &G);
    void set_tstp(int tstp, const green_func_tstp &G);
    void set_tstp_zero(int tstp);
    // give tstp to G
		void get_tstp(int tstp, green_func_tstp &G) const;


		/*Multiplications*/
		void smul(int tstp, cplx weight);
    // Only for calculating self-energy in the hubbard model
		void right_multiply(int tstp, const function &ft, cplx weight=1.);
		void left_multiply(int tstp, const function &ft, cplx weight=1.);
		void right_multiply(int tstp, cplx *f0, cplx *ft, cplx weight=1.);
		void left_multiply(int tstp, cplx *f0, cplx *ft, cplx weight=1.);



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

double distance_norm2(int tstp, green_func &G1, green_func &G2); 

class tti_green_func{
	public:	
		/*construct, destruct*/
		tti_green_func();
		~tti_green_func();
		tti_green_func(int nt, int ntau, int size1, int sig);
		tti_green_func(const tti_green_func &g);
		tti_green_func &operator=(const tti_green_func &g);
		void resize(int nt, int ntau, int size1);
		
		/*get sizes*/
		int element_size(void) const {return element_size_;}
		int size1(void) const {return size1_;}
		int ntau(void) const {return ntau_;}
		int nt(void) const {return nt_;}
		int sig(void) const {return sig_;}

		/*get pointers*/
		cplx *lesptr(int i) const ;
		cplx *retptr(int i) const ;
		cplx *tvptr(int i, int j) const ;
		cplx *matptr(int i) const ;

		/*fill matricies*/
		#define read_element                 \
		        {                            \
		          int r,s;                   \
		          M.resize(size1_,size1_);   \
		          for(r=0;r<size1_;r++)      \
		            for(s=0;s<size1_;s++)    \
		              M(r,s)=x[r*size1_+s];  \
		        }
		
		#define read_element_minus_conj      \
		        {                            \
		          int r,s;                   \
		          cplx w;                    \
		          M.resize(size1_,size1_);   \
		          for(r=0;r<size1_;r++)      \
		            for(s=0;s<size1_;s++){   \
		              w=x[r*size1_+s];       \
		              M(s,r)=std::complex<double>(-w.real(),w.imag());  \
		            }                                              \
		        } 
		
		template <class Matrix>
		void get_les(int i, Matrix &M) const {
		  assert(i<=nt_ && i>=-nt_);
		  cplx *x;
		  if(i<=0){
		    x=lesptr(-i);
		    read_element
		  }else{
		    x=lesptr(i);
		    read_element_minus_conj
		  }
		}
		
		
		template <class Matrix>
		void get_ret(int i, Matrix &M) const {
		  assert(i<=nt_ && i>=-nt_);
		  cplx *x;
		  if(i>=0){
		    x=retptr(i);
		    read_element
		  }else{
		    x=retptr(-i);
		    read_element_minus_conj
		  }
		}
		
		
		template <class Matrix>
		void get_adv(int i, Matrix &M) const {
		  assert(i<=nt_ && i>-=nt_);
		  get_ret(-i,M);
		  M.adjointInPlace();
		}
		
		
		template <class Matrix>
		void get_tv(int i, int j, Matrix &M) const {
		  assert(i<=nt_ && j<=ntau_);
		  cplx *x=tvptr(i,j);
		  read_element
		}
		
		
		template <class Matrix>
		void get_vt(int i, int j, Matrix &M) const {
		  assert(i<=ntau_ && j<=nt_);
		  cplx *x= tvptr(j,ntau_-i);
		  read_element_minus_conj if (sig_==-1) M=-M;
		}
		
		
		template <class Matrix>
		void get_mat(int j, Matrix &M) const {
		  assert(j<=ntau_);
		  cplx *x=matptr(j);
		  read_element
		}
		
		
		template <class Matrix>
		void get_grt(int i, Matrix &M) const {
		  assert(i<=nt_ && i>=-nt_);
		  Matrix M1;
		  get_ret(i,M);
		  get_les(i,M1);
		  M+=M1;
		}

		template <class Matrix>
		void get_dm(int i, Matrix &M) const {
 		  assert(M.rows()==size1_ && M.cols()==size1_ && i>=-1 && i<=nt_);
		  if(i==-1){
		    get_mat(ntau_,M);
		    M*=-1.;
		  } else {
		    get_les(0,M);
		    M*=cplx(0.0,1.0*sig_);
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

		#undef read_element
		#undef read_element_minus_conj

    #define set_element          \
      {                          \
        int r,s;                 \
        for(r=0;r<size1_;r++)    \
          for(s=0;s<size1_;s++)  \
            x[r*size1_+s]=M(r,s);\
      }
      template <class Matrix>
      void set_ret(int i, Matrix &M) {
        assert(i <= nt_ && i >= 0);
        cplx *x = retptr(i);
        set_element
      }

      template <class Matrix>
      void set_les(int i, Matrix &M) {
        assert(i <= 0 && i >= -nt_);
        cplx *x = lesptr(i);
        set_element
      }

      template <class Matrix>
      void set_tv(int i, int j, Matrix &M) {
        assert(i <= nt_ && j <= ntau_);
        cplx *x = tvptr(i, j);
        set_element
      }

      template <class Matrix>
      void set_mat(int j, Matrix &M) {
        assert(j <= ntau_);
        cplx *x = matptr(j);
        set_element
      }
    #undef set_element

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
