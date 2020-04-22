#ifndef FUNCTION_DECL
#define FUNCTION_DECL

#include <complex>
#include <Eigen/Eigen>
#include <cassert>
#include <iostream>

#define cplx std::complex<double>
#define cdmatrix Eigen::MatrixXcd
namespace NEdyson{


class function{
	public:
		function();
		~function();
		function(int nt, int size1);
		function(const function &f);
		function &operator=(const function &f);
		
		int element_size(void) {return element_size_;}
		int size1(void) {return size1_;}
		int nt(void) {return nt_;}
		cplx *ptr(int t) {return data_ + (t+1)*element_size_;}

		void resize(int nt, int size1);
		void set_Zero(void);
		void set_value(int tstp, cdmatrix M);
		void set_value(int tstp, cplx x);
		void set_value(int tstp, cplx *x);
		void set_constant(cdmatrix M);
		void set_constant(cplx x);
		void set_constant(cplx *x);
		
		void get_value(int tstp, cdmatrix &M);
		void get_value(int tstp, cplx &x);
		void get_value(int tstp, cplx *x);

		void printf(void);

	private:
		cplx *data_;
		int nt_;
		int size1_;
		int element_size_;

};
}//Namespace


#endif //FUNCTION_DECL
