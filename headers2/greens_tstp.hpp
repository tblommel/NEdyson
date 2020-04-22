#ifndef GREEN_FUNC_TSTP_DECL
#define GREEN_FUNC_TSTP_DECL

#include <complex>
#include "function.hpp"
#include "elementops.hpp"

#define CFUNC NEdyson::function

namespace NEdyson{

class green_func_tstp{
	private:
		cplx *data_;
		int tstp_;
		int ntau_;
		int size1_;
		int element_size_;
		int sig_;

	public:
		green_func_tstp();
		~green_func_tstp();
		green_func_tstp(int tstp, int ntau, int size1, int sig);
		green_func_tstp(const green_func_tstp &g);
		green_func_tstp & operator=(const green_func_tstp &g);

		int size1(void) const {return size1_;}
		int tstp(void) const {return tstp_;}
		int ntau(void) const {return ntau_;}
		int element_size(void) const {return element_size_;}
		int sig(void) const {return sig_;}

		inline cplx *retptr(int i) {return data_ + i*element_size_;}
		inline cplx *tvptr(int i) {return data_ + (tstp_+1+i)*element_size_;}
		inline cplx *lesptr(int i) {return data_ + (tstp_+2+ntau_+i)*element_size_;}
		inline cplx *matptr(int i) {return data_ + i* element_size_;}

    void left_multiply(cplx *f0, cplx *ft, cplx weight = 1.0);
    void right_multiply(cplx *f0, cplx *ft, cplx weight = 1.0);
    void left_multiply(CFUNC &ft, cplx weight = 1.0);
    void right_multiply(CFUNC &ft, cplx weight = 1.0);
		
		void smul(cplx weight);

};// class
} // namespace
#endif
