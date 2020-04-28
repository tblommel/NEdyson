#ifndef GREEN_FUNC_TSTP_IMPL
#define GREEN_FUNC_TSTP_IMPL

#include "greens_tstp.h"

namespace NEdyson{

green_func_tstp::green_func_tstp(){
	data_=0;
	ntau_=0;
	tstp_=0;
	size1_=0;
	element_size_=0;
	sig_=-1;
}

green_func_tstp::~green_func_tstp(){
	delete[] data_;
}


green_func_tstp::green_func_tstp(int tstp, int ntau, int size1, int sig){
	assert(size1>0 && tstp >=-1 && ntau>=0 && sig*sig==1);
	int len = ((tstp+1)*2+(ntau+1))*size1*size1;
	if(len==0) data_=0;
	else{
		data_ = new cplx[len];
		memset(data_,0,sizeof(cplx)*len);
	}
	size1_=size1;
	element_size_=size1*size1;
	tstp_=tstp;
	ntau_=ntau;
	sig_=sig;
}


green_func_tstp::green_func_tstp(const green_func_tstp &g){
	tstp_=g.tstp_;
	ntau_=g.ntau_;
	size1_=g.size1_;
	element_size_=g.element_size_;
	int len = ((tstp_+1)*2+(ntau_+1))*element_size_;
	if(len==0) data_=0;
	else{
		data_ = new cplx[len];
		memcpy(data_,g.data_,sizeof(cplx)*len);
	}
	sig_=g.sig_;
}


green_func_tstp &green_func_tstp::operator=(const green_func_tstp &g){
	if(this == &g) return *this;
	int len1 = ((tstp_+1)*2+(ntau_+1))*element_size_;
	int len2 = ((g.tstp_+1)*2+(g.ntau_+1))*g.element_size_;
	if(len1 != len2){
		delete[] data_;
		if(len2==0) data_=0;
		else data_ = new cplx[len2];
	}
	size1_=g.size1_;
	element_size_=size1_*size1_;
	tstp_=g.tstp_;
	ntau_=g.ntau_;
	sig_=g.sig_;
	return *this;
}

void green_func_tstp::left_multiply(cplx *f0, cplx *ft, cplx weight){
	int m;
	cplx *x0, *xtmp, *ftmp;
	xtmp = new cplx[element_size_];
	if(tstp_==-1){
		// Matsubara
		x0 = data_;
		for(m=0;m<=ntau_;m++){
			element_mult(size1_,xtmp,f0,x0);
			element_smul(size1_,xtmp,weight);
			element_set(size1_,x0,xtmp);
			x0+=element_size_;
		}
	}
	else{
		// Multiply the Retarded Component
		ftmp = ft + tstp_*element_size_;
		x0=data_;
		for(m=0;m<=tstp_;m++){
			element_mult(size1_,xtmp,ftmp,x0);
			element_smul(size1_,xtmp,weight);
			element_set(size1_,x0,xtmp);
			x0+=element_size_;
		}
		// Multiply the TV Component
		for(m=0;m<=ntau_;m++){
			element_mult(size1_,xtmp,ftmp,x0);
			element_smul(size1_,xtmp,weight);
			element_set(size1_,x0,xtmp);
			x0+=element_size_;
		}
		// Multiply the Lesser Component
		for(m=0;m<=tstp_;m++){
			element_mult(size1_,xtmp,ft+m*element_size_,x0);
			element_smul(size1_,xtmp,weight);
			element_set(size1_,x0,xtmp);
			x0+=element_size_;
		}
	}
	delete[] xtmp;
}

void green_func_tstp::left_multiply(const CFUNC &ft, cplx weight){
	this->left_multiply(ft.ptr(-1),ft.ptr(0),weight);
}


void green_func_tstp::right_multiply(cplx *f0, cplx *ft, cplx weight){
	int m;
	cplx *x0, *xtmp, *ftmp;
	xtmp = new cplx[element_size_];
	if(tstp_==-1){
		x0 = data_;
		// Multiply the Matsubara
		for(m=0;m<=ntau_;m++){
			element_mult(size1_,xtmp,x0,f0);
			element_smul(size1_,xtmp,weight);
			element_set(size1_,x0,xtmp);
			x0+=element_size_;
		}
	}
	else{
		x0=data_;
		ftmp = ft;
		// Multiply the Retarded
		for(m=0;m<=tstp_;m++){
			element_mult(size1_,xtmp,x0,ftmp);
			element_smul(size1_,xtmp,weight);
			element_set(size1_,x0,xtmp);
			x0+=element_size_;
			ftmp+=element_size_;
		}
		// Multiply the TV
		for(m=0;m<=ntau_;m++){
			element_mult(size1_,xtmp,x0,f0);
			element_smul(size1_,xtmp,weight);
			element_set(size1_,x0,xtmp);
			x0+=element_size_;
		}
		// Multiply the Lesser
		ftmp = ft + tstp_*element_size_;
		for(m=0;m<=tstp_;m++){
			element_mult(size1_,xtmp,x0,ftmp);
			element_smul(size1_,xtmp,weight);
			element_set(size1_,x0,xtmp);
			x0+=element_size_;
		}
	}
	delete[] xtmp;
}


void green_func_tstp::right_multiply(const CFUNC &ft, cplx weight){
	this->right_multiply(ft.ptr(-1),ft.ptr(0),weight);
}

void green_func_tstp::smul(cplx weight){
	int len = (ntau_+1 + (tstp_+1)*2)*element_size_;
	int i;
	for(i=0;i<len;i++){
		data_[i]*=weight;
	}
}


} // Namespace
#endif
