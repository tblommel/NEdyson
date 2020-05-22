#ifndef FUNCTION_IMPL
#define FUNCTION_IMPL
	
#include "function.h"
namespace NEdyson{


function::function(){
	data_=0;
	nt_=-2;
	size1_=0;
	element_size_=0;
}


function::~function(){
	delete[] data_;
}


function::function(int nt, int size1){
	assert(size1>0 && nt >=-1);
	size1_=size1;
	nt_=nt;
	element_size_=size1*size1;
	data_=new cplx[(nt+2)*element_size_];
}


function::function(const function &f){
	nt_=f.nt_;
	size1_=f.size1_;
	element_size_=size1_*size1_;
	int len=(nt_+2)*element_size_;
	if(len>0){
		data_=new cplx[len];
		memcpy(data_,f.data_, sizeof(cplx)*len);
	}
	else data_=0;
}


function &function::operator=(const function &f){
	if(this == &f) return *this;
	delete[] data_;
	nt_=f.nt_;
	size1_=f.size1_;
	element_size_=size1_*size1_;
	int len=(nt_+2)*element_size_;
	if(len>0){
		data_=new cplx[len];
		memcpy(data_,f.data_, sizeof(cplx)*len);
	}
	else data_=0;
	return *this;
}


void function::resize(int nt, int size1){
	assert(nt>=-1 && size1>0);
	delete[] data_;
	int len=(nt+2)*size1*size1;
	if(len == 0) data_=0;
	else data_ = new cplx[len];

	size1_=size1;
	element_size_=size1*size1;
	nt_=nt;
}


void function::set_Zero(void){
	int len=element_size_*(nt_+2);
	if(len>0) memset(data_,0,sizeof(cplx)*len);
}



void function::set_value(int tstp, ZMatrix M){
	assert(M.rows()==M.cols() && M.rows()==size1_ && tstp<=nt_ && tstp >= -1);
	cplx *ft=ptr(tstp);
	ZMatrixMap(ft,size1_,size1_)=M;
}



void function::set_value(int tstp, cplx x){
	assert(tstp<=nt_ && tstp >=-1 && size1_==1);
	*ptr(tstp)=x;
}


void function::set_value(int tstp, cplx *x){
	assert(tstp<=nt_ && tstp >=-1);
	memcpy(ptr(tstp),x,sizeof(cplx)*element_size_);	
}



void function::set_constant(ZMatrix M){
	assert(M.rows()==M.cols() && M.rows()==size1_);
	for(int i=-1;i<=nt_;i++) set_value(i,M);
}


void function::set_constant(cplx x){
  assert(size1_==1);
	for(int i=-1;i<=nt_;i++) set_value(i,x);
}


void function::set_constant(cplx *x){
	for(int i=-1;i<=nt_;i++) set_value(i,x);
}


void function::get_value(int tstp, ZMatrix &M) const {
	assert(M.rows()==M.cols() && M.rows()==size1_ && tstp<=nt_ && tstp >= -1);
	cplx *ft=ptr(tstp);
	M=ZMatrixMap(ft,size1_,size1_);
}


void function::get_value(int tstp, cplx &x) const {
	assert(tstp<=nt_ && tstp >=-1 && size1_==1);
	x=*ptr(tstp);
}


void function::get_value(int tstp, cplx *x) const {
	assert(tstp<=nt_ && tstp >=-1);
	memcpy(x,ptr(tstp),sizeof(cplx)*element_size_);	
}


void function::printf(void) const {
	cplx *x=0;
	for(int i=-1;i<=nt_;i++){
		x=ptr(i);
		for(int j=0;j<element_size_;j++){
			std::cout<<x[j]<<std::endl;
		}
		std::cout<<std::endl;
	}
}




} //namespace func
#endif// FUNCTION_IMPL
