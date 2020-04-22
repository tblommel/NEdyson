#ifndef SPECTRAL_IMPL
#define SPECTRAL_IMPL

#define PI 3.1415926535897932384626433832795028841971693

#include "spectral.hpp"
namespace NEdyson{

spectral::spectral(){
	nt_=0;
	nw_=0;
	size1_=0;
	es_=0;
	wmax_=0;
	dw_=0;
	dt_=0;
	A_=0;
}

spectral::~spectral(){
	delete[] A_;
}

spectral::spectral(int nt, int nw, int size1, double wmax, double dt){
	nt_ = nt;
	nw_ = nw;
	size1_ = size1;
	es_ = size1*size1;
	wmax_ = wmax;
	dw_ = 2*wmax/(nw-1);
	dt_ = dt;
	A_ = new double[(nt_+1)*nw_*es_];
}

spectral::spectral(const spectral &A){
	nt_ = A.nt_;
	nw_ = A.nw_;
	size1_ = A.size1_;
	es_ = A.es_;
	wmax_ = A.wmax_;
	dw_ = A.dw_;
	dt_ = A.dt_;
	memcpy(A_,A.A_,sizeof(double)*(nt_+1)*nw_*es_);
}

spectral& spectral::operator=(const spectral &A){
	if(this == &A) return *this;
	wmax_ = A.wmax_;
	dw_ = A.dw_;
	dt_=A.dt_;
	if(nt_ != A.nt_ || nw_ != A.nw_ || es_ != A.es_){
		delete[] A_;
		A_ = new double[(nt_+1)*nw_*es_];	
	}
	memcpy(A_,A.A_,sizeof(double)*(nt_+1)*nw_*es_);
	return *this;
}

// Get pointer
double *spectral::ptr(int t, int w){
	return A_ + (t*nw_ + w)*es_;
}

// Input/Output
void spectral::print_to_file(std::string file, int precision){
	int i,j,l;
	double *x = NULL;
	std::ofstream out;
	std::string actfile;
	actfile = "";
	actfile += file;
	actfile += "_A.dat";
	out.open(actfile);
	out.precision(precision);
	
	out << nt_ << " " << nw_ << " " << size1_ << " " << wmax_ << " " << dw_ << " " << dt_ << std::endl;

	for(i=0;i<=nt_;i++){
		for(j=0;j<nw_;j++){
			x = ptr(i,j);
			for(l=0;l<es_;l++) out << x[l] <<" ";
			out<<std::endl;
		}
	}
	out.close();
}

void spectral::read_from_file(const char *file){
	int i,j,l;
	double *x=NULL;

	std::ifstream in;
	std::string actfile;
	actfile = "";
	actfile += file;
	actfile += "_A.dat";
	in.open(actfile);

	if(!(in >> nt_ >> nw_ >> size1_ >> wmax_ >> dw_ >> dt_)){
		std::cerr << "Error in reading spectral function from file " << actfile << std::endl;
		abort();
	}

	delete[] A_;
	A_ = new double[(nt_+1)*nw_*es_];

	for(i =0;i<=nt_;i++){
		for(j=0;j<nw_;j++){
			x = ptr(i,j);
			for(l=0;l<es_;l++){
				if(!(in>>x[l])){
					std::cerr<<"Error in reading spectral function from file " << actfile<< std::endl;
				}
			}
		}
	}
	in.close();
}

void spectral::AfromG(green_func &G, int nw, double wmax, double dt){
	nt_ = 2*G.nt();
	nw_ = nw;
	size1_ = G.size1();
	es_ = size1_*size1_;
	wmax_ = wmax;
	dw_ = 2*wmax/(nw-1);
	dt_ = dt;
	delete[] A_;
	A_ = new double[(nt_+1)*nw_*es_];


	int ita, iw, itr;
	cplx *tmp;
	cplx weight;
	
	tmp = new cplx[es_];
	cplx arg, cplxi = cplx(0.,1.);
	for(ita=0;ita<=nt_;ita++){
		int max1 = ita;
		int max2 = nt_-ita;	
		int max = max1>=max2?max2:max1;
		for(itr=ita%2;itr<=max;itr+=2){
			element_set(size1_,tmp,G.retptr((ita+itr)/2,(ita-itr)/2));
			for(iw=0;iw<nw;iw++){
				arg = (iw-(nw-1)/2)*dw_*itr*dt_*cplxi;
				weight = std::exp(arg);
				element_imag_incr(size1_,ptr(ita,iw),weight,tmp);
			}
		}
	}

	double a = -2*dt_/PI;
	for(ita=0;ita<=nt_;ita++){
		for(iw=0;iw<nw;iw++){
			element_smul(size1_,ptr(ita,iw),a);
		}
	}

	delete[] tmp;
}


void spectral::AfromG(green_func &G, int nw, double wmax, double dt, cplx *extdata, int nfit, int ntp){
	nt_ = 2*G.nt();
	nw_ = nw;
	size1_ = G.size1();
	es_ = size1_*size1_;
	wmax_ = wmax;
	dw_ = 2*wmax/(nw-1);
	dt_ = dt;
	delete[] A_;
	A_ = new double[(nt_+1)*nw_*es_];
	memset(A_,0,sizeof(double)*(nt_+1)*nw_*es_);

	int ita, iw, itr;
	cplx *tmp;
	cplx weight;
	
	tmp = new cplx[es_];
	cplx arg, cplxi = cplx(0.,1.);
	for(ita=0;ita<=nt_;ita++){
		int max1 = ita;
		int max2 = nt_-ita;	
		int max = max1>=max2?max2:max1;
		for(itr=ita%2;itr<=max;itr+=2){
			element_set(size1_,tmp,G.retptr((ita+itr)/2,(ita-itr)/2));
			for(iw=0;iw<nw;iw++){
				arg = (iw-(nw-1)/2)*dw_*itr*dt_*cplxi;
				weight = std::exp(arg);
				element_imag_incr(size1_,ptr(ita,iw),weight,tmp);
			}
		}
	}
	
	int t, tp, Gnt=G.nt(), tpmax=Gnt-nfit-ntp;
	for(tp=0;tp<tpmax;tp++){
		for(t=0;t<Gnt-tp;t++){
			element_set(size1_,tmp,extdata+tp*Gnt*es_+t*es_);
			for(iw=0;iw<nw;iw++){
				arg = (iw-(nw-1)/2)*dw_*(Gnt+1+t-tp)*dt_*cplxi;
				weight = std::exp(arg);
				element_imag_incr(size1_,ptr(Gnt+1+t+tp,iw),weight,tmp);
			}
		}
	}

	double a = -2*dt_/PI;
	for(ita=0;ita<=nt_;ita++){
		for(iw=0;iw<nw;iw++){
			element_smul(size1_,ptr(ita,iw),a);
		}
	}

	delete[] tmp;
}






} // namespace
	
#endif
