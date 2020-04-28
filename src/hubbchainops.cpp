#ifndef HUBB_CHAIN_IMPL
#define HUBB_CHAIN_IMPL

#include "hubbchainops.h"

namespace Hubb{

cdmatrix EWMult(cdmatrix X1, cdmatrix X2, cdmatrix X3){
	return X1.cwiseProduct(X2.cwiseProduct(X3));
}

// C_{c1,c2}(tau) = -A_{a1,a2}(tau)*B_{b2,b1}(-tau)
//                = -sigb* A_{a1,a2}(tau)*B_{b2,b1}(beta-tau)
void get_bubble_1_mat(cplx* cmat, int sc, int c1, int c2, const cplx *amat, int sa, int a1, int a2, const cplx *bmat, int sb, int b1, int b2, int sigb, int ntau){
	int m;
	int a12 = a1*sa+a2, sa2=sa*sa;
	int b21 = b2*sb+b1, sb2=sb*sb;
	int c12 = c1*sc+c2, sc2=sc*sc;
	double sig = -1.0*sigb;
	for(m=0;m<=ntau;m++) cmat[c12 + m*sc2]=sig*amat[a12+m*sa2]*bmat[b21+(ntau-m)*sb2];
}


void get_bubble_1_tstp(int tstp, cplx *cret, cplx *ctv, cplx *cles, int sc, int c1, int c2,
																 const cplx *aret, const cplx *atv, const cplx *ales, int sa, int a1, int a2,
																 const cplx *bret, const cplx *btv, const cplx *bles, int sb, int b1, int b2, int bsig, int ntau){
	int m;
//	int a12 = a1*sa+a2, a21 = a21*sa+a1, sa2 = sa*sa;
	int b12 = b1*sb+b2, b21 = b21*sb+b1, sb2 = sb*sb;
//	int c12 = c1*sc+c2, sc2 = sc*sc;
	int count12 = a1*sa+a2;
	int count21 = a2*sa+a1;
	int size2 = sa*sa;
	cplx tmp;
	cplx I=cplx(0.0,1.0);

	// First do tv
	// C^{tv}_{i,j}(t,tau) = i*A^{tv}_{ij}(t,tau)*B^{vt}_{ji}(tau,t)
	for(m=0;m<=ntau;m++){
		tmp = -(double)bsig*I*std::conj(btv[b12 + (ntau-m)*sb2]);
		ctv[count12]=tmp*atv[count12];
		count12 += size2;
	}

	// Next do Ret and Les
	count12 = a1*sa+a2;
	// C^R_{i,j}(t,t') = i*A^R_{i,j}(t,t')B^<_{j,i}(t',t)+i*A^<_{i,j}(t,t')B^A(t',t)
	for(m=0;m<=tstp;m++){
		cret[count12] = I*(aret[count12]*bles[count21]-std::conj(ales[count21])*std::conj(bret[count12]));
		cles[count12] = I*(ales[count12]*(bret[count21]-std::conj(bles[count12])));
		count12+=size2;
		count21+=size2;
	}
}


void get_bubble_2_mat(cplx *smat, int ss,int s1, int s2, const cplx *gmat, int sg, int g1, int g2, const cplx *pmat,int sp,int p1,int p2,int ntau){
  int m;
	int es = ss*ss;
	int count12 = s1*ss+s2;
	for(m=0;m<=ntau;m++){
		smat[count12]=-gmat[count12]*pmat[count12];
		count12+=es;
	}
}


void get_bubble_2_tstp(int tstp, cplx *sret, cplx *stv, cplx *sles, int ss, int s1, int s2,
																 const cplx *gret, const cplx *gtv, const cplx *gles, int gs, int g1, int g2,
																 const cplx *pret, const cplx *ptv, const cplx *ples, int ps, int p1, int p2, int ntau){
	int m;
	int es = ss*ss;
	int count12 = s1*ss+s2;
	int count21 = s2*ss+s1;
	cplx I = cplx(0.0,1.0);
	// First do TV
	for(m=0;m<=ntau;m++){
		stv[count12] = I*gtv[count12]*ptv[count12];
		count12+=es;
	}

	// Next do Ret and Les
	count12 = s1*ss+s2;
	for(m=0;m<=tstp;m++){
		sret[count12]=I*(gret[count12]*pret[count12]-std::conj(gles[count21])*pret[count12]-gret[count12]*std::conj(ples[count21]));
		sles[count12]=I*(ples[count12]*gles[count12]);
		count12+=es;
		count21+=es;
	}
}


void Bubble2(int tstp, GREEN &Sigma, int s1, int s2, const GREEN &G, int g1, int g2, const GREEN_TSTP &Pol, int p1, int p2){
	int ntau = Pol.ntau();
	if(tstp == -1) get_bubble_2_mat(Sigma.matptr(0),Sigma.size1(),s1,s2,G.matptr(0),G.size1(),g1,g2,Pol.matptr(0),Pol.size1(),p1,p2,ntau);
	else get_bubble_2_tstp(tstp,Sigma.retptr(tstp,0),Sigma.tvptr(tstp,0),Sigma.lesptr(0,tstp),Sigma.size1(),s1,s2,
															G.retptr(tstp,0),G.tvptr(tstp,0),G.lesptr(0,tstp),G.size1(),g1,g2,
															Pol.retptr(0),Pol.tvptr(0),Pol.lesptr(0),Pol.size1(),p1,p2,ntau);

}


void Bubble2(int tstp, GREEN &Sigma, int s1, int s2, const GREEN &G, int g1, int g2, const GREEN &Pol, int p1, int p2){
	int ntau = G.ntau();
	if(tstp == -1) get_bubble_2_mat(Sigma.matptr(0),Sigma.size1(),s1,s2,G.matptr(0),G.size1(),g1,g2,Pol.matptr(0),Pol.size1(),p1,p2,ntau);
	else get_bubble_2_tstp(tstp,Sigma.retptr(tstp,0),Sigma.tvptr(tstp,0),Sigma.lesptr(0,tstp),Sigma.size1(),s1,s2,
															G.retptr(tstp,0),G.tvptr(tstp,0),G.lesptr(0,tstp),G.size1(),g1,g2,
															Pol.retptr(tstp,0),Pol.tvptr(tstp,0),Pol.lesptr(0,tstp),Pol.size1(),p1,p2,ntau);
}



void Bubble1(int tstp, GREEN_TSTP &Pol, int p1, int p2, const GREEN &A, int a1, int a2, const GREEN &B, int b1, int b2){
	int ntau = Pol.ntau();
	if(tstp == -1) get_bubble_1_mat(Pol.matptr(0),Pol.size1(),p1,p2,A.matptr(0),A.size1(),a1,a2,B.matptr(0),B.size1(),b1,b2,B.sig(), ntau);
	else get_bubble_1_tstp(tstp, Pol.retptr(0),Pol.tvptr(0),Pol.lesptr(0),Pol.size1(),p1,p2,
															 A.retptr(tstp,0),A.tvptr(tstp,0),A.lesptr(0,tstp),A.size1(),a1,a2,
															 B.retptr(tstp,0),B.tvptr(tstp,0),B.lesptr(0,tstp),B.size1(),b1,b2, B.sig(),ntau);
}


void Bubble1(int tstp, GREEN &Pol, int p1, int p2, const GREEN_TSTP &A, int a1, int a2, const GREEN &B, int b1, int b2){
	int ntau = Pol.ntau();
	if(tstp == -1) get_bubble_1_mat(Pol.matptr(0),Pol.size1(),p1,p2,A.matptr(0),A.size1(),a1,a2,B.matptr(0),B.size1(),b1,b2,B.sig(), ntau);
	else get_bubble_1_tstp(tstp, Pol.retptr(tstp,0),Pol.tvptr(tstp,0),Pol.lesptr(0,tstp),Pol.size1(),p1,p2,
															 A.retptr(0),A.tvptr(0),A.lesptr(0),A.size1(),a1,a2,
															 B.retptr(tstp,0),B.tvptr(tstp,0),B.lesptr(0,tstp),B.size1(),b1,b2, B.sig(),ntau);
}


// P_ij(t,t') = -iG_ij(t,t')G_ji(t',t)
void Polarization(int tstp, const GREEN &G, GREEN_TSTP &Pol){
	int nsites=G.size1(),i,j;
	for(i=0;i<nsites;i++){
		for(j=0;j<nsites;j++){
			Bubble1(tstp,Pol,i,j,G,i,j,G,i,j);
		}
	}
	Pol.smul(-1.0);
}


// HMF_{ij}(t) = h0_{ij} + delta_{ij} U(i,i)*rho(i,i)
void Ham_MF(int tstp, const GREEN &G, const CFUNC &Ut, const cdmatrix &h0, CFUNC &hmf){
	int nsites=G.size1();
	cdmatrix rho(nsites,nsites), U(nsites,nsites), hmft(nsites,nsites);
	hmft=h0;
	G.get_dm(tstp,rho);
	Ut.get_value(tstp,U);
	for(int i=0;i<nsites;i++) hmft(i,i) += U(i,i)*rho(i,i);

	hmf.set_value(tstp,hmft);
}


void Sigma_2B(int tstp, const GREEN &G, const CFUNC &Ut, GREEN &Sigma){
	int nsites = G.size1();
	int ntau = G.ntau(),i,j;
	GREEN_TSTP Pol(tstp, ntau, nsites,1);
	
	Polarization(tstp, G, Pol);
	Pol.right_multiply(Ut);
	Pol.left_multiply(Ut);

	for(i=0;i<nsites;i++){
		for(j=0;j<nsites;j++){
			Bubble2(tstp, Sigma,i,j,G,i,j,Pol,i,j);
		}
	}
}


void GenTPP(int tstp, double dt, double beta, const GREEN &G, GREEN &Phi, const CFUNC &Ut, GREEN &PhixU, GREEN &UxPhi, GREEN &TPP, const INTEG &I){
	int size1=G.size1(),n1,n2;

	for(n1=0;n1<size1;n1++){
		for(n2=0;n2<size1;n2++){
			Bubble2(tstp,Phi,n1,n2,G,n1,n2,G,n1,n2);
		}
	}
	Phi.smul(tstp,-1);
	
	UxPhi.set_tstp(tstp,Phi);
	PhixU.set_tstp(tstp,Phi);

	UxPhi.left_multiply(tstp,Ut);
	PhixU.right_multiply(tstp,Ut);

	if(tstp==-1)	NEdyson::vie2_mat_fourier(TPP, PhixU, UxPhi, Phi, beta);
	else NEdyson::vie2_timestep(tstp, TPP, PhixU, UxPhi, Phi, beta, dt, I);
}


void GenTPP(double dt, double beta, const GREEN &G, GREEN &Phi, const CFUNC &Ut, GREEN &PhixU, GREEN &UxPhi, GREEN &TPP, const INTEG &I){
	int size1=G.size1(),n1,n2,k=I.k();
	
	for(int tstp=0;tstp<=k;tstp++){
		for(n1=0;n1<size1;n1++){
			for(n2=0;n2<size1;n2++){
				Bubble2(tstp,Phi,n1,n2,G,n1,n2,G,n1,n2);
			}
		}
		Phi.smul(tstp,-1);

		UxPhi.set_tstp(tstp,Phi);
		PhixU.set_tstp(tstp,Phi);

		UxPhi.left_multiply(tstp,Ut);
		PhixU.right_multiply(tstp,Ut);
	}
	NEdyson::vie2_start(TPP, PhixU, UxPhi, Phi, beta, dt, I);	

}


void Sigma_TPP(int tstp, const GREEN &G, const CFUNC &Ut, const GREEN &TPP, GREEN &Sigma){
	int ntau = G.ntau();
	int size1 = G.size1(),n1,n2;

	GREEN_TSTP UTU(tstp, ntau, size1, 1);

	TPP.get_tstp(tstp,UTU);
	UTU.right_multiply(Ut);
	UTU.left_multiply(Ut);

	for(n1=0;n1<size1;n1++){
		for(n2=0;n2<size1;n2++){
			Bubble1(tstp,Sigma,n1,n2,UTU,n1,n2,G,n1,n2);
		}
	}
}
}

#endif
