#ifndef VIE2_IMPL
#define VIE2_IMPL

#include "vie2.h"
#include "mat_utils.h"
#include <chrono>

#define PI 3.1415926535897932384626433832795028841971693


namespace NEdyson{

// This solves G^M(\tau) + \int_0^\beta d\tau' F^M(\tau-\tau')G^M(\tau') = Q^M(\tau)
// via fourier iteration
void vie2_mat_fourier(GREEN &G, const GREEN &F, const GREEN &Fcc, const GREEN &Q, double beta, int pcf){
	cplx *fdft, *fiomn, *qiomn, *qiomn1, *qdft, *qmasy, *z1, *zinv, *one;
	cplx expfac, *gmat, iomn;
	int ntau = G.ntau(), m, r, p, l, sig=G.sig(), es=G.element_size(), matsub_one, size1=G.size1();
	double dtau = beta/ntau;
	matsub_one = (sig==-1) ? 1.0 : 0.0;

  assert(G.ntau()==F.ntau());
  assert(G.ntau()==Fcc.ntau());
  assert(G.ntau()==Q.ntau());
  assert(G.size1()==F.size1());
  assert(G.size1()==Fcc.size1());
  assert(G.size1()==Q.size1());
	
	if(ntau%2==1){
		std::cout<<"must have ntau even. aborting..."<<std::endl;
		abort();
	}
	
	fdft = new cplx[(ntau+1)*es];
	fiomn = new cplx[es];
	qiomn = new cplx[es];
	qdft = new cplx[(ntau+1)*es];
	qmasy= new cplx[es];
	z1 = new cplx[es];
	zinv = new cplx[es];
	one = new cplx[es];
	gmat = new cplx[(ntau+1)*es];
	
	element_iden(size1,one);
	memset(gmat,0,sizeof(cplx)*(ntau+1)*es);

	// Calculate the offset for the tail
	// qmasy = -Q(0)+sig*Q(\beta)
	element_set(size1,qmasy,Q.matptr(0));
	element_smul(size1,qmasy,-1.);
	element_incr(size1,qmasy,(cplx)sig,Q.matptr(ntau));

	// Set first order tail to deal with discontinuity
	set_first_order_tail(gmat,qmasy,beta,es,ntau,sig,size1);

	// Get the discrete FT of Q and F
	matsubara_dft(fdft, F, sig);
	matsubara_dft(qdft, Q, sig);
	int m2 = ntau/2;
	for(p=-pcf;p<=pcf;p++){
		int h=0,l=0;
		if(sig==1&&p==0) h=1;
		if(sig==1&&p>0){
			h=1;
			l=-1;
		}
		for(m=-m2-l;m<=m2-1+h;m++){
			iomn=cplx(0,get_omega(m+p*ntau,beta,sig));
			// Cubically correct the dft
			matsubara_ft(fiomn,m+p*ntau,F,fdft,sig,beta);
			matsubara_ft(qiomn,m+p*ntau,Q,qdft,sig,beta);
			
			// z1 = 1/beta*( (1+F)^{-1}Q-1/iomn*qmasy )
			element_set(size1,z1,fiomn);
			element_incr(size1,z1,one);
			element_inverse(size1,zinv,z1);
			element_mult(size1,z1,zinv,qiomn);
			if(!(sig==1&&m+p*ntau==0)) element_incr(size1,z1,-1./iomn,qmasy);
			element_smul(size1,z1,1./beta);

			// add to gmat
			for(r=0;r<=ntau;r++){
				double arg = -r*beta/ntau*get_omega(m+p*ntau,beta,sig);
				expfac = cplx(cos(arg),sin(arg));
				for(l=0;l<es;l++) gmat[r*es+l]+=z1[l]*expfac;
			}
		}
	}

	for(r=0;r<=ntau;r++) element_set(size1,G.matptr(r),gmat+r*es);


	delete[] fdft;
	delete[] fiomn;
	delete[] qiomn;
	delete[] qdft;
	delete[] qmasy;
	delete[] z1;
	delete[] zinv;
	delete[] one;
	delete[] gmat;
}


// Solves (1+F)*G = Q retarded component on the (k+1)x(k+1) square
// \sum_{l=m+1}^k M_{nl}G^R_{lm} = Q^R_{nm}-\sum_{l=0}^m M_{nl}G^R{lm}
// M_{nl} = \delta_{nl}+dtI_{mnl}F^R{nl}
void vie2_ret_start(GREEN &G, const GREEN &F, const GREEN &Fcc, const GREEN &Q, const INTEG &I, double dt){
  assert(G.size1()==F.size1());
  assert(G.size1()==Fcc.size1());
  assert(G.size1()==Q.size1());
  assert(G.nt()>=I.k());
  assert(Q.nt()>=I.k());
  assert(F.nt()>=I.k());
  assert(Fcc.nt()>=I.k());

	int k=I.k(), size1=G.size1(), es=G.element_size(), m,n,l,i;
	cplx *M, *X, *Y, *iden, *tmp, *stmp;
	double dnl,weight;

	M = new cplx[k*k*es];
	X = new cplx[k*es];
	Y = new cplx[k*es];
	iden = new cplx[es];
	tmp = new cplx[es];
	stmp = new cplx[es];

	element_iden(size1,iden);

	// Initial Condition
	for(i=0;i<=k;i++){
		element_set(size1,G.retptr(i,i),Q.retptr(i,i));
	}

	// Iterate through t'=m
	for(m=0;m<k;m++){
		memset(M,0,sizeof(cplx)*k*k*es);
		memset(Y,0,sizeof(cplx)*k*es);
		// Solve for t=n=m+1...k
		for(n=m+1;n<=k;n++){
			for(l=0;l<=k;l++){
				if(l<=m){ // We know these G, so put them into Y
					element_set(size1, tmp, F.retptr(n,l));
					element_smul(size1,tmp, dt*I.poly_integ(m,n,l));
					element_conj(size1,stmp, G.retptr(m,l));
					element_incr(size1,Y+(n-m-1)*es,tmp,stmp);
				}
				else{ // Place these into M
					if(n>=l){ // Use F
						element_set(size1,tmp,F.retptr(n,l));
					}
					else{ // Use Fcc
						element_set(size1,tmp,Fcc.retptr(l,n));
						element_conj(size1,tmp);
						element_smul(size1,tmp,-1);
					}
					dnl=(n==l)?(1.):(0.);
					weight = dt*I.poly_integ(m,n,l);
					for(i=0;i<es;i++){
						M[(n-m-1)*(k-m)*es+(i/size1)*(k-m)*size1+(l-m-1)*size1+i%size1] = dnl*iden[i]+weight*tmp[i];
					}
				}
			}
			// Put Q into Y
			element_incr(size1,Y+(n-m-1)*es,Q.retptr(n,m));
		}
		// Solve MX=Y
		element_linsolve_left((k-m)*size1,(k-m)*size1,size1,M,X,Y);

		//Put result into G
		for(n=0;n<k-m;n++){
			element_set(size1,G.retptr(n+m+1,m),X+n*es);
		}
	}

	delete[] M;
	delete[] X;
	delete[] Y;
	delete[] iden;
	delete[] tmp;
	delete[] stmp;
}



// Solves G^TV(t,\tau) + \int_0^\t dt~ F^R(t,t~)G^TV(t~,\tau) 
// = Q^TV(t,\tau) - \int_0^\beta d\tau' F^TV(t,\tau') G^M(\tau'-\tau)
// for t = 0...k, \tau = 0...ntau
// at each m we have MX=Y
// M_nl = \delta_nl = dt*w_nl*F^R_nl
// Y_n = Q_nm - C_2^TV[F,G](n,m) - C_3^TV[F,G](n,m) - M_n0 G^TV_0m
void vie2_tv_start(GREEN &G, const GREEN &F, const GREEN &Fcc, const GREEN &Q, const INTEG &I, double dt, double beta){
  assert(G.size1()==F.size1());
  assert(G.size1()==Fcc.size1());
  assert(G.size1()==Q.size1());
  assert(G.nt()>=I.k());
  assert(Q.nt()>=I.k());
  assert(F.nt()>=I.k());
  assert(Fcc.nt()>=I.k());
  assert(G.ntau()==F.ntau());
  assert(G.ntau()==Fcc.ntau());
  assert(G.ntau()==Q.ntau());
	
  int k=I.k(), ntau=G.ntau(), size1=G.size1(), es=G.element_size(),m,n,l,i;
	int k1=k+1;
	cplx *M, *Y, *X, *iden, *tmp, *stmp, weight, cplxi=cplx(0.,1.);
	double dnl;

	M = new cplx[k1*k1*es];
	Y = new cplx[k1*es];
	X = new cplx[k1*es];
	iden = new cplx[es];
	tmp = new cplx[es];
	stmp = new cplx[es];
	
	element_iden(size1,iden);

	// Do the integrals and store them in Gtv(n,m)
	for(m=0;m<=ntau;m++){
		for(n=0;n<=k;n++){
			CTV2(I, F, G, n, m, beta, tmp);
			CTV3(I, F, G, n, m, beta, stmp);
			for(i=0;i<es;i++) G.tvptr(n,m)[i] = -tmp[i] - stmp[i];
		}
	}


	// At each m, get values for n=1...k
	for(m=0;m<=ntau;m++){
		memset(M,0,sizeof(cplx)*k1*k1*es);
		memset(Y,0,sizeof(cplx)*k1*es);
		
		// Set up the kxk linear problem
		for(n=0;n<=k;n++){
			for(l=0;l<=k;l++){
				dnl = n==l?1.:0.;
				weight = I.gregory_weights(n,l)*dt;
				if(n>=l){ // use F
					element_set(size1,tmp,F.retptr(n,l));
				}
				else{ // use Fcc
					element_set(size1,stmp,Fcc.retptr(l,n));
					element_conj(size1,tmp,stmp);
					weight*=-1.;
				}
				for(i=0;i<es;i++){
					M[(n)*k1*es+(i/size1)*k1*size1+(l)*size1+i%size1] = dnl*iden[i]+weight*tmp[i];
				}
			}
		
			// add the integrals
			element_incr(size1, Y+(n)*es, G.tvptr(n,m));
			
			// add Q
			element_incr(size1, Y+(n)*es, Q.tvptr(n,m));
		}

		// Solve MX=Y
		element_linsolve_left(k1*size1,k1*size1,size1,M,X,Y);

		// Put result into G
		for(l=0;l<=k;l++){
      element_set(size1, G.tvptr(l,m), X+l*es);
    }
	}

	delete[] M;
	delete[] Y;
	delete[] X;
	delete[] iden;
	delete[] tmp;
	delete[] stmp;

}


void vie2_les_timestep(int n, GREEN &G, const GREEN &F, const GREEN &Fcc, const GREEN &Q, double beta, double dt, const INTEG &I){
  assert(G.size1()==F.size1());
  assert(G.size1()==Fcc.size1());
  assert(G.size1()==Q.size1());
  assert(G.nt()>=n);
  assert(Q.nt()>=n);
  assert(F.nt()>=n);
  assert(Fcc.nt()>=I.k());
  assert(G.ntau()==F.ntau());
  assert(G.ntau()==Fcc.ntau());
  assert(G.ntau()==Q.ntau());
	
  int k=I.k(), ntau=G.ntau(), sig=G.sig(), m,l,i, k1=k+1, size1=G.size1(), es=size1*size1;
	double dtau=beta/ntau, dml;
	cplx *M, *X, *Y, *tmp, *stmp, weight, cplxi, *iden, *Fptr;
	int num = (n>=k)?n:k;
	
	M = new cplx[k1*k1*es];
	X = new cplx[(num+1)*es];
	Y = new cplx[(num+1)*es];
	tmp = new cplx[es];
	stmp = new cplx[es];
	iden = new cplx[es];
	cplxi = cplx(0.,1.);
	element_iden(size1,iden);

	// Integrals go into Y
	memset(Y,0,sizeof(cplx)*(num+1)*es);
	Cles2_tstp(I,F,Fcc,G,G,n,dt,Y);
	Cles3_tstp(I,F,Fcc,G,G,n,beta,Y);
	for(i=0;i<es*(num+1);i++) Y[i]*=-1;

	// Set up k1xk1 linear problem
	for(m=0;m<=k;m++){
		for(l=0;l<=k;l++){
			dml=m==l?1:0;
			weight=dt*I.gregory_weights(m,l);
			if(m>=l){ // use F
				element_set(size1,tmp,F.retptr(m,l));
			}
			else{ // use Fcc
				element_conj(size1,tmp,Fcc.retptr(l,m));
				weight *= -1.;
			}
			for(i=0;i<es;i++) M[m*es*k1+(i/size1)*k1*size1+l*size1+i%size1] = dml*iden[i] + weight*tmp[i];
		}

		// Put Q into Y
		if(n>=m) element_incr(size1,Y+m*es,Q.lesptr(m,n));
		else{
			element_conj(size1,tmp,Q.lesptr(n,m));
			element_incr(size1,Y+m*es,-1,tmp);
		}
	}


	// Solve MX=Y
	element_linsolve_left(k1*size1,k1*size1,size1,M,X,Y);

	
	// Timestepping
	for(m=k+1;m<=n;m++){
		// Put Q into Y
		element_incr(size1, Y+m*es, Q.lesptr(m,n));
	
		// Put Ret Les integral into Y
		for(l=0;l<m;l++){
			element_incr(size1,Y+m*es, -dt*I.gregory_weights(m,l), F.retptr(m,l), X+l*es);
		}

		// set up M
		weight=dt*I.gregory_weights(m,m);
		Fptr = F.retptr(m,m);
		for(i=0;i<es;i++) M[i] = iden[i]+weight*Fptr[i];

		// Solve MX=Y
		element_linsolve_left(size1,size1,size1,M,X+m*es,Y+m*es);
	}


	// Put into G
	for(i=0;i<=n;i++) element_set(size1,G.lesptr(i,n),X+i*es);



	delete[] M;
	delete[] X;
	delete[] Y;
	delete[] tmp;
	delete[] stmp;
	delete[] iden;
}



void vie2_les_start(GREEN &G, const GREEN &F, const GREEN &Fcc, const GREEN &Q, const INTEG &I, double dt, double beta){
  assert(G.size1()==F.size1());
  assert(G.size1()==Fcc.size1());
  assert(G.size1()==Q.size1());
  assert(G.nt()>=I.k());
  assert(Q.nt()>=I.k());
  assert(F.nt()>=I.k());
  assert(Fcc.nt()>=I.k());
  assert(G.ntau()==F.ntau());
  assert(G.ntau()==Fcc.ntau());
  assert(G.ntau()==Q.ntau());

  int k=I.k();
	for(int i=0;i<=k;i++){
		vie2_les_timestep(i,G,F,Fcc,Q,beta,dt,I);
	}
}


void vie2_start(GREEN &TPP, const GREEN &PhixU, const GREEN &UxPhi, const GREEN &Phi, double beta, double dt, const INTEG &I){
  assert(TPP.size1()==PhixU.size1());
  assert(TPP.size1()==UxPhi.size1());
  assert(TPP.size1()==Phi.size1());
  assert(TPP.nt()>=I.k());
  assert(PhixU.nt()>=I.k());
  assert(UxPhi.nt()>=I.k());
  assert(Phi.nt()>=I.k());
  assert(TPP.ntau()==Phi.ntau());
  assert(TPP.ntau()==UxPhi.ntau());
  assert(TPP.ntau()==PhixU.ntau());
	
  vie2_ret_start(TPP, PhixU, UxPhi, Phi, I, dt);
	vie2_tv_start( TPP, PhixU, UxPhi, Phi, I, dt, beta);
	vie2_les_start(TPP, PhixU, UxPhi, Phi, I, dt, beta);
}


// same as dyson ret start
// solves G^R(n,n-t') + \int_0^{t'} G^R(n,n-t~) F^R(n-t~,n-t') = Q^R(n,n-t')
// for t' = 0...n
void vie2_ret_timestep(int n, GREEN &G, const GREEN &F, const GREEN &Fcc, const GREEN &Q, const INTEG &I, double dt){
  assert(G.size1()==F.size1());
  assert(G.size1()==Fcc.size1());
  assert(G.size1()==Q.size1());
  assert(n>I.k());
  assert(G.nt()>=n);
  assert(Q.nt()>=n);
  assert(F.nt()>=n);
  assert(Fcc.nt()>=I.k());
	
	int k=I.k(), size1=G.size1(), es=size1*size1, m,i,l,j;
	cplx *M, *Y, *X, *tmp, *stmp, *iden, cplxi=cplx(0.,1.);
	double dml, weight;

	M = new cplx[k*k*es];
	Y = new cplx[(n+1)*es];
	X = new cplx[k*es];
	tmp = new cplx[es];
	stmp = new cplx[es];
	iden = new cplx[es];
	element_iden(size1, iden);
	// Initial value
	element_set(size1,G.retptr(n,n),Q.retptr(n,n));
	
	memset(M,0,k*k*es*sizeof(cplx));
	memset(Y,0,(n+1)*es*sizeof(cplx));
	
	// Start up routine
	for(m=1;m<=k;m++){
		element_set(size1,tmp,Q.retptr(n,n-m));
		element_incr(size1,tmp,-dt*I.gregory_weights(m,0),G.retptr(n,n),F.retptr(n,n-m));
		for(i=0;i<es;i++) Y[(i/size1)*size1*k+(m-1)*size1+i%size1]=tmp[i];

		for(l=1;l<=k;l++){
			dml=m==l?1.:0.;
			weight = dt*I.gregory_weights(m,l);
			if(n-l>=n-m) element_set(size1,tmp,F.retptr(n-l,n-m));
			else{
				element_conj(size1,tmp,Fcc.retptr(n-m,n-l));
				weight*=-1;
			}
			for(i=0;i<es;i++) M[(l-1)*k*es+(i/size1)*size1*k+(m-1)*size1+i%size1] = dml*iden[i]+weight*tmp[i];
		}
	}
	
	element_linsolve_right(size1,k*size1,k*size1,X,M,Y);

	for(m=1;m<=k;m++){
		for(i=0;i<es;i++) G.retptr(n,n-m)[i] = X[(i/size1)*size1*k+(m-1)*size1+i%size1];
	}

	// Do the integrals
	for(m=k+1;m<=n;m++){
		for(l=0;l<=k;l++){
			element_incr(size1,Y+m*es,I.gregory_weights(m,l),G.retptr(n,n-l),F.retptr(n-l,n-m));
		}
	}

	// Do the rest of the timesteps
	cplx *yptr, *gptr;
	for(m=k+1;m<=n;m++){
		weight=dt*I.gregory_weights(m,m);
		// Set up X[1+dt*w_{mm}*F^R_{n-m n-m}] = Q^R_{n,n-m}-dt*\sum_l=0^{m-1} G^R_{n n-l} F^R_{n-l,n-m}
		element_smul(size1,Y+m*es,-dt);
		element_incr(size1,Y+m*es,Q.retptr(n,n-m));
		element_set(size1,M,iden);
		element_incr(size1,M,weight,F.retptr(n-m,n-m));

		// Solve
		element_linsolve_right(size1,size1,size1,G.retptr(n,n-m),M,Y+m*es);
		
		// Add into remaining integrals
		gptr=G.retptr(n,n-m);
		for(l=m+1;l<=n;l++) element_incr(size1,Y+l*es,I.gregory_weights(l,m),gptr,F.retptr(n-m,n-l));
	}


	delete[] M;
	delete[] Y;
	delete[] X;
	delete[] tmp;
	delete[] stmp;
	delete[] iden;
}


// Same as dyson tv stepping routine
void vie2_tv_timestep(int n, GREEN &G, const GREEN &F, const GREEN &Fcc, const GREEN &Phi, const INTEG &I, double dt, double beta){
  assert(G.size1()==F.size1());
  assert(G.size1()==Fcc.size1());
  assert(G.size1()==Phi.size1());
  assert(n>I.k());
  assert(G.nt()>=n);
  assert(Phi.nt()>=n);
  assert(F.nt()>=n);
  assert(Fcc.nt()>=I.k());
  assert(G.ntau()==F.ntau());
  assert(G.ntau()==Fcc.ntau());
  assert(G.ntau()==Phi.ntau());

	int ntau=G.ntau(), k=I.k(), m,i, size1=G.size1(), es=size1*size1;
	cplx *M, *Y, *iden, weight, *tmp, *stmp, *Gptr, *Pptr;
	double dtau=beta/ntau;

	M = new cplx[es];
	Y = new cplx[(ntau+1)*es];
	iden = new cplx[es];
	tmp = new cplx[es];
	stmp = new cplx[es];
	element_iden(size1,iden);
/*
	// First do the FTV GM integral.  Put into GTV(tstp, m)
	for(m=0;m<=ntau;m++){
		CTV2(I,F,G,n,m,beta,tmp);
		CTV3(I,F,G,n,m,beta,stmp);
		Gptr=G.tvptr(n,m);
		for(i=0;i<es;i++) Gptr[i] = tmp[i]+stmp[i];
	}

	// Next do the FR GTV integral. Adds to GTV(tstp,m) for all m
	CTV1(I,F,Fcc,G,n,dt);
*/
  Ctv_tstp(n, G, F, Fcc, G, G, I, beta, dt);
  
	weight = dt*I.omega(0);
	// Do the solving
	for(m=0;m<=ntau;m++){
		// Y = Phi - integrals
		Gptr = G.tvptr(n,m);
		Pptr = Phi.tvptr(n,m);
		for(i=0;i<es;i++){
			Y[i] = Pptr[i]-Gptr[i];
		}

		// M = I + dt w_{nn} F^R_{nn}
		Pptr = F.retptr(n,n);
		for(i=0;i<es;i++){
			M[i] = iden[i] + weight*Pptr[i];
		}

		// Solve MX=Y
		element_linsolve_left(size1,size1,size1,M,Gptr,Y);
	}


	delete[] M;
	delete[] Y;
	delete[] iden;
	delete[] tmp;
	delete[] stmp;
}



void vie2_timestep(int tstp, GREEN &TPP, const GREEN &PhixU, const GREEN &UxPhi, const GREEN &Phi, double beta, double dt, const INTEG &I){
  assert(TPP.size1()==PhixU.size1());
  assert(TPP.size1()==UxPhi.size1());
  assert(TPP.size1()==Phi.size1());
  assert(tstp>I.k());
  assert(TPP.nt()>=tstp);
  assert(Phi.nt()>=tstp);
  assert(PhixU.nt()>=tstp);
  assert(UxPhi.nt()>=I.k());
  assert(TPP.ntau()==PhixU.ntau());
  assert(TPP.ntau()==UxPhi.ntau());
  assert(TPP.ntau()==Phi.ntau());

  vie2_ret_timestep(tstp, TPP, UxPhi, PhixU, Phi, I, dt);
	vie2_tv_timestep(tstp, TPP, PhixU, UxPhi, Phi, I, dt, beta);
	vie2_les_timestep(tstp, TPP, PhixU, UxPhi, Phi, beta, dt, I);
}



}//namespace
#endif
