#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "hydro_math.h"


const ihf ihf_default = { .calculated = 0};
const mf mf_default = { .calculated = 0, .part1=0, .part2=0, .part3=0};

Hydrodynamic_Function::Hydrodynamic_Function(){

};
Hydrodynamic_Function::~Hydrodynamic_Function(){

};
void Hydrodynamic_Function::alloc_mem() {
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
	ihf_all = (ihf **) malloc( sizeof(ihf*) * kn*xn );
	hf_all = (double *) malloc( sizeof(double) *kn );
	mf_all = (mf *) malloc( sizeof(mf) *xn );

	int i =0 ;
	while ( i < (kn*xn) ) {
		ihf_all[i] = NULL;	
		i++;
	};
	

	i=0; 
	while ( i < kn ) {
		hf_all[i] = 0;	
		i++;

	};

	i=0; 
	while ( i < xn ) {
		mf_all[i] = mf_default; 
		i++;
	};

	
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
	for (int xi = 0; xi < xn; xi++) {
		for (int ki = 0; ki < kn; ki++) {
			i = xi*ki;
			if ( ihf_all[i] == NULL) {
				ihf_all[i] = (ihf *) malloc( sizeof(ihf) ) ;
				ihf_all[i]->calculated = 0 ;
			}
			
		}
	}
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
}
void Hydrodynamic_Function::init(int _kn, int _xn, int _N,double _r1, double  _k1, double _a,  double _V,  double* _g000, double* _h000) {
	if ( ihf_all != NULL) {
		puts("call init double!! Error occur!!");
	}
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
	kn=_kn;
	xn=_xn;
	r1=_r1;
	k1=_k1;
	g000=_g000;
	h000=_h000;
	a=_a;
	N=_N;
	V=_V;
	alloc_mem();
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
}


void Hydrodynamic_Function::prepare_calculation() {
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
	double kr1= r1*k1;
	double kr_i;

	int i, previ;
	ihf ihf_0 = { .calculated = 1, 
		.x= 0,.x2= 0, .x3=0,
		.invx3= 0, .cosx=1, .cosx_x= 0,
		.sinx = 0, .sinx_3_minus_x2 = 0,
		.h1 = 1./3.,  .h2= 0 };

	copy_ihf( ihf_all[0], &ihf_0);
	set_ihf(  ihf_all[1], kr1);

	
	
	double prevcos, prevsin, cosi, sini, cos1, sin1;
	pihf pihf_i, pihf_previ;
	for (int xi = 1; xi < xn; xi++) {
		int ki=0; 
		{
			i = xi*ki; 
		}
		ki = 1;
		{
			i = xi*ki; 
			pihf_i = ihf_all[i] ;
			cos1 = pihf_i->cosx;
			sin1 = pihf_i->sinx;
//			printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
			printf("i=%3d,k=%3d,ik=%5d : %.5lf    %.5lf\n", xi,ki,i,cos1, sin1);
		}
		for (int ki = 2; ki < kn; ki++) {
			i = xi*ki; 
			pihf_i = ihf_all[i] ;
			if ( pihf_i != NULL  && pihf_i->calculated == 0 ) {
/* 				previ = xi*(ki-1);
 * 				pihf_previ = ihf_all[previ] ;
 * 				prevcos = pihf_previ->cosx; 
 * 				prevsin = pihf_previ->sinx; 
 * 				cosi = prevcos*cos1 - prevsin*sin1;
 * 				sini = prevsin*cos1 + sin1*prevcos;
 * 				kr_i = i * kr1;
 * 				set_ihf_sinx_cosx(pihf_i , kr_i, sini,cosi);
 */
				kr_i = i * kr1;
				set_ihf(pihf_i, kr_i);
				cosi = pihf_i->cosx;
				sini = pihf_i->sinx;
				printf("i=%3d,k=%3d,ik=%5d : %.3lf    %.3lf\n", xi,ki,i,cosi, sini);
			}
		} //for int ki 
	} // for int xi 


	double rho0 = N/V;
	double int_dr = r1;
	double Cpart1 = 6.*M_PI*rho0* int_dr  , Cpart2 =4.*M_PI*rho0*int_dr;
	double Spart = 4./3.*M_PI* rho0 * int_dr;;

	D_s = 1;
	

	for (int xi = 1; xi < xn; xi++) {
		double ri = xi*r1;
		double r2    =  ri*ri;
		double ainvr = a/ ri;
		double ainvr2 = ainvr*ainvr;
		double ainvr3 = ainvr2*ainvr;
		double ainvr4 = ainvr3*ainvr;
		double ainvr6 = ainvr4*ainvr2;
		double ainvr7 = ainvr6*ainvr;
		double Ac_star  = -ainvr3 + 75./4. * ainvr7;
		double Bc_star  = .5*ainvr3;
		double As = -15./4. * ainvr4 + 11./2. * ainvr6;
		double Bs = -17./16. * ainvr6;

		mf_all[xi].part1 = Cpart1* h000[xi] * r2  * ainvr;
		mf_all[xi].part2 = Cpart2* g000[xi] * r2  * (Ac_star + 2.*Bc_star);
		mf_all[xi].part3 = Cpart2* g000[xi] * r2  * Ac_star;

		D_s += Spart* g000[xi]* r2 * ( As+ 2.*Bs); 



	}

	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
}


void Hydrodynamic_Function::free_mem() {
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
	int i =0 ;
	int i_end = kn*xn;

	while ( i < i_end ) {
		if ( ihf_all[i] != NULL) {
			free(ihf_all[i] );
			ihf_all[i] = NULL;
		}
		i++;
	};
	free(ihf_all);
	free(hf_all);
	free(mf_all);
	ihf_all = NULL;
	hf_all= NULL;
	mf_all=NULL;

}


void set_ihf_sinx_cosx(ihf* pihf,double x, double sinx, double cosx) {

	if (x == 0 ) {
		puts ("Error x==0");
		exit(1);
	}
	pihf->calculated =1;

	pihf->x = x;
	pihf->x2 = x*x;
	pihf->x3 = pihf->x2*x;
	pihf->invx3 = 1./pihf->x3;
	pihf->cosx = cosx;
	pihf->cosx_x = x*cosx;
	pihf->sinx = sinx;
	pihf->sinx_3_minus_x2 = pihf->sinx *( 3.- pihf->x2);

	pihf->h1 = pihf->invx3 * ( sinx - pihf->cosx_x ) ;
	pihf->h2 = pihf->invx3 * (3. * pihf->cosx_x - pihf->sinx_3_minus_x2 ); 
}

void set_ihf(ihf* pihf,double x) {

	if (x == 0 ) {
		puts ("Error x==0");
		exit(1);
	}
	pihf->calculated =1;

	pihf->x = x;
	pihf->x2 = x*x;
	pihf->x3 = pihf->x2*x;
	pihf->invx3 = 1./pihf->x3;
	pihf->cosx = cos(x);
	pihf->cosx_x = pihf->x* pihf->cosx;
	pihf->sinx = sin(x);
	pihf->sinx_3_minus_x2 = pihf->sinx *( 3.- pihf->x2);

	pihf->h1 = pihf->invx3 * ( pihf->sinx - pihf->cosx_x ) ;
	pihf->h2 = pihf->invx3 * (3. * pihf->cosx_x - pihf->sinx_3_minus_x2 ); 
}
void copy_ihf(ihf* lhs,const ihf* rhs) {

	if ( rhs->calculated == 0) {
		puts( "Error rhs not calculated");
		exit(1);
	}
	*lhs = *rhs;
}

void Hydrodynamic_Function::end() {
	free_mem();
}


int Hydrodynamic_Function::run()
{
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);


	double part1, part2;
	int i ;
	double h1, h2;
	prepare_calculation();
	pihf pihf_i;
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
	for (int ki = 1; ki < kn; ki++) {
		part1 = 0 , part2 = 0 ;
		for (int xi = 1; xi < xn; xi++) {
			i = xi*ki; 
			pihf_i = ihf_all[i] ;
			h1 = pihf_i->h1;
			h2 = pihf_i->h2;
			part1 += mf_all[xi].part1 * ( 2.*h1+ h2);
			part2 += h1*mf_all[xi].part2 + h2*mf_all[xi].part3;
		}
		hf_all[ki] = D_s+part1+part2;
	}
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
	return 0;
}
int Hydrodynamic_Function::run(int argc, char const* argv[])
{
	printf("%s() at %s::%05d\n", __FUNCTION__, __FILE__, __LINE__);
	double part1, part2;
	int i ;
	double h1, h2;
	prepare_calculation();
	pihf pihf_i;
	for (int ki = 1; ki < kn; ki++) {
		part1 = 0 , part2 = 0 ;
		for (int xi = 1; xi < xn; xi++) {
			i = xi*ki; 
			pihf_i = ihf_all[i] ;
			h1 = pihf_i->h1;
			h2 = pihf_i->h2;
			part1 += mf_all[xi].part1 * ( 2.*h1+ h2);
			part2 += h1*mf_all[xi].part2 + h2*mf_all[xi].part3;
		}
		hf_all[ki] = D_s+part1+part2;
	}
	return 0;
}


