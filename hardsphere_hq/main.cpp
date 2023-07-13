/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2022년 04월 06일 13시 53분 17초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  KIM Hyeok (), ekh0324@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include "hydro_math.h"
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <cstdio>
#include <array>

int main(int argc, char const* argv[])
{
	Hydrodynamic_Function hf;
	double L = 64;
	int kn = 200;
	int xn = 300;
	int N  = 96;
	double r1 = L/300.;
	double k1 = 2.*M_PI/L;
	double  a = 10.061591832087*0.8/2.;
	double V= pow(L,3);
	double g000[300];
	double h000[300];

	for( int i=0; i<300; i++) {
		if( i* r1 < 2.*a) 
			g000[i] = 0;
		else
			g000[i] = 1;
		h000[i] =g000[i] -1;
	}


  hf.init(kn,xn,N,r1,k1,a,V,g000,h000);
	hf.run();

 const	double * hf_all = hf.getHq();
 const	double * ssf_all = hf.getSSF();
 for (int i =0; i<100; i++) {
	 double ki= k1*i;
	 double ssf = ssf_all[i]; 
	 double hf = hf_all[i];
	 double dq = hf/ssf;
	 printf("%f %f %f %f \n" , ki,ssf_all[i], hf_all[i], dq);
 }


	


	return 0;
}
