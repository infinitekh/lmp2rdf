/*!
 *    \file  kh_math_fourier.cpp
 *   \brief  
 *
 *  Implementation of kh_math_fourier.h
 *
 *  \author  KIM Hyeok (kh), ekh0324@gmail.com
 *
 *  \internal
 *       Created:  2017년 09월 12일
 *      Revision:  none
 *      Compiler:  gcc
 *  Organization:  Konkuk University
 *     Copyright:  Copyright (c) 2017, KIM Hyeok
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */


#include <numeric>
#include <cmath>
#include "kh_math_fourier.h"
real trapz_iso3D_forward_with_l2_n ( real* a, real* b, real q, int n) {
	real sum=0,dr,qri; real a1, aN, b1,bN,f1,fN,ai,bi,fi;
	real norm = 4*M_PI;

	aN= a[n-1]; bN= b[n-1];
	a1= a[0]; b1= b[0];

	dr = a[2] - a[1];
	if( q== 0) {
		return 0;
	}
	else{
		real qr1 = a1*q, qrN = aN*q;
		f1= a1*a1 *( (3./(qr1*qr1)-1)* sin(qr1)/qr1 - 3.*cos(qr1)/(qr1*qr1)     )* b1;
		fN= aN*aN *( (3./(qrN*qrN)-1)* sin(qrN)/qrN - 3.*cos(qrN)/(qrN*qrN)     )* bN;
		sum = .5*(f1+fN);
		for (int i=1 ;i<n-1; i++) {
			ai= a[i]; bi= b[i];
			qri = q* ai;
			fi= ai*ai *( (3./(qri*qri)-1)* sin(qri)/qri - 3.*cos(qri)/(qri*qri)     )* bi;
			sum += fi;
		}
		return -norm*sum*dr;
	}

}
real trapz_iso3D_forward_with_l2( vector<real> a, vector<real> b,real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation fourier transform for l=2 iso3D function   
	 *   -  4pi \int_0^infty dr h_112(r) * j_2(q r)  *r^2
	 *   j_2 (x) = (3/x^2 - 1) sin(x)/x   - 3cos(x)/x^2
	 *   j_2 (0) == 0
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr,qr; real a1, aN, b1,bN,f1,fN,ai,bi,fi;
	vector<real>::iterator a_iter = a.begin();
	vector<real>::reverse_iterator a_riter = a.rbegin();
	real norm = 4*M_PI;
	vector<real>::iterator iter = b.begin();
	vector<real>::reverse_iterator riter = b.rbegin();
	vector<real>::iterator end = b.end();

	b1= *iter; bN= *riter;
	a1= *a_iter; aN= *a_riter;

	iter++;a_iter++; end--;

	dr = *a_iter - a1;
	if( q== 0) {
		return 0;
	}
	else{
		real qr1 = a1*q, qrN = aN*q;
		f1= a1*a1 *( (3./(qr1*qr1)-1)* sin(qr1)/qr1 - 3.*cos(qr1)/(qr1*qr1)     )* b1;
		fN= aN*aN *( (3./(qrN*qrN)-1)* sin(qrN)/qrN - 3.*cos(qrN)/(qrN*qrN)     )* bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++,a_iter++) {
			ai= *a_iter; bi= *iter;
			qr = q* ai;
			fi= ai*ai *( (3./(qr*qr)-1)* sin(qr)/qr - 3.*cos(qr)/(qr*qr)     )* bi;
			sum += fi;
		}
		return -norm*sum*dr;
	}
}

real trapz_iso3D_backward( map<real,real> map, real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation isotropic function fourier transform  F[h(r)] (q)  
	 *  \int_0^infty dr h(r) * sin(qr)/(qr) * 4*pi*r^2
	 *  if q =0     ->  \int_0^infty dr h(r) * 4*pi*r^2           
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi, a2;
/* 	map<real, real >::iterator iter = map.begin();
 * 	map<real, real >::reverse_iterator riter = map.rbegin();
 * 	map<real, real >::iterator end = map.end();
 */
	auto iter = map.begin();
	auto riter = map.rbegin();
	auto end = map.end();
	real norm = 4*M_PI;

	a1= iter->first;  aN= riter->first;
	b1= iter->second; bN= riter->second;

	iter++;riter++; end--;
	a2 = iter->first;
	dr = a2 - a1;

	if( q== 0) {
		f1= norm* a1 * a1 * b1;
		fN=  norm* aN * aN * bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++) {
			ai= iter->first; bi= iter->second;
			fi=  norm* ai*ai* bi;
			sum += fi;
		}
		return sum*dr/pow(2.*M_PI,3);
	}
	else{
		f1= norm* a1 *(sin(q*a1)/q)* b1;
		fN=  norm* aN *(sin(q*aN)/q)* bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++) {
			ai= iter->first; bi= iter->second;
			fi=  norm* ai*(sin(ai*q)/q) * bi;
			sum += fi;
		}
		return sum*dr/pow(2.*M_PI,3);
	}
}
real trapz_iso3D_forward( map<real,real> map, real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation isotropic function fourier transform  F[h(r)] (q)  
	 *  \int_0^infty dr h(r) * sin(qr)/(qr) * 4*pi*r^2
	 *  if q =0     ->  \int_0^infty dr h(r) * 4*pi*r^2           
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi, a2;
/* 	map<real,real>::iterator iter = map.begin();
 * 	map<real,real>::reverse_iterator riter = map.rbegin();
 * 	map<real,real>::iterator end = map.end();
 */
	auto iter = map.begin();
	auto riter = map.rbegin();
	auto end = map.end();
	real norm = 4*M_PI;

	a1= iter->first;  aN= riter->first;
	b1= iter->second; bN= riter->second;

	iter++;riter++; end--;
	a2 = iter->first;
	dr = a2 - a1;

	if( q== 0) {
		f1= norm* a1 * a1 * b1;
		fN=  norm* aN * aN * bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++) {
			ai= iter->first; bi= iter->second;
			fi=  norm* ai*ai* bi;
			sum += fi;
		}
		return sum*dr;
	}
	else{
		f1= norm* a1 *(sin(q*a1)/q)* b1;
		fN=  norm* aN *(sin(q*aN)/q)* bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++) {
			ai= iter->first; bi= iter->second;
			fi=  norm* ai*(sin(ai*q)/q) * bi;
			sum += fi;
		}
		return sum*dr;
	}
}
real trapz_iso3D_forward_n ( real* a, real* b, real q, int n) {
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi;
	real norm = 4*M_PI;

	a1= a[0]; aN= a[n-1];
	b1= b[0]; bN= b[n-1];

	dr = a[2]-a[1];

	if( q== 0) {
		f1= norm* a1 * a1 * b1;
		fN=  norm* aN * aN * bN;
		sum = .5*(f1+fN);
		for (int i=1 ;i<n-1; i++) {
			ai= a[i]; bi= b[i];
			fi=  norm* ai*ai* bi;
			sum += fi;
		}
		return sum*dr;
	}
	else{
		f1= norm* a1 *(sin(q*a1)/q)* b1;
		fN=  norm* aN *(sin(q*aN)/q)* bN;
		sum = .5*(f1+fN);
		for (int i=1 ;i<n-1; i++) {
			ai= a[i]; bi= b[i];
			fi=  norm* ai*(sin(ai*q)/q) * bi;
			sum += fi;
		}
		return sum*dr;
	}

}
real trapz_iso3D_forward( vector<real> a, vector<real> b,real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation isotropic function fourier transform  F[h(r)] (q)  
	 *  \int_0^infty dr h(r) * sin(qr)/(qr) * 4*pi*r^2
	 *  if q =0     ->  \int_0^infty dr h(r) * 4*pi*r^2           
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi;
	vector<real>::iterator a_iter = a.begin();
	vector<real>::reverse_iterator a_riter = a.rbegin();
	real norm = 4*M_PI;
	vector<real>::iterator iter = b.begin();
	vector<real>::reverse_iterator riter = b.rbegin();
	vector<real>::iterator end = b.end();

	b1= *iter; bN= *riter;
	a1= *a_iter; aN= *a_riter;

	iter++;a_iter++; end--;

	dr = *a_iter - a1;

	if( q== 0) {
		f1= norm* a1 * a1 * b1;
		fN=  norm* aN * aN * bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++,a_iter++) {
			ai= *a_iter; bi= *iter;
			fi=  norm* ai*ai* bi;
			sum += fi;
		}
		return sum*dr;
	}
	else{
		f1= norm* a1 *(sin(q*a1)/q)* b1;
		fN=  norm* aN *(sin(q*aN)/q)* bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++,a_iter++) {
			ai= *a_iter; bi= *iter;
			fi=  norm* ai*(sin(ai*q)/q) * bi;
			sum += fi;
		}
		return sum*dr;
	}
}
real trapz_iso3D_backward_n( real *  a, real *  b,real q, int n) {
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi;
	real norm = 4*M_PI;

	b1= b[0]; bN= b[n-1];
	a1= a[0]; aN= a[n-1];


	dr = a[2]-a[1];

	if( q== 0) {
		f1= norm* a1 * a1 * b1;
		fN=  norm* aN * aN * bN;
		sum = .5*(f1+fN);
		for (int i=1 ;i<n-1; i++) {
			ai= a[i]; bi= b[i];
			fi=  norm* ai*ai* bi;
			sum += fi;
		}
		return sum*dr/pow(2.*M_PI,3);
	}
	else{
		f1= norm* a1 *(sin(q*a1)/q)* b1;
		fN=  norm* aN *(sin(q*aN)/q)* bN;
		sum = .5*(f1+fN);
		for (int i=1 ;i<n-1; i++) {
			ai= a[i]; bi= b[i];
			fi=  norm* ai*(sin(ai*q)/q) * bi;
			sum += fi;
		}
		return sum*dr/pow(2.*M_PI,3);
	}
}
real trapz_iso3D_backward( vector<real> a, vector<real> b,real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation isotropic function fourier transform  F[h(r)] (q)  
	 *  \int_0^infty dr h(r) * sin(qr)/(qr) * 4*pi*r^2
	 *  if q =0     ->  \int_0^infty dr h(r) * 4*pi*r^2           
	 *  backward   q-> r    r -> q   and   -> frac 2 pi  
	 *  fourier transform convention : r -> k (1)   k -> r ( 1/2pi)  
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi;
	vector<real>::iterator a_iter = a.begin();
	vector<real>::reverse_iterator a_riter = a.rbegin();
	real norm = 4*M_PI;
	vector<real>::iterator iter = b.begin();
	vector<real>::reverse_iterator riter = b.rbegin();
	vector<real>::iterator end = b.end();

	b1= *iter; bN= *riter;
	a1= *a_iter; aN= *a_riter;

	iter++;a_iter++; end--;

	dr = *a_iter - a1;

	if( q== 0) {
		f1= norm* a1 * a1 * b1;
		fN=  norm* aN * aN * bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++,a_iter++) {
			ai= *a_iter; bi= *iter;
			fi=  norm* ai*ai* bi;
			sum += fi;
		}
		return sum*dr/pow(2.*M_PI,3);
	}
	else{
		f1= norm* a1 *(sin(q*a1)/q)* b1;
		fN=  norm* aN *(sin(q*aN)/q)* bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++,a_iter++) {
			ai= *a_iter; bi= *iter;
			fi=  norm* ai*(sin(ai*q)/q) * bi;
			sum += fi;
		}
		return sum*dr/pow(2.*M_PI,3);
	}
}
real trapz_iso3D( vector<real> a, vector<real> b,real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation isotropic function fourier transform  F[h(r)] (q)  
	 *  \int_0^infty dr h(r) * sin(qr)/(qr) * 4*pi*r^2
	 *  if q =0     ->  \int_0^infty dr h(r) * 4*pi*r^2           
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi;

	vector<real>::iterator a_iter = a.begin();
	vector<real>::reverse_iterator a_riter = a.rbegin();

	real norm = 4*M_PI;
	vector<real>::iterator iter = b.begin();
	vector<real>::reverse_iterator riter = b.rbegin();
	vector<real>::iterator end = b.end();

	b1= *iter; bN= *riter;
	a1= *a_iter; aN= *a_riter;

	iter++;a_iter++; end--;

	dr = *a_iter - a1;

	if( q== 0) {
		f1= norm* a1 * a1 * b1;
		fN=  norm* aN * aN * bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++,a_iter++) {
			ai= *a_iter; bi= *iter;
			fi=  norm* ai*ai* bi;
			sum += fi;
		}
		return sum*dr;
	}
	else{
		f1= norm* a1 *(sin(q*a1)/q)* b1;
		fN=  norm* aN *(sin(q*aN)/q)* bN;
		sum = .5*(f1+fN);
		for ( ;iter!=end; iter++,a_iter++) {
			ai= *a_iter; bi= *iter;
			fi=  norm* ai*(sin(ai*q)/q) * bi;
			sum += fi;
		}
		return sum*dr;
	}
}
real trapzoidal( vector<real> b, real dx) {
	real sum=0; real f1, fN;
	f1 = *(b.begin()); fN= *(b.end());
//	for (  real y : b) sum += y;
	sum = accumulate( b.begin(), b.end(), 0);
	sum -=.5* (f1+fN);
		
	return sum*dx;
}
// simson 3/8
real simson(vector<real> b, real dx) {
	real sum=0; real f1, f2,f3,fNm2,fNm1, fN;
	vector<real>::iterator begin = b.begin();
	vector<real>::reverse_iterator rbegin = b.rbegin();
	f1 = *begin; begin++;
	f2 = *begin; begin++;
	f3 = *begin;
	fN = *rbegin; rbegin++;
	fNm1=*rbegin; rbegin++;
	fNm2=*rbegin;

//	for (  real y : b) sum += y;
	sum = accumulate( b.begin(), b.end(), 0);
	sum -= 5./8.* (f1+fN);
	sum += 1./6.* (f2+fNm1);
	sum -= 1./24.*(f3+fNm2);
		
	return sum*dx;
}
