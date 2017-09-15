/*!
 *    \file  kh_math_fourier.h
 *   \brief  KIM Hyeok custom; math header
 *
 *  This is support to user that fourier transform for varius obital number correlation functions.
 *  But current version is only support l=2 l=0;
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
#pragma once
#include <vector>
#include <map>
#include <iterator>
typedef long int bigint;
typedef double real;
using namespace std;
real simson(vector<real> b, real dx) ;
real trapz_iso3D( vector<real> a, vector<real> b,real q) ;
real trapz_iso3D_backward( map<real,real> map, real q) ;
real trapz_iso3D_backward( vector<real> a, vector<real> b,real q) ;
real trapz_iso3D_backward_n( real *  a, real *  b,real q, int n) ;
real trapz_iso3D_forward( map<real,real> map, real q) ;
real trapz_iso3D_forward( vector<real> a, vector<real> b,real q) ;
real trapz_iso3D_forward_n( real* a, real* b,real q, int n) ;
real trapz_iso3D_forward_with_l2( vector<real> a, vector<real> b,real q) ;
real trapz_iso3D_forward_with_l2_n( real* a, real* b,real q,int n) ;
real trapzoidal( vector<real> b, real dx) ;
