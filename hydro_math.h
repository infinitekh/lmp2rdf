/*
 * =====================================================================================
 *
 *       Filename:  hydro_math.h
 *
 *    Description: jkj 
 *
 *        Version:  1.0
 *        Created:  2022년 03월 24일 23시 27분 47초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  KIM Hyeok (), ekh0324@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#pragma once


typedef struct mf_s {
	int calculated   :1 ;
	double part1 , part2, part3;
} mf;
typedef mf* pmf;

typedef struct ihf_s {
	int calculated   :1 ;
	double x,x2,x3;
	double invx3,cosx, cosx_x;
	double sinx,sinx_3_minus_x2;
	double h1,h2;
} ihf;
typedef ihf* pihf;

void copy_ihf(ihf* lhs,const ihf* rhs);
void set_ihf(ihf* pihf,double x) ;
void set_ihf_sinx_cosx(ihf* pihf,double x, double sinx, double cosx) ;



class Hydrodynamic_Function {
	public:
		Hydrodynamic_Function();
		~Hydrodynamic_Function();
		int run( int, char const* argv []);
		int run( );
		void end( );
		void init ( int,int, int, double, double, double, double, double*, double*);
		double getD_s() { return D_s;};
		double* getHq() { return hf_all;};
	protected:
		void alloc_mem();
		void free_mem();
		void prepare_calculation();
	private:
		ihf** ihf_all =NULL;
		double* hf_all =NULL;
		mf* mf_all = NULL;
		double D_s;

		int kn, xn;
		double r1, k1;

		double* g000;
		double* h000;
		double a;
		int N;
		double V;
};


