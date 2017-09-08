#include "makeRDF.h"
#include <unistd.h>
#include <memory>
#include <vector>
#include <numeric>
#define Dz     (2)
#define halfDz     (Dz/2.0)
using namespace std;

real trapz_iso3D_forward_with_l2( std::vector<real> a, std::vector<real> b,real q) {
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

real trapz_iso3D_backward( std::map<real,real> map, real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation isotropic function fourier transform  F[h(r)] (q)  
	 *  \int_0^infty dr h(r) * sin(qr)/(qr) * 4*pi*r^2
	 *  if q =0     ->  \int_0^infty dr h(r) * 4*pi*r^2           
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi, a2;
	std::map<real,real>::iterator iter = map.begin();
	std::map<real,real>::reverse_iterator riter = map.rbegin();
	std::map<real,real>::iterator end = map.end();
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
real trapz_iso3D_forward( std::map<real,real> map, real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation isotropic function fourier transform  F[h(r)] (q)  
	 *  \int_0^infty dr h(r) * sin(qr)/(qr) * 4*pi*r^2
	 *  if q =0     ->  \int_0^infty dr h(r) * 4*pi*r^2           
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi, a2;
	std::map<real,real>::iterator iter = map.begin();
	std::map<real,real>::reverse_iterator riter = map.rbegin();
	std::map<real,real>::iterator end = map.end();
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
real trapz_iso3D_forward( std::vector<real> a, std::vector<real> b,real q) {
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
real trapz_iso3D_backward( std::vector<real> a, std::vector<real> b,real q) {
	/*-----------------------------------------------------------------------------
	 *  calulation isotropic function fourier transform  F[h(r)] (q)  
	 *  \int_0^infty dr h(r) * sin(qr)/(qr) * 4*pi*r^2
	 *  if q =0     ->  \int_0^infty dr h(r) * 4*pi*r^2           
	 *  backward   q-> r    r -> q   and   -> frac 2 pi  
	 *  fourier transform convention : r -> k (1)   k -> r ( 1/2pi)  
	 *-----------------------------------------------------------------------------*/
	real sum=0,dr; real a1, aN, b1,bN,f1,fN,ai,bi,fi;
	vector<real>::iterator a_iter = a.begin();
	std::vector<real>::reverse_iterator a_riter = a.rbegin();
	real norm = 4*M_PI;
	vector<real>::iterator iter = b.begin();
	std::vector<real>::reverse_iterator riter = b.rbegin();
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
	sum = std::accumulate( b.begin(), b.end(), 0);
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
	sum = std::accumulate( b.begin(), b.end(), 0);
	sum -= 5./8.* (f1+fN);
	sum += 1./6.* (f2+fNm1);
	sum -= 1./24.*(f3+fNm2);
		
	return sum*dx;
}
void makeRDF::calcSSF () {
	if(flag_SSF_from_g)
		calcSSF_from_g ();
	else {
	  int bin,n,i,iii,jjj,kkk;	
		real qx,qy,qz,q1;	
		real x,y,z;
		int id;
		real* qlist = new real [2*maxbinq+1];
		real* hist = new real [maxbinq+2];
		real* S_qr       = new real [maxbinq+2];
		real* S_qi       = new real [maxbinq+2];
		FILE* fp_out;
		fp_out = fopen("SSF_direct.info", "w");

		dq = (4.*M_PI / maxbinq);
		for ( int i = 0; i < 2*maxbinq+1; i += 1 )
			qlist[i] = -2.*M_PI+ dq * i;

		puts("");
		int files=0;
		for (list<Snapshot*>::iterator iter =snaplist.begin() ; iter != snaplist.end();  iter++) {
			maxAtom = (*iter)->n_atoms;
			atoms   = (*iter)->atoms;
			box     = &(*iter)->box;
			box_x = box->xhigh-box->xlow; 
			box_y = box->yhigh-box->ylow; 
			box_z = box->zhigh-box->zlow;
			files++;
			for ( id =0; id<maxAtom; id++){
			  ppi = &atoms[id];
				printf("\r%5dth-Snapshots(%ld) %4dth-atoms(%d)", files,snaplist.size(), id+1,maxAtom);
				fflush(stdout);
				x = ppi->x; y = ppi->y; z = ppi->z;
				for ( iii = 0; iii < 2*maxbinq+1; iii += 1 ) {
					qx = qlist[iii];		
					for ( jjj = 0; jjj < 2*maxbinq+1; jjj += 1 ) {
						qy = qlist[jjj];
						for ( kkk = 0; kkk < 2*maxbinq+1; kkk += 1 ) {
							qz = qlist[kkk];
							q1 = sqrt(qx*qx+qy*qy+qz*qz);
							bin = int(ceil(q1/dq));
							if (bin <= maxbinq) {
								hist[bin] += 1;
								S_qr[bin] += cos( qx* x+qy*y+qz* z);
								S_qi[bin] += sin( qx* x+qy*y+qz* z);
							}

						}
					}
				}
			}
		}
		
		real norm_q = 4.0 * M_PI * dq * maxAtom * snaplist.size();

		real qqq, value,value2;
		for ( i = 0; i < maxbinq+1; i += 1 ){
			qqq = ( i-0.5) * dq;
			if (hist[i] ==0) {
				value =0;
				value2=0;
			}
			else {
				value = sqrt(S_qr[i]*S_qr[i] + S_qi[i]*S_qi[i])/ norm_q/ ((qqq*qqq)+(dq*dq/12.0) ) ;
				value2= sqrt(S_qr[i]*S_qr[i] + S_qi[i]*S_qi[i])/ maxAtom / hist[i]; 
			}
			
			fprintf(fp_out, "%le ", qqq);
			fprintf(fp_out, "%le ", value);
			fprintf(fp_out, "%le ", value2);
			fprintf(fp_out, "\n");

		}

		fclose (fp_out);
		delete [] hist;
		delete [] S_qr;
		delete [] S_qi;
	}
}

void makeRDF::calcSSF_from_g() {
	int size = h000.size();
	real r, h0_r, c_r, h_mod, h112_r;
	real q,h_k, c_k, h112_k;
	real S_q;
	std::vector<real> x,y, y2;
	FILE* fp_S = fopen("SSF.info","w");
	fprintf(fp_S,"##q S_k h_k c_k h112_k\n");
	for(std::map<real,real>::iterator  iter=h000.begin(); iter != h000.end(); iter++) {
		 r = iter->first;
		 h0_r = iter->second;
		 x.push_back ( r);
		 y.push_back ( h0_r);
		 y2.push_back ( h112[r] );
	}

/* 	for( int i=0; i<size; i++) {
 * 		 r = h000[i].first;
 * 		 h0_r = h000[i].second;
 * 		 h112_r = h112[i].second;
 * 		 x.push_back ( r);
 * 		 y.push_back ( h0_r);
 * 		 y2.push_back ( h112_r);
 * 	}
 */
	real Vol = box_x*box_y*box_z;
	real phi = maxAtom / Vol;
	for(int i=0; i<maxbinq ; i++) {
		q = i * dq;
		h_k = trapz_iso3D_forward( x,y,q);
		h112_k = trapz_iso3D_forward_with_l2( x,y,q);
		S_q = 1.+ phi*h_k;
		c_k = h_k / S_q;
    S.insert( std::pair<real,real>(q,S_q ) );
    h_q.insert( std::pair<real,real>(q, h_k ) );
    c_q.insert( std::pair<real,real>(q, c_k ) );

		fprintf(fp_S,"%lf %lf %lf %lf %lf\n", q, S_q, h_k, c_k, h112_k);
	}
	fclose(fp_S);
	/*-----------------------------------------------------------------------------
	 * we get c_r from c_k using Backward Fourier transform
	 *-----------------------------------------------------------------------------*/


	FILE* fp_c_r = fopen("DCF.info","w");
	fprintf(fp_c_r,"##r c_r h_mod\n");

	for(int i=0; i<maxbin ; i++) {
		r = i * var_r;
		c_r = trapz_iso3D_backward( c_q,r);
		h_mod = trapz_iso3D_backward( h_q,r);

		fprintf(fp_c_r,"%lf %lf %lf \n", r, c_r, h_mod );
	}

	fclose(fp_c_r);
	
}




void makeRDF::calcRDF() {
	if (flag_anisotropy )
		calcRDF_anisotropy();
	else 
		calcRDF_isotropy();
}
void makeRDF::calcRDF_isotropy () {
	bigint i, ii, jj;
	int nsnap =0;

	int bin, rhobin, zibin, zjbin, maxallbin, allbin;
	real xij,yij,zij,zi,zj, r1, rho1, z1, 
			 abszi,abszj,var_z=Dz,var_r = r_cut/maxbin;
	real hbox_x,hbox_y,hbox_z;
	
	maxallbin = maxbinz*maxbinz*maxbin;
	bigint* hist000 = new bigint[maxbin+1];
	bigint* histszz1 = new bigint[maxallbin+1];
	bigint* histszz2 = new bigint[maxallbin+1];

	// for slab
	calcP1z (1);
	
	for(i=0; i<=maxallbin; i++) {
		histszz1[i] = 0;
		histszz2[i] = 0;
	}
	for(i=0; i<=maxbin; i++) {
		hist000[i] = 0;
	}
	for (list<Snapshot*>::iterator iter =snaplist.begin() ; iter != snaplist.end();  iter++) {
		nsnap ++;
		maxAtom = (*iter)->n_atoms;
		atoms   = (*iter)->atoms;
		box     = &(*iter)->box;
		box_x = box->xhigh-box->xlow; 
		box_y = box->yhigh-box->ylow; 
		box_z = box->zhigh-box->zlow;
		hbox_x= box_x/2.;
		hbox_y= box_y/2.;
		hbox_z= box_z/2.;
		for ( ii =0; ii<maxAtom; ii++){
			ppi = &atoms[ii];
			zi = ppi->z;
			for ( jj =0; jj<maxAtom; jj++) { 
				ppj = &atoms[jj];
				zj = ppj->z;
				if ( ppi->id > ppj->id){
					xij = ppi->x - ppj->x;
					yij = ppi->y - ppj->y;
					zij = zi - zj ;
					if ( periodicity[0]){
						while ( xij <-hbox_x) xij+=box_x;
						while ( xij > hbox_x) xij-=box_x;
					}
					if ( periodicity[1]){
						while ( yij <-hbox_y) yij+=box_y;
						while ( yij > hbox_y) yij-=box_y;
					}
					if ( periodicity[2]){
						while ( zij <-hbox_z) zij+=box_z;
						while ( zij > hbox_z) zij-=box_z;
					}
					r1 = sqrt(xij*xij+yij*yij+zij*zij);
					rho1 = sqrt(xij*xij+yij*yij);
					z1 = abs(zij); abszi=abs(zi); abszj=abs(zj);
					bin = static_cast<int>( r1/var_r);
					rhobin = static_cast<int>( rho1/var_r);
					zibin = static_cast<int>( floor(abszi/var_z));
					zjbin = static_cast<int>( floor(abszj/var_z));
					if (bin <=maxbin)
						hist000[bin] +=1;
					if ( rhobin <=maxbin && zibin < maxbinz && zjbin < maxbinz) {
						allbin = maxbin*maxbinz*zibin + zjbin*maxbin+ rhobin;
						if( zi* zj >0 ) histszz1[allbin] +=1 ;
						else histszz2[allbin] +=1 ;
					}
				}
			}
		}
	}
	real Vol = box_x*box_y*box_z;
	real Area= box_x*box_y;
	real phi = maxAtom / Vol;
	real val000, valszz1, valszz2;
	real norm000    = 2.0* M_PI *var_r * phi * nsnap* maxAtom;
	real normcyl000 = M_PI * Sqr(var_z) * Sqr(var_r)* Area  *  nsnap;
	real rrr,rr, zz1, zz2;
	FILE* fp_rdf = fopen("rdf.info","w");
	FILE* fp_rdfcyl = fopen("rdfcyl.info","w");
	if ( (periodicity[2])) 
		fputs("##3d periodicity\n", fp_rdf);
	else
		fputs("##2d periodicity and z axis bounded\n", fp_rdf);
	for (i=1; i<=maxbin; i++){
		rrr = (i + 0.5) * var_r;
		rr  = (i - 0.5) * var_r;
		val000 = hist000[i]/norm000/(rr*rr+var_r*var_r*12.);
		g000.insert( std::pair<real,real>( rrr,val000 ) );
		h000.insert( std::pair<real,real>( rrr,val000-1. ) );
		h110.insert( std::pair<real,real>( rrr,0.0 ) );
		h112.insert( std::pair<real,real>( rrr,0.0 ) );
		h220.insert( std::pair<real,real>( rrr,0.0 ) );
		fprintf(fp_rdf,"%lf %lf\n", rrr,val000);
	}
	i=0; rhobin =0; zibin=0; zjbin=0;
	while ( i < maxallbin) {
		rrr = (rhobin + 0.5) * var_r;
		zz1 = (zibin + 0.5) * var_z;
		zz2 = (zjbin + 0.5) * var_z;
		
		valszz1 = histszz1[i]/(normcyl000*(rhobin+.5) * vP1z[zibin]* vP1z[zjbin]);
		valszz2 = histszz2[i]/(normcyl000*(rhobin+.5) * vP1z[zibin]* vP1z[zjbin]);
		fprintf(fp_rdfcyl,"%lf %lf %lf %lf %lf\n", rrr, zz1, zz2,valszz1, valszz2);

		i++; rhobin ++;
		if (rhobin == maxbin) { rhobin=0; zjbin++ ;}
		if ( zjbin == maxbinz) { zjbin=0; zibin++ ;}
	}
	fclose(fp_rdf);
	fclose(fp_rdfcyl);
}
void makeRDF::calcRDF_anisotropy () {
	int i, ii, jj;
	int nsnap =0;

	int bin, rhobin;
	real xij,yij,zij, r1, rho1, z1, var_z=Dz,var_r = r_cut/maxbin;
	real hbox_x,hbox_y,hbox_z;


	bigint* hist000 = new bigint[maxbin+1];
	double* hist110 = new double[maxbin+1];
	double* hist112 = new double[maxbin+1];
	double* hist220 = new double[maxbin+1];
	bigint* histcyl000 = new bigint[(maxbin+1)];
	double* histcyl110 = new double[(maxbin+1)];
	double* histcyl112 = new double[(maxbin+1)];
	double* histcyl220 = new double[(maxbin+1)];
	for(i=0; i<=maxbin; i++) {
		hist000[i] = 0;
		hist110[i] = 0.;
		hist112[i] = 0.;
		hist220[i] = 0.;
		histcyl000[i] = 0;
		histcyl110[i] = 0;
		histcyl112[i] = 0;
		histcyl220[i] = 0;
	}

	real si_dot_sj, si_dot_rij, sj_dot_rij, si[3],sj[3];

	for (list<Snapshot*>::iterator iter =snaplist.begin() ; iter != snaplist.end();  iter++) {
		nsnap ++;
		maxAtom = (*iter)->n_atoms;
		atoms   = (*iter)->atoms;
		box     = &(*iter)->box;
		box_x = box->xhigh-box->xlow; 
		box_y = box->yhigh-box->ylow; 
		box_z = box->zhigh-box->zlow;
		hbox_x= box_x/2.;
		hbox_y= box_y/2.;
		hbox_z= box_z/2.;
		for ( ii =0; ii<maxAtom; ii++){
			ppi = &atoms[ii];
			si[0]=ppi->mux;si[1]=ppi->muy;si[2]=ppi->muz;

			for ( jj =0; jj<maxAtom; jj++) { 
				ppj = &atoms[jj];
				sj[0]=ppj->mux;sj[1]=ppj->muy;sj[2]=ppj->muz;
				if ( ppi->id > ppj->id){

					xij = ppi->x - ppj->x;
					yij = ppi->y - ppj->y;
					zij = ppi->z - ppj->z;

					if ( periodicity[0]){
						while ( xij <-hbox_x) xij+=box_x;
						while ( xij > hbox_x) xij-=box_x;
					}
					if ( periodicity[1]){
						while ( yij <-hbox_y) yij+=box_y;
						while ( yij > hbox_y) yij-=box_y;
					}
					if ( periodicity[2]){
						while ( zij <-hbox_z) zij+=box_z;
						while ( zij > hbox_z) zij-=box_z;
					}

					r1 = sqrt(xij*xij+yij*yij+zij*zij);
					rho1 = sqrt(xij*xij+yij*yij);
					z1 = abs(zij);
					bin = static_cast<int>( r1/var_r);
					rhobin = static_cast<int>( rho1/var_r);
					real rij[3] = {xij,yij,zij};

					if (bin <=maxbin){
						si_dot_sj = VDOT3(si,sj);
						si_dot_rij= VDOT3(si,rij)/r1;
						sj_dot_rij= VDOT3(sj,rij)/r1;
						hist000[bin] +=1;
						hist110[bin] += si_dot_sj;
						hist112[bin] += 3.0*si_dot_rij*sj_dot_rij-si_dot_sj;
						hist220[bin] += 3.0*(si_dot_sj*si_dot_sj) -1.;
					}
					if ( rhobin <=maxbin && z1 < halfDz) {
						si_dot_sj = VDOT3(si,sj);
						si_dot_rij= VDOT3(si,rij)/r1;
						sj_dot_rij= VDOT3(sj,rij)/r1;
						histcyl000[rhobin] +=1 ;
						histcyl110[rhobin] +=si_dot_sj;
//						printf("%f %f %f\n", si_dot_sj, si_dot_rij, histcyl110[rhobin] );
						histcyl112[rhobin] +=3.0*si_dot_rij*sj_dot_rij-si_dot_sj;
						histcyl220[rhobin] +=3.0*(si_dot_sj*si_dot_sj) -1.;
					}
				}
			}
		}
	}
	real Vol = box_x*box_y*box_z;
	real phi = maxAtom / Vol;
	real val000,val110,val112,val220;
	real norm000 = 1./(2.0* M_PI *var_r * phi * nsnap* maxAtom);
	real norm110 = norm000*3.0;
	real norm112 = norm000*(2./3.);
	real norm220 = norm000*2.5;

	real normcyl000 = 1./(M_PI * var_r * var_z * phi * nsnap* maxAtom);
	real normcyl110 = normcyl000*3.0;
	real normcyl112 = normcyl000*(2./3.);
	real normcyl220 = normcyl000*2.5;
	real rrr,rr, invr2, invr, zzz;
	
	char filename1[100] ="rdf00.info" ;
	char filename2[100] ="rdfcyl00.info";
	int nfile = 0;
	while( 0 == (access(filename1,F_OK) + access(filename2,F_OK) )) {
		nfile++;
		sprintf(filename1, "rdf%02d.info",nfile);
		sprintf(filename2, "rdfcyl%02d.info",nfile);
	}



	FILE* fp_rdf = fopen(filename1,"w");
	FILE* fp_rdfcyl = fopen(filename2,"w");
/* 	if ( (periodicity[2])) 
 * 		fputs("##3d periodicity\n", fp_rdf);
 * 	else
 * 		fputs("##2d periodicity and z axis bounded\n", fp_rdf);
 */
	for (i=1; i<=maxbin; i++){
		rrr = (i + 0.5) * var_r;
		rr  = (i - 0.5) * var_r;
		invr2= 1./(rr*rr+var_r*var_r/12.);
		val000 = double(hist000[i])*norm000*invr2;
		val110 = hist110[i]*norm110*invr2;
		val112 = hist112[i]*norm112*invr2;
		val220 = hist220[i]*norm220*invr2;
		g000.insert( std::pair<real,real>( rrr,val000 ) );
		h000.insert( std::pair<real,real>( rrr,val000-1. ) );
		h110.insert( std::pair<real,real>( rrr,val110 ) );
		h112.insert( std::pair<real,real>( rrr,val112 ) );
		h220.insert( std::pair<real,real>( rrr,val220 ) );
		fprintf(fp_rdf,"%lf %lf %lf %lf %lf\n", rrr,val000,val110,val112,val220);
	}
	for(i=1; i<=(maxbin); i++) {
		rhobin = i;
		rrr = (rhobin + 0.5) * var_r;
		rr  = rrr- var_r;
		invr = 1/rr;

		val000 = double(histcyl000[i])*normcyl000*invr;
		val110 = histcyl110[i]*normcyl110*invr;
		val112 = histcyl112[i]*normcyl112*invr;
		val220 = histcyl220[i]*normcyl220*invr;
		fprintf(fp_rdfcyl,"%.5e %.5e %.5e %.5e %.5e %.5e \n", rrr, var_z,val000,val110,val112,val220);
	}
	fclose(fp_rdf);
	fclose(fp_rdfcyl);
}
void makeRDF::calcRDF_cotype (int itype) {
}
void makeRDF::calcRDF_inter_type (int itype, int jtype, T_RDF rdftype) {
	int i, ii, jj;
	int nsnap =0;
	int ltype=   itype<jtype?itype:jtype;
	int gtype=   itype>jtype?itype:jtype;
	int n_itype=0, n_jtype=0;	

	int bin, rhobin, zibin, zjbin, maxallbin, allbin;
	real xij,yij,zij,zi,zj, r1, rho1, z1 , invr2,
			 abszi,abszj,var_z=Dz,var_r = r_cut/maxbin;
	real hbox_x,hbox_y,hbox_z;

	typedef bigint *   p_bigint;
	typedef double *  p_double;
	p_bigint hist000, histszz1, histszz2;
	p_double hist110, hist112, hist220;
	real si_dot_sj, si_dot_rij, sj_dot_rij, si[3],sj[3];

	maxallbin = maxbinz*maxbinz*maxbin;
	hist000 = new bigint[maxbin+1];
	histszz1 = new bigint[maxallbin+1];
	histszz2 = new bigint[maxallbin+1];
	if ( rdftype == ANISO ) {
		hist110 = new double[maxbin+1];
		hist112 = new double[maxbin+1];
		hist220 = new double[maxbin+1];
	}
	int n_types_2body = (*snaplist.begin())->calc_n_type();
	n_types_2body = (n_types_2body* (n_types_2body+1)) /2;

	// for slab
	calcP1z ( itype) ;
	vP1zi.resize( vP1z.size() );
	vP1zj.resize( vP1z.size() );
	printf("%ld  \n", vP1z.size() );
	printf("%ld  \n", vP1zi.size() );
	printf("%ld  \n", vP1zj.size() );
	copy(vP1z.begin(), vP1z.end(), vP1zi.begin());
	if ( itype!= jtype) calcP1z(jtype);
	copy(vP1z.begin(), vP1z.end(), vP1zj.begin());

	
	for(i=0; i<=maxallbin; i++) {
		histszz1[i] = 0;
		histszz2[i] = 0;
	}
	for(i=0; i<=maxbin; i++) {
		hist000[i] = 0;
	}
	if (rdftype == ANISO) {
		for(i=0; i<=maxbin; i++) {
			hist110[i] = 0.; 
			hist112[i] = 0.; 
			hist220[i] = 0.; 
		}
	}
	list<Snapshot*>::iterator iter = snaplist.begin();
	atoms = (*iter)->atoms;
	maxAtom = (*iter)->n_atoms;
	for ( ii =0; ii<maxAtom; ii++){
		ppi = &atoms[ii];

		if ( ppi->type == itype) n_itype++;
		if ( ppi->type == jtype) n_jtype++;
	}
	for (; iter != snaplist.end();  iter++) {
		nsnap ++;
		maxAtom = (*iter)->n_atoms; atoms   = (*iter)->atoms;
		box     = &(*iter)->box; box_x = box->xhigh-box->xlow; 
		box_y = box->yhigh-box->ylow; 
		box_z = box->zhigh-box->zlow;
		hbox_x= box_x/2.;
		hbox_y= box_y/2.;
		hbox_z= box_z/2.;
		for ( ii =0; ii<maxAtom; ii++){
			ppi = &atoms[ii];
			zi = ppi->z-hbox_z;
			si[0]=ppi->mux;si[1]=ppi->muy;si[2]=ppi->muz;
			if ( ppi->type != itype) continue;
			for ( jj =0; jj<maxAtom; jj++) { 
				ppj = &atoms[jj];
				zj = ppj->z - hbox_z;
				sj[0]=ppj->mux;sj[1]=ppj->muy;sj[2]=ppj->muz;
				if ( ppj->type != jtype) continue;
				{
					xij = ppi->x - ppj->x;
					yij = ppi->y - ppj->y;
					zij = ppi->z - ppj->z;
					if ( periodicity[0]){
						while ( xij <-hbox_x) xij+=box_x;
						while ( xij > hbox_x) xij-=box_x;
					}
					if ( periodicity[1]){
						while ( yij <-hbox_y) yij+=box_y;
						while ( yij > hbox_y) yij-=box_y;
					}
					if ( periodicity[2]){
						while ( zij <-hbox_z) zij+=box_z;
						while ( zij > hbox_z) zij-=box_z;
					}
					r1 = sqrt(xij*xij+yij*yij+zij*zij);
					rho1 = sqrt(xij*xij+yij*yij);
					z1 = abs(zij); abszi=abs(zi); abszj=abs(zj);
					bin    = static_cast<int>( r1/var_r);
					rhobin = static_cast<int>( rho1/var_r);
					zibin  = static_cast<int>( floor(abszi/var_z));
					zjbin  = static_cast<int>( floor(abszj/var_z));

					if (bin <=maxbin){ 
						hist000[bin] +=1;
						if (rdftype == ANISO) {
							real rij[3] = {xij,yij,zij};

							si_dot_sj = VDOT3(si,sj);
							si_dot_rij= VDOT3(si,rij);
							sj_dot_rij= VDOT3(sj,rij);
							hist110[bin] += si_dot_sj;
							hist112[bin] += 3.0*si_dot_rij*sj_dot_rij/(r1*r1)-si_dot_sj;
							hist220[bin] += 3.0*(si_dot_sj*si_dot_sj) -1.;
						}
					}
					if ( rhobin <=maxbin && zibin < maxbinz && zjbin < maxbinz) {
						allbin = maxbin*maxbinz*zibin + zjbin*maxbin+ rhobin;
						if( (signbit(zi) xor signbit(zj)) ==0 ) histszz1[allbin] +=1 ;
						else histszz2[allbin] +=1 ;
					}
				}
			}
		}
	}
	real Vol = box_x*box_y*box_z;
	real Area= box_x*box_y;
	real phi_i = n_itype / Vol;
	real phi_j = n_jtype / Vol;
	real val000, val110, val112, val220,valszz1, valszz2;
	real norm000    = 1./ (4.0* M_PI *var_r * phi_i * phi_j * nsnap* Vol);
	real norm110 = norm000*3.0;
	real norm112 = norm000*(2./3.);
	real norm220 = norm000*2.5;

	real normcyl000 = 2.* 2.0* M_PI * Sqr(var_z) * Sqr(var_r)  *  nsnap* Area;
	//  11 -1-1    1 -1   -1 1
	real rrr,rr, zz1, zz2;
	char rdffilename[100] ;
	sprintf(rdffilename, "rdf00_ij%d%d.info",itype, jtype);
	if (rdftype == ANISO) {
		sprintf(rdffilename, "rdf11_ij%d%d.info",itype, jtype);
	}	
	FILE* fp_rdf = fopen(rdffilename,"w");
	sprintf(rdffilename, "rdfcyl_ij%d%d.info",itype, jtype);
	FILE* fp_rdfcyl = fopen(rdffilename,"w");

	if ( (periodicity[2])) 
		fputs("##3d periodicity\n", fp_rdf);
	else
		fputs("##2d periodicity and z axis bounded\n", fp_rdf);
	if (rdftype == ANISO) {
		for (i=1; i<=maxbin; i++){
			rrr = (i + 0.5) * var_r;
			rr  = (i - 0.5) * var_r;
			invr2= 1./(rr*rr+var_r*var_r/12.);
			val000 = double(hist000[i])*norm000*invr2;
			val110 = hist110[i]*norm110*invr2;
			val112 = hist112[i]*norm112*invr2;
			val220 = hist220[i]*norm220*invr2;
			g000.insert( std::pair<real,real>( rrr,val000 ) );
			h000.insert( std::pair<real,real>( rrr,val000-1. ) );
			h110.insert( std::pair<real,real>( rrr,val110 ) );
			h112.insert( std::pair<real,real>( rrr,val112 ) );
			h220.insert( std::pair<real,real>( rrr,val220 ) );
			fprintf(fp_rdf,"%lf %lf %lf %lf %lf\n", rrr,val000,val110,val112,val220);
		}
	} else {
		for (i=1; i<=maxbin; i++){
			rrr = (i + 0.5) * var_r;
			rr  = (i - 0.5) * var_r;
			invr2= 1./(rr*rr+var_r*var_r/12.);
			val000 = hist000[i]*norm000*(rr*rr+var_r*var_r*12.);
			g000.insert( std::pair<real,real>( rrr,val000 ) );
			h000.insert( std::pair<real,real>( rrr,val000-1. ) );
			fprintf(fp_rdf,"%lf %lf\n", rrr,val000);
		}
	}
	i=0; rhobin =0; zibin=0; zjbin=0;
	while ( i < maxallbin) {
		rrr = (rhobin + 0.5) * var_r;
		zz1 = (zibin + 0.5) * var_z;
		zz2 = (zjbin + 0.5) * var_z;
		
		valszz1 = histszz1[i]/(normcyl000*(rhobin+.5) * vP1zi[zibin]* vP1zj[zjbin]);
		valszz2 = histszz2[i]/(normcyl000*(rhobin+.5) * vP1zi[zibin]* vP1zj[zjbin]);
		fprintf(fp_rdfcyl,"%lf %lf %lf %lf %lf\n", rrr, zz1, zz2,valszz1, valszz2);

		i++; rhobin ++;
		if (rhobin == maxbin) { rhobin=0; zjbin++ ;}
		if ( zjbin == maxbinz) { zjbin=0; zibin++ ;}
	}

	delete [] hist000, hist110, hist112,hist220, histszz1,histszz2;
	fclose(fp_rdf);
	fclose(fp_rdfcyl);
}
// P1(z) = A*P(x,y,z)
void makeRDF::calcP1z( int type) {
	int i, ii, jj;
	int nsnap =0; int bin;
	real z1, zi, var_z=Dz;

	if ( periodicity[2])
		fprintf(stderr, "WARNING : z-direction is not periodic.\n");

	bigint* numSum = new bigint[maxbinz+1];
	double* szSum = new double[maxbinz+1];
	for(i=0; i<=maxbinz; i++) {
		numSum[i] = szSum[i] =0.;
	}

	real sz,zz, val, val1, hbox_z;

	list<Snapshot*>::iterator iter = snaplist.begin();
	box     = &(*iter)->box;
	box_x = box->xhigh-box->xlow; 
	box_y = box->yhigh-box->ylow; 
	box_z = box->zhigh-box->zlow;
	hbox_z = box_z/2.;

	for ( ; iter != snaplist.end();  iter++) {
		nsnap ++;
		maxAtom = (*iter)->n_atoms;
		atoms   = (*iter)->atoms;
		for ( ii =0; ii<maxAtom; ii++){
			ppi = &atoms[ii];
			if( ppi->type != type) continue;
			sz=ppi->muz;
			zi = ( ppi->z - hbox_z );
			z1 = abs( zi );

			bin = static_cast<int>( floor(z1/var_z));

			if (bin <maxbinz){
				numSum[bin] +=1;
				szSum[bin] += zi*sz;
			}
		}
	}
	real Area= box_x*box_y;
	real norm = 1./(2.*Dz*nsnap*Area);
	
	char filename1[100] ="P1z00.info" ;
	int nfile = 0;
	while( 0 == (access(filename1,F_OK)  )) {
		nfile++;
		sprintf(filename1, "P1z%02d.info",nfile);
	}

	FILE* fp = fopen(filename1,"w");
	vP1z.clear(); vP1z.resize(maxbinz);
	vPmuz1z.clear(); vPmuz1z.resize(maxbinz);
	vz.clear(); vz.resize(maxbinz);

	for (i=0; i<maxbinz; i++){
		zz = (i+0.5) * var_z;
		val  = double(numSum[i])*norm;
		val1 = double(szSum[i])*norm;
		fprintf(fp,"%lf %lf %lf\n", zz, val, val1);
		vz[i]   = zz;
		vP1z[i] = val;
		vPmuz1z[i] = val;
	}
	fclose(fp);
	delete numSum; delete szSum;
}

void makeRDF::calcP1s(int type) {
	puts("calcP1s() function is not completed.");
	exit (1);
	int i, ii, jj;
	int nsnap =0; int bin;
	real s1, zi, var_s=Dz;

	if (periodicity[0]&&periodicity[1] &&!periodicity[2])
		return ;

	bigint* numSum = new bigint[maxbins+1];
	double* DivMuSum = new double[maxbins+1];
	for(i=0; i<=maxbins; i++) {
		numSum[i] = DivMuSum[i] =0.;
	}

	real mus, val, val1;

	for (list<Snapshot*>::iterator iter =snaplist.begin() ; iter != snaplist.end();  iter++) {
		nsnap ++;
		maxAtom = (*iter)->n_atoms;
		atoms   = (*iter)->atoms;
		for ( ii =0; ii<maxAtom; ii++){
			ppi = &atoms[ii];
			mus= (ppi->mux * ppi->x + ppi->muy*ppi->y);
			s1= sqrt(ppi->x * ppi->x + ppi->y*ppi->y);
			mus /= s1;
			bin = static_cast<int>( floor(s1/var_s));

			if (bin <=maxbins){
				numSum[bin] +=1;
				DivMuSum[bin] += mus;
			}
		}
	}
	real norm = 1./(2.*M_PI*var_s*var_s*nsnap);
	
	char filename1[100] ="P1s00.info" ;
	int nfile = 0;
	while( 0 == (access(filename1,F_OK)  )) {
		nfile++;
		sprintf(filename1, "P1s%02d.info",nfile);
	}

	FILE* fp = fopen(filename1,"w");

	vP1s.clear(); vP1s.resize(maxbins+1);
	vPmus1s.clear(); vPmus1s.resize(maxbins+1);
	vs.clear(); vs.resize(maxbins+1);

	for (i=0; i<=maxbins; i++){
		mus = (i+0.5) * var_s;
		val  = double(numSum[i])*norm/ (i+0.5);
		val1 = double(DivMuSum[i])*norm/ (i+0.5);
		fprintf(fp,"%lf %lf %lf\n", mus, val, val1);
		vs[i]   = mus;
		vP1s[i] = val;
		vPmus1s[i] = val;
	}
	fclose(fp);
	delete numSum; delete DivMuSum;

}
void makeRDF::EvalSpacetimeCorr () {
	real b, c, c0, c1, c2,  s, s1, s2, w;
	int j, k, m, n, nb, nc, ni, nv;
	int i;
	real x,y,z,sx,sy,sz;
	real pos[3];
	real hbox_x=box_x/2.;
	real hbox_y=box_y/2.;
	real hbox_z=box_z/2.;

	// 8 means 
  // sum_s_x cos(kr), sum s_x sin(kr)	
  // sum_s_y cos(kr), sum s_y sin(kr)	
  // sum_s_z cos(kr), sum s_z sin(kr)	
	// sum_cos(kr), sum_sin(kr)
	// nFuncorr ?  nFuncorr*var_k*r
	error("Call EvalSpacetimeCorr ()");
	for (j =0; j < n_valST; j ++) valST[j] =0.;

	for ( i =0; i<maxAtom; i++){
		ppi = &atoms[i];
		x = ppi->x;y=ppi->y;z=ppi->z;
		sx = ppi->mux;sy=ppi->muy;sz=ppi->muz;
		//printf("%p x=%f y=%f z=%f sx=%f sy=%f sz=%f\n",ppi,x,y,z,sx,sy,sz);
		if ( periodicity[0]){
			while ( x <-hbox_x) x+=box_x;
			while ( x > hbox_x) x-=box_x;
		}
		if ( periodicity[1]){
			while ( y <-hbox_y) y+=box_y;
			while ( y > hbox_y) y-=box_y;
		}
		if ( periodicity[2]){
			while ( z <-hbox_z) z+=box_z;
			while ( z > hbox_z) z-=box_z;
		}
		pos[0]=x;pos[1]=y;pos[2]=z;
		j =0;
	//var_k -> var_k
		// build sin cos table
		for ( k=0; k<3; k++) { // k_[x,y,z]
			for ( m=0; m<nFunCorr; m++) {
				if (m == 0) { //k1_[x,y,z]
					b = var_k *  pos[k];
					c = cos (b);
					s = sin (b);
					c0= c;
				}else if(m==1) { //k2_[x,y,z]
					c1 = c;
					s1 = s;
					c= 2.* c0 * c1 - 1.;//cos2x=2*cos^2 x -1
					s = 2. * c0 * s1;   //sin2x=2*sinx*cosx
				}else { // ... 
					c2=c1;
					s2=s1;
					c1=c;
					s1=s;
					c = 2. * c0 * c1 - c2;
					s = 2. * c0 * s1 - s2;
				}
				valST[j ++] += sx * c;
				valST[j ++] += sx * s;
				valST[j ++] += sy * c;
				valST[j ++] += sy * s;
				valST[j ++] += sz * c;
				valST[j ++] += sz * s;
				valST[j ++] += 1. * c;
				valST[j ++] += 1. * s;
			}
		}
	} // loop all atoms

	for ( nb =0; nb < nBuffCorr; nb ++) { //t0=step1,step2..
		if ( tBuf[nb].count ==0) { // count==0 reset initial time value
			for ( j =0; j < n_valST; j++)
				tBuf[nb].orgST[j] = valST[j]; 
		}

		if ( tBuf[nb].count >=0) {
			for ( j=0; j<3 * nFunCorr; j++) 
				tBuf[nb].acfST[j][tBuf[nb].count] =0.; 

			j=0;
			for ( k=0; k <3; k++ ) {
				//k=0 (k_x,0,0) line sum 
				//k=1 (k_y,0,0) line sum
				for (m =0 ; m < nFunCorr; m ++) {
					for ( nc = 0 ; nc < 4; nc ++) {
						nv = 3*m +2;
						if (nc<3) {
							//w = Sqr (var_k * (m + 1));
							w = 1.;
							-- nv; //nc 0 1 2 spin, 3 density
							if (nc == k) -- nv;
							else w *= .5;
						} else w=1.;
						tBuf[nb].acfST[nv][tBuf[nb].count] +=
							 w * (valST[j] * tBuf[nb].orgST[j] +
									valST[j +1] * tBuf[nb].orgST[j+1]);
						j += 2;
					}
				}
			}
		} 
/****************************************************
acfST[0] += w *(                          nc=0 k=0
acfST[1] +=w/2*( w... part none           nc=1
acfST[1] +=w/2*( 0 1 1 2 -> all 1         nc=2
acfST[2] +=1.0*(                          nc=3
 0- longitudinal 1-transverse 2 density
	*****************************************************/
		++ tBuf[nb].count;
	}
	AccumSpacetimeCorr ();
}

void makeRDF::AccumSpacetimeCorr () {
	int j, n, nb;

	error("Call AccumSpacetimeCorr()");
	for( nb =0; nb < nBuffCorr; nb ++) {
		if (tBuf[nb].count == nValCorr) {
			tBuf[nb].count =0;

			for (j=0;j<3*nFunCorr;j++) 
				for (n=0; n< nValCorr; n ++) 
					avAcfST[j][n] += tBuf[nb].acfST[j][n];
			
			++ countCorrAv;
			if ( countCorrAv == limitCorrAv) {
				for (j=0;j<3*nFunCorr;j++) 
					for (n=0; n< nValCorr; n ++) 
						avAcfST[j][n] /= 3. * maxAtom * limitCorrAv;
				
				PrintSpacetimeCorr ();
				ZeroSpacetimeCorr ();
			}
		}
	}
}
void makeRDF::ZeroSpacetimeCorr () {
	error("Call ZeroSpacetimeCorr ()");
	int j, n;
	countCorrAv =0;
	for (j =0; j<3*nFunCorr;j++) 
		for (n=0; n< nValCorr; n ++) 
			avAcfST[j][n] = 0.;
}


void makeRDF::calcISF () {
	for (list<Snapshot*>::iterator iter =snaplist.begin() ; iter != snaplist.end();  iter++) {
		atoms   = (*iter)->atoms;
		EvalSpacetimeCorr ();
	}
}

void makeRDF::InitSpacetimeCorr(){
	int nb;

	n_valST = 24*nFunCorr;
	valST = new real [n_valST];

	avAcfST = new real* [3*nFunCorr];
	avAcfST[0] = new real [ 3*nFunCorr*nValCorr];
	for (int k=1; k<3*nFunCorr; k++) avAcfST[k] = avAcfST[k-1] +nValCorr;

	tBuf = new TBuf[nBuffCorr];
	for(int nb=0; nb<nBuffCorr; nb ++) {
		tBuf[nb].orgST = new real[n_valST];
		tBuf[nb].acfST = new real* [3*nFunCorr]; // max nv = 3*nFunCorr
		tBuf[nb].acfST[0] = new real [ 3*nFunCorr*nValCorr];
		for (int k=1; k<3*nFunCorr; k++) 
			tBuf[nb].acfST[k] = tBuf[nb].acfST[k-1] +nValCorr;
	}
	for (nb =0; nb < nBuffCorr; nb ++)
		tBuf[nb].count = - nb * nValCorr/ nBuffCorr;
	ZeroSpacetimeCorr();
}
char* makeRDF::filename_template= "Corr_%s.out%ld";
makeRDF::makeRDF(list<Snapshot*> &_sl) { 
	// C++99 not allow non-const static member init
	flag_anisotropy=0;
	flag_SSF_from_g=0;
	nFunCorr = 4;
	limitCorrAv=200;
	nBuffCorr=10;   // nValCorr = number x nBuffCorr
	nValCorr=500; // # of average sets
	step,stepCorr=1;
	flag_anisotropy=0;
	flag_SSF_from_g=0;
	deltaT = 0.0001;
	maxbin = 200;
	maxbinz = 5;
	maxbins = 10;
	maxbinq = 100;
	


	snaplist = _sl;
	list<Snapshot*>::iterator firstsnap = (snaplist.begin());
	firstsnap++;
	box     = &((*firstsnap)->box);
	double large_r = 1000, min_r, max_r;
	box_x = box->xhigh-box->xlow; 
	box_y = box->yhigh-box->ylow; 
	box_z = box->zhigh-box->zlow;
	periodicity = box->pbc;
	maxSnap = snaplist.size();
	maxAtom = (*firstsnap)->n_atoms;
	step    = (*firstsnap)->timestep ; firstsnap ++;
	stepCorr= abs( (*firstsnap)->timestep - step);
	printf("step Corr = %d %p %p\n", stepCorr,&(*firstsnap),&(*snaplist.begin()));
	if(periodicity[0]||periodicity[1]||periodicity[2]) {
		min_r = large_r;
		if( periodicity[0])  min_r = box_x < min_r ? box_x:min_r;
		if( periodicity[1])  min_r = box_y < min_r ? box_y:min_r;
		if( periodicity[2])  min_r = box_z < min_r ? box_z:min_r;
	} else {
		min_r= (box_x<box_y)?((box_x<box_z)?box_x:box_z):((box_y<box_z)?box_y:box_z);
	}
	r_cut= min_r/2.;
	max_r= (box_x>box_y)?((box_x>box_z)?box_x:box_z):((box_y>box_z)?box_y:box_z);
	q_cut= M_PI / max_r;
	dq   = 2.*M_PI/ maxbinq;
	var_r = r_cut/maxbin;
	var_k = (4. /nFunCorr )*M_PI;
	InitSpacetimeCorr ();
};
void makeRDF::PrintSpacetimeCorr (){
	static int filenumber = 0;
	char *header[] = {"longi", "trans", "density"};
	sprintf(filename0,filename_template,header[0],filenumber);
	sprintf(filename1,filename_template,header[1],filenumber);
	sprintf(filename2,filename_template,header[2],filenumber);
	fp[0] = fopen( filename0, "w");
	fp[1] = fopen( filename1, "w");
	fp[2] = fopen( filename2, "w");

	real tVal;
	int j, k, n;
	real qVal;
	//		fprintf (fp, "space-time corr\n");
	for (k=0; k<3; k++) {
		fprintf (fp[k], "%s", header[k]);
		for (j =0; j < nFunCorr; j ++ ) {
			qVal = (j+1)* var_k;
			fprintf(fp[k], " %8.4f", qVal);
		}
		fprintf(fp[k],"\n");

		for (n = 0; n < nValCorr; n ++) {
			tVal = n * stepCorr * deltaT;
			//	fprintf (stdout, "%7.3f %d %d %f\n", tVal,n,stepCorr,deltaT);
			fprintf (fp[k], "%7.3f", tVal);
			for (j =0; j < nFunCorr; j ++ )
				fprintf(fp[k], " %8.4f", avAcfST[3*j+k][n]);
			fprintf(fp[k],"\n");
		}
	}
	filenumber ++;
	fclose(fp[0]);
	fclose(fp[1]);
	fclose(fp[2]);
}
