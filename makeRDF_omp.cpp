#include "makeRDF.h"
#include <unistd.h>
#include <memory>
#include <vector>
#include <algorithm> // std::fill
#include "kh_math_fourier.h"
#include <omp.h>
#define Dz     (2)
#define halfDz     (Dz/2.0)
#define R_CUT   120
#include "hydro_math.h"

void makeRDF::calcSSF (T_RDF rdftype) {
	if(flag_SSF_from_g)
		calcSSF_from_g (rdftype);
	else {
	  int bin,i,iii,jjj,kkk;	
		real qx,qy,qz,q1;	
		real x,y,z;
		int id;
		real* qlist = new real [maxbinq];
		real* hist  = new real [maxbinq+2];
		real* S_qr  = new real [maxbinq+2];
		real* S_qi  = new real [maxbinq+2];
		FILE* fp_out;
		fp_out = fopen("SSF_direct.info", "w");

		dq = 2*M_PI/64;
		for ( int i = 0; i <= maxbinq; i += 1 ){
			qlist[i] =  dq * i;
		}

		puts("");
		int files=0;
		for (citer iter =snapbegin ; iter != snapend;  iter++) {
			maxAtom     = (*iter)->n_atoms;
			atom* atoms = (*iter)->atoms;
			box         = &(*iter)->box;
			box_x       = box->xhigh-box->xlow;
			box_y       = box->yhigh-box->ylow;
			box_z       = box->zhigh-box->zlow;
			files++;
			for ( id =0; id<maxAtom; id++){
			  atom* ppi = &atoms[id];
				printf("\r%5dth-Snapshots(%ld) %4dth-atoms(%d)", files,snaplist.size(), id+1,maxAtom);
				fflush(stdout);
				x = ppi->x; y = ppi->y; z = ppi->z;
				for ( iii = 0; iii <= maxbinq; iii += 1 ) {
					qx = qlist[iii];		
					qy = qlist[iii];		
					qz = qlist[iii];		
					bin = iii;
					hist[bin] += 1;
					S_qr[bin] += cos( qx* x);
					S_qi[bin] += sin( qx* x);

					hist[bin] += 1;
					S_qr[bin] += cos( qy*y);
					S_qi[bin] += sin( qy*y);

					hist[bin] += 1;
					S_qr[bin] += cos( qz* z);
					S_qi[bin] += sin( qz* z);



					
				}
			}
		}
		
		real norm_q = 4.0 * M_PI * dq * maxAtom * snaplist.size();

		real qqq, value,value2;
		for ( i = 1; i <= maxbinq; i += 1 ){
			qqq = qlist[i];
			if (hist[i] ==0.) {
				value =0;
				value2=0;
			}
			else {
/* 				value = sqrt(S_qr[i]*S_qr[i] + S_qi[i]*S_qi[i])/ norm_q/ ((qqq*qqq)+(dq*dq/12.0) ) ;
 * 				value2= sqrt(S_qr[i]*S_qr[i] + S_qi[i]*S_qi[i])/ maxAtom / hist[i]; 
 */
				value = (S_qr[i]*S_qr[i]+ S_qi[i]*S_qi[i]) / maxAtom/snaplist.size();
			}
			
			fprintf(fp_out, "%le ", qqq);
			fprintf(fp_out, "%le ", value);
//			fprintf(fp_out, "%le ", value2);
			fprintf(fp_out, "\n");

		}

		fclose (fp_out);
		delete [] hist;
		delete [] S_qr;
		delete [] S_qi;
	}
}

void makeRDF::calcSSF_from_g(T_RDF rdftype) {
	real Vol = box_x*box_y*box_z;
	real phi = maxAtom / Vol;

#pragma omp parallel for 
	for(int i=0; i<maxbinq ; i++) {
		real q         = i * dq;
		ca_q_radius[i] = q;
		real h000_k    = trapz_iso3D_forward_n( ca_radius,ca_h000,q,maxbin);
		ca_h000_q[i]   = h000_k ;
		if (rdftype == ANISO) {
			ca_h110_q[i] = trapz_iso3D_forward_n( ca_radius,ca_h110,q,maxbin);
			ca_h220_q[i] = trapz_iso3D_forward_n( ca_radius,ca_h220,q,maxbin);
			ca_h112_q[i] = trapz_iso3D_forward_with_l2_n( ca_radius,ca_h112,q,maxbin);
		}
		real S_q  = 1.+ phi*h000_k;
		ca_S_q[i] = S_q;
		ca_c_q[i] = h000_k / S_q;
	}
	char ssf_filename[100];
	sprintf(ssf_filename,"SSF.info%s", ext);
	FILE* fp_S = fopen(ssf_filename,"w");
	fprintf(fp_S,"##q S_k h_k c_k h112_k\n");
	for(int i=0; i<maxbinq ; i++) {
		fprintf(fp_S,"%lf %lf %lf %lf %lf\n", ca_q_radius[i], 
				ca_S_q[i], ca_h000_q[i], ca_c_q[i], ca_h112_q[i]);
	}
	fclose(fp_S);
	/*-----------------------------------------------------------------------------
	 * we get c_r from c_k using Backward Fourier transform
	 *-----------------------------------------------------------------------------*/

	
#pragma omp parallel for 
	for(int i=0; i<maxbin ; i++) {
		real r_     = ca_radius[i];
		ca_c_r[i]   = trapz_iso3D_backward_n(ca_q_radius, ca_c_q,
				r_,maxbinq);
		ca_h_mod[i] = trapz_iso3D_backward_n(ca_q_radius, ca_h000_q,
				r_,maxbinq);
	}

	FILE* fp_c_r = fopen("DCF.info","w");
	fprintf(fp_c_r,"##r c_r h_mod\n");
	for(int i=0; i<maxbin ; i++) {
		fprintf(fp_c_r,"%lf %lf %lf \n", ca_radius[i], 
				ca_c_r[i], ca_h_mod[i] );
	}
	fclose(fp_c_r);
	
}

void makeRDF::calcRDF() {
	if (flag_anisotropy )
		calcRDF_anisotropy();
	else 
		calcRDF_isotropy();
}
void makeRDF::calcRDF_anisotropy(){
	calcRDF_inter_type (1,1, ANISO);
}
void makeRDF::calcRDF_isotropy(){
	calcRDF_inter_type (1,1, ISO);
}
void makeRDF::calcRDF_anisotropy(int i, int j){
	calcRDF_inter_type (i,j, ANISO);
}
void makeRDF::calcRDF_isotropy(int i, int j){
	calcRDF_inter_type (i,j, ISO);
}

void makeRDF::calcRDF_inter_type (int itype, int jtype, T_RDF rdftype) {
	int ltype=   itype<jtype?itype:jtype;
	int gtype=   itype>jtype?itype:jtype;
	bool flag_equaltype;
	if (itype == jtype) flag_equaltype=true;
	else flag_equaltype=false;

	int  maxallbin, allbin;
	maxallbin = maxbinz*maxbinz*maxbin;
	real abszi,abszj,var_z=Dz,var_r = r_cut/maxbin;

	typedef long *   p_long;
	typedef double *  p_double;
	real  si[3],sj[3];

/* 	long*	hist000    = new long[rbin_t];
 * 	double*	hist110    = new double[rbin_t];
 * 	double*	hist112    = new double[rbin_t];
 * 	double*	hist220    = new double[rbin_t];
 */
	vector<long>	hist000(rbin_t);
	vector<double>	hist110(rbin_t);
	vector<double>	hist112(rbin_t);
	vector<double>	hist220(rbin_t);
	/* 	histcyl000 = new long[(rbin_t)];
	 * 	histcyl110 = new double[(rbin_t)];
	 * 	histcyl112 = new double[(rbin_t)];
	 * 	histcyl220 = new double[(rbin_t)];
	 */
	std::fill(hist000.begin(),hist000.end(), 0l);
	std::fill(hist110.begin(),hist110.end(), 0.0);
	std::fill(hist112.begin(),hist112.end(), 0.0);
	std::fill(hist220.begin(),hist220.end(), 0.0);


	

/* 	omp_lock_t writelock;
 * 	omp_init_lock(&writelock);
 */
	int progress_index = 0 ;
	int progress= 0;
	int progress_index_max = snaplist.size();
#pragma omp parallel for 
	for (int nsnap=0 ; nsnap< snaplist.size();  nsnap++) {
		atom* atoms   = snaplist[nsnap]->atoms;

/* 		long * local_hist000 = new long[maxbin+5];
 * 		real * local_hist110 =   new double[maxbin+5];
 * 		real * local_hist112 =   new double[maxbin+5];
 * 		real * local_hist220 =   new double[maxbin+5];
 */
		vector<long> local_hist000 (maxbin+5);
		vector<real> local_hist110 (maxbin+5);
		vector<real> local_hist112 (maxbin+5);
		vector<real> local_hist220 (maxbin+5);
		for(int i=0; i<=maxbin; i++) {
			std::fill(local_hist000.begin(),local_hist000.end(), 0l);
		}
		if (rdftype == ANISO) {
			std::fill(local_hist110.begin(),local_hist110.end(), 0.0);
			std::fill(local_hist112.begin(),local_hist112.end(), 0.0);
			std::fill(local_hist220.begin(),local_hist220.end(), 0.0);
		}
		int Box_replica= round(r_cut / box_x * 2.) ;
		for ( int  ii =0; ii<maxAtom; ii++){
			atom* ppi = &atoms[ii];
			if ( ppi->type != itype) continue;
			real si[3]= {ppi->mux,ppi->muy,ppi->muz};
			//			zi = ppi->z-hbox_z;
			for (int  jj =0; jj<maxAtom; jj++) { 
				if ( ii==jj ) //continue ;                // self term out
				{
					// LOOK AT ME ... 
					continue;
					if (Box_replica == 0 ) continue;
					for (int xxx=-Box_replica; xxx <=Box_replica ; xxx++) {
						real xij =  xxx * box_x;
						if (abs(xij) > r_cut ) continue;
						for (int yyy=-Box_replica; yyy <=Box_replica ; yyy++) {
							real yij =  yyy * box_y;
							if (abs(yij) > r_cut ) continue;
							for (int zzz=-Box_replica; zzz <=Box_replica ; zzz++) {
								if ( (xxx == 0) && (yyy == 0) && (zzz == 0) ) continue; // self term out
								real zij =  zzz * box_z;
								if (abs(zij) > r_cut ) continue;

									real r1 = sqrt(xij*xij + yij*yij + zij*zij);

									if (r1 <r_cut){ 
										int bin    = floor( r1/var_r);
										local_hist000[bin] +=1;
										if (rdftype == ANISO) {
											real rij[3] = {xij,yij,zij};
											real sj[3]= {si[0], si[1], si[2] };

											real si_dot_sj      = VDOT3(si,sj);
											real si_dot_rij     = VDOT3(si,rij);
											real sj_dot_rij     = VDOT3(sj,rij);
											local_hist110[bin] += si_dot_sj;
											local_hist112[bin] += 3.0*si_dot_rij*sj_dot_rij/(r1*r1)-si_dot_sj;
											local_hist220[bin] += 3.0*(si_dot_sj*si_dot_sj) -1.;
										}
									}
							}
						}
					}
				}
				else {
					atom* ppj = &atoms[jj];
					if ( ppj->type != jtype) continue;
					//				zj = ppj->z - hbox_z;

					real xij_orig = ppi->x - ppj->x;
					if ( periodicity[0]){
						while ( xij_orig <-hbox_x) xij_orig+=box_x;
						while ( xij_orig > hbox_x) xij_orig-=box_x;
					}

					real yij_orig = ppi->y - ppj->y;
					if ( periodicity[1]){
						while ( yij_orig <-hbox_y) yij_orig+=box_y;
						while ( yij_orig > hbox_y) yij_orig-=box_y;
					}

					real zij_orig = ppi->z - ppj->z;
					if ( periodicity[2]){
						while ( zij_orig <-hbox_z) zij_orig+=box_z;
						while ( zij_orig > hbox_z) zij_orig-=box_z;
					}

					for (int xxx=-Box_replica; xxx <=Box_replica ; xxx++) {
						real xij =  xxx*box_x + xij_orig;
						if (abs(xij) > r_cut ) continue;
						for (int yyy=-Box_replica; yyy <=Box_replica ; yyy++) {
							real yij =  yyy*box_y + yij_orig;
							if (abs(yij) > r_cut ) continue;
							for (int zzz=-Box_replica; zzz <=Box_replica ; zzz++) {
								real zij =  zzz*box_z + zij_orig;
								if (abs(zij) > r_cut ) continue;
								real r1 = sqrt(xij*xij + yij*yij + zij*zij);

								if (r1 <r_cut){ 
									int bin    = floor( r1/var_r);
									local_hist000[bin] +=1;
									if (rdftype == ANISO) {
										real rij[3] = {xij,yij,zij};
										real sj[3]= {ppj->mux,ppj->muy,ppj->muz};

										real si_dot_sj      = VDOT3(si,sj);
										real si_dot_rij     = VDOT3(si,rij);
										real sj_dot_rij     = VDOT3(sj,rij);
										local_hist110[bin] += si_dot_sj;
										local_hist112[bin] += 3.0*si_dot_rij*sj_dot_rij/(r1*r1)-si_dot_sj;
										local_hist220[bin] += 3.0*(si_dot_sj*si_dot_sj) -1.;
									}
								}
							}
						}
					}

				} 

			}// jj < maxAtom
		}// ii < maxAtom
#pragma omp critical 
		{
			for( int i=0; i<=maxbin; i++) {
				hist000[i] +=	local_hist000[i] ;
			}
			if (rdftype == ANISO) {
				for(int i=0; i<=maxbin; i++) {
					hist110[i] +=local_hist110[i] ;
					hist112[i] +=local_hist112[i]; 
					hist220[i] +=local_hist220[i] ;
				}
			}
			progress_index ++;
			int new_progress = 1000.0*progress_index/progress_index_max;
			if (  new_progress != progress) {
				progress = new_progress;
				fprintf(stderr, "\r %4.1f%%(", progress*.1);
			}
		}  // omp parallel critical


	} /// pragma omp for

//	printf("thread num : %d\n", omp_get_thread_num());
//	omp_destroy_lock(&writelock);
//
	int n_itype=0, n_jtype=0;	
/* 	for (int  ii =0; ii<maxAtom; ii++){
 * 		atom *ppi = &first_atoms[ii];
 * 
 * 		if ( ppi->type == itype) n_itype++;
 * 		if ( ppi->type == jtype) n_jtype++;
 * 	}
 */
	n_itype = (*snapbegin)->get_n_ptls(itype);
	n_jtype = (*snapbegin)->get_n_ptls(jtype);

	real val000, val110, val112, val220,valszz1, valszz2;
	int nsnap    = snaplist.size();

	real Vol     = box_x*box_y*box_z;
	real Area    = box_x*box_y;
	real phi_i,phi_j;
	if (flag_equaltype) {
		phi_i   = n_itype / Vol;
		phi_j   = (n_jtype-1) / Vol;
	}
	else {
		phi_i   = n_itype / Vol;
		phi_j   = n_jtype / Vol;
	}
	real norm000 = 1./ (4.0* M_PI *var_r * phi_i * phi_j * nsnap* Vol);
	real norm110 = norm000*3.0;
	real norm112 = norm000*(2./3.);
	real norm220 = norm000*2.5;

	real normcyl000 = 2.* 2.0* M_PI * Sqr(var_z) * Sqr(var_r)  *  nsnap* Area;
	//  11 -1-1    1 -1   -1 1

	ca_g000[0] = 0;
	ca_h000[0] = -1;
	ca_h110[0] = 0;
	ca_h112[0] = 0;
	ca_h220[0] = 0;
	real drdr3 = var_r*var_r/3.0;
	if (rdftype == ANISO) {
#pragma omp parallel for
		for (int i=0; i<=maxbin; i++){
			real rrr     = (i + 0.5) * var_r;
			real rr      = (i ) * var_r;
			real invr2   = 1./(rr*rr+ rr*var_r + drdr3 );
			ca_radius[i] = rrr;
			ca_g000[i]   = double(hist000[i])*norm000*invr2;
			ca_h000[i]   = ca_g000[i]-1.0;
			ca_h110[i]   = hist110[i]*norm110*invr2;
			ca_h112[i]   = hist112[i]*norm112*invr2;
			ca_h220[i]   = hist220[i]*norm220*invr2;
		}
	} else {
#pragma omp parallel for
		for (int i=0; i<=maxbin; i++){
			real rrr     = (i + 0.5) * var_r;
			real rr      = (i ) * var_r;
			real invr2   = 1./(rr*rr+ rr*var_r + drdr3 );
			ca_radius[i] = rrr;
			ca_g000[i]   = hist000[i]*norm000*invr2;
			ca_h000[i]   = ca_g000[i]-1.0;
		}
	}
	printRDF( itype,jtype, rdftype);
}
void makeRDF::printRDF (int itype, int jtype, T_RDF rdftype )
{
	char rdffilename[100] ; 
	sprintf(rdffilename, "rdf00_ij%d%d.info%s",itype, jtype, ext);
	if (rdftype == ANISO) {
		sprintf(rdffilename, "rdf11_ij%d%d.info%s",itype, jtype, ext);
	}	
	FILE* fp_rdf = fopen(rdffilename,"w");
/* 	sprintf(rdffilename, "rdfcyl_ij%d%d.info",itype, jtype);
 * 	FILE* fp_rdfcyl = fopen(rdffilename,"w");
 */

	fprintf(stderr, "makeRDF::printRDF %d - %d (type : %d)\n", itype,jtype,(int)rdftype);
	if ( (periodicity[2])) 
		fputs("##3d periodicity\n", fp_rdf);
	else
		fputs("##2d periodicity and z axis bounded\n", fp_rdf);
	if (rdftype == ANISO) {
		for (int i=0; i<maxbin; i++){
			fprintf(fp_rdf,"%lf %lf %lf %lf %lf\n", ca_radius[i],ca_g000[i],ca_h110[i],ca_h112[i],ca_h220[i]);
		}
	} else {
		for (int i=0; i<=maxbin; i++){
			fprintf(fp_rdf,"%lf %lf\n", ca_radius[i],ca_g000[i]);
		}
	}

	fclose(fp_rdf);
/* 	fclose(fp_rdfcyl);
 */
}
// P1(z) = A*P(x,y,z)
void makeRDF::calcP1z( int type) {
	int i, ii, jj;
	int nsnap =0; int bin;
	real z1, zi, var_z=Dz;

	if ( periodicity[2])
		fprintf(stderr, "WARNING : z-direction is not periodic.\n");

	long* numSum = new long[maxbinz+1];
	double* szSum = new double[maxbinz+1];
	for(i=0; i<=maxbinz; i++) {
		numSum[i] = szSum[i] =0.;
	}

	real sz,zz, val, val1, hbox_z;


	for (	citer iter = snapbegin ; iter != snapend;  iter++) {
		nsnap ++;
		maxAtom     = (*iter)->n_atoms;
		atom* atoms = (*iter)->atoms;
		for ( ii =0; ii<maxAtom; ii++){
			atom* ppi = &atoms[ii];
			if( ppi->type != type) continue;
			sz  = ppi->muz;
			zi  = ( ppi->z - hbox_z );
			z1  = abs( zi );

			bin = static_cast<int>( floor(z1/var_z));

			if (bin <maxbinz){
				numSum[bin] += 1;
				szSum[bin]  += zi*sz;
			}
		}
	}
	real Area = box_x*box_y;
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
		zz         = (i+0.5) * var_z;
		val        = double(numSum[i])*norm;
		val1       = double(szSum[i])*norm;
		fprintf(fp,"%lf %lf %lf\n", zz, val, val1);
		vz[i]      = zz;
		vP1z[i]    = val;
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

	long* numSum = new long[maxbins+1];
	double* DivMuSum = new double[maxbins+1];
	for(i=0; i<=maxbins; i++) {
		numSum[i] = DivMuSum[i] =0.;
	}

	real mus, val, val1;

	for (citer iter =snapbegin ; iter != snapend;  iter++) {
		nsnap ++;
		maxAtom     = (*iter)->n_atoms;
		atom* atoms = (*iter)->atoms;
		for ( ii = 0; ii<maxAtom; ii++){
			atom* ppi  = &atoms[ii];
			mus        = (ppi->mux * ppi->x + ppi->muy*ppi->y);
			s1         = sqrt(ppi->x * ppi->x + ppi->y*ppi->y);
			mus       /= s1;
			bin        = static_cast<int>( floor(s1/var_s));

			if (bin         <= maxbins){
				numSum[bin]   += 1;
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
		mus        = (i+0.5) * var_s;
		val        = double(numSum[i])*norm/ (i+0.5);
		val1       = double(DivMuSum[i])*norm/ (i+0.5);
		fprintf(fp,"%lf %lf %lf\n", mus, val, val1);
		vs[i]      = mus;
		vP1s[i]    = val;
		vPmus1s[i] = val;
	}
	fclose(fp);
	delete numSum; delete DivMuSum;

}

void makeRDF::freeMem() 
{

	delete [] 	ca_radius ;
	delete [] 	ca_g000 ;
	delete [] 	ca_h110 ;
	delete [] 	ca_h112 ;
	delete [] 	ca_h220 ;
	/* delete []	ca_gcyl000 ;
		 delete []  * 	ca_hcyl110 ;
		 delete []  * 	ca_hcyl112 ;
		 delete []  * 	ca_hcyl220 ;
		 */
	/*!
	 *  \brief Reserving memory for  q_space calcuation 
	 */
	delete [] 	ca_q_radius ;
	delete [] 	ca_c_q ;
	delete [] 	ca_h000_q ;
	delete [] 	ca_h110_q ;
	delete [] 	ca_h112_q ;
	delete [] 	ca_h220_q ;
	delete [] 	ca_S_q ;

	delete [] 	ca_c_r ;
	delete [] 	ca_h_mod ;

}
void makeRDF::AllocMem() 
{
	rbin_t = maxbin+5;

	ca_radius = new double[rbin_t];
	ca_g000   = new double[rbin_t];
	ca_h000   = new double[rbin_t];
	ca_h110   = new double[rbin_t];
	ca_h112   = new double[rbin_t];
	ca_h220   = new double[rbin_t];
/* 	ca_gcyl000 = new double[rbin_t];
 * 	ca_hcyl110 = new double[rbin_t];
 * 	ca_hcyl112 = new double[rbin_t];
 * 	ca_hcyl220 = new double[rbin_t];
 */
	for(int i=0; i<rbin_t; i++) {
		ca_radius[i] = 0.;
		ca_g000[i]   = 0 ;
		ca_h000[i]   = 0 ;
		ca_h110[i]   = 0 ;
		ca_h112[i]   = 0 ;
		ca_h220[i]   = 0 ;
	}
	/*!
	 *  \brief Reserving memory for  q_space calcuation 
	 */
	ca_q_radius = new double [maxbinq];
	ca_c_q      = new double[maxbinq];
	ca_h000_q   = new double[maxbinq];
	ca_h110_q   = new double[maxbinq];
	ca_h112_q   = new double[maxbinq];
	ca_h220_q   = new double[maxbinq];
	ca_S_q      = new double[maxbinq];
	
	ca_c_r   = new double[rbin_t];
	ca_h_mod = new double[rbin_t];

}

std::string  makeRDF::filename_template= "Corr_%s.out%ld";
makeRDF::makeRDF(vector<Snapshot*> &_sl) :
		snaplist(_sl)
{ 
	// C++99 not allow non-const static member init
//	InitSpacetimeCorr ();
	StartPreProcess ();
	AllocMem();
};
void makeRDF::StartPreProcess() {
	// Get information from first snapshot.
	snapbegin = (snaplist.begin());
	snapend = (snaplist.end());

	first_atoms     = (*snapbegin)->atoms;
	step            = (*snapbegin)->timestep ;

	maxAtom         = (*snapbegin)->n_atoms;
	maxSnap         = snaplist.size();

	box             = &((*snapbegin)->box);
	box_x           = box->xhigh-box->xlow;
	box_y           = box->yhigh-box->ylow;
	box_z           = box->zhigh-box->zlow;
	periodicity     = box->pbc;

	double large_r = 1000, min_L, max_L;
	hbox_x= box_x/2.; hbox_y= box_y/2.; hbox_z= box_z/2.;
	// Get time step between first two snapshots.
	snapbegin++;
	stepCorr= abs( (*snapbegin)->timestep - step);

	printf("step Corr = %d %p\n", stepCorr,&(*snapbegin));
	if(periodicity[0]||periodicity[1]||periodicity[2]) {
		min_L = large_r;
		if( periodicity[0])  
			min_L = box_x < min_L ? box_x:min_L;
		if( periodicity[1])  
			min_L = box_y < min_L ? box_y:min_L;
		if( periodicity[2])  
			min_L = box_z < min_L ? box_z:min_L;
	} else {
		min_L= (box_x<box_y)?((box_x<box_z)?box_x:box_z):((box_y<box_z)?box_y:box_z);
	}
//	r_cut= min(32.0,min_L/2.);
	r_cut =  R_CUT;
	max_L = (box_x>box_y)?((box_x>box_z)?box_x:box_z):((box_y>box_z)?box_y:box_z);
	q_cut = M_PI / max_L;
	dq    = 2.*M_PI/ maxbinq;
	dq    = 2.*M_PI/ max_L;
	var_r = r_cut/maxbin;
	var_k = (4. /nFunCorr )*M_PI;

}
void makeRDF::StartMainProcess() {
	var_r = r_cut/maxbin;
	AllocMem();
}
void makeRDF::hydro_init() {
		real vol  = box->getVol();
		real a = 4.024635;

		hydro_function = new Hydrodynamic_Function();
	  hydro_function->init(maxbinq,maxbin,maxAtom, var_r,var_k,a,vol,ca_g000,ca_h000);
}
void makeRDF::hydro_run(){
	hydro_function->run();
};
void makeRDF::hydro_print() 
{
/* 	char mobilityfilename[100];
 * 	sprintf(mobilityfilename, "rdf_mobility.info%s", ext);
 * 	FILE* fp_mobility = fopen(mobilityfilename,"w");
 * 	fprintf(fp_mobility,"radius g000 h000 Ds=%lf", hydro_function->getD_s());
 * 	for (int i=0; i<=maxbin; i++){
 * 		fprintf(fp_mobility,"%lf %lf %lf\n", ca_radius[i],ca_g000[i], ca_h000[i]);
 * 	}
 * 	fclose(fp_mobility);
 */


	char mobility2filename[100];
	sprintf(mobility2filename, "rdf_mobility2.info%s", ext);
	FILE* fp_mobility2 = fopen(mobility2filename,"w");

	fprintf(fp_mobility2,"#q S_k h000 cq  c112 h(q) Ds=%lf\n", hydro_function->getD_s());
	double q_back = -1.;
	for (int i=0; i<=maxbinq; i++){
		if (q_back < ca_q_radius[i]) {
			fprintf(fp_mobility2,"%lf %lf %lf %lf %lf %lf\n", ca_q_radius[i], 
					ca_S_q[i], ca_h000_q[i], ca_c_q[i], ca_h112_q[i], hydro_function->getHq()[i]);
			q_back = ca_q_radius[i];
		}
	}
	fclose(fp_mobility2);
};
void makeRDF::hydro_end(){
	hydro_function->end();
}
