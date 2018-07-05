/*!
 *    \file  makeRDF_spacetime.cpp
 *   \brief  
 *
 *  <+DETAILED+>
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

typedef std::vector<Snapshot*>::const_iterator citer;

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
void makeRDF::ZeroSpacetimeCorr () {
	error("Call ZeroSpacetimeCorr ()");
	int j, n;
	countCorrAv =0;
	for (j =0; j<3*nFunCorr;j++) 
		for (n=0; n< nValCorr; n ++) 
            avAcfST[j][n] = 0.;
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
void makeRDF::calcISF () {
	for (vector<Snapshot*>::iterator iter =snaplist.begin() ; iter != snaplist.end();  iter++) {
		atoms   = (*iter)->atoms;
		EvalSpacetimeCorr ();
	}
}
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
