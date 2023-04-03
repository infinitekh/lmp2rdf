/*
 * =====================================================================================
 *
 *       Filename:  a.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2023년 03월 27일 01시 31분 28초
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <string>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <array>
#include <vector>
#include <random>
using std::string;
using std::stringstream;
using std::cout;
using std::endl;
using std::array;
using std::vector;
typedef struct real3 {
	double x,y,z;
}real3;

#define CUT (10.)

int main(int argc, char* argv[])
{
	char*  ldfile ( argv[1]);
	char*  srdfile (argv[2]);
	 std::random_device rd;

  // random_device 를 통해 난수 생성 엔진을 초기화 한다.
  std::mt19937 gen(rd());

  // 0 부터 99 까지 균등하게 나타나는 난수열을 생성하기 위해 균등 분포 정의.

	FILE* ldfp = fopen(ldfile,"r");
	FILE* srdfp = fopen(srdfile,"rw");

	FILE* ofp = fopen("output", "w");

	char* line1 = NULL;
	char* line2 = NULL;
	size_t len1 = 0 ;
	size_t len2 = 0 ;

	getline(&line1, &len1, ldfp);
	getline(&line2, &len2, srdfp);
	puts(line1);
	fputs(line2, ofp);

	//blank line
	getline(&line1, &len1, ldfp);
	getline(&line2, &len2, srdfp);

	fputs(line2, ofp);
	// atoms/
	getline(&line1, &len1, ldfp);
	getline(&line2, &len2, srdfp);

	fputs(line2, ofp);
	stringstream stream1, stream2;

	stream1.str(line1);
	stream2.str(line2);

	int natom1; 
	int natom2;

	stream1 >> natom1;
	stream2 >> natom2;

	if ( natom1>natom2) {
		std::cout << "Error natom1>natom2"<< std::endl;
		exit(1);
	}

	getline(&line1, &len1, ldfp);
	getline(&line2, &len2, srdfp);
	stream1.str(line1);
	stream2.str(line2);

	fputs(line2, ofp);

	int int1, int2;

	stream1 >> int1;
	stream2 >> int2;
	if ( int1>=int2) {
		puts("Error atomstypes 1>=2");
		exit(1);
	}
	std::array<double,6> d1;
	std::array<double,6> d2;
	double dt1 ,dt2;
	//
	// blank
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);

	for (int i=0; i<6; i=i+2){
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);
/* 		printf("<<%s\n" ,line1);
 * 		printf(">>%s\n" ,line2);
 * 		stream1.str(line1);
 * 		stream2.str(line2);
 * 		stream1 >> dt1; d1[i] = dt1;
 * 		stream2 >> dt2; d2[i] = dt2;
 * 		i=i+1;
 * 		stream1 >> dt1; d1[i] = dt1;
 * 		stream2 >> dt2; d2[i] = dt2;
 */

		sscanf(line1, "%le %le ", &d1[i],&d1[i+1]);
		sscanf(line2, "%le %le ", &d2[i],&d2[i+1]);

	}

	if ( d1!= d2) {
		puts("Error box type different");
		for ( auto i : d1) cout << i;
		for ( auto i : d2) cout << i;
		exit(1);
	}

		// blank
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);

		//Masses
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);

		//blank
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);

		//Mass1
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);
		if (string(line1) != string(line2) ) {
			puts("Error Mass1 different");
			exit(1);
		}

		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);

		//blank
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);


		//Atoms
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);
		if (string(line1) != string(line2) ) {
			puts("Error not Atoms # hybrid");
			exit(1);
		}
		//blank
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);

		auto atom_index=0;
		vector<real3> pos;
		
		for (atom_index  = 1 ; atom_index <= natom1;atom_index ++) { 
			getline(&line1, &len1, ldfp);
			getline(&line2, &len2, srdfp);
			fputs(line1, ofp);
			double x,y,z;
			int typenum;
			int num_index2;
			char textend[1000];
			sscanf(line1,"%d %d %le %le %le %s", &num_index2, &typenum, &x,&y,&z,textend);
			real3 p1 = {x,y,z};
			pos.push_back(p1);
		}
		std::uniform_real_distribution<double> ranf(0,64.0);

		for (     ; atom_index <= natom2;atom_index ++) { 
			getline(&line2, &len2, srdfp);

			double x,y,z;
			int typenum;
			int num_index2;
			char textend[1000];
			sscanf(line2,"%d %d %le %le %le %[\001-\377]", &num_index2, &typenum, &x,&y,&z,textend);
			real3 p2 = {x,y,z};
			bool pos_ok=true;
			auto retry = 0 ;
			while (pos_ok) {
				pos_ok=false;
				for(auto i : pos) {
					double dx = i.x - p2.x; 
					double dy = i.y - p2.y; 
					double dz = i.z - p2.z; 
					while ( dx<-32) dx+=64;
					while ( dy<-32) dy+=64;
					while ( dz<-32) dz+=64;
					while ( dx>32) dx-=64;
					while ( dy>32) dy-=64;
					while ( dz>32) dz-=64;
					
					if ( abs(dx) < CUT && abs(dy) < CUT && abs(dz) <CUT && (dx*dx+dy*dy+dz*dz)< CUT*CUT  ) {
						p2. x = ranf(gen);
						p2. y = ranf(gen);
						p2. z = ranf(gen);
						retry ++;
						pos_ok =true;
						break;
					}
				}
			}
			x= p2.x;
			y=p2.y;
			z=p2.z;

			fprintf(ofp,    "%d %d %.16le %.16le %.16le %s",    num_index2, typenum, x,y,z,textend);
			fprintf(stdout, "%d %d %le %le %le %s retry %d\n", num_index2, typenum, x,y,z,textend,retry);


//			fputs(line2, ofp);
		}
		//blank
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);
		//Velocities
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);
		if (string(line1) != string(line2) ) {
			puts("Error not Atoms # hybrid");
			exit(1);
		}
		//blank
		getline(&line1, &len1, ldfp);
		getline(&line2, &len2, srdfp);
		fputs(line2, ofp);


		for (atom_index  = 1 ; atom_index <= natom1;atom_index ++) { 
			getline(&line1, &len1, ldfp);
			getline(&line2, &len2, srdfp);
			fputs(line1, ofp);
		}
		for (     ; atom_index <= natom2;atom_index ++) { 
			getline(&line2, &len2, srdfp);
			fputs(line2, ofp);
		}

	return 0;
}
