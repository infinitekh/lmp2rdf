#include<cstdio>
#include<cstring>
#include<cstdlib>
#include <map>
#include<list>
#include<cmath>
#include"snapshot.h"
#include"makeRDF.h"


using namespace std;
typedef std::list<Snapshot*> t_snaplist;
typedef Snapshot* p_snapshot;


int main(int argc, char** argv) {
	char filename[100];
	int typei, typej;
	FILE* fp;
	if(argc <2) {
		perror("#run inputfilename");
		return 1;
	}
	if ( argc == 3 && isdigit(argv[2][0]) ) {
		 typei = atoi(argv[2]);
	}
	else if (argc ==4 && isdigit(argv[3][0] ))  {
		 typei = atoi(argv[2]);
		 typej = atoi(argv[3]);
	}
	


	strcpy( filename,argv[1]);
	fp = fopen( filename ,"r");
	t_snaplist snaplist;
	p_snapshot snap;
	
	while(1)
	{
		snap =	read_dump(fp);

		if (snap==NULL) {
			puts("read dump end");
			break;
		}
		snaplist.push_back(snap);
	}

	if (snaplist.size() <5){
		fprintf(stderr, "The # of snap is too small(%d<5)!!\n", snaplist.size());
		return 23;
	}
	
	makeRDF makerRdf(snaplist);
	makerRdf.flag_anisotropy =0;
	
	puts("calc RDF begin");
//	makerRdf.calcRDF();
	puts("calc RDF intertype");
	makerRdf.calcRDF_inter_type(1,1, makerRdf.ANISO);
	makerRdf.calcRDF_inter_type(1,2);
	puts("calc RDF end");
	makerRdf.calcSSF_from_g ();
/* 	makerRdf.calcSSF();
 * 	makerRdf.calcISF();
 */

}


