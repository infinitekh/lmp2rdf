#include<cstdio>
#include<cstring>
#include<cstdlib>
#include <map>
#include<list>
#include<cmath>
#include <cstdlib>
#include"snapshot.h"
#include"makeRDF.h"
#include <cerrno>


using namespace std;
typedef std::list<Snapshot*> t_snaplist;
typedef Snapshot* p_snapshot;

void PrintHelp ( char *pName);

int main(int argc, char** argv) {
	char* filename;
	int typei=1, typej=1;
	FILE* input;
	bool aniso=false;
	int n =1;
	if(-- argc <1 || ! strcmp (argv[1], "-h")) PrintHelp(argv[0]);
	
	while (-- argc >= 0) {
		if (! strcmp (argv[n], "-aniso")) aniso =true;
		else if (! strcmp (argv[n], "-i")) {typei = atoi (argv[n+1])  ; n++;argc--;}
		else if (! strcmp (argv[n], "-j")) {typej = atoi (argv[n+1]); ; n++;argc--;}
		else {
			filename = argv[n];
			break;
		}
		++ n;
	}

	
	if (argc>0) PrintHelp (filename);

	printf("typei %d typej %d  argc %d\n", typei, typej, argc); 


	if(!strcmp(filename,"-")) {
		input = stdin;
	} else {
		input = fopen(filename,"r");
		if (NULL == input) {
			fprintf(stderr, "Unable to open '%s': %s\n",
					filename, strerror(errno));
			exit(EXIT_FAILURE);
		}   
	}

	input = fopen( filename ,"r");
	t_snaplist snaplist;
	p_snapshot snap;

	while(1)
	{
		snap =	read_dump(input);

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
	makerRdf.calcRDF_inter_type(1,2, makerRdf.ISO);
	puts("calc RDF end");
	makerRdf.calcSSF_from_g ();
/* 	makerRdf.calcSSF();
 * 	makerRdf.calcISF();
 */

}

void PrintHelp ( char *pName)
{
	printf ("Usage: %s [-aniso  -i itype -j jtype "
			" input-file \n"
			" if you want to use stdin, you should used -  \n"
			, pName);
//	exit(0);
}
