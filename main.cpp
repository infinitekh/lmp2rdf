#include<cstdio>
#include<cstring>
#include<cstdlib>
#include <map>
#include<cmath>
#include <cstdlib>
#include"snapshot.h"
#include"makeRDF.h"
#include "kh_math_fourier.h"
#include <cerrno>
#include <getopt.h>


using namespace std;
typedef std::vector<Snapshot*> t_snaplist;
typedef Snapshot* p_snapshot;
/* 
			{"help",no_argument, 0, 'h'},
			{"verbose", no_argument, 0, &verbose_flag, 1},
			{"aniso", no_argument,0, 'I'},
			{"iso", no_argument,0, 'i'},
			{"maxbin", required_arguement,0, 'b'},
			{"rcut", required_arguement,0, 'c'},
			{0,0,0,0},
 */
void PrintHelp ( char *pName)
{
	printf ("Usage: %s  \n"
			" if you want to use stdin, you should used /dev/stdin  \n"
			"--aniso (I) :       get non isotropic terms h112 etc. \n"
			"--iso   (i) : 			 get only isotropic terms g000 h000 \n"
			"--maxbin (b) <maxbin(int)>  set Maxbin \n"
			"--rcut (c) <rcut(double)>   set rcut   \n"
			"--ext (x) .com .bat .1 .2  \n"
			"default : iso maxbin(200) rcut(30) \n" 
			, pName);
	exit(0);
}

int verbose_flag=0;
int main(int argc, char** argv) {
	char filename[200];
	int typei=1, typej=1;
	FILE* input;
	bool aniso=false;
	int n =1;
	int maxbin=200;
	double r_cut = 30;
	char ext[30]= "";
	
	int opt;
	int option_index=0;

	while (1) {
		static struct option long_options[] = 
		{ 
			{"help",no_argument, 0, 'h'},
			{"help",no_argument, 0, 'h'},
			{"verbose", no_argument,  &verbose_flag, 1},
			{"aniso", no_argument,0, 'I'},
			{"iso", no_argument,0, 'i'},
			{"maxbin", required_argument,0, 'b'},
			{"ext", required_argument,0, 'x'},
			{"rcut", required_argument,0, 'c'},
			{0,0,0,0},
		};
		opt = getopt_long (argc,argv, "hIib:c:x:",
				long_options, &option_index);
		if ( opt == -1) break;

		switch(opt) {
			case 0: break;
			case 'i': 
							aniso = false; 
							puts("isotropic on");
							break;
			case 'I': 
							aniso = true; 
							puts("anisotropic on");
							break;
			case 'x':  
							strcpy(ext, optarg);
							printf("maxbin %d\n",maxbin);
							break;
			case 'b': 
							maxbin = atoi(optarg);
							printf("maxbin %d\n",maxbin);
							break;
			case 'c': 
							r_cut = atof(optarg);
							printf("cut-off length %.2g\n", r_cut);
							break;
			case 'h':
			case '?': 
							PrintHelp(argv[0]);
							break;
			default: 
							printf("pass through -%c \n", opt);
							break;
		}
	}



	printf("typei %d typej %d  argc %d\n", typei, typej, argc); 
	for (int i =0;  i!= argc; ++i) {
		printf("%s ", argv[i]);
	}
	puts("");


	//	for (int i = optind; i!= argc; ++i){
	//		strcpy(filename, argv
	//	}

	if (optind == argc ) {
		PrintHelp(argv[0]);
	}
	t_snaplist snaplist;
	for( int opt_num = optind;   opt_num < argc; opt_num++)  {
		strcpy( filename, argv[opt_num]);
		printf("filename is %s \n", filename);

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

		//	input = fopen( filename ,"r");

		while(1)
		{
			p_snapshot snap =	read_dump(input);

			if (snap==NULL) {
				puts("read dump end");
				break;
			}
			snaplist.push_back(snap);
		}
	}

	if (snaplist.size() <5){
		fprintf(stderr, "The # of snap is too small(%ld<5)!!\n", 
				snaplist.size());
		return 23;
	}

	makeRDF makerRdf(snaplist);
	makerRdf.flag_anisotropy =0;
	makerRdf.maxbin = maxbin;
	makerRdf.r_cut = r_cut;
	strcpy(	makerRdf.ext , ext);

	makerRdf.StartMainProcess ();
	puts("calc RDF begin");
	//	makerRdf.calcRDF();
	puts("calc RDF intertype");
	makerRdf.calcRDF_inter_type(1,1, makerRdf.ANISO);
//	makerRdf.calcRDF_inter_type(1,2, makerRdf.ISO);
	puts("calc RDF end");
	makerRdf.calcSSF_from_g ();
/* 	makerRdf.calcSSF();
 * 	makerRdf.calcISF();
 */
	exit(0);

}

