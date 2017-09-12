#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include <map>
#include<list>
#include<cmath>
#include <cstdlib>
#include"snapshot.h"
#include"makeRDF.h"
#include <cerrno>
#include <mpi.h>


using namespace std;
typedef std::vector<Snapshot*> t_snaplist;
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
//		else if (! strcmp (argv[n], "-i")) {typei = atoi (argv[n+1])  ; n++;argc--;}
//		else if (! strcmp (argv[n], "-j")) {typej = atoi (argv[n+1]); ; n++;argc--;}
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
	MPI_Init(NULL, NULL);

	

	int world_size;
	int world_rank;
	int ret;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int array_size[world_size];
	int array_bytes[world_size];
	int local_array_size;
	int array_min[world_size];
	int array_min_bytes[world_size];
	int array_max[world_size];
	char * global_data;
	char * local_data;
	int buffer_size;
	int data_bytes;


	t_snaplist snaplist, local_snaplist;
	int snapshot_bytes;
	int n_atoms;
	int n_atoms_bytes;
	p_snapshot snap;
	Snapshot* pp_snap[world_size];
	if(world_rank == 0) {

		input = fopen( filename ,"r");

		while(1)
		{
			snap =	read_dump(input);

			if (ret!=SUCCESS) {
				puts("read dump end");
				break;
			}

			snaplist.push_back(snap);
		}

		int snapsize= snaplist.size();

		if (snapsize <5){
			fprintf(stderr, "The # of snap is too small(%ld<5)!!\n", snaplist.size());
			return 23;
		}
		fclose(input);

		snapshot_bytes = snaplist.at(0)->bytes;
		n_atoms = snaplist.at(0)->n_atoms;
		n_atoms_bytes = sizeof(atom)*n_atoms;
		data_bytes = snapshot_bytes*snapsize;
		buffer_size = 0;

		for (int i=0; i<world_size; i++) {
			array_min[i]= i/world_size*snapsize;
			array_max[i]= (i+1)/world_size*snapsize;
			array_size[i] = array_max-array_min;
			buffer_size = max(array_size[i], buffer_size);
			array_bytes[i] = array_size[i]*snapshot_bytes;
			array_min_bytes[i] = array_min[i]*snapshot_bytes;
		}

		buffer_size = buffer_size * snapshot_bytes;
		global_data = (char *)malloc( data_bytes );
		
		int write_bytes =0;
		for (int i = 0; i<snapsize ; i++)
		{
			memcpy(global_data+write_bytes,snaplist.at(i), sizeof(Snapshot));
			write_bytes += sizeof(Snapshot);
			memcpy(global_data+write_bytes,snaplist.at(i)->atoms, n_atoms_bytes);
			write_bytes += n_atoms_bytes;

		}
		assert(write_bytes == data_bytes);
	}
	MPI_Bcast( &buffer_size , 1, MPI_INT, 0,MPI_COMM_WORLD);
	MPI_Bcast( &snapshot_bytes , 1, MPI_INT, 0,MPI_COMM_WORLD);

	local_data = (char*) malloc( buffer_size);
	
	MPI_Scatter(array_size, 1, MPI_INTEGER, 
			&local_array_size, 1, MPI_INTEGER, 0,
			MPI_COMM_WORLD);
	MPI_Scatter(array_size, 1, MPI_INTEGER, 
			&local_array_size, 1, MPI_INTEGER, 0,
			MPI_COMM_WORLD);
	MPI_Scatterv(global_data, array_bytes , array_min_bytes ,MPI_CHAR, 
			local_data, buffer_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	local_snaplist.resize(local_array_size);
	for (int j=0; j<local_array_size; j++) {
		local_snaplist.at(j) = (Snapshot*) (local_data + (j*snapshot_bytes) );
		local_snaplist.at(j)->atoms = (atom*) (local_data + (j*snapshot_bytes)+ sizeof(Snapshot) );
	}


	makeRDF makerRdf(local_snaplist);
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

	MPI_Finalize();	
}

void PrintHelp ( char *pName)
{
/* 	printf ("Usage: %s [-aniso  "
 * 			" input-file \n"
 * 			" if you want to use stdin, you should used -  \n"
 * 			, pName);
 */
	printf ("Usage: %s  \n"
			" if you want to use stdin, you should used /dev/stdin  \n"
			, pName);
/* 	printf ("Usage: %s [-aniso  -i itype -j jtype "
 * 			" input-file \n"
 * 			" if you want to use stdin, you should used -  \n"
 * 			, pName);
 */
//	exit(0);
}
