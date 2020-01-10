
#include "snapshot.h"

ENUM_DUMP get_dump( FILE* fp, Snapshot* snap) {
	const char s_timestep[] = "ITEM: TIMESTEP";
	const char s_n_atoms[] = "ITEM: NUMBER OF ATOMS";
	const char s_box_bounds[] = "ITEM: BOX BOUNDS pp pp pp";
	const char s_box_bounds_z[] = "ITEM: BOX BOUNDS pp pp ff";
	const char s_box_bounds_xy[] = "ITEM: BOX BOUNDS ff ff pp";
	const char s_box_bounds_xyz[] = "ITEM: BOX BOUNDS ff ff ff";
	const char s_atoms_rm[]    = "ITEM: ATOMS id type xu yu zu mux muy muz";
	const char s_atoms_rvm[]    = "ITEM: ATOMS id type xu yu zu vx vy vz mux muy muz";
	const char s_atoms_pos_r[]    = "ITEM: ATOMS id type xu yu zu";
	const int i_timestep = strlen(s_timestep);
	const int i_n_atoms = strlen(s_n_atoms);
	const int i_box_bounds = strlen(s_box_bounds);
	const int i_atoms_rm = strlen(s_atoms_rm);
	const int i_atoms_pos_r = strlen(s_atoms_pos_r);
	const int i_atoms_rvm = strlen(s_atoms_rvm);
	const char delimeter[] = " ";


	long timestep;
	int n_atoms;
	real xlow,xhigh,ylow,yhigh,zlow,zhigh;
	int id,type;
	real xu,yu,zu;
	real vx,vy,vz; 
	real mux,muy,muz;
	int i;
	struct atom* p_atom;

	read_lines(1,fp);
	if( strncmp(s_timestep,line,i_timestep) !=0) {
		error("not ITEM: TIMESTEP");
		return ERR1;
	}

	read_lines(1,fp);
	timestep = atol(line);
	fprintf(stderr,"atol(line) = %ld\n"
			"atoi(line) = %d\n"
			, timestep,atoi(line));

	read_lines(1,fp);
	if( strncmp(s_n_atoms,line,i_n_atoms) !=0) {
		return ERR1;
	}
	read_lines(1,fp);
	n_atoms = atoi(line);
	bool pbc[3];
	read_lines(1,fp);
	if( strncmp(s_box_bounds,line,i_box_bounds) ==0) {
		pbc[0] = true; pbc[1]=true; pbc[2]=true;
	}
	else if( strncmp(s_box_bounds_z,line,i_box_bounds) ==0) {
		pbc[0] = true; pbc[1]=true; pbc[2]=false;
	}
	else if( strncmp(s_box_bounds_xy,line,i_box_bounds) ==0) {
		pbc[0] = false; pbc[1]=false; pbc[2]=true;
	}
	else if( strncmp(s_box_bounds_xyz,line,i_box_bounds) ==0) {
		pbc[0] = false; pbc[1]=false; pbc[2]=false;
	}
	else {
		(error("not ITEM: BOX BOUNDS pp pp pp(ff)"));
		return ERR1;
	}

	read_lines(1,fp); 
	xlow = atof(strtok(line,delimeter));
	xhigh = atof(strtok(NULL,delimeter));
	read_lines(1,fp); 
	ylow = atof(strtok(line,delimeter));
	yhigh = atof(strtok(NULL,delimeter));
	read_lines(1,fp); 
	zlow = atof(strtok(line,delimeter));
	zhigh = atof(strtok(NULL,delimeter));

	
	struct box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, pbc[0], pbc[1], pbc[2]};
	snap->init(timestep,n_atoms);
	snap->setBox(box);


	error(s_timestep);
//	fprintf(stderr,"%ld\n", timestep);
	error(s_n_atoms);
//	fprintf(stderr,"%d\n", snap->n_atoms);
	error(s_box_bounds);
//	fprintf(stderr,"%f %f\n", snap->box.xlow,snap->box.xhigh);
//	fprintf(stderr,"%f %f\n", snap->box.ylow,snap->box.yhigh);
//	fprintf(stderr,"%f %f\n", snap->box.zlow,snap->box.zhigh);

	read_lines(1,fp);

	if( strncmp(s_atoms_rvm,line,i_atoms_rvm) ==0) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			id =atoi(strtok(line, delimeter));

			type = atoi(strtok(NULL,delimeter));
			xu   = atof(strtok(NULL,delimeter));
			yu   = atof(strtok(NULL,delimeter));
			zu   = atof(strtok(NULL,delimeter));
			vx   = atof(strtok(NULL,delimeter));
			vy   = atof(strtok(NULL,delimeter));
			vz   = atof(strtok(NULL,delimeter));
			mux   = atof(strtok(NULL,delimeter));
			muy   = atof(strtok(NULL,delimeter));
			muz   = atof(strtok(NULL,delimeter));

			p_atom = &(*snap).atoms[i];
			p_atom->init(id,type,xu,yu,zu,mux,muy,muz);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}else if( strncmp(s_atoms_rm,line,i_atoms_rm) ==0) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			id =atoi(strtok(line, delimeter));


			type = atoi(strtok(NULL,delimeter));
			xu   = atof(strtok(NULL,delimeter));
			yu   = atof(strtok(NULL,delimeter));
			zu   = atof(strtok(NULL,delimeter));
			mux   = atof(strtok(NULL,delimeter));
			muy   = atof(strtok(NULL,delimeter));
			muz   = atof(strtok(NULL,delimeter));

			p_atom = &(*snap).atoms[i];
			p_atom->init(id,type,xu,yu,zu,mux,muy,muz);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else if (strncmp(s_atoms_pos_r,line,i_atoms_pos_r) ==0 ) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			id =atoi(strtok(line, delimeter));
			type = atoi(strtok(NULL,delimeter));
			xu   = atof(strtok(NULL,delimeter));
			yu   = atof(strtok(NULL,delimeter));
			zu   = atof(strtok(NULL,delimeter));
			mux   = 0.;
			muy   = 0.;
			muz   = 0.;

			p_atom = &(*snap).atoms[i];
			p_atom->init(id,type,xu,yu,zu,mux,muy,muz);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}

	}
	else {
		error(line);
		error("not ITEM: ATOMS id type xu yu zu mux muy muz"
				"\n nor ITEM: ATOMS id type xu yu zu");
		return ERR1;
	}


	return SUCCESS;
}
Snapshot* read_dump( FILE* fp) {
	struct Snapshot* snap = new struct Snapshot;
	get_dump(fp, snap);
	return snap;
}
void* error( const char * string ) {
	fputs(string, stderr);
	fputs("\n", stderr);

	return NULL;
	//	exit(1);
}
void read_lines(int n,FILE* fp)  // from lammps reader_native.cpp
{
	char *eof;
	for (int i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
	if (eof == NULL) error("Unexpected end of dump file");
}


int compare_atom_type(const void * first, const void * second){
	if (   ((atom*) first) ->type < ((atom*) second )->type )
		return 1;
	else if (   ((atom*) first) ->type > ((atom*) second )->type )
		return -1;
	else return 0;
}

