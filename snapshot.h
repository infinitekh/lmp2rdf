#ifndef __SNAPSHOT__ 
#define __SNAPSHOT__ 

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include <map>
#include<list>
#include<cmath>
#include "kh_math_fourier.h"
#define MAXLINE 256


static char line[MAXLINE];
struct atom {
	int id,type;
	real x,y,z;
	real mux,muy,muz;
	real mu1;
	atom& operator = (const atom & p){
		id = p.id;
		type=p.type;
		x=p.x;y=p.y;z=p.z;
		mux=p.mux;muy=p.muy;muz=p.muz;
	}
	atom() {}
	bool operator< (const atom & rhs) {
		return type < rhs.type;
	}
	atom(int id, int type, 
			real x, real y, real z,
			real mux,real muy, real muz) {
		init( id,  type, x,  y,  z, mux, muy,  muz);
	}
	void init(int id, int type, 
			real x, real y, real z,
			real mux,real muy, real muz) {
		this->id=id; this->type=type;
		this->x=x;this->y=y;this->z=z;
		mu1 = sqrt(mux*mux+muy*muy+muz*muz);
		if (mu1>0.) { 
			this->mux=mux/mu1;this->muy=muy/mu1;this->muz=muz/mu1;
		}
	}
};
struct box3 { 
	real xlow,xhigh;
	real ylow,yhigh;
	real zlow,zhigh;
	bool pbc[3];
};

typedef struct Snapshot {
	long timestep;
	int    n_atoms;
	int    n_types;
	long  bytes;
	struct atom* atoms;
	struct box3 box;
/* 	Snapshot(){
 * 	}
 */
	Snapshot(long _ts,int _n_atoms){
		init(_ts,_n_atoms);
	}
/* 	Snapshot(const  Snaptshot &T) {
 * 		timestep    = T.timestep;
 * 		n_atoms     = T.n_atoms;
 * 		box         = T.box;
 * 		atoms       = T.atoms;
 * 		n_types     = T.n_types;
 * 	}
 */
	void init(long _ts,int _n_atoms){
		this->timestep = _ts;
		n_atoms = _n_atoms;
		atoms = new struct atom[n_atoms];
		bytes = sizeof(Snapshot) + sizeof(atom) * n_atoms;
	}
	void setBox(real xl,real xh, real yl,real yh, real zl,real zh) {
		box.xlow=xl;box.xhigh=xh;
		box.ylow=yl;box.yhigh=yh;
		box.zlow=zl;box.zhigh=zh;
	}
	void setBox(box3 otherbox){
		memcpy(&box,&otherbox,sizeof(otherbox));
	}
	int get_n_ptls (int type){
		int _n_types = calc_n_type () ;
		if( type<_n_types) return 0;
		
		int n_ptls=0;
		for( int i =0; i<n_atoms; i++ ) {
			if (type == atoms[i].type) 
				n_ptls ++ ;
		}
		return n_ptls;
	}

	int calc_n_type (){
		if(flag_calc_n_type) return n_types;

		std::list<int> unique_type_list;
		for( int i =0; i<n_atoms; i++ ) {
			unique_type_list.push_back(atoms[i].type);
		}
		unique_type_list.sort();
		unique_type_list.unique();

		n_types = unique_type_list.size();
		flag_calc_n_type = true;
		return n_types;
	}
	
	~Snapshot() {

		delete [] atoms;
	}
	protected:
	bool flag_calc_n_type=false;
} Snapshot;


int compare_atom_type(const void *, const void *);
struct Snapshot* read_dump(FILE*);
typedef enum { 
	SUCCESS, ERR1,ERR2,ERR3,ERR4
} ENUM_DUMP;
ENUM_DUMP get_dump(FILE*,Snapshot*);
void read_lines(int n,FILE*);
void* error( const char *);


typedef std::list<Snapshot*> l_snapshot;

#endif
