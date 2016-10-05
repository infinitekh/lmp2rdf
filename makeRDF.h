#ifndef  MAKERDF_H_INC
#define  MAKERDF_H_INC

#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include "snapshot.h"
using namespace std;
struct s_xyz {
	real x, y ,z;
};
typedef struct {
	real **acfST, *orgST;
	int count;
} TBuf;

#define VDOT3(a,b) ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] )
#define Sqr(rrrr)   ( (rrrr)*(rrrr))


class makeRDF {
public:
	makeRDF(std::list<Snapshot*> &_sl) ;
	~makeRDF() { 
		
		delete [] tBuf;
		delete [] avAcfST[0];
		delete [] avAcfST;
		delete [] valST;
	};
	typedef enum { ISO, ANISO } T_RDF;
	bool flag_anisotropy;
	bool flag_SSF_from_g;

	void calcRDF();
	void calcRDF_inter_type (int i, int j, T_RDF rdftype );
	void calcP1z(int type);
	void calcP1s(int type);
	void calcSSF();
	void calcSSF_from_g();
	void calcISF();
	void EvalSpacetimeCorr();
	void InitSpacetimeCorr();

	bool* periodicity;
	// intput value
	int nFunCorr ;
	int limitCorrAv;
	int nBuffCorr;   // nValCorr = number x nBuffCorr
	int nValCorr; // # of average sets
	int step,stepCorr;
	int countCorrAv;
	static char* filename_template;
	char filename0[100];
	char filename1[100];
	char filename2[100];
	FILE* fp[3];
	// Avf
	TBuf *tBuf;
	real **avAcfST, *valST;
	int  n_valST;
	real deltaT ;
	// Avf function
	

	void PrintSpacetimeCorr ();
	void ZeroSpacetimeCorr ();

	map<real, real>& get_g000();	
	map<real, real>& get_h000();	
	map<real, real>& get_h110();	
	map<real, real>& get_h112();	
	map<real, real>& get_h220();	
	map<real, real>& get_S();	
	map<real, real>& get_F();	
	map<real, real>& get_M();	
	map<real, real>& get_M_L();	
	map<real, real>& get_M_T();	

	real box_x,box_y,box_z;
	box3* box;
	int maxSnap ;
	real r_cut,q_cut;
	int maxbin ;
	int maxbinz ;
	int maxbins ;
	int maxbinq ;
	real dq;
	int maxAtom;
private:
	vector<real> vP1s;
	vector<real> vPmus1s;
	vector<real> vP1z;
	vector<real> vP1zi;
	vector<real> vP1zj;
	vector<real> vPmuz1z;
	vector<real> vz,vs;
	map<real, real> g000;	
	map<real, real> h000;	
	map<real, real> h110;	
	map<real, real> h112;	
	map<real, real> h220;	
	map<real, real> S;	
/* 	map<real, real> S_qr;	
 * 	map<real, real> S_qi;	
 * 	map<real, real> qlist;	
 */

	map<real, real> F;	
	map<real, real> M;	
	map<real, real> M_L;	
	map<real, real> M_T;	
	void calcRDF_isotropy ();
	void calcRDF_cotype (int i);
	void calcRDF_anisotropy ();
	void AccumSpacetimeCorr ();

	void calcIntegrated_h000();
	real var_k, var_r;
	std::list<Snapshot*> snaplist;
	atom *atoms,*ppi,*ppj;
};

#endif   /* ----- #ifndef MAKERDF_H_INC  ----- */
