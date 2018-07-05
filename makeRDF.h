#ifndef  MAKERDF_H_INC
#define  MAKERDF_H_INC

#include <map>
//#include <list>
#include <vector>
#include <algorithm>
#include "snapshot.h"
#include "kh_math_fourier.h"
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
	makeRDF(std::vector<Snapshot*> &_sl) ;
	~makeRDF() { 
		freeMem();
		
/* 		delete [] tBuf;
 * 		delete [] avAcfST[0];
 * 		delete [] avAcfST;
 * 		delete [] valST;
 */
	};
	typedef enum { ISO, ANISO } T_RDF;
	bool flag_anisotropy= false;
	bool flag_SSF_from_g= false;

	void calcRDF();
	void calcRDF_inter_type (int i, int j, T_RDF rdftype );
	void printRDF (int i, int j, T_RDF rdftype );

	void calcP1z(int type);
	void calcP1s(int type);
	void calcSSF(T_RDF rdftype = ISO);
	void calcSSF_from_g(T_RDF rdftype = ISO);
	void calcISF();
	void EvalSpacetimeCorr();
	void InitSpacetimeCorr();
	void StartMainProcess();
	void StartPreProcess();

	bool* periodicity;
	// intput value
	int nFunCorr =4 ;
	int limitCorrAv =200;
	int nBuffCorr = 10;   // nValCorr = number x nBuffCorr
	int nValCorr =500; // # of average sets
	int step=1,stepCorr=1;
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
	real deltaT = 0.0001 ;
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
	real hbox_x,hbox_y,hbox_z;
	box3* box;
	int maxSnap ;
	real r_cut,q_cut;
	int maxbin =200;
	int maxbinz = 5;
	int maxbins = 5;
	int maxbinq = 100;
	real dq;
	int maxAtom;
	void set_maxbin(int _maxbin){
		maxbin = _maxbin;
		var_r = r_cut/maxbin;
	}

	size_t  rbin_t;
/* 	long* hist000;
 * 	double* hist110;
 * 	double* hist112;
 * 	double* hist220;
 * 
 * 	long* histcyl000;
 * 	double* histcyl110;
 * 	double* histcyl112;
 * 	double* histcyl220;
 */

	double* ca_radius;
	double* ca_g000;
	double* ca_h000;
	double* ca_h110;
	double* ca_h112;
	double* ca_h220;
	double* ca_gcyl000;
	double* ca_hcyl110;
	double* ca_hcyl112;
	double* ca_hcyl220;


	double* ca_c_q;
	double* ca_q_radius;
	double* ca_S_q;
	double* ca_h000_q;
	double* ca_h110_q;
	double* ca_h112_q;
	double* ca_h220_q;

	double* ca_c_r;
	double* ca_h_mod;
private:
	vector<real> vP1s;
	vector<real> vP1zj;
	vector<real> vz;
	vector<real> vs;
	vector<real> vPmus1s;
	vector<real> vPmuz1z;
	vector<real> vP1z;
	vector<real> vP1zi;

/* 	map<real, real> S_qr;	
 * 	map<real, real> S_qi;	
 * 	map<real, real> qlist;	
 */

	void AllocMem();
	void freeMem();
	void calcRDF_anisotropy_snap (Snapshot* snap);

	void calcRDF_isotropy ();
	void calcRDF_isotropy (int i, int j);
	void calcRDF_anisotropy ();
	void calcRDF_anisotropy (int i, int j);
	void AccumSpacetimeCorr ();

	void calcIntegrated_h000();
	real var_k, var_r;
	const std::vector<Snapshot*> snaplist;
	atom *first_atoms;
};

#endif   /* ----- #ifndef MAKERDF_H_INC  ----- */
