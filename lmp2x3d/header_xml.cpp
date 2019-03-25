/*
 * =====================================================================================
 *
 *       Filename:  header_xml.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2015년 07월 13일 19시 00분 37초
 *       Revision:  1.1
 *       Compiler:  g++
 *
 *         Author:  Dr.ing. KIM Hyeok (Authorref??), ekh0324@gmail.com
 *        Company:  Konkuk University
 *
 * =====================================================================================
 */

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include <map>
#include<list>
#include<vector>
#include<cmath>
#include"snapshot.h"
#define NCell  0

typedef struct { real x,y,z; } real3d;
using namespace std;
real shift_periodic_position(const real x, real min,real max)
{
	real value=x;
	real width = max - min;
	real ret = 0 ;
	if (width <0 ) return x;
	if (value >= max) {
		while (value >=max)  {
			value-=width;  ret -=width;
		}
	}
	else {
		while (value <min)  {
			value+=width; ret += width;
		}
	}
	return ret;
}
real len3(real* s3)
{
	real x = s3[0],y =s3[1],z=s3[2];
	return sqrt(x*x+y*y+z*z);
}
void angle2rot4(real* s4,real* s3)
{
	real x = s3[0],y =s3[1],z=s3[2];
	real len_s = sqrt(x*x+y*y+z*z);
	s4[0] = z/len_s;
	s4[1] = 0;
	s4[2] = -x/len_s;
	s4[3] = acos(y/len_s);
}
#define BIG_RADIUS (4.3/2.)
#define BIG_HEIGHT 4.3
#define FRAMERATE  30
#define SMALL_RADIUS 0.1
double big_radius = BIG_RADIUS;
class PrintX3D {
	public:
		PrintX3D (std::list<Snapshot*> *_sl);


		~PrintX3D(); 
		string pathKey;
		string* colPath;
		string* colOriPath;

		real box_x,box_y,box_z;
		box3* box;
		std::list<Snapshot* > snaplist;
		
		real3d * firstsnap_shift;
		atom *atoms;
		int maxSnap ;
		int maxAtom;
		bool pbc[3]={true,true,false};
		const string first_text = R"long(<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE X3D PUBLIC "ISO//Web3D//DTD X3D 3.2//EN" "../../schemas/x3d-3.2.dtd">
<X3D profile='Interchange' version='3.1' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation='http://www.web3d.org/specifications/x3d-3.1.xsd'>
<head>
<meta content='magneticColloid.x3d' name='title'/>
<meta content='Example for Group node' name='description'/>
<meta content='1.1' name='version'/>
</head>
<!-- create an almost white background for printing -->
<Scene DEF='scene'>
<Viewpoint
  position="50.829875946044922 17.521451950073242 125.79796600341797"
	  orientation="-0.030014935880899429 0.99576050043106079 0.086928412318229675 0.30863013863563538"
/>
<Background skyColor='1 1 1' />
<TimeSensor DEF='Clock' cycleInterval='%f' loop='true'/>
	)long";
    const string group_begin=  R"(<Group DEF="%s">)";
    const string group_end=  R"(</Group>)";
		const string format_pack = R"(<Transform translation="%.3f %.3f %.3f">
<Group USE="%s"/>
</Transform>
			)";

		const string last_text  = R"(</Scene> 
			</X3D> )";

		const string format_path = R"( <PositionInterpolator DEF='ColPath%d' key='%s' keyValue='%s'/> )";
		const string format_oriPath= R"( <OrientationInterpolator DEF='ColOriPath%d' key='%s' keyValue='%s'/> )";



		const string bounding_box_z = R"(<Shape>
<Appearance>
<Material diffuseColor='0 0 0'/>
</Appearance> 
<IndexedFaceSet coordIndex='0 1 2 3 0 -1 4 5 6 7 4 -1 0 3 2 1 0 -1 4 7 6 5 4 -1 '>
<Coordinate point='%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f'/>
</IndexedFaceSet> </Shape> )"; 
const string format_space = R"(<Shape>
<Appearance>
<Material emissiveColor='0 0 0'/>
</Appearance> 
<IndexedLineSet coordIndex='0 1 2 3 0 -1 4 5 6 7 4 -1 0 4 -1 1 5 -1 2 6 -1 3 7 -1'>
<Coordinate point='%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f 
%.2f %.2f %.2f'/>
</IndexedLineSet> </Shape> 
			)"; 

			const string format_event= R"(<ROUTE fromNode='Clock' fromField='fraction_changed' toNode='ColPath%d'
toField='set_fraction'/>
<ROUTE fromNode='ColPath%d' fromField='value_changed' toNode='Col%d'
toField='set_translation'/> 
<ROUTE fromNode='Clock' fromField='fraction_changed' toNode='ColOriPath%d'
toField='set_fraction'/>
<ROUTE fromNode='ColOriPath%d' fromField='value_changed' toNode='Col%d'
toField='set_rotation'/> 
			)";

/*		const string format_dipole  =R"( <Transform DEF='Col%d' translation="%.3f %.3f %.3f" 
rotation="%.3f %.3f %.3f %.3f">
<Shape>
<Appearance>
<Material/>
<PixelTexture image='1 4 4 0xFF0000ff 0xFF0000ff 0x0000FFFF 0x0000FFFF' repeatS='FALSE' repeatT='FALSE' />
</Appearance>
<Sphere radius="%.3f"/>
</Shape>
<Transform translation='0 0 0.5'>
<Shape>
<Text string='"%d"' > <FontStyle justify='"MIDDLE" "MIDDLE"' size='%f'/> </Text>
<Appearance>
<Material diffuseColor='0 0 0'   />
</Appearance>
</Shape>
</Transform>
</Transform>
			)";*/
		const string format_dipole  =R"(<Transform DEF='Col%d' translation="%.3f %.3f %.3f" 
rotation="%.3f %.3f %.3f %.3f" scale="0.8 1.2 0.8">
<Shape>
<Appearance>
<Material/>
<PixelTexture image='1 2 4 0xFF0000ff  0x0000FFFF' repeatS='FALSE' repeatT='FALSE' />
</Appearance>
<Sphere radius="%.3f"/>
</Shape>
</Transform>
			)";

		const string format_big   =   R"( <Transform DEF='Big%d' translation="%.3f %.3f %.3f" >
<Shape> <Appearance>
<Material diffuseColor='0 1 0' />
</Appearance>
<Sphere radius="%.3f"/>
</Shape> </Transform>
			)";
		const string format_bigPath= R"( <PositionInterpolator DEF='BigPath%d' key='%s' keyValue='%s'/>)";
		const string format_bigEvent= R"( <ROUTE fromNode='Clock' fromField='fraction_changed' toNode='BigPath%d'
toField='set_fraction'/>
<ROUTE fromNode='BigPath%d' fromField='value_changed' toNode='Big%d'
toField='set_translation'/> 
			)";

		const string format_small   =   R"( <Transform DEF='Small%d' translation="%.3f %.3f %.3f" >
<Shape> <Appearance>
<Material diffuseColor='0 1 0' />
</Appearance>
<Sphere radius="%.3f"/>
</Shape> </Transform>
			)";
		const string format_smallPath= R"( <PositionInterpolator DEF='SmallPath%d' key='%s' keyValue='%s'/>)";
		const string format_smallEvent= R"( <ROUTE fromNode='Clock' fromField='fraction_changed' toNode='SmallPath%d'
toField='set_fraction'/>
<ROUTE fromNode='SmallPath%d' fromField='value_changed' toNode='Small%d'
toField='set_translation'/> 
			)";

		void draw_bounding_box_z(FILE* fp) {
			//coordIndex='0 1 2 3 0 -1 4 5 6 7 4 -1 0 4 -1 1 5 -1 2 6 -1 3 7 -1'
			fprintf(fp,bounding_box_z.c_str(), 
					box->xlow, box->ylow, box->zlow, //0
					box->xhigh,box->ylow, box->zlow, //1
					box->xhigh,box->yhigh,box->zlow, //2
					box->xlow, box->yhigh,box->zlow, //3
					box->xlow, box->ylow, box->zhigh,//4
					box->xhigh,box->ylow, box->zhigh,//5
					box->xhigh,box->yhigh,box->zhigh,//7
					box->xlow, box->yhigh,box->zhigh //6
					);
		}
		void draw_space(FILE* fp) {
			//coordIndex='0 1 2 3 0 -1 4 5 6 7 4 -1 0 4 -1 1 5 -1 2 6 -1 3 7 -1'
			fprintf(fp,format_space.c_str(), 
					box->xlow, box->ylow, box->zlow, //0
					box->xlow, box->ylow, box->zhigh,//1
					box->xhigh,box->ylow, box->zhigh,//2
					box->xhigh,box->ylow, box->zlow, //3
					box->xlow,box->yhigh,box->zlow, //4
					box->xlow,box->yhigh,box->zhigh,//5
					box->xhigh,box->yhigh,box->zhigh,//6
					box->xhigh,box->yhigh,box->zlow );//7
		}
		void draw_dipole(FILE* fp,int id,real* v3, real* s4)
		{ 
		//	fprintf(fp,format_dipole.c_str(),id, v3[0],v3[1],v3[2],s4[0],s4[1],s4[2],s4[3] , BIG_RADIUS,id,BIG_RADIUS);
		fprintf(fp,format_dipole.c_str(),id, v3[0],v3[1],v3[2],s4[0],s4[1],s4[2],s4[3]  ,big_radius );

		}

		void draw_path(FILE* fp, int id) {
			fprintf (fp,format_path.c_str(),
					id, pathKey.c_str(), colPath[id].c_str());
			fprintf (fp,format_oriPath.c_str(),
					id, pathKey.c_str(), colOriPath[id].c_str());
		}
		void draw_event(FILE* fp,int id){
			fprintf (fp,format_event.c_str(), id,id,id,id,id,id);
		}

		void draw_big(FILE* fp,
				int id,real rx,real ry,real rz){
			fprintf (fp,format_big.c_str(),
					id,rx,ry,rz,big_radius);
		}
		void draw_bigPath(FILE* fp,int id){
			fprintf (fp,format_bigPath.c_str(),
					id,pathKey.c_str(), colPath[id].c_str());
		}
		void draw_bigEvent(FILE* fp,int id){
			fprintf (fp,format_bigEvent.c_str(), id,id,id);
		}

		void draw_small(FILE* fp,
				int id,real rx,real ry,real rz){
			fprintf (fp,format_small.c_str(),
					id,rx,ry,rz,SMALL_RADIUS);
		}
		void draw_smallPath(FILE* fp,int id){
			fprintf (fp,format_smallPath.c_str(),
					id,pathKey.c_str(), colPath[id].c_str());
		}
		void draw_smallEvent(FILE* fp,int id){
			fprintf (fp,format_smallEvent.c_str(), id,id,id);
		}
		void draw_images(FILE* fp) {
			fprintf(fp, group_begin.c_str(), "IMAGE_XYZ");
hello:
			{
				int imin,imax,jmin,jmax,kmin,kmax;
				if(pbc[0]) {
					imin = -NCell; imax=NCell;
				}
				else {
					imin =0 ; imax =0;
				}
				if(pbc[1]) {
					jmin = -NCell; jmax=NCell;
				}
				else {
					jmin =0 ; jmax =0;
				}
				if(pbc[2]) {
					kmin = -NCell; kmax=NCell;
				}
				else {
					kmin =0 ; kmax =0;
				}

				for (int i=imin; i<=imax; i++) {
					for(int j=jmin; j<=jmax; j++) {
						for(int k=kmin; k<=kmax; k++) {
							if( (i==0)&&(j==0) &&(k==0)){
								fprintf(stderr, "%d %d %d  i j k\n",i,j,k);
								continue;
							}
							fprintf(fp, format_pack.c_str()
									,i*box_x, j*box_y,k*box_z
									, "Reality" );
						}
					}
				}
			}
			fputs(group_end.c_str(), fp);
		}  // end void draw_images

		void drawX3D(FILE* fp) {
			real x,y,z,s3[3],s4[4], x3[3];
			int type, id; bool bigFlag =false;
			char temp[300];
			fprintf(fp,first_text.c_str(), 1.*maxSnap/FRAMERATE);
			fprintf(fp, group_begin.c_str(), "Reality");
			draw_space(fp);
			
			if(!(pbc[2]))
				draw_bounding_box_z(fp);
			// auto snap -> auto s
			auto snap = snaplist.begin();
			for (int n=0; n<maxAtom; n++){
					atoms = &((*snap)->atoms[n]);
					firstsnap_shift[n].x=shift_periodic_position(atoms->x,0,box_x);
					firstsnap_shift[n].y=shift_periodic_position(atoms->y,0,box_y);
					firstsnap_shift[n].z=shift_periodic_position(atoms->z,0,box_z);
			}

			for (auto s=snaplist.begin(); s!=snaplist.end();s++){
				for (int n=0; n<maxAtom; n++){
					atoms = &((*s)->atoms[n]);
					real3d * shift = &(firstsnap_shift[n]);
					id=atoms->id;
					x=atoms->x+ shift->x;
					y=atoms->y+ shift->y;
					z=atoms->z+ shift->z;
					type = atoms->type;
					if ( type==1) {
						s3[0]=atoms->mux;
						s3[1]=atoms->muy;
						s3[2]=atoms->muz;
						angle2rot4(s4,s3);
						sprintf(temp, "%f %f %f %f\n", s4[0],s4[1],s4[2],s4[3]); 
						colOriPath[id].append(temp);
					}

					if( s == snaplist.begin()) {
						if ( type ==2 )
							draw_small(fp,id,x,y,z);
						if (type ==1  ) {
							x3[0]=x;x3[1]=y;x3[2]=z;
							if(bigFlag||len3(s3)<0.1) {
								draw_big(fp,id,x,y,z);
								bigFlag=true;
							}
							else {
								draw_dipole(fp,id, x3,s4); 
							}
						}
					}
					sprintf(temp, "%f %f %f\n", x,y,z); 
					colPath[id].append(temp);
				}

			}
			for (int n=0; n<maxAtom; n++){
				atoms = &((*snap)->atoms[n]);
				id=atoms->id;
				x=atoms->x;
				y=atoms->y;
				z=atoms->z;
			}
			fputs(group_end.c_str(), fp);
			// Path printing
			for (int i =1; i<=maxSnap; i++) {
				sprintf(temp, "%f ", 1.*i/(1.*maxSnap));
				pathKey.append(temp); 
				if (i%4==0 ) pathKey.append("\n");
			}

			for (int n=0; n<maxAtom; n++) {
				atoms = &((*snap)->atoms[n]);
				type = atoms->type;
				id=atoms->id;
				if (type == 1) {
					if (bigFlag) {
						draw_bigPath(fp,id);
						draw_bigEvent(fp,id);
					}
					else {
						draw_path(fp,id);
						draw_event(fp,id);
					}
				}
				if (type== 2) {
					draw_smallPath(fp,id);
					draw_smallEvent(fp,id);
				}
			}
			// path printing end
			draw_images(fp);

			fputs(last_text.c_str(),fp);
		}
};


int main(int argc, char** argv) {
	char filename[100];
	if(argc <2) {
		perror("#run inputfilename\n ENVIRON RADIUS is big radius.");
		return 1;
	}
	char* env = getenv("RADIUS");

	if (env!= NULL) {
		big_radius = atof(env);
		if( big_radius<0.0 || big_radius > 1e3)
			big_radius = BIG_RADIUS;
	}
	
	strcpy( filename,argv[1]);
	FILE* fp = fopen( filename ,"r");
	std::list<Snapshot*> snaplist;
	Snapshot* snap;

	while(1)
	{
		snap =	read_dump(fp);

		if (snap==NULL) {
			perror("read dump end");
			break;
		}
		snaplist.push_back(snap);
	}

	if (snaplist.size() <1){
		perror("The # of snap is too small(<1)!!");
		return 23;
	}

	PrintX3D a(&snaplist);
	a.drawX3D(stdout);
	/* 	makeRDF makerRdf(snaplist);
	 * 	makerRdf.flag_anisotropy =1;
	 * 	makerRdf.calcRDF();
	 * 			puts("calc RDF end");
	 * 	makerRdf.calcSSF();
	 */

}







	PrintX3D::PrintX3D
(std::list<Snapshot*> *_sl)
{
	snaplist = *_sl ;
	auto firstsnap = (snaplist.begin());
	box     = &(*firstsnap)->box;
	box_x = box->xhigh-box->xlow; 
	box_y = box->yhigh-box->ylow; 
	box_z = box->zhigh-box->zlow;
	memcpy( pbc, box->pbc, 3);
	
	maxSnap = snaplist.size();
	maxAtom = (*firstsnap)->n_atoms;
	
	firstsnap_shift = new real3d[maxAtom];
	colPath = new string[maxAtom+1];
	colOriPath = new string[maxAtom+1];
}

PrintX3D::~PrintX3D()
{
	delete [] firstsnap_shift;
	delete [] colPath;
	delete [] colOriPath;
} 

