#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define NODE_HAVE_PARTICLE 12
#define EPSILON (0.0)
#define MinNumMem 20
#define MaxLinkedParticles 1000000

#define minstarmass 1.e6

#define INDXTYPE long long
#ifdef XYZDBL
#	define POSTYPE double
#else
#	define POSTYPE float
#endif

enum {TYPE_TREE = 0,TYPE_PTL = 1, TYPE_STAR = 1, TYPE_DM=2, TYPE_GAS=3, TYPE_SINK=4, TYPE_AGN=4};


typedef struct Box{
	dptype x,y,z,width;
} Box;

typedef struct HaloQ{
    size_t np,npstar,npgas,npdm,npsink;
    POSTYPE x,y,z;
    float mass, mstar,mgas,mdm,msink;
    float vx,vy,vz;
}HaloQ;

/*      */
#define GENERAL_TYPE unsigned int type
#define GENERALTPtlPOINTER GENERAL_TYPE;void *sibling
/*      */
typedef struct TYPE {
	GENERAL_TYPE;
} TYPE;
typedef struct GENERAL_TPtl_POINTER {
	GENERALTPtlPOINTER;
} GENERAL_TPtl_POINTER;
typedef struct TPtlStruct{
	GENERALTPtlPOINTER;
	dptype x,y,z;
	dptype mass;
	int indx;
} TPtlStruct;

typedef struct TStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	dptype L;
	dptype dist2;
	dptype x0,y0,z0;
	dptype mass;
	int Nparticle;
	dptype monox,monoy,monoz;
	dptype mrr;
	dptype quad[6];
} TStruct;

typedef struct BeginEndTree{
	TStruct *start;
	void *EndTree;
	void *EndPtl;
} BeginEndTree;
/*      */
typedef struct DMParticle{
	dptype x,y,z;
	dptype vx,vy,vz;
} DMParticle;


typedef struct FoFTStruct{
	unsigned int type;
	void *sibling;
	POSTYPE L;
	POSTYPE dist2;
	POSTYPE dist;
	POSTYPE x0,y0,z0;
	POSTYPE monox,monoy,monoz;
	int Nparticle;
	void *daughter;
	float maxlink02;
} FoFTStruct;
typedef struct FoFBeginEndTree{
	FoFTStruct *start;
	void *EndTree;
	void *EndPtl;
} FoFBeginEndTree;


/* for Fof */
typedef struct FoFTPtlStruct{
	unsigned int type;
//	void *sibling;
//	enum boolean included;
	dptype x,y,z,vx,vy,vz,mass;
//	float link02;
/*
	union {
		DmType dm;
		StarType star;
		SinkType sink;
		GasType gas;
	}p;
	size_t haloindx;
	int indx;
*/
} FoFTPtlStruct;

typedef struct FoFTPtlStruct particle;
/*
typedef particle FoFTPtlStruct;
*/


typedef struct SimpleBasicParticleType{
	unsigned int type;
	dptype mass,x,y,z;
	dptype vx,vy,vz;
	dptype link02;
	idtype indx;
} SimpleBasicParticleType;
typedef struct pforce {
	dptype x,y,z;
} pforce;

#define EraseFromTree(optr,ptr,nptr) {\
	switch(((TYPE*)optr)->type) {\
		case TYPE_TREE:\
			if(((FoFTStruct*)optr)->daughter == ptr) \
				((FoFTStruct*)optr)->daughter == nptr;\
			else ((FoFTStruct*)optr)->sibling = nptr;\
			break;\
		default :\
			((FoFTPtlStruct*)ptr)->sibling = nptr;\
	}\
}\
	
void treeforce(particle*,float, TStruct *,TPtlStruct *,pforce *);
void Make_Tree(TStruct *,TPtlStruct *,int ,Box );
float treeplumpotential(particle*,float, TStruct *,TPtlStruct *);
BeginEndTree divide_node(TStruct *,TStruct *, TPtlStruct *, Box ,TStruct *); 

typedef struct HaloBound{
	size_t nmem;
	POSTYPE zmin,zmax;
	int boundflag;
	FoFTPtlStruct *sibling;
} HaloBound;



#define EraseFromTree(optr,ptr,nptr) {\
	switch(((TYPE*)optr)->type) {\
		case TYPE_TREE:\
			if(((FoFTStruct*)optr)->daughter == ptr) \
				((FoFTStruct*)optr)->daughter == nptr;\
			else ((FoFTStruct*)optr)->sibling = nptr;\
			break;\
		default :\
			((FoFTPtlStruct*)ptr)->sibling = nptr;\
	}\
}\
	
FoFBeginEndTree FoF_divide_node(FoFTStruct *,FoFTStruct *, FoFTPtlStruct *, 
		Box ,FoFTStruct *); 
void FoF_Make_Tree(FoFTStruct *,FoFTPtlStruct *,size_t ,Box );
int new_fof_link(particle*,POSTYPE, FoFTStruct *, FoFTPtlStruct *,particle *);
size_t pnew_fof_link(particle*,POSTYPE, FoFTStruct *, FoFTPtlStruct *,particle *,
		size_t nhalo, POSTYPE,POSTYPE,POSTYPE);
