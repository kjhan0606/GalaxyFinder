
typedef struct Center{
	dptype x,y,z,R;
}Center;
typedef struct GalCenter{
	idtype haloid, subgalid;
	Center gal, dmhalo,gas;
}GalCenter;


#define NO_CENTER_POS -999999
#define MIN_SHIFT 0.005
#define rSTEP 0.142857142857142857 /* ==1/7 */
#define RSTEP 0.142857142857142857

void find_centers(GalCenter *, FoFTPtlStruct *, lint);
