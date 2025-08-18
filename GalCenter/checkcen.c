#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>

typedef double dptype;
typedef long idtype;



typedef struct Center{
	dptype x,y,z,R;
}Center;
typedef struct GalCenter{
	idtype haloid, subgalid;
	Center gal, dmhalo,gas;
}GalCenter;


int main(int argc, char **argv){
	FILE *fp = fopen(argv[1],"r");



	GalCenter galcen[1000000];


	int np;

	while( (np=fread(galcen, sizeof(GalCenter), 1000000, fp))>0){
		int i;
		for(i=0;i<np;i++){
			if(( galcen[i].gal.x > -999 && galcen[i].gal.x <0) || galcen[i].gal.x > 717.225){
				printf("Please check this %g %g %g\n", galcen[i].gal.x, galcen[i].dmhalo.x, galcen[i].gas.x);
			}
		}
	}



	fclose(fp);

	return 0;
}
