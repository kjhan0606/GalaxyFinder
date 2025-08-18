/* 
icc -g -o check check.c -lm  -DINDEX -DVarPM   -DXYZDBL  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL 
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
//#include<mpi.h>
#include "Memory.h"
#include "header.h"
#include "ramses.h"
#include "tree.h"
#include "defs.h"
//#include "hfind.h"
#include "galcenter.h"

/*
#define dptype double
#define idtype long
*/



#define NGAL 1000000

GalCenter gal[NGAL];

int main(int argc, char **argv){

	int i,j,k;
	int mp,np;
	char infile[190];

	sprintf(infile,"./FoF_Data/FoF.%.5d/GALFIND.CENTER.%.5d",atoi(argv[1]),atoi(argv[1]));

	FILE *fp = fopen(infile,"r");


	while( (np = fread(gal, sizeof(GalCenter), NGAL,fp))>0){
		for(i=0;i<np;i++){
			dptype xxx = gal[i].gal.x*gal[i].gal.y*gal[i].gal.z;
			xxx *= gal[i].dmhalo.x *gal[i].dmhalo.y *gal[i].dmhalo.z;
			xxx *= gal[i].gas.x *gal[i].gas.y *gal[i].gas.z;
			if(!isfinite(xxx)){
				printf("error in the data\n");
			}
		}
	}
}

