/*
 
icc -g -o how2readfof how2readfof.c -DINDEX -DVarPM   -DXYZDBL  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL  -lm
 
*/
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
/*
#include "header.h"
#include "tree.h"
#include "defs.h"
*/
#include "ramses.h"
#include "fof.h"





void FREAD(HaloQ *haloq, FILE *fp, DmType *dm, GasType *gas,
		SinkType *sink, StarType *star){
    fread(dm, sizeof(DmType), haloq->npdm,fp);
    fread(gas, sizeof(GasType), haloq->npgas,fp);
    fread(sink, sizeof(SinkType), haloq->npsink,fp);
    fread(star, sizeof(StarType), haloq->npstar,fp);
}

int main(int argc, char *argv[]) {

	int nstep = atoi(argv[1]);
	char rheader[190];
	char dir[190];

	sprintf(dir,"./FoF_Data/FoF.%.5d/",nstep);
	sprintf(rheader,"%s/FoF_halo_cat.%.5d",dir,nstep);
	FILE *rhfp = fopen(rheader,"r");
	float bias,astep,anow,npower;
	char rfilename[190];
	sprintf(rfilename,"%sFoF_member_particle.%.5d",dir,nstep); 
	FILE *fp = fopen(rfilename,"r"); 
	{ 
		float bias; 
		float astep,anow; 
		int npower; 
		float size, hubble, omep,omepb,omeplam,amax;
		fread(&size,sizeof(float),1,rhfp); 
		fread(&hubble,sizeof(float),1,rhfp); 
		fread(&omep,sizeof(float),1,rhfp); 
		fread(&omepb,sizeof(float),1,rhfp); 
		fread(&omeplam,sizeof(float),1,rhfp); 
		fread(&amax,sizeof(float),1,rhfp); 
		fread(&anow,sizeof(float),1,rhfp); 
	}



	HaloQ haloq;
	while(fread(&haloq,sizeof(HaloQ), 1, rhfp) == 1){

		DmType *dm = (DmType*)malloc(sizeof(DmType)*haloq.npdm);
		GasType *gas = (GasType*)malloc(sizeof(GasType)*haloq.npgas);
		StarType *star = (StarType*)malloc(sizeof(StarType)*haloq.npstar);
		SinkType *sink = (SinkType*)malloc(sizeof(SinkType)*haloq.npsink);
		FREAD(&haloq,fp, dm, gas, sink,star);

		/* Do something if you want here */
		free(dm);free(gas);free(star);free(sink);
	}


}
