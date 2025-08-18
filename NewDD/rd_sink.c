/* NAME: RD_PaRT
 * PURPOSE: This proceduer reads particles from a RAMSES PART file. And it is rewritten from the IDL version of rd_part.pro.
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "ramses.h"
#include "Memory.h"



int rd_sink(RamsesType *ram, char *infile){
	FILE *fp = fopen(infile,"r");
	int chip,i,j,k;
	int npartp,mstar;


	F77read(&(ram->nsink), sizeof(int), 1, fp);
	F77read(&(ram->nindsink), sizeof(int), 1, fp);



	npartp = ram->nsink;
	dptype *xbuff = (dptype*)Malloc(sizeof(dptype)*npartp,PPTR(xbuff));
	int *idbuff = (int*)Malloc(sizeof(int)*npartp,PPTR(idbuff));
	ram->sink = (SinkType*)Malloc(sizeof(SinkType)*npartp,PPTR(ram->sink));


	GetPart(idbuff,sizeof(int), npartp, fp, ram,sink,id);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,mass);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,x);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,y);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,z);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,vx);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,vy);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,vz);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,tbirth);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,dMsmbh);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,dMBH_coarse);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,dMEd_coarse);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,Esave);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,Jx);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,Jy);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,Jz);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,Sx);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,Sy);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,Sz);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,Smag);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,eps);

	/*
	for(i=0;i<ram->ndim*2+1;i++){
		for(j=0;j<ram->nlevelmax-ram->levelmin+1;j++){
			GetPart(xbuff,sizeof(dptype), npartp, fp, ram,sink,stat[i][j]);
		}
	}
	*/

	SinkType *sink = ram->sink;
	for(i=0;i<npartp;i++){
		sink[i].mass *= ram->scale_m;
		sink[i].x *= ram->mpcscale_l;
		sink[i].y *= ram->mpcscale_l;
		sink[i].z *= ram->mpcscale_l;
		sink[i].vx *= ram->kmscale_v;
		sink[i].vy *= ram->kmscale_v;
		sink[i].vz *= ram->kmscale_v;
		sink[i].dMsmbh *= ram->scale_m;
		sink[i].dMBH_coarse *= ram->scale_m;
		sink[i].dMEd_coarse *= ram->scale_m;
		sink[i].Jx *= ram->scale_m *ram->scale_l/kpc * ram->kmscale_v;
		sink[i].Jy *= ram->scale_m *ram->scale_l/kpc * ram->kmscale_v;
		sink[i].Jz *= ram->scale_m *ram->scale_l/kpc * ram->kmscale_v;
		sink[i].Smag *= ram->scale_m *ram->scale_l/kpc * ram->kmscale_v;
	}
	fclose(fp);
	Free(idbuff);
	Free(xbuff);
	return npartp;
}
