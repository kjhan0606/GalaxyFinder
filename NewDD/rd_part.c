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




int rd_part(RamsesType *ram, char *infile){
	FILE *fp = fopen(infile,"r");
	size_t i;
	int chip,nchem;
	int npartp,mstar;
	F77read(&(ram->ncpu), sizeof(int), 1, fp);
	F77read(&(ram->ndim), sizeof(int), 1, fp);
	F77read(&(ram->npart), sizeof(int), 1, fp);
	npartp = ram->npart;
	F77read(&(ram->localseed), sizeof(dptype), IRandNumSize, fp);
	/*
	F77read(&(ram->nstar_tot), sizeof(long), 1, fp);
	*/
	
 	F77read(&(ram->nstar_tot), sizeof(int), 1, fp);
	
	F77read(&(ram->mstar_tot), sizeof(dptype), 1, fp);
	F77read(&(ram->mstar_lost), sizeof(dptype), 1, fp);
	F77read(&(ram->nsink), sizeof(int), 1, fp);
	dptype *xbuff = (dptype*)Malloc(sizeof(dptype)*npartp,PPTR(xbuff));
#ifdef NCHEM
	nchem=NCHEM;
#endif
	ram->particle = (PmType*)Malloc(sizeof(PmType)*npartp,PPTR(ram->particle));


	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,x);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,y);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,z);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,vx);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,vy);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,vz);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,mass);

	idtype *idbuff = (idtype*)Malloc(sizeof(idtype)*npartp, PPTR(idbuff));
	GetPart(idbuff,sizeof(idtype), npartp, fp, ram,particle,id);
	int *ibuff = (int*)Malloc(sizeof(int)*npartp,PPTR(ibuff));
	GetPart(ibuff,sizeof(int), npartp, fp, ram,particle,levelp);

	familytype *bytbuff = (familytype*)Malloc(sizeof(familytype)*npartp, PPTR(bytbuff));
	GetPart(bytbuff,sizeof(familytype), npartp, fp, ram,particle,family);
	GetPart(bytbuff,sizeof(familytype), npartp, fp, ram,particle,tag);

#ifdef OUTPUT_PARTICLE_POTENTIAL
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,potent);
#endif

#ifndef NBODY
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,tp);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,zp);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,mass0);
#ifdef NCHEM
	for(i=0;i<nchem;i++){
		GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,chem[i]);
	}
#endif
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,birth_d);
	GetPart(ibuff,sizeof(int), npartp, fp, ram,particle,partp);
#endif
	fclose(fp);
	Free(ibuff);
	Free(idbuff);
	Free(bytbuff);
	Free(xbuff);
	PmType *part = ram->particle;
	int ipartp = 0;
	long nstar = 0;
	for(i=0;i<npartp;i++){
		if(part[i].family <= 2) {
			part[ipartp] = part[i];
			part[ipartp].mass = part[i].mass* ram->scale_m;
			part[ipartp].x = part[i].x * ram->mpcscale_l;
			part[ipartp].y = part[i].y * ram->mpcscale_l;
			part[ipartp].z = part[i].z * ram->mpcscale_l;
			part[ipartp].vx = part[i].vx * ram->kmscale_v;
			part[ipartp].vy = part[i].vy * ram->kmscale_v;
			part[ipartp].vz = part[i].vz * ram->kmscale_v;
#ifndef NBODY
			part[ipartp].mass0 = part[i].mass0* ram->scale_m;
			if(part[ipartp].family == 2) nstar ++;
#endif
			ipartp ++;
		}
	}
	ram->npart = (npartp= ipartp);
	part = (ram->particle = (PmType*)Realloc(ram->particle, sizeof(PmType)*npartp));

	printf("Total star and dm particles are %d %d from total np=  %d\n", nstar, npartp-nstar, npartp);
	return npartp;
}
