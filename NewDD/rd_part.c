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

static dptype periodic_position(dptype x, dptype boxlen){
	/* Float-backed particle structs can round a position just below one to
	 * exactly boxlen after unit conversion.  Keep every coordinate in the
	 * half-open periodic interval required by SplitDump and opFoF. */
	if(x >= boxlen || x < 0) {
		x = fmod(x, boxlen);
		if(x < 0) x += boxlen;
	}
	return x;
}




int rd_part(RamsesType *ram, char *infile){
	FILE *fp = fopen(infile,"r");
	size_t i;
	int chip,nchem;
	int npartp,mstar;
	F77read(&(ram->ncpu), sizeof(int), 1, fp);
	F77read(&(ram->ndim), sizeof(int), 1, fp);
	F77read(&(ram->npart), sizeof(int), 1, fp);
	npartp = ram->npart;
	/* RAMSES writes four local seeds, or eight when MC_tracer also writes
	 * tracer_seed, in one unformatted record.  Accept either record length. */
	{
		int rec_bytes = 0, rec_end = 0;
		int seedbuf[2*IRandNumSize];
		fread(&rec_bytes, sizeof(int), 1, fp);
		if(rec_bytes != (int)(sizeof(int)*IRandNumSize) &&
		   rec_bytes != (int)(2*sizeof(int)*IRandNumSize)) {
			fprintf(stderr, "Unexpected particle seed record length %d\n", rec_bytes);
			exit(99);
		}
		fread(seedbuf, sizeof(int), rec_bytes/sizeof(int), fp);
		fread(&rec_end, sizeof(int), 1, fp);
		if(rec_end != rec_bytes) { fprintf(stderr, "Bad particle seed record\n"); exit(99); }
		memcpy(ram->localseed, seedbuf, sizeof(ram->localseed));
	}
	/*
	F77read(&(ram->nstar_tot), sizeof(long), 1, fp);
	*/
	
	{
		long nstar_file = 0;
		F77read(&nstar_file, sizeof(long), 1, fp);
		ram->nstar_tot = (int)nstar_file;
	}
	
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
	/* Particle IDs are the first standard record after mass.  Preserve them
	 * even in the DM-only fast path: they are needed to trace a z=0 halo back
	 * to its Lagrangian cells for the zoom mask. */
	{
		int rec_bytes = 0, rec_end = 0;
		fread(&rec_bytes, sizeof(int), 1, fp);
		if(rec_bytes == (int)(sizeof(idtype)*npartp)) {
			fread(idbuff, sizeof(idtype), npartp, fp);
		} else if(rec_bytes == (int)(sizeof(int)*npartp)) {
			int *id32 = (int*)Malloc(sizeof(int)*npartp, PPTR(id32));
			fread(id32, sizeof(int), npartp, fp);
			for(i=0;i<(size_t)npartp;i++) idbuff[i] = (idtype)id32[i];
			Free(id32);
		} else if(rec_bytes == (int)(sizeof(long long)*npartp)) {
			long long *id64 = (long long*)Malloc(sizeof(long long)*npartp, PPTR(id64));
			fread(id64, sizeof(long long), npartp, fp);
			for(i=0;i<(size_t)npartp;i++) idbuff[i] = (idtype)id64[i];
			Free(id64);
		} else {
			fprintf(stderr, "Unexpected particle ID record length %d\n", rec_bytes);
			exit(99);
		}
		fread(&rec_end, sizeof(int), 1, fp);
		if(rec_end != rec_bytes) { fprintf(stderr, "Bad particle ID record\n"); exit(99); }
	}
	if(getenv("NEWDD_DM_ONLY") != NULL) {
		/* The lagRamses DM-only production snapshot may carry additional
		 * legacy/type records after the standard particle ID with ABI-dependent
		 * widths.  Preserve positions, velocities, masses, and IDs; classify
		 * every retained particle as DM and stop before optional records. */
		for(i=0;i<(size_t)npartp;i++) {
			ram->particle[i].id = idbuff[i];
			ram->particle[i].levelp = 1;
			ram->particle[i].family = 1;
			ram->particle[i].tag = 0;
			ram->particle[i].mass *= ram->scale_m;
			ram->particle[i].x *= ram->mpcscale_l;
			ram->particle[i].y *= ram->mpcscale_l;
			ram->particle[i].z *= ram->mpcscale_l;
			ram->particle[i].x = periodic_position(ram->particle[i].x, ram->boxlen_ini);
			ram->particle[i].y = periodic_position(ram->particle[i].y, ram->boxlen_ini);
			ram->particle[i].z = periodic_position(ram->particle[i].z, ram->boxlen_ini);
			ram->particle[i].vx *= ram->kmscale_v;
			ram->particle[i].vy *= ram->kmscale_v;
			ram->particle[i].vz *= ram->kmscale_v;
		}
		fclose(fp);
		Free(idbuff);
		Free(xbuff);
		ram->npart = npartp;
		return npartp;
	}

	int *ibuff = (int*)Malloc(sizeof(int)*npartp,PPTR(ibuff));
	GetPart(ibuff,sizeof(int), npartp, fp, ram,particle,levelp);

	familytype *bytbuff = (familytype*)Malloc(sizeof(familytype)*npartp, PPTR(bytbuff));
	/* family/tag are byte records in the lagRamses dump.  Read the record
	 * markers explicitly rather than relying on the legacy macro's type-size
	 * assumptions. */
	for(int ibyte_field=0; ibyte_field<2; ibyte_field++) {
		int rec_bytes=0, rec_end=0;
		fread(&rec_bytes,sizeof(int),1,fp);
		if(rec_bytes != npartp) {
			fprintf(stderr,"Unexpected family/tag record length %d\n",rec_bytes);
			exit(99);
		}
		fread(bytbuff,sizeof(familytype),npartp,fp);
		fread(&rec_end,sizeof(int),1,fp);
		if(rec_end != rec_bytes) { fprintf(stderr,"Bad family/tag record\n"); exit(99); }
		for(i=0;i<(size_t)npartp;i++) {
			if(ibyte_field==0) ram->particle[i].family=bytbuff[i];
			else ram->particle[i].tag=bytbuff[i];
		}
	}

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
	/* Legacy birth-time field remains on disk but is not retained in PmType. */
	if(npartp>0) F77read(xbuff,sizeof(dptype),npartp,fp);
	/*
	 * The current PmType intentionally omits the legacy RAMSES `partp`
	 * member (see ramses.h: the on-disk field is still present in old
	 * snapshots).  Consume the record to keep the Fortran-unformatted file
	 * position aligned, but do not write into a non-existent struct member.
	 */
	if(npartp>0) F77read(ibuff,sizeof(int),npartp,fp);
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

	printf("Total star and dm particles are %ld %ld from total np=  %d\n", nstar, npartp-nstar, npartp);
	return npartp;
}
