/* NAME: RD_PaRT
 * PURPOSE: This proceduer reads particles from a RAMSES PART file. And it is rewritten from the IDL version of rd_part.pro.
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>

#include "ramses.h"
#include "Memory.h"

#define MAX(a,b) ( (a)>(b) ? (a): (b) )


#define TWICE(a) ( (a)*(a))
#define TRIPLE(a) ( (a)*(a)*(a))

void GetHydroQ(RamsesType *ram, int *ind_leaf, int ncache, HydroCellType *cell){
	int ngridmax = ram->ngridmax;
	int ndim = ram->ndim;
	HydroType *hydro=&(ram->hydro);
	dptype *uold = (hydro->uold);
//	dptype ekk;
	dptype *temp;
	dptype smallr = ram->smallr;
	int idim,i;
	int imetal = ram->imetal;
	for(i=0;i<ncache;i++){
	//	ekk = 0;
		dptype Hden = MAX(uold[ind_leaf[i]-1],smallr);
		cell[i].metallicity = uold[ind_leaf[i]-1+ngridmax*imetal]/Hden/0.02;
//		for(idim=0;idim<ndim;idim++){
//			ekk += 0.5*TWICE(uold[ind_leaf[i]-1+ngridmax*idim])/Hden;
//		}
		dptype Temp = uold[ind_leaf[i]-1+ngridmax*(ndim+1)];
		dptype err=0;
#if NENER>0
		int irad;
		for(irad=0;irad<nener;irad++){
			err += uold[ind_leaf[i]-1 + ngridmax*(inener+irad)];
		}
#endif
		dptype emag=0;
#ifdef SOLVERmhd
		for(idim=0;idim<ndim;idim++){
			emag += 0.125L*TWICE(uold[ind_leaf[i]-1+ngridmax*(idim+ndim+1)] +
					uold[ind_leaf[i]-1+ngridmax*(idim+nvar)]);
		}
#endif
//		Temp = ram->gamma*(Temp - ekk - err -emag);
		Temp = Temp/Hden*ram->scale_T2;
		cell[i].temp = Temp;
		cell[i].den = Hden;
	}
}

void units(RamsesType *ram){
	dptype scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2, scale_m;
	int cosmo = ram->cosmo;
  /*-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------*/

  // scale_d converts mass density from user units into g/cc
  scale_d = ram->unit_d;
  if(cosmo) scale_d = ram->omega_m * rhoc *TWICE(ram->H0/100.L) / TRIPLE(ram->aexp);

  // scale_t converts time from user units into seconds
  scale_t = ram->unit_t;
  if(cosmo) scale_t = TWICE(ram->aexp) / (ram->H0*1e5L/Mpc);


  // scale_l converts distance from user units into cm
  scale_l = ram->unit_l;
  if(cosmo) scale_l = ram->aexp * ram->boxlen_ini * Mpc / (ram->H0/100L);


  // scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t;

  // scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v*scale_v;

  // scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d;
  scale_m =scale_d*scale_l*scale_l*scale_l/Msun * ram->H0/100.L;

  ram->scale_d = scale_d;
  ram->scale_t = scale_t;
  ram->scale_l = scale_l;
  ram->scale_v = scale_v;
  ram->scale_T2 = scale_T2;
  ram->scale_nH = scale_nH;
  ram->scale_m = scale_m;

  ram->mpcscale_l = ram->boxlen_ini; /* in cMpc/h */
  ram->kmscale_v = scale_v /1.E5L; /* in km/second */
  ram->scale_Gyr = scale_t /3600.L/24.L/365.2425L/1.E9L;
}



void rd_info(RamsesType *header, char *infile){
	FILE *fp = fopen(infile,"r");
	int chip,i,j,k,istep;
	int npartp,mstar;
	char in1[190], in2[190], in[190];

	header->smallr = 1.e-10L;
	header->cosmo = YES;
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->ncpu = atoi(in);
//	fscanf(fp, "%s %s %s\n", in1, in2, in); header->nthr = atoi(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->ndim = atoi(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->levelmin = atoi(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->nlevelmax = atoi(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->ngridmax = atol(in);
	fscanf(fp, "%s %s\n", in1, in); header->nstep_coarse = atoi(in);
	fscanf(fp,"\n");
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->boxlen = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->time = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->aexp = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->H0 = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->omega_m = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->omega_l = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->omega_k = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->omega_b = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->unit_l = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->unit_d = atof(in);
	fscanf(fp, "%s %s %s\n", in1, in2, in); header->unit_t = atof(in);
	fscanf(fp,"\n");
	fscanf(fp, "ordering type=%s\n", header->ordering);
	fclose(fp);
	header->twotondim = 1;
	for(i=0;i<header->ndim;i++) header->twotondim *= 2;

	sscanf(infile,"./output_%5d/info_%d.txt", &istep, &(header->nrestart));
	header->nrestart_quad = 0;
	header->overload = 1;
#ifndef NENER
	int nener=0;
#else
	int nener=NENER;
#endif
	header->nener = nener;

//#ifndef NVAR
	header->nvar=header->ndim+2+nener;
//#else
//	header->nvar=NVAR;
//#endif

//#ifdef NCHEM
//	header->nvar=header->nvar+NCHEM
//#endif

	for(i=0;i<512;i++) header->gamma_rad[i] = 1.33333333334L;
	header->gamma = 1.4L;
	header->neq_chem = NO;
#ifdef RT
	header->rt = YES;
#else
	header->rt = NO;
#endif 
	header->imetal = header->nener + header->ndim + 2;
	header->ichem = header->imetal + 1;
	header->idust = header->ichem + 1;
#ifdef NCHEM
	header->idust = header->imetal + 1 +NCHEM;
#endif

}


void cleanup_mesh(RamsesType *header, int simple_boundary){
	MeshType *mesh = &(header->mesh);
#ifndef NBODY
	Free(mesh->headl);
	Free(mesh->taill);
	Free(mesh->numbl);
	Free(mesh->numbtot);
	Free(mesh->father);
	Free(mesh->nbor);

	if(simple_boundary==YES){
		Free(mesh->headb);
		Free(mesh->tailb);
		Free(mesh->numbb);
	}

	Free(mesh->flag1);
	Free(mesh->flag2);
	Free(mesh->son);
	Free(mesh->cpu_map);
	Free(mesh->cpu_map2);


	Free(mesh->xg);
	Free(mesh->prev);
	Free(mesh->next);
	Free(mesh->bound_key);
#endif
	Free(mesh->ilevel);
	Free(mesh->dtold);
	Free(mesh->dtnew);
}

void cleanup_ramses(RamsesType *header){
	Free(header->particle);
	
	/*
	Free(header->hcell);
	*/
#ifndef NBODY
	Free(header->gas);
//	Free(header->hydro.phi);
//	Free(header->hydro.fxyz);
	Free(header->hydro.uold);
#endif
	cleanup_mesh(header, NO);
}
