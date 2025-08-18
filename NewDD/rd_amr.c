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



int rd_amr(RamsesType *ram, char *infile, int simple_boundary){
	FILE *fp = fopen(infile,"r");
	int ilevel,ibound;
	int chip,i;
	int npartp,mstar,ng;
	int *next = ram->mesh.next;
	int ncpu2,ndim2, nx2,ny2,nz2,nlevelmax2,ngridmax2,nboundary2,noutput2;
	int nlevelmax,ncpu;
	int twondim;
	int ndim;
	size_t ncell;
	MeshType *mesh;

	mesh = &(ram->mesh);

	ndim = ram->ndim;
	ncpu = ram->ncpu;
	nlevelmax = ram->nlevelmax;

	ram->ndomain = ram->ncpu*ram->overload;
	twondim = 2*ndim;




	printf("Opening %s\n", infile);
	F77read(&(ncpu2), sizeof(int), 1, fp);
	F77read(&(ndim2), sizeof(int), 1, fp);
	{
		int _nchip;
		fread(&_nchip, sizeof(int), 1, fp);
		fread(&(nx2), sizeof(int), 1, fp);
		fread(&(ny2), sizeof(int), 1, fp);
		fread(&(nz2), sizeof(int), 1, fp);
		fread(&_nchip, sizeof(int), 1, fp);
	}
	ram->nx = nx2;
	ram->ny = ny2;
	ram->nz = nz2;

    ram->ncoarse = nx2*ny2*nz2;

	ncell = ram->ncoarse + ram->twotondim * ram->ngridmax;

#ifndef NBODY

	mesh->flag1 = (int*)Malloc(sizeof(int)*(1+ncell),PPTR(mesh->flag1));
	mesh->flag2 = (int*)Malloc(sizeof(int)*(1+ncell),PPTR(mesh->flag2));
	mesh->son = (int*)Malloc(sizeof(int)*ncell,PPTR(mesh->son));
	mesh->cpu_map = (int*)Malloc(sizeof(int)*ncell,PPTR(mesh->cpu_map));
	mesh->cpu_map2 = (int*)Malloc(sizeof(int)*ncell,PPTR(mesh->cpu_map2));
#endif

	ram->nbinodes = 1;
	for(i=0;i<ndim;i++) ram->nbinodes *= 2;
	ram->nbinodes -= 1;

#ifndef NBODY
	for(i=0;i<ncell;i++){
		mesh->son[i] = 0;
		mesh->flag1[i] = 0;
		mesh->flag2[i] = 0;
	}
#endif


	F77read(&(nlevelmax2), sizeof(int), 1, fp);
	F77read(&(ngridmax2), sizeof(int), 1, fp);
	F77read(&(nboundary2), sizeof(int), 1, fp);
	/*         */
	ram->nboundary = nboundary2;
	F77read(&(ram->ngrid_current), sizeof(int), 1, fp);
	F77read(&(ram->boxlen), sizeof(dptype), 1, fp);
	{
		int _nchip;
		fread(&_nchip, sizeof(int), 1, fp);
		fread(&(noutput2), sizeof(int), 1, fp);
		fread(&(ram->iout), sizeof(int), 1, fp);
		fread(&(ram->ifout), sizeof(int), 1, fp);
		fread(&_nchip, sizeof(int), 1, fp);
	}
	if(1){
		long offset = sizeof(int)*4 + sizeof(dptype)*noutput2*2;
		fseek(fp, offset, SEEK_CUR);
	}
	else {
		/*
		Fread(&(ram->tout), sizeof(dptype), ram->noutput, fp);
		Fread(&(ram->aout), sizeof(dptype), ram->noutput, fp);
		*/
	}

	mesh->ilevel = (int*)Malloc(sizeof(int)*ram->ngridmax,PPTR(mesh->ilevel));
	for(i=0;i<ram->ngridmax;i++) mesh->ilevel[i] = 0;

	F77read(&(ram->t), sizeof(dptype), 1, fp);
	mesh->dtold = (dptype*)Malloc(sizeof(dptype)*nlevelmax2,PPTR(mesh->dtold));
	mesh->dtnew = (dptype*)Malloc(sizeof(dptype)*nlevelmax2,PPTR(mesh->dtnew));
	F77read((mesh->dtold), sizeof(dptype), nlevelmax2, fp);
	F77read((mesh->dtnew), sizeof(dptype), nlevelmax2, fp);
	{
		int _nchip;
		fread(&_nchip, sizeof(int), 1, fp);
		fread(&(ram->nstep), sizeof(int), 1, fp);
		fread(&(ram->nstep_coarse), sizeof(int), 1, fp);
		fread(&_nchip, sizeof(int), 1, fp);
	}
	{
		int _nchip;
		fread(&_nchip, sizeof(int), 1, fp);
		fread(&(ram->constant), sizeof(dptype), 1, fp);
		fread(&(ram->mass_tot_0), sizeof(dptype), 1, fp);
		fread(&(ram->rho_tot), sizeof(dptype), 1, fp);
		fread(&_nchip, sizeof(int), 1, fp);
	}
	{
		int _nchip;
		fread(&_nchip, sizeof(int), 1, fp);
		fread(&(ram->omega_m), sizeof(dptype), 1, fp);
		fread(&(ram->omega_l), sizeof(dptype), 1, fp);
		fread(&(ram->omega_k), sizeof(dptype), 1, fp);
		fread(&(ram->omega_b), sizeof(dptype), 1, fp);
		fread(&(ram->omega_h0), sizeof(dptype), 1, fp);
		fread(&(ram->aexp_ini), sizeof(dptype), 1, fp);
		fread(&(ram->boxlen_ini), sizeof(dptype), 1, fp);
		fread(&_nchip, sizeof(int), 1, fp);
	}
#ifndef NBODY	
	{
		int _nchip;
		fread(&_nchip, sizeof(int), 1, fp);
		fread(&(ram->aexp), sizeof(dptype), 1, fp);
		fread(&(ram->hexp), sizeof(dptype), 1, fp);
		fread(&(ram->aexp_old), sizeof(dptype), 1, fp);
		fread(&(ram->epot_tot_int), sizeof(dptype), 1, fp);
		fread(&(ram->epot_tot_old), sizeof(dptype), 1, fp);
		fread(&_nchip, sizeof(int), 1, fp);
	}
	F77read(&(ram->mass_sph), sizeof(dptype), 1, fp);
	mesh->headl = (int*)Malloc(sizeof(int)*ram->ncpu*nlevelmax2,PPTR(mesh->headl));
	mesh->taill = (int*)Malloc(sizeof(int)*ram->ncpu*nlevelmax2,PPTR(mesh->taill));
	mesh->numbl = (int*)Malloc(sizeof(int)*ram->ncpu*nlevelmax2,PPTR(mesh->numbl));
	mesh->numbtot = (int*)Malloc(sizeof(int)*10*nlevelmax2,PPTR(mesh->numbtot));
	F77read((mesh->headl), sizeof(int), ram->ncpu*nlevelmax2, fp);
	F77read((mesh->taill), sizeof(int), ram->ncpu*nlevelmax2, fp);
	F77read((mesh->numbl), sizeof(int), ram->ncpu*nlevelmax2, fp);
	F77read((mesh->numbtot), sizeof(int), 10*nlevelmax2, fp);




	if(simple_boundary == YES) {
		mesh->headb = (int*)Malloc(sizeof(int)*ram->nboundary*nlevelmax2,PPTR(mesh->headb));
		mesh->tailb = (int*)Malloc(sizeof(int)*ram->nboundary*nlevelmax2,PPTR(mesh->tailb));
		mesh->numbb = (int*)Malloc(sizeof(int)*ram->nboundary*nlevelmax2,PPTR(mesh->numbb));
		F77read((mesh->headb), sizeof(int), ram->nboundary*nlevelmax2, fp);
		F77read((mesh->tailb), sizeof(int), ram->nboundary*nlevelmax2, fp);
		F77read((mesh->numbb), sizeof(int), ram->nboundary*nlevelmax2, fp);
	}

	mesh->xg = (dptype*)Malloc(sizeof(dptype)*ram->ngridmax*ram->ndim,PPTR(mesh->xg));
	mesh->father = (int*)Malloc(sizeof(int)*ram->ngridmax,PPTR(mesh->father));
	mesh->nbor = (int*)Malloc(sizeof(int)*ram->ngridmax * twondim,PPTR(mesh->nbor));
	mesh->next = (int*)Malloc(sizeof(int)*ram->ngridmax,PPTR(mesh->next));
	mesh->prev = (int*)Malloc(sizeof(int)*ram->ngridmax,PPTR(mesh->prev));



	{
		int _nchip;
		fread(&_nchip, sizeof(int), 1, fp);
		fread(&(mesh->headf), sizeof(int), 1, fp);
		fread(&(mesh->tailf), sizeof(int), 1, fp);
		fread(&(mesh->numbf), sizeof(int), 1, fp);
		fread(&(ram->used_mem), sizeof(int), 1, fp);
		fread(&(ram->used_mem_tot), sizeof(int), 1, fp);
		fread(&_nchip, sizeof(int), 1, fp);
	}

	mesh->headf = ram->ngrid_current;
	mesh->tailf = ram->ngridmax-1;
	mesh->numbf = ram->ngridmax-ram->ngrid_current;
	mesh->prev[mesh->headf-1] = 0;
	mesh->next[mesh->tailf-1] = 0;


	char ordering2[129];
	F77read((ordering2),sizeof(char), 128, fp);
	if(strncmp(ordering2, ram->ordering, 5)!=0){
		ERRORPRINT("Ordering is uncompatible %s :: %s\n", ordering2, ram->ordering);
		exit(999);
	}


	mesh->bound_key = (qdptype*)Malloc(sizeof(qdptype)*(ram->ndomain+1),PPTR(mesh->bound_key));
	if(ram->nrestart_quad == ram->nrestart){
		mesh->bound_key_restart = (dptype*)Malloc(sizeof(dptype)*(ram->ndomain+1),PPTR(mesh->bound_key));
	}


	char testaa[190];
	sprintf(testaa,"bisection");
	if(strncmp(ram->ordering, testaa, 5)==0){
		long offset = sizeof(int)*5*2 + 
			sizeof(dptype)*ram->nbinodes + 
			sizeof(int)*ram->nbinodes*2 + 
			sizeof(int)*ram->nbinodes + 
			sizeof(dptype)*ram->ndim*ram->ncpu*2;
		fseek(fp, offset, SEEK_CUR);
	}
	else {
#ifdef QUADHILBERT
		if(ram->nrestart_quad == ram->nrestart){
			F77read(mesh->bound_key_restart, sizeof(dptype), (ram->ndomain+1), fp);
			for(i=0;i<ram->ndomain+1;i++) mesh->bound_key[i] = mesh->bound_key_restart[i];
			Free(mesh->bound_key_restart);
		}
		else{
			F77read(mesh->bound_key, sizeof(qdptype), (ram->ndomain+1), fp);
		}
#else 
		F77read(mesh->bound_key, sizeof(qdptype), (ram->ndomain+1), fp);
#endif
	}
	F77read((mesh->son), sizeof(int), ram->ncoarse, fp);
	F77read((mesh->flag1), sizeof(int), ram->ncoarse, fp);
	/*
	F77read((mesh->cpu_map), sizeof(int), ram->ncoarse, fp);
	*/
	{
		int chip;
		fread(&chip, sizeof(int), 1,fp);
		fread(mesh->cpu_map, sizeof(int), ram->ncoarse, fp);
		fread(&chip, sizeof(int), 1,fp);

	}

	for(ilevel=0;ilevel<ram->nlevelmax;ilevel++){
		for(ibound=0;ibound<ram->nboundary+ram->ncpu;ibound++){
			int ncache;
			if(ibound < ram->ncpu){
				ncache = mesh->numbl[ibound + ram->ncpu*ilevel];
			}
			else {
				ncache = mesh->numbb[ibound - ram->ncpu + ram->nboundary*ilevel];
			}
			if(ncache >0 ) {
				int *ind_grid = (int*)Malloc(sizeof(int)*ncache, PPTR(ind_grid));
				int *pos = (int*)Malloc(sizeof(int)*ncache, PPTR(pos));
				int *grid = (int*)Malloc(sizeof(int)*ncache, PPTR(grid));
				dptype *xdp = (dptype*)Malloc(sizeof(dptype)*ncache, PPTR(xdp));
				dptype *xxg = (dptype*)Malloc(sizeof(dptype)*ncache, PPTR(xxg));
				int *iig = (int*)Malloc(sizeof(int)*ncache, PPTR(iig));
				F77read(ind_grid, sizeof(int), ncache, fp);
//				for(i=0;i<ncache;i++) ind_grid[i] = ind_grid[i];
				F77read(iig, sizeof(int), ncache, fp);
				for(i=0;i<ncache;i++) mesh->next[ind_grid[i]-1] = iig[i];
				F77read(iig, sizeof(int), ncache, fp);
				for(i=0;i<ncache;i++) mesh->prev[ind_grid[i]-1] = iig[i];
				for(i=0;i<ncache;i++) mesh->ilevel[ind_grid[i]-1] = ilevel+1;
				int idim;
				for(idim=0;idim<ram->ndim;idim++) {
					F77read(xxg, sizeof(dptype), ncache, fp);
					for(i=0;i<ncache;i++) mesh->xg[ind_grid[i]-1+ram->ngridmax*idim] = xxg[i];
				}
				F77read(iig, sizeof(int), ncache, fp);
				if(ram->ngridmax != ngridmax2 && ilevel > 0) {
					for(i=0;i<ncache;i++) {
						pos[i] = (iig[i]-ram->ncoarse-1)/ngridmax2;
						grid[i] = (iig[i]-ram->ncoarse-pos[i])/ngridmax2;
						iig[i] = ram->ncoarse+pos[i]*ram->ngridmax + grid[i];
					}
				}
				for(i=0;i<ncache;i++) mesh->father[ind_grid[i]-1] = iig[i];
				int ind;
				for(ind=0;ind<twondim;ind++){
					F77read(iig, sizeof(int), ncache, fp);
					if(ram->ngridmax != ngridmax2 && ilevel >0){
						for(i=0;i<ncache;i++) {
							pos[i] = (iig[i]-ram->ncoarse-1)/ngridmax2;
							grid[i] = (iig[i]-ram->ncoarse-pos[i])/ngridmax2;
							iig[i] = ram->ncoarse+pos[i]*ram->ngridmax + grid[i];
						}
					}
					for(i=0;i<ncache;i++) mesh->nbor[ind_grid[i]-1 + ram->ngridmax*ind] = iig[i];
				}
				for(ind=0;ind<ram->twotondim;ind++){
					int iskip = ram->ncoarse + ind*ram->ngridmax;
					F77read(iig, sizeof(int), ncache, fp);
					for(i=0;i<ncache;i++) {
						mesh->son[ind_grid[i]-1 + iskip] = iig[i];
					}
				}
				for(ind=0;ind<ram->twotondim;ind++){
					int iskip = ram->ncoarse + ind*ram->ngridmax;
					F77read(iig, sizeof(int), ncache, fp);
					for(i=0;i<ncache;i++) mesh->cpu_map[ind_grid[i]-1 + iskip] = iig[i];
				}
				for(ind=0;ind<ram->twotondim;ind++){
					int iskip = ram->ncoarse + ind*ram->ngridmax;
					F77read(iig, sizeof(int), ncache, fp);
					for(i=0;i<ncache;i++) mesh->flag1[ind_grid[i]-1 + iskip] = iig[i];
				}
				Free(iig);
				Free(xxg);
				Free(xdp);
				Free(grid);
				Free(pos);
				Free(ind_grid);
			}
		}
	}

#endif	
	fclose(fp);
	units(ram);

	return 0;
}
