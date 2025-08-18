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

#define MAX(a,b) ( (a) >(b) ? (a):(b) )
#define MIN(a,b) ( (a) <(b) ? (a):(b) )



int rd_hydro(RamsesType *ram, char *infile){
	FILE *fp = fopen(infile,"r");
	int ilevel,ibound;
	int chip,i;
	int npartp,mstar,ng;
	int *next;
	size_t ncell;
	HydroType *hydro;
	MeshType *mesh;
	int ncpu2, nvar2, ndim2, nlevelmax2, nboundary2;
	int nvar, nener;
	int twotondim;
	int ngridmax, ncoarse,ndim;
	dptype gamma2,gamma,smallr;
	dptype *uold, *unew;
	int rt, neq_chem;

	twotondim = ram->twotondim;
	ngridmax = ram->ngridmax;
	ncoarse = ram->ncoarse;
	ndim = ram->ndim;
	nener = ram->nener;
//	nvar = ram->nvar;
	gamma = ram->gamma;
	rt = ram->rt;
	neq_chem = ram->neq_chem;
	smallr = ram->smallr;


	mesh = &(ram->mesh);
	hydro = &(ram->hydro);
	next = mesh->next;

	ncell = ram->ncoarse + ram->twotondim*ram->ngridmax;
	/*
	unew = hydro->unew = (dptype*)Malloc(sizeof(dptype)*ncell*nvar);
	*/


//	if(ram->nrestart >0){
		F77read(&ncpu2, sizeof(int), 1, fp);
		F77read(&nvar2, sizeof(int), 1, fp);
		F77read(&ndim2, sizeof(int), 1, fp);
		F77read(&nlevelmax2, sizeof(int), 1, fp);
		F77read(&nboundary2, sizeof(int), 1, fp);
		F77read(&gamma2, sizeof(dptype), 1, fp);
	
	ram->nvar=nvar2;
	nvar=nvar2;

	uold = hydro->uold = (dptype*)Malloc(sizeof(dptype)*ncell*nvar,PPTR(hydro->uold));

//	}
	if( !(neq_chem || rt) && nvar2 != ram->nvar){
		DEBUGPRINT0("FILE hydro.tmp is not compatible\n");
		DEBUGPRINT("Found     = %d\n", nvar2);
		DEBUGPRINT("Expected  = %d\n", ram->nvar);
		exit(999);
	}
	
	for(ilevel=0;ilevel<nlevelmax2;ilevel++){
		for(ibound=0;ibound<ram->nboundary+ram->ncpu;ibound++){
			int ncache,istart;
			if(ibound < ram->ncpu){
				ncache = mesh->numbl[ibound + ram->ncpu*ilevel];
				istart = mesh->headl[ibound + ram->ncpu*ilevel];
			}
			else {
				ncache = mesh->numbb[ibound - ram->ncpu + ram->nboundary*ilevel];
				istart = mesh->headb[ibound - ram->ncpu + ram->nboundary*ilevel];
			}
			int ilevel2,numbl2;
			F77read(&ilevel2, sizeof(int), 1, fp);
			F77read(&numbl2, sizeof(int), 1, fp);
			if(numbl2 != ncache){
				DEBUGPRINT0("File hydro.tmp is not compatible\n");
				DEBUGPRINT("Found    = %d for level %d \n", numbl2, ilevel2);
				DEBUGPRINT("Expected = %d for level %d \n", ncache, ilevel+1);
			}
			if(ncache >0 ) {
				int *ind_grid = (int*)Malloc(sizeof(int)*ncache,PPTR(ind_grid));
				dptype *xx = (dptype*)Malloc(sizeof(dptype)*ncache,PPTR(xx));

				int igrid = istart;
				for(i=0;i<ncache;i++){
					ind_grid[i] = igrid;
					igrid = next[igrid-1];
				}
				int ind;
				for(ind=0;ind<twotondim;ind++){
					size_t iskip = ncoarse + ind*ngridmax;
					int ivar;
					for(ivar=1;ivar<=ndim+1;ivar++){
						F77read(xx, sizeof(dptype), ncache, fp);
						/* Density */
						if(ivar == 1) {
							for(i=0;i<ncache;i++){
								uold[ind_grid[i]-1+iskip] = xx[i];
							}
						}
						/* Velocity */
						else if(ivar >= 2 && ivar <= ndim+1){
							for(i=0;i<ncache;i++){
								uold[ind_grid[i]-1+iskip+ncell*(ivar-1)] = 	xx[i]; //* MAX(uold[ind_grid[i]-1+iskip],smallr);
							}
						}
					}
#if NENER >0
					for(ivar=ndim+3;ivar<=ndim+2+nener; ivar++){
						F77read(xx, sizeof(dptype), ncache, fp);
						for(i=0;i<ncache;i++){
							uold[ind_grid[i]-1+iskip+ncell*(ivar-1)] = xx[i];//(ram->gamma_rad[ivar-1-ndim-2]-1.L)
						}
					}
#endif
					F77read(xx, sizeof(dptype), ncache, fp);

					for(i=0;i<ncache;i++){
/*
						xx[i] = xx[i]/(gamma-1.);
						if(uold[ind_grid[i]-1+iskip] > 0) {
							xx[i] = xx[i]+0.5L*pow(uold[ind_grid[i]-1+iskip + ncell*1],2)/
								MAX(uold[ind_grid[i]-1+iskip],smallr);
#if NDIM>1
							xx[i] = xx[i]+0.5L*pow(uold[ind_grid[i]-1+iskip + ncell*2],2)/
								MAX(uold[ind_grid[i]-1+iskip],smallr);
#endif
#if NDIM>2
							xx[i] = xx[i]+0.5L*pow(uold[ind_grid[i]-1+iskip + ncell*3],2)/
								MAX(uold[ind_grid[i]-1+iskip],smallr);
#endif
#if NENER>0
							for(irad=0;irad<nener;irad++){
								xx[i] = xx[i] + uold[ind_grid[i]-1+iskip + ncell*(ndim+1+irad)];
							}
#endif
						}
						else{
							xx[i] = 0;
						}
						// total fluid energy */

						uold[ind_grid[i]-1+iskip+ncell*(ndim+1)] = xx[i];
					}

//#if NVAR>NDIM+2+NENER
					for(ivar=ndim+3+nener;ivar<=nvar;ivar++){
						F77read(xx, sizeof(dptype), ncache, fp);
						/* passive scalars */
						for(i=0;i<ncache;i++){
							uold[ind_grid[i]-1+iskip + ncell*(ivar-1)] = xx[i]; //* MAX(uold[ind_grid[i]-1+iskip], smallr);
						}
					}
//#endif
				}
				Free(xx);
				Free(ind_grid); 
			}
		}
	}
	fclose(fp);
}
 
