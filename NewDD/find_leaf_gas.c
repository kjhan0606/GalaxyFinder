#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "ramses.h"
#include "Memory.h"

#define MAX(a,b) ( (a)>(b) ? (a): (b) )



int find_leaf_gas(RamsesType *ram, int icpu, char *outfile){

	int i,j,k,ig,nchem,ndust;

	MeshType *mesh = &(ram->mesh);
	int *numbl = mesh->numbl;
	int *headl = mesh->headl;
	HydroType *hydro = &(ram->hydro);
	int ldim, ncpu, nlevelmax;
	size_t ncell;
	int ngridmax;
	dptype dx;
	int ndim;
	int ncoarse = ram->ncoarse;
	int *next = (mesh->next);

	ndim = ram->ndim;
	ncpu = ram->ncpu;
	nlevelmax = ram->nlevelmax;
	ngridmax = ram->ngridmax;
	dptype boxlen_ini = ram->boxlen_ini;
	dptype hubin = ram->H0 / 100.L;
	int nboundary = ram->nboundary;
	int ilevel,istart;
	int *numbb = mesh->numbb;
	int *headb = mesh->headb;
	int igrid, ind;
	int twotondim = ram->twotondim;
	size_t iskip;
	int *son = mesh->son;
	float CPI[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
	float *cpi;

	ncell = ncoarse + twotondim * ngridmax;


	dptype Lbox = boxlen_ini;
	dptype Lbox05 = Lbox/2;
	dptype hh = ram->H0/100.L;
	dptype smallr = ram->smallr;
	int imetal = ram->imetal;
	int ichem = ram->ichem;
	int idust = ram->idust;

	dptype *uold = hydro->uold;
//	dptype *phi = hydro->phi;
//	dptype *fxyz = hydro->fxyz;




	int ntot = 0;
	int ncache = 512;
//	int *ind_leaf = (int*)Malloc(sizeof(int)*ncache, PPTR(ind_leaf));
	int ileaf=0;


	int nleaf = 0;
	int ibound;
	for(ibound=1;ibound <= nboundary+ncpu;ibound++){
		if(ibound == icpu){
			for(ilevel=1;ilevel<=nlevelmax;ilevel++){
				if(ibound <=ncpu) {
					ncache = numbl[ibound -1 + (ilevel-1)*ncpu]; 
					istart = headl[ibound -1 + (ilevel-1)*ncpu]; 
				}
				else {
					ncache = numbb[ibound -1 + (ilevel-1)*nboundary]; 
					istart = headb[ibound -1 + (ilevel-1)*nboundary]; 
				}
				if(ncache>0) {
					int *ind_grid = (int*)Malloc(sizeof(int)*ncache, PPTR(ind_grid));
					igrid = istart;
					for(i=0;i<ncache;i++){
						ind_grid[i] = igrid;
						igrid = next[igrid-1];
					}
					for(ind=0;ind<twotondim;ind++){
						iskip = ncoarse + ind*ngridmax;
						for(i=0;i<ncache;i++){
							if(son[ind_grid[i]-1+iskip] ==0){
								nleaf ++;
							}
						}
					}
					Free(ind_grid);
					ntot+=ncache;
				}
			}
		}
	}
	printf("Total nleaf= %d\n", nleaf);
	GasType *gas;
	gas = (ram->gas = (GasType*)Malloc(sizeof(GasType)*nleaf, PPTR(ram->gas)));

#ifdef NCHEM
	nchem=NCHEM;
//	gas->chem = (float*)Malloc(sizeof(float)*NCHEM,PPTR(gas->chem));
#endif
#ifdef NDUST
	ndust=NDUST;
//	gas->chem = (float*)Malloc(sizeof(float)*NCHEM,PPTR(gas->chem));
#endif

	ileaf = 0;
	for(ibound=1;ibound <= nboundary+ncpu;ibound++){
		if(ibound == icpu){
			for(ilevel=1;ilevel<=nlevelmax;ilevel++){
				dptype dx = pow(0.5, ilevel);
				dptype cellsize = dx*Lbox; 
				dptype volcell = pow(cellsize*(Mpc/hh) * ram->aexp,3.L);
				dptype den2massilevel = ram->scale_d *volcell/(Msun/hh);
				if(ibound <=ncpu) {
					ncache = numbl[ibound -1 + (ilevel-1)*ncpu]; 
					istart = headl[ibound -1 + (ilevel-1)*ncpu]; 
				}
				else {
					ncache = numbb[ibound -1 -ram->ncpu + (ilevel-1)*nboundary]; 
					istart = headb[ibound -1 -ram->ncpu + (ilevel-1)*nboundary]; 
				}
				if(ncache>0) {
					int *ind_grid = (int*)Malloc(sizeof(int)*ncache, PPTR(ind_grid));
					igrid = istart;
					nleaf = 0;
					for(i=0;i<ncache;i++){
						ind_grid[i] = igrid;
						igrid = next[igrid-1];
					}
					for(ind=0;ind<twotondim;ind++){
						iskip = ncoarse + ind*ngridmax;
						for(i=0;i<ncache;i++){
							if(son[ind_grid[i]-1+iskip] ==0){
								nleaf ++;
							}
						}
					}
					if(nleaf >0) {
						j = 0;
						for(ind=0;ind<twotondim;ind++){
							iskip = ncoarse + ind*ngridmax;
							cpi = CPI[ind];
							for(i=0;i<ncache;i++){
								if(son[ind_grid[i]-1+iskip]==0){
									gas[ileaf].den = uold[ind_grid[i]-1+iskip];
									gas[ileaf].x  = (mesh->xg[ind_grid[i]-1+ngridmax*0] + (cpi[0]-0.5)*dx)* Lbox;
									gas[ileaf].y  = (mesh->xg[ind_grid[i]-1+ngridmax*1] + (cpi[1]-0.5)*dx)* Lbox;
									gas[ileaf].z  = (mesh->xg[ind_grid[i]-1+ngridmax*2] + (cpi[2]-0.5)*dx)* Lbox;
									gas[ileaf].vx = uold[ind_grid[i]-1+iskip+1*ncell];//MAX(gas[ileaf].den,smallr);
									gas[ileaf].vy = uold[ind_grid[i]-1+iskip+2*ncell];//MAX(gas[ileaf].den,smallr);
									gas[ileaf].vz = uold[ind_grid[i]-1+iskip+3*ncell];//MAX(gas[ileaf].den,smallr);
//									gas[ileaf].potent = phi[ind_grid[i]-1+iskip];
//									gas[ileaf].fx = fxyz[ind_grid[i]-1+iskip+0*ncell];
//									gas[ileaf].fy = fxyz[ind_grid[i]-1+iskip+1*ncell];
//									gas[ileaf].fz = fxyz[ind_grid[i]-1+iskip+2*ncell];

									gas[ileaf].cellsize = cellsize; /* in cMpc/h */
									gas[ileaf].metallicity = uold[ind_grid[i]-1+iskip + ncell*imetal];//gas[ileaf].den;i
#ifdef NCHEM
									for (k=0;k<nchem;k++){
										gas[ileaf].chem[k]=uold[ind_grid[i]-1+iskip + ncell*(ichem+k)];
								//		printf("Mfrac chem= %g\n", gas[ileaf].chem[k]);
									}
#endif
#ifdef NDUST
									for (k=0;k<ndust;k++){
										gas[ileaf].dust[k]=uold[ind_grid[i]-1+iskip + ncell*(idust+k)];
								//		printf("Mfrac ion= %g\n", gas[ileaf].ion[k]);
									}
#endif


			//						dptype ekk = 0;
				//					ekk += 0.5*gas[ileaf].den*gas[ileaf].vx*gas[ileaf].vx;
				//					ekk += 0.5*gas[ileaf].den*gas[ileaf].vy*gas[ileaf].vy;
				//					ekk += 0.5*gas[ileaf].den*gas[ileaf].vz*gas[ileaf].vz;

//									dptype Temp = uold[ind_grid[i]-1+iskip+(ndim+1)*ncell];
									dptype err = 0;
/*
#if NENER>0
									int irad;
									for(irad=0;irad<nener;irad++){ 
										err += uold[ind_grid[i]-1 + iskip + ncell*(inener+irad)]; 
									}
#endif
									dptype emag = 0;
#ifdef SOLVERmhd 
									for(idim=0;idim<ndim;idim++){ 
										emag += 0.125L*TWICE(uold[ind_grid[i]-1+iskip+ncell*(idim+ndim+2)] + 
												uold[ind_grid[i]-1+iskip+ncell*(idim+nvar)]); 
									}
#endif
*/
									gas[ileaf].temp = (uold[ind_grid[i]-1+iskip+ncell*(ndim+1)])/gas[ileaf].den;
//									gas[ileaf].temp = (ram->gamma-1)*(uold[ind_grid[i]-1+iskip+ncell*(ndim+1)]-ekk-err-emag);
									gas[ileaf].temp *= ram->scale_T2; /* in K/mu */
									gas[ileaf].vx *= ram->kmscale_v; /* in km/sec */
									gas[ileaf].vy *= ram->kmscale_v; /* in km/sec */
									gas[ileaf].vz *= ram->kmscale_v; /* in km/sec */

						//			gas[ileaf].H  = uold[ind_grid[i]-1+iskip+6*ncell]/MAX(gas[ileaf].den,smallr);
						//			gas[ileaf].O  = uold[ind_grid[i]-1+iskip+7*ncell]/MAX(gas[ileaf].den,smallr);
						//			gas[ileaf].Fe = uold[ind_grid[i]-1+iskip+8*ncell]/MAX(gas[ileaf].den,smallr);

						//			gas[ileaf].H  = 1-Y-gas[ileaf].metallicity;

						//			gas[ileaf].ilevel = ilevel;
									gas[ileaf].mass = gas[ileaf].den * den2massilevel;
						//			gas[ileaf].indx = ncell*(icpu-1) + ind_grid[i]+iskip;

									ileaf++;
								}
							}
						}
					}
					Free(ind_grid);
				}
			}
		}
	}
	ram->nleafcell = ileaf;
	ram->ngas = ileaf;
	printf("Total number of leaf cells are %d from %d\n", ileaf,ntot);
	if(0){
		FILE *wp = fopen(outfile,"w");
		for(i=0;i<ileaf;i++){
			/*
			if(hcell[i].metallicity > 1.E-5 && hcell[i].temp < 1.E3)
			*/
			if(gas[i].temp < 1.E3)
			{
				fprintf(wp, "%d %g %g \n",i, gas[i].metallicity, gas[i].temp);
		//		fprintf(wp, "%d %g %g %g %g %g\n",i, gas[i].metallicity, gas[i].temp, gas[i].H, gas[i].O, gas[i].Fe);
			}
		}
		fclose(wp);
	}
}
