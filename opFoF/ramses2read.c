/* 
 * icc -o how2read how2read.c -L -lmyram -g -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL
 *
 * To Run,  ./how2read [THE FILENAME YOU WANT TO READ]  [STAR/SINK/HCELL/DM]
 * */

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<sys/stat.h>
#include<math.h>
#include<mpi.h>
#include "ramses.h"
#include "fof.h"
#include "params.h"
#include "pmheader.h"

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))


/*
inline idtype GetHCellID(HydroCellType *hcell, idtype nx){
	idtype id;
	idtype ix,iy,iz;
	dptype x,y,z;

	x = (hcell->x)/simpar.boxlen_ini;
	y = (hcell->y)/simpar.boxlen_ini;
	z = (hcell->z)/simpar.boxlen_ini;
	ix = x * nx;
	iy = y * nx;
	iz = z * nx;
	id = ix + nx * (iy + nx*iz);
	id += IDOFFSET;
	return id;
}
*/


dptype getlink02(dptype mass){
	dptype link02;
	link02 = 0.2*pow(mass/2.7755e11L/simpar.omega_m, 0.33333333333333333333L);
	return link02;
}

size_t read_ramses_data(FoFTPtlStruct **Bp, size_t np, char *file1, char *type, double xoffset){
	size_t i,j,k;

	DmType *dm;
	StarType *star;
	SinkType *sink;
	HydroCellType *hcell;
	GasType *gas;
	struct stat st;
	size_t size = 0,mp;

	int ierr = stat(file1, &st);
	if(ierr == 0) size = st.st_size;
	size_t dsize;
	void *aa = NULL;

	if(strcmp(type,"STAR")==0){
		dsize = sizeof(StarType);
		star = NULL;
	}
	else if(strcmp(type,"SINK")==0){
		dsize = sizeof(SinkType);
		sink = NULL;
	}
	else if(strcmp(type,"GAS")==0){
		dsize = sizeof(GasType);
		gas = NULL;
	}
	else if(strcmp(type,"DM")==0){
		dsize = sizeof(DmType);
		dm = NULL;
	}
	else {
		fprintf(stderr,"Oooooooooops. Wrong size in file and type\n");
		MPI_Abort(MPI_COMM_WORLD, 98);
		return 0;
	}
	if(ierr == 0){
		aa = malloc(size > 0 ? size : 1);
		if(aa == NULL){
			fprintf(stderr,"Error allocating %ld bytes for %s\n",size,file1);
			MPI_Abort(MPI_COMM_WORLD, 97);
			return 0;
		}
		if(strcmp(type,"STAR")==0) star = (StarType*)aa;
		else if(strcmp(type,"SINK")==0) sink = (SinkType*)aa;
		else if(strcmp(type,"GAS")==0) gas = (GasType*)aa;
		else dm = (DmType*)aa;
	}

	if(ierr == 0 && size%dsize !=0){
		fprintf(stderr,"File size is not a multiple of the %s record ABI: %s\n",
				type,file1);
		MPI_Abort(MPI_COMM_WORLD, 95);
		return 0;
	}

	int nid,myid;
	MPI_Status status;
	int itag=1;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	int isend,iget;
	isend = iget = 1;
	int WGroupSize = WGROUPSIZE;
	int src = myid -1;
	int tgt = myid + 1;

	if(RANKINGROUP(myid,WGroupSize) !=0) MPI_Recv(&i,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	if(ierr ==0){
		FILE *fp = NULL;
		fp = fopen(file1,"rb");
		if(fp == NULL){
			fprintf(stderr,"Error opening %s\n",file1);
			MPI_Abort(MPI_COMM_WORLD, 94);
			return 0;
		}
		if(fread(aa, sizeof(char), size, fp) != size){
			fprintf(stderr,"Short read from %s\n",file1);
			MPI_Abort(MPI_COMM_WORLD, 96);
		}
		mp = size/dsize;
		fclose(fp);
	}
	else {
		mp = 0;
	}
	printf("P%d has read %s with np= %ld\n", myid, file1, mp);fflush(stdout);
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) MPI_Send(&i,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);

	FoFTPtlStruct *p, *bp;
	bp = *Bp;
	bp = (FoFTPtlStruct *)realloc(bp,sizeof(FoFTPtlStruct)*(np+mp));
	*Bp = bp;
	p = bp + np;

	idtype nx;
	nx = 1;
	for(i=0;i<simpar.nlevelmax;i++) nx = nx * 2;

	dptype yzmin,yzmax;
	yzmin = 0.445312;
	yzmax = 0.554688;
	yzmin = 0L;
	yzmax = 1L;
	yzmin = yzmin * simpar.boxlen_ini;
	yzmax = yzmax * simpar.boxlen_ini;

	if(strcmp(type,"STAR")==0){
		for(i=0;i<mp;i++){
			if(star[i].y >= yzmin && star[i].y < yzmax && star[i].z>=yzmin && star[i].z<yzmax){
				star[i].x += xoffset;
				p->type = TYPE_STAR;
				p->z = star[i].x;
				p->y = star[i].y;
				p->x = star[i].z;
				/*
				p->rv[2] = star[i].vx;
				p->rv[1] = star[i].vy;
				p->rv[0] = star[i].vz;
				p->indx = star[i].id;
				p->mass = star[i].mass;
				*/
				p->p.star = star[i];

				p->link02 = getlink02(star[i].mass);
				p++;
			}
		}
	}
	else if(strcmp(type,"SINK")==0){
		for(i=0;i<mp;i++){
			if(sink[i].y >= yzmin && sink[i].y < yzmax && sink[i].z>=yzmin && sink[i].z<yzmax){
				sink[i].x += xoffset;
				p->type = TYPE_SINK;
				p->z = sink[i].x;
				p->y = sink[i].y;
				p->x = sink[i].z;
				/*
				p->rv[2] = sink[i].vx;
				p->rv[1] = sink[i].vy;
				p->rv[0] = sink[i].vz;
				p->mass = sink[i].mass;
				p->indx = -sink[i].id;
				*/
				p->p.sink = sink[i];
				p->link02 = getlink02(sink[i].mass);
				p++;
			}
		}
	}
	else if(strcmp(type,"GAS")==0){
		for(i=0;i<mp;i++){
			if(gas[i].y >= yzmin && gas[i].y < yzmax && gas[i].z>=yzmin && gas[i].z<yzmax){
				gas[i].x += xoffset;
				p->type = TYPE_GAS;
				p->z = gas[i].x;
				p->y = gas[i].y;
				p->x = gas[i].z;
				/*
				p->rv[2] = gas[i].vx;
				p->rv[1] = gas[i].vy;
				p->rv[0] = gas[i].vz;
				p->mass = gas[i].mass;
				p->indx = gas[i].indx;
				*/
				p->p.gas = gas[i];
				p->link02 = getlink02(gas[i].mass);
				p++;
			}
		}
	}
	else if(strcmp(type,"DM")==0){
		for(i=0;i<mp;i++){
			if(dm[i].y >= yzmin && dm[i].y < yzmax && dm[i].z>=yzmin && dm[i].z<yzmax){
				dm[i].x += xoffset;
				p->type = TYPE_DM;
				p->z = dm[i].x;
				p->y = dm[i].y;
				p->x = dm[i].z;
				/*
				p->rv[2] = dm[i].vx;
				p->rv[1] = dm[i].vy;
				p->rv[0] = dm[i].vz;
				p->mass = dm[i].mass;
				p->indx = dm[i].id;
				*/
				p->p.dm = dm[i];
				p->link02 = getlink02(dm[i].mass);
				p++;
			}
		}
	}
	free(aa);
	dptype mass = 1.E11;
	printf("P has example link02= %g Mpc/h for a mass of %g\n", getlink02(mass),mass);
	mp = p-(bp+np);
	*Bp = (FoFTPtlStruct *)realloc(*Bp,sizeof(FoFTPtlStruct)*(np+mp));
	return mp;
}
