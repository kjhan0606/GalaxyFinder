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
#include "fof.h"
#include "ramses.h"
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


inline dptype getlink02(dptype mass){
	dptype link02;
	link02 = 0.2*pow(mass/2.7755e11L/simpar.omega_m, 0.33333333333333333333L);
	return link02;
}

size_t read_ramses_data(FoFTPtlStruct **Bp, size_t np, char *file1, char *type){
	size_t i,j,k;

	DmType *dm;
	StarType *star;
	SinkType *sink;
	HydroCellType *hcell;
	GasType *gas;
	struct stat st;
	size_t size,mp;

	stat(file1, &st);
	size = st.st_size;
	size_t dsize;
	void *aa;

	if(strcmp(type,"STAR")==0){
		dsize = sizeof(StarType);
		aa = (void*) star = (StarType*)malloc(sizeof(char)*size);
	}
	else if(strcmp(type,"SINK")==0){
		dsize = sizeof(SinkType);
		aa = (void*)sink = (SinkType*)malloc(sizeof(char)*size);
	}
	else if(strcmp(type,"GAS")==0){
		dsize = sizeof(GasType);
		aa = (void*)gas = (GasType*)malloc(sizeof(char)*size);
	}
	else if(strcmp(type,"DM")==0){
		dsize = sizeof(DmType);
		aa = (void*)dm = (DmType*)malloc(sizeof(char)*size);
	}
	else {
		fprintf(stderr,"Oooooooooops. Wrong size in file and type\n");
	}

	if(size%dsize !=0){
		printf("Error in file size,,,, %s\n", file1);
		exit(99);
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
	FILE *fp = fopen(file1,"r");
	fread(aa, sizeof(char), size, fp);
	mp = size/dsize;
	fclose(fp);
	printf("P%d has read %s with np= %ld\n", myid, file1, mp);fflush(stdout);
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) MPI_Send(&i,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);

	FoFTPtlStruct *p, *bp;
	bp = *Bp;
	bp = (FoFTPtlStruct *)realloc(bp,sizeof(FoFTPtlStruct)*(np+mp));
	*Bp = bp;
	p = bp + np;
	dptype Mpc = 3.0857E24;
	dptype Msun = 2E33;

	idtype nx;
	nx = 1;
	for(i=0;i<simpar.nlevelmax;i++) nx = nx * 2;

	dptype yzmin,yzmax;
	yzmin = 0.445312;
	yzmax = 0.554688;
	yzmin = yzmin * simpar.boxlen_ini;
	yzmax = yzmax * simpar.boxlen_ini;

	if(strcmp(type,"STAR")==0){
		for(i=0;i<mp;i++){
			if(star[i].y >= yzmin && star[i].y < yzmax && star[i].z>=yzmin && star[i].z<yzmax){
				p->type = TYPE_STAR;
				p->r[2] = star[i].x;
				p->r[1] = star[i].y;
				p->r[0] = star[i].z;
				p->indx = star[i].id;
				p->mass = star[i].mass;
				p->link02 = getlink02(star[i].mass);
				p++;
			}
		}
	}
	else if(strcmp(type,"SINK")==0){
		for(i=0;i<mp;i++){
			if(sink[i].y >= yzmin && sink[i].y < yzmax && sink[i].z>=yzmin && sink[i].z<yzmax){
				p->type = TYPE_SINK;
				p->r[2] = sink[i].x;
				p->r[1] = sink[i].y;
				p->r[0] = sink[i].z;
				p->mass = sink[i].mass;
				p->indx = -sink[i].id;
				p->link02 = getlink02(sink[i].mass);
				p++;
			}
		}
	}
	else if(strcmp(type,"GAS")==0){
		for(i=0;i<mp;i++){
			if(gas[i].y >= yzmin && gas[i].y < yzmax && gas[i].z>=yzmin && gas[i].z<yzmax){
				p->type = TYPE_GAS;
				p->r[2] = gas[i].x;
				p->r[1] = gas[i].y;
				p->r[0] = gas[i].z;
				p->mass = gas[i].mass;
				p->indx = gas[i].indx;
				p->link02 = getlink02(gas[i].mass);
				p++;
			}
		}
	}
	else if(strcmp(type,"DM")==0){
		for(i=0;i<mp;i++){
			if(dm[i].y >= yzmin && dm[i].y < yzmax && dm[i].z>=yzmin && dm[i].z<yzmax){
				p->type = TYPE_DM;
				p->r[2] = dm[i].x;
				p->r[1] = dm[i].y;
				p->r[0] = dm[i].z;
				p->mass = dm[i].mass;
				p->indx = dm[i].id;
				p->link02 = getlink02(dm[i].mass);
				p++;
			}
		}
	}
	dptype mass = 1.E11;
	printf("P has example link02= %g Mpc/h for a mass of %g\n", getlink02(mass),mass);
	free(aa);
	mp = p-(bp+np);
	*Bp = (FoFTPtlStruct *)realloc(*Bp,sizeof(FoFTPtlStruct)*(np+mp));
	return mp;
}
