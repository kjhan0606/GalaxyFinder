#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "ramses.h"
#include "params.h"
void write_head(FILE *wp, RamsesType wsimpar){
	int isize,ncnt = 0;
	FILE_HEADER(fprintf,wp,,wsimpar);
}
RamsesType read_head(FILE *fp){
	RamsesType rsimpar;
	int isize,ncnt=0;
	char line[MAX_LINE_LENGTH];
	/* Loop until it reaches P_Closing */
	while(strcmp(P_Closing,fgets(line, MAX_LINE_LENGTH,fp)) != 0){
		int ncnt=0;
		/* check whether it is a comment */
		if(line[0] != '#'){
			if(strstr(line,"define") != NULL) FILE_HEADER(sscanf,line,&,rsimpar)
			if(strstr(line,"=") != NULL && strstr(line,"define")!=NULL && ncnt ==0) {
				fprintf(stderr,"Warning: the following parameter line is unknown:\n %s\n\n\n",line);
			}
		}
	}
	/* Now change the stepnum to this step */
	return rsimpar;
}
void write_default_sim_parameter_file(FILE *wp, RamsesType wsimpar){
	int isize,ncnt = 0;
	DEFAULT_PARAMS(fprintf,wp,,wsimpar);
}

void write_sim_parameter_file(FILE *wp, RamsesType wsimpar){
	int isize,ncnt = 0;
	PARAMS(fprintf,wp,,wsimpar);
}
RamsesType  read_sim_parameter_file(FILE *fp){
	RamsesType rsimpar;
	int isize,ncnt=0;
	char line[MAX_LINE_LENGTH];
	while(strcmp(P_Closing,fgets(line, MAX_LINE_LENGTH,fp)) != 0){ /* Loop until it reaches P_Closing */
		int ncnt=0;
		if(line[0] != '#'){ /* check whether it is a comment */
			if(strstr(line,"define") != NULL) PARAMS(sscanf,line,&,rsimpar)
			if(strstr(line,"=") != NULL && strstr(line,"define")!=NULL && ncnt ==0) {
				fprintf(stderr,"Warning: the following parameter line is unknown:\n %s\n\n\n",line);
			}
		}
	}
	return rsimpar;
}
void mk_default_param(RamsesType *defsim, char *cosmology){
	double zi;
	/* Simulation Parameters for WMAP 5-year */
	if(strcmp(cosmology,"WMAP3")==0){
		defsim->omega_m = 0.238;
		defsim->omega_l = 0.762;
		defsim->omega_b = 0.042;
		defsim->H0 = 0.732;
	}
	else {
		defsim->omega_m = 0.26;
		defsim->omega_l = 0.74;
		defsim->omega_b = 0.044;
		defsim->H0 = 0.72;
	}
	defsim->boxlen_ini = 1024.;
	defsim->amax = 1.;
	defsim->aexp = 1.;
	defsim->nx = defsim->ny = defsim->nz = 1024;
}

/*
#include "mpi.h"
void determine_mpi_misc_param(RamsesType *simpar){
	int nid,myid;
	float zi;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	simpar->nid = nid;
	simpar->myid = myid;
#ifdef XYZDBL
	simpar->xyzshiftflag = 1;
#else
	simpar->xyzshiftflag = 0;
#endif
	zi = simpar->amax-1;
	simpar->omei = simpar->omega_m*pow(1+zi,3)/(simpar->omega_m*pow(1+zi,3) +
			simpar->omega_l+(1-simpar->omega_m-simpar->omega_l)*pow(1+zi,2));
	simpar->mx = simpar->nx/simpar->nspace;
	simpar->my = simpar->ny/simpar->nspace;
	simpar->mz = simpar->nz/simpar->nspace;
	simpar->mxmy = simpar->mx * simpar->my;
	simpar->lnx = simpar->nx;
	simpar->lny = simpar->ny;
	simpar->lnz = simpar->nz;
	simpar->sphere_radius = 4;
	simpar->particle_radius = 4;
	simpar->rth = 8/(simpar->boxsize/simpar->nx);
	simpar->zinit = simpar->amax -1;
}
*/
