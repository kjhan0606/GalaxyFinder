#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<limits.h>
#include <mpi.h>
#include "../ramses.h"
#define DEFINE_SIM_PARA
#include "pmheader.h"
#undef DEFINE_SIM_PARA

#include "fof.h"
#include "distributed.h"
#include "Time.h"
#define MIN(a,b) a > b? b: a
#define MAX(a,b) a > b? a: b
#define pow2(a) ((a)*(a))
#define pow3(a) ((a)*(a)*(a))


/* Number of file to be read */
#define MAXNFILE 100000
#define NLEN 190


particle p;
size_t readparticle(FoFTPtlStruct **,size_t ,int ,int ,int ,char *,int);
size_t read_ramses_data(FoFTPtlStruct **, size_t , char *, char *, double);
POSTYPE fof_link = 0.2;
POSTYPE lx,ly,lz;

int pflag;

	float size,hubble,npower,omep,omepb,omeplam,bias,smooth;
	int nx,ny,nz,nspace;
	float ntree,ntree1,theta;
	float zinit,amax,astep,anow,a;
	int iseed;
	float vscale,rscale;
	FILE *hfp,*pfp;
	char halofile[190],memparticlefile[190];
void xhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].x < ra[j+1].x) j++;
            if (rra.x < ra[j].x) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}
void yhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].y < ra[j+1].y) j++;
            if (rra.y < ra[j].y) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}
void zhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].z < ra[j+1].z) j++;
            if (rra.z < ra[j].z) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}

/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */

HaloQ haloproperty(particle *member,size_t nmem){
	HaloQ halo;
	particle *pp;
	size_t i,j,k;
	double cx,cy,cz;
	double vx,vy,vz;
	POSTYPE xmin,xmax,ymin,ymax,zmin,zmax;
	POSTYPE x1,x2,y1,y2,z1,z2;
	halo.np=nmem;
	cx=cy=cz=vx=vy=vz=0;
	pp=member;
	xmax = ymax = zmax = -1.E25;
	xmin = ymin = zmin = +1.E25;
	halo.npstar = halo.npgas = halo.npdm = halo.npsink = 0;
	for(i=0;i<nmem;i++){
		if(pp->type == TYPE_GAS) halo.npgas ++;
		else if(pp->type == TYPE_DM) halo.npdm ++;
		else if(pp->type == TYPE_SINK) halo.npsink ++;
		else if(pp->type == TYPE_STAR) halo.npstar ++;
		xmax = MAX(xmax,(pp->x));
		xmin = MIN(xmin,(pp->x));
		ymax = MAX(ymax,(pp->y));
		ymin = MIN(ymin,(pp->y));
		zmax = MAX(zmax,(pp->z));
		zmin = MIN(zmin,(pp->z));
		pp++;
	}
	if(pflag ==1) {
		if(xmin <= fof_link && xmax >=lx-fof_link){
			xhpsort(nmem,member-1);
			pp=member+1;
			for(i=1;i<nmem;i++){
				if(((pp->x)-(pp-1)->x) > 1.5*fof_link) {
					break;
				}
				pp++;
			}
			for(j=i;j<nmem;j++){
				member[j].x -= lx;
			}
		}
		if(ymin <= fof_link && ymax >=ly-fof_link){
                yhpsort(nmem,member-1);
                pp=member+1;
                for(i=1;i<nmem;i++){
                        if(((pp->y)-(pp-1)->y) > 1.5*fof_link) {
                                break;
                        }
						pp++;
                }
                for(j=i;j<nmem;j++){
                        member[j].y -= ly;
                }
     	}
		if(zmin <= fof_link && zmax >=lz-fof_link){
                zhpsort(nmem,member-1);
                pp=member+1;
                for(i=1;i<nmem;i++){
                        if(((pp->z)-(pp-1)->z) > 1.5*fof_link) {
                                break;
                        }
						pp++;
                }
                for(j=i;j<nmem;j++){
                        member[j].z -= lz;
                }
     	}
	}


	pp=member;
	halo.mass= 0;
	halo.mstar = halo.mdm = halo.msink = halo.mgas = 0;
	for(i=0;i<nmem;i++){
		if(pp->type == TYPE_GAS) {
			pp->p.gas.z = pp->x;
			pp->p.gas.y = pp->y;
			pp->p.gas.x = pp->z;
			halo.mgas += pp->p.gas.mass;
			cx+=(pp->z) * pp->p.gas.mass;
			cy+=(pp->y) * pp->p.gas.mass;
			cz+=(pp->x) * pp->p.gas.mass;
			vx += pp->p.gas.vx * pp->p.gas.mass;
			vy += pp->p.gas.vy * pp->p.gas.mass;
			vz += pp->p.gas.vz * pp->p.gas.mass;
			halo.mass += pp->p.gas.mass;
		}
		else if(pp->type == TYPE_DM) {
			pp->p.dm.z = pp->x;
			pp->p.dm.y = pp->y;
			pp->p.dm.x = pp->z;
			halo.mdm += pp->p.dm.mass;
			cx+=(pp->z) * pp->p.dm.mass;
			cy+=(pp->y) * pp->p.dm.mass;
			cz+=(pp->x) * pp->p.dm.mass;
			vx += pp->p.dm.vx * pp->p.dm.mass;
			vy += pp->p.dm.vy * pp->p.dm.mass;
			vz += pp->p.dm.vz * pp->p.dm.mass;
			halo.mass += pp->p.dm.mass;
		}
		else if(pp->type == TYPE_STAR) {
			pp->p.star.z = pp->x;
			pp->p.star.y = pp->y;
			pp->p.star.x = pp->z;
			halo.mstar += pp->p.star.mass;
			cx+=(pp->z) * pp->p.star.mass;
			cy+=(pp->y) * pp->p.star.mass;
			cz+=(pp->x) * pp->p.star.mass;
			vx += pp->p.star.vx * pp->p.star.mass;
			vy += pp->p.star.vy * pp->p.star.mass;
			vz += pp->p.star.vz * pp->p.star.mass;
			halo.mass += pp->p.star.mass;
		}
		else if(pp->type == TYPE_SINK) {
			pp->p.sink.z = pp->x;
			pp->p.sink.y = pp->y;
			pp->p.sink.x = pp->z;
			halo.msink += pp->p.sink.mass;
			cx+=(pp->z) * pp->p.sink.mass;
			cy+=(pp->y) * pp->p.sink.mass;
			cz+=(pp->x) * pp->p.sink.mass;
			vx += pp->p.sink.vx * pp->p.sink.mass;
			vy += pp->p.sink.vy * pp->p.sink.mass;
			vz += pp->p.sink.vz * pp->p.sink.mass;
			halo.mass += pp->p.sink.mass;
		}
		pp++;
	}
	/* position in h^-1 Mpc */
	halo.x = cx/halo.mass *rscale;
	halo.y = cy/halo.mass *rscale;
	halo.z = cz/halo.mass *rscale;
	/* velocity in km/sec */
	halo.vx = vx/halo.mass*vscale;
	halo.vy = vy/halo.mass*vscale;
	halo.vz = vz/halo.mass*vscale;
	return halo;
}
int myid,nid,mid;
int main(int argc,char *argv[]){
	double std,mean;
	int ntmp;
	size_t num,ii;
	POSTYPE tmpx,tmpy,tmpz,dist2;
	POSTYPE fplmf,ptlmass;
	size_t i,j,k;
	int si,snp;
	int nowfile;
	int N,M;
	size_t N3;
	size_t np,addhere;
	int nstep;
	FoFTPtlStruct *ptl;
	FoFTPtlStruct *ptr;
	particle *linked;
	int nfof;
	FoFBeginEndTree beginend;
	Box box;
	FoFTStruct *TREE;
	long long ntreemax = 9000000L;
	float wtime;
	int nfile;
	FILE *fp;
	size_t nhalo;
	POSTYPE zmin,zmax,zminlocal;
	size_t npadd,npwrite,npread;
	MPI_Status status;
	int initfile,finalfile;
	char infile[190], infolder[190], outfolder[190], garfolder[190];


	char infilegas[MAXNFILE*NLEN], infiledm[MAXNFILE*NLEN], infilestar[MAXNFILE*NLEN], infilesink[MAXNFILE*NLEN];
	double xoffset[MAXNFILE];

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);



	if(argc != 3 && argc !=4) {
		if(myid==0){
			fprintf(stderr,"Error in # of arguments\n");
			fprintf(stderr,"%%fof fileleadingheadername nstep nfiles\n");
			fprintf(stderr,"or\n");
			fprintf(stderr,"%%fof lightcone nstart nfinal nstepfiles\n");
		}
		MPI_Finalize();
		exit(199);
	}
	else if(argc ==4){
		char inlist[190];
		int nstart, nfinal;
		int nstepfile;
		if(myid==0){
			nfile = 0;
			/*
			printf("please input NewDD.xxx: nstart, nfinal\n");
			scanf("%d %d",&nstart, &nfinal);
			printf("please input nstepfile per step\n");
			scanf("%d",&nstepfile);
			*/
			nstart = atoi(argv[1]);
			nfinal = atoi(argv[2]);
			nstepfile = atoi(argv[3]);
			for(i=nfinal;i>=nstart;i--){
				struct stat buff;
				sprintf(infolder,"./FoF_Data/NewDD.%.5d/",i);
				for(j=0;j<nstepfile;j++){
					char type[100];
					sprintf(type,"DM");
					sprintf(infiledm+nfile*NLEN,"%sSN.%.5d.%s.%.5d.dat",infolder,i,type,j);
					int ierr = stat(infiledm+nfile*NLEN,&buff);
					if( ierr ==0 &&  buff.st_size >0){
						if(nfile==0) sprintf(infile,"%sSN.%.5d.%.5d.info",infolder,i,j);
						sprintf(type,"GAS");
						sprintf(infilegas+nfile*NLEN,"%sSN.%.5d.%s.%.5d.dat",infolder,i,type,j);
						sprintf(type,"SINK");
						sprintf(infilesink+nfile*NLEN,"%sSN.%.5d.%s.%.5d.dat",infolder,i,type,j);
						sprintf(type,"STAR");
						sprintf(infilestar+nfile*NLEN,"%sSN.%.5d.%s.%.5d.dat",infolder,i,type,j);
						xoffset[nfile] = (13-(i-100))* 717.229039925849L;
						nfile ++;
					}
				}
				printf("Now found nfile = %d for %d step direction with xoffset= %g\n",nfile,i,xoffset[nfile-1]);
			}
			FILE *ffp = fopen(infile,"r");
			printf("reading %s\n",infile);fflush(stdout);
			RamsesType read_head(FILE *);
			simpar = read_head(ffp);
			if(sizeof(simpar) == simpar.ramses_sizeof){
				fread(&simpar, sizeof(RamsesType), 1, ffp);
			}
			else {
				fprintf(stderr,"Warning: different size of RamsesType.\n");
				fprintf(stderr,"Warning: We cannot read the full precision parameter values from the binary format.\n");
				fprintf(stderr,"So, we degrade the precision from double to float in each parameters.\n");
			}
			fclose(ffp);
		}
		MPI_Bcast(&simpar,sizeof(SimParameters),MPI_BYTE,0,MPI_COMM_WORLD);
		MPI_Bcast(&nfile,sizeof(int),MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(infiledm,sizeof(char)*MAXNFILE*NLEN,MPI_CHAR,0,MPI_COMM_WORLD);
		MPI_Bcast(infilegas,sizeof(char)*MAXNFILE*NLEN,MPI_CHAR,0,MPI_COMM_WORLD);
		MPI_Bcast(infilesink,sizeof(char)*MAXNFILE*NLEN,MPI_CHAR,0,MPI_COMM_WORLD);
		MPI_Bcast(infilestar,sizeof(char)*MAXNFILE*NLEN,MPI_CHAR,0,MPI_COMM_WORLD);
		MPI_Bcast(xoffset,MAXNFILE,MPI_DOUBLE,0,MPI_COMM_WORLD);
		if(nfile%nid !=0){
			printf("nfile= %d is not divisible by nid= %d\n",nfile,nid);
			MPI_Finalize();
		}
		nstep = 0;
		sprintf(outfolder,"./FoF_Data/FoF.%.5d/",nstep);
		sprintf(garfolder,"./FoF_Garbage/Garb.%.5d/",nstep);
		mkfolder(outfolder);
		mkfolder(garfolder);
//		lx = ly = lz = simpar.boxlen_ini*10;
		lx = ly = lz = 15000L;
		pflag = 0;
		box.x = box.y = box.z = 0.;
		box.width = lx;
		sprintf(halofile,"%s/FoF_halo_cat.%.5d",outfolder,nstep);
		sprintf(memparticlefile,"%s/FoF_member_particle.%.5d",outfolder,nstep);
	}
	else if(argc == 3){
		FILE *ffp;
		double r2kineticfact;

		RamsesType read_head(FILE *);
		nstep = atoi(argv[1]);

		sprintf(infolder,"./FoF_Data/NewDD.%.5d/",nstep);
		sprintf(outfolder,"./FoF_Data/FoF.%.5d/",nstep);
		sprintf(garfolder,"./FoF_Garbage/Garb.%.5d/",nstep);
		if(myid==0){
			mkfolder(outfolder);
			mkfolder(garfolder);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		sprintf(infile,"%sSN.%.5d.%.5d.info",infolder,nstep,myid);
		if(myid==0){
			ffp = fopen(infile,"r");
			printf("reading %s\n",infile);fflush(stdout);
			simpar = read_head(ffp);
			if(sizeof(simpar) == simpar.ramses_sizeof){
				fread(&simpar, sizeof(RamsesType), 1, ffp);
			}
			else {
				fprintf(stderr,"Warning: different size of RamsesType.\n");
				fprintf(stderr,"Warning: We cannot read the full precision parameter values from the binary format.\n");
				fprintf(stderr,"So, we degrade the precision from double to float in each parameters.\n");
			}
			fclose(ffp);
		}
		pflag = 1;

		MPI_Bcast(&simpar,sizeof(SimParameters),MPI_BYTE,0,MPI_COMM_WORLD);
		nfile = atoi(argv[2]);

		char type[100]; 
	 	sprintf(type,"DM"); 
		for(i=0;i<nfile;i++) sprintf(infiledm+i*NLEN,"%sSN.%.5d.%s.%.5d.dat",infolder,nstep,type,i);
	 	sprintf(type,"GAS"); 
		for(i=0;i<nfile;i++) sprintf(infilegas+i*NLEN,"%sSN.%.5d.%s.%.5d.dat",infolder,nstep,type,i);
	 	sprintf(type,"SINK"); 
		for(i=0;i<nfile;i++) sprintf(infilesink+i*NLEN,"%sSN.%.5d.%s.%.5d.dat",infolder,nstep,type,i);
	 	sprintf(type,"STAR"); 
		for(i=0;i<nfile;i++) sprintf(infilestar+i*NLEN,"%sSN.%.5d.%s.%.5d.dat",infolder,nstep,type,i);
		for(i=0;i<nfile;i++) xoffset[i] = 0.L;
		lx = simpar.boxlen_ini;
		ly = simpar.boxlen_ini;
		lz = simpar.boxlen_ini;
		box.x = box.y = box.z = 0.;
		box.width = lx;
		printf("P%d has a periodic box of Lx=%g Ly=%g Lz=%g\n",myid,lx,ly,lz);
		sprintf(halofile,"%s/FoF_halo_cat.%.5d",outfolder,nstep);
		sprintf(memparticlefile,"%s/FoF_member_particle.%.5d",outfolder,nstep);
	}

	if(1){ 
		char hostname[HOST_NAME_MAX+1];
		gethostname(hostname, HOST_NAME_MAX+1);
		pid_t pid = getpid(); 
		printf("P%d has pid %d in hostname %s\n", myid, pid, hostname); 
	} 


	if(0){
		int kkk = 1;
		while(kkk) {
			kkk = 2;
		}
	}
	size = simpar.boxlen_ini;
	hubble = simpar.H0;
	omep = simpar.omega_m;
	omepb = simpar.omega_b;
	omeplam = simpar.omega_l;
	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	anow = simpar.aexp;


    ny=nz=nx;
	a = anow;
	amax = simpar.amax;
	N3 = (size_t)(nx)*(size_t)(ny)*(size_t)(nz)/nfile;
	N3 = N3*0.5;
	rscale = 1;
	amax = 1;
	double HSUB = sqrt(omep*pow3(amax/a)+omeplam+(1.-omep-omeplam)*pow2(amax/a));
	vscale = 1;
	fof_link = 0.2; /* the maximum fof link distance (input) in Mpc/h */

	if(myid==0){
		printf("P%d size = %g hubble = %g\n",myid,size,hubble);
		printf("P%d  omep = %g omeplam = %g smooth=%g\n",
				myid,omep,omeplam,smooth);
		printf("P%d  anow = %g \n",myid,anow);
		printf("P%d nx = %d ny= %d nz= %d \n",myid,nx,ny,nz);
		printf("P%d rscale = %g vscale= %g\n",myid,rscale,vscale);
	}
	if((ptl = (FoFTPtlStruct *) malloc(sizeof(FoFTPtlStruct)*100)) == NULL){
		fprintf(stderr,"Error allocating ptl %ld\n",N3);
		exit(99);
	}
	M = 1;
	/*
	if((p = (particle *) malloc(sizeof(particle)*M)) == NULL){
                fprintf(stderr,"Error allocating p\n");
                exit(99);
        }
		*/
	//linked = (particle *) malloc(sizeof(particle)*MaxLinkedParticles);
	linked = NULL;
	if((TREE = (FoFTStruct *) malloc(sizeof(FoFTStruct)*ntreemax)) == NULL){
                fprintf(stderr,"Error allocating TREE\n");
                exit(99);
        }
	wtime = WALLCLOCK();
	(void)WALLCLOCK();
	if(myid==0)
	{
		hfp=fopen(halofile,"w");
		pfp=fopen(memparticlefile,"w");
		fwrite(&size,sizeof(float),1,hfp);
		fwrite(&hubble,sizeof(float),1,hfp);
		fwrite(&omep,sizeof(float),1,hfp);
		fwrite(&omepb,sizeof(float),1,hfp);
		fwrite(&omeplam,sizeof(float),1,hfp);
		fwrite(&amax,sizeof(float),1,hfp);
		fwrite(&anow,sizeof(float),1,hfp);
		fclose(hfp);
		fclose(pfp);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	np = 0;
	npwrite = 0;
#ifdef USE_MASTER
	mid = nid -1;
#error Not yet implmented
#else
	mid = nid;
#endif
	/*
	 * Each input slab has exactly one owning rank.  Local FoF components that
	 * do not touch a rank face are final.  Components that touch either face
	 * exchange only face-band particles with adjacent ranks.  Component labels
	 * then propagate to convergence, so a halo may span any number of domains
	 * without gathering all boundary particles on rank 0.
	 */
	initfile = nfile*myid/mid;
	finalfile = nfile*(myid+1)/mid;
	printf("P%d has file ranges %d %d\n", myid, initfile, finalfile);
	for(nowfile=initfile;nowfile<finalfile;nowfile++){
		sprintf(infile,"SN.%.5d.%.5d",nstep,nowfile);
		char type[100]; 
//		sprintf(infile,"%sHR5.%.5d.%s.%.5d.dat",infolder,nstep,type,nowfile);
		sprintf(type,"DM"); 
		npadd = read_ramses_data(&ptl,np,infiledm+nowfile*NLEN,type, xoffset[nowfile]); np += npadd;
//		sprintf(type,"STAR"); sprintf(infile,"%sHR5.%.5d.%s.%.5d.dat",infolder,nstep,type,nowfile);
		sprintf(type,"STAR"); 
		npadd = read_ramses_data(&ptl,np,infilestar+nowfile*NLEN,type, xoffset[nowfile]); np += npadd;
//		sprintf(type,"SINK"); sprintf(infile,"%sHR5.%.5d.%s.%.5d.dat",infolder,nstep,type,nowfile);
		sprintf(type,"SINK"); 
		npadd = read_ramses_data(&ptl,np,infilesink+nowfile*NLEN,type, xoffset[nowfile]); np += npadd;
//		sprintf(type,"GAS"); sprintf(infile,"%sHR5.%.5d.%s.%.5d.dat",infolder,nstep,type,nowfile);
		sprintf(type,"GAS"); 
		npadd = read_ramses_data(&ptl,np,infilegas+nowfile*NLEN,type, xoffset[nowfile]); np += npadd;
		printf("P%d has np= %ld : %s \n", myid, np, infile);
		
	}
	if(np == 0){
		fprintf(stderr,"P%d owns no input particles\n",myid);
		MPI_Abort(MPI_COMM_WORLD,91);
	}

	ptl = (FoFTPtlStruct *)realloc(ptl,sizeof(FoFTPtlStruct)*np);
	printf("P%d: loaded %ld particles for its complete rank domain\n",myid,np);
	if(np/NODE_HAVE_PARTICLE*10 > ntreemax){
		ntreemax = np/NODE_HAVE_PARTICLE*10;
		TREE = (FoFTStruct *)realloc(TREE,sizeof(FoFTStruct)*ntreemax);
	}
	for(i=0;i<np;i++){
		ptl[i].sibling = (i+1 < np) ? &ptl[i+1] : NULL;
		ptl[i].haloindx = (size_t)-1;
	}
	linked = (particle*)malloc(sizeof(particle)*np);
	FoF_Make_Tree(TREE,ptl,np,box);
	nhalo = 0;
	for(i=0;i<np;i++){
		if(ptl[i].included == NO){
			p.x = ptl[i].x;
			p.y = ptl[i].y;
			p.z = ptl[i].z;
			p.link02 = ptl[i].link02;
			if(pflag==0) num=new_fof_link(&p,fof_link,TREE,ptl,linked,nhalo);
			else num=pnew_fof_link(&p,fof_link,TREE,ptl,linked,nhalo,lx,ly,lz);
			nhalo ++;
		}
	}
	free(linked);
	printf("P%d has %ld local components from %ld particles\n",myid,nhalo,np);

	{
		double local_max_link = 0.0;
		double boundary_link = 0.0;
		double domain_low = lx*(double)initfile/(double)nfile;
		double domain_high = lx*(double)finalfile/(double)nfile;
		for(i=0;i<np;i++) local_max_link = MAX(local_max_link,ptl[i].link02);
		MPI_Allreduce(&local_max_link,&boundary_link,1,MPI_DOUBLE,MPI_MAX,
				MPI_COMM_WORLD);
		if(myid==0) printf("Boundary linking width = %.9g Mpc/h\n",boundary_link);
		DistributedMergeAndWrite(ptl,np,nhalo,&TREE,&ntreemax,box,
				boundary_link,domain_low,domain_high,fof_link,pflag,lx,ly,lz,
				halofile,memparticlefile);
	}
	free(ptl);
	free(TREE);
	MPI_Finalize();
	return 0;
}
