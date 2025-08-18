#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include <mpi.h>
#include "../ramses.h"
#define DEFINE_SIM_PARA
#include "pmheader.h"
#undef DEFINE_SIM_PARA

#include "fof.h"
#include "Time.h"
#define MIN(a,b) a > b? b: a
#define MAX(a,b) a > b? a: b
#define pow2(a) ((a)*(a))
#define pow3(a) ((a)*(a)*(a))
particle p;
size_t readparticle(FoFTPtlStruct **,size_t ,int ,int ,int ,char *,int);
size_t read_ramses_data(FoFTPtlStruct **, size_t , char *, char *);
POSTYPE fof_link = 0.2;
POSTYPE lx,ly,lz;

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
	long long ntreemax = 1000000L;
	float wtime;
	int nfile;
	FILE *fp;
	size_t nhalo;
	POSTYPE zmin,zmax,zminlocal;
	size_t npadd,npwrite,npread;
	MPI_Status status;
	int initfile,finalfile;
	char infile[190], infolder[190], outfolder[190], garfolder[190];

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);



	if(argc != 4) {
		if(myid==0){
			fprintf(stderr,"Error in # of arguments\n");
			fprintf(stderr,"%%fof fileleadingheadername nstep nfiles\n");
		}
		MPI_Finalize();
		exit(199);
	}
	else {
		FILE *ffp;
		double r2kineticfact,HSUB;

		{
			FILE *ffp;
			RamsesType read_head(FILE *);
			nstep = atoi(argv[2]);

			sprintf(infolder,"./FoF_Data/NewDD.%.5d/",nstep);
			sprintf(outfolder,"./FoF_Data/FoF.%.5d/",nstep);
			sprintf(garfolder,"./FoF_Garbage/Garb.%.5d/",nstep);
			mkfolder(outfolder);
			mkfolder(garfolder);

			sprintf(infile,"%s%s.%.5d.%.5d.info",infolder,argv[1],nstep,myid);
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
			MPI_Bcast(&simpar,sizeof(SimParameters),MPI_BYTE,0,MPI_COMM_WORLD);
			nfile = atoi(argv[3]);
			size = simpar.boxlen_ini;
			hubble = simpar.H0;
			omep = simpar.omega_m;
			omepb = simpar.omega_b;
			omeplam = simpar.omega_l;
			nx = simpar.nx;
			ny = simpar.ny;
			nz = simpar.nz;
			anow = simpar.aexp;
		}

	if(0){
		int kkk = 1;
		while(kkk) {
			kkk = 2;
		}
	}


        ny=nz=nx;
		a = anow;
		amax = simpar.amax;
		N3 = (size_t)(nx)*(size_t)(ny)*(size_t)(nz)/nfile;
		N3 = N3*0.5;
		/*
		rscale = size/nx;
		*/
		rscale = 1;
		amax = 1;

		HSUB = sqrt(omep*pow3(amax/a)+omeplam+(1.-omep-omeplam)*pow2(amax/a));
		/*
		vscale = size/amax*a*a*100.*HSUB;
		*/
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
	sprintf(halofile,"%s/FoF_halo_cat.%.5d",outfolder,nstep);
	sprintf(memparticlefile,"%s/FoF_member_particle.%.5d",outfolder,nstep);
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
	lx = simpar.boxlen_ini;
	ly = simpar.boxlen_ini;
	lz = simpar.boxlen_ini;
	box.x = box.y = box.z = 0.;
	box.width = lx;

	printf("P%d has a periodic box of Lx=%g Ly=%g Lz=%g\n",myid,lx,ly,lz);
	np = 0;
	npwrite = 0;
#ifdef USE_MASTER
	mid = nid -1;
#error Not yet implmented
#else
	mid = nid;
#endif
	initfile = nfile/mid*myid;
	finalfile = nfile/mid*(myid+1);
	for(nowfile=initfile;nowfile<finalfile;nowfile++){
		sprintf(infile,"%s.%.5d%.5d",argv[1],nstep,nowfile);
		/*
		npadd = readparticle(&ptl,np,nstep,nowfile,nz,infile,N3);
		*/
		char type[100]; 
		sprintf(type,"DM"); sprintf(infile,"%sHR5.%.5d.%s.%.5d.dat",infolder,nstep,type,nowfile);
		npadd = read_ramses_data(&ptl,np,infile,type); np += npadd;
		sprintf(type,"STAR"); sprintf(infile,"%sHR5.%.5d.%s.%.5d.dat",infolder,nstep,type,nowfile);
		npadd = read_ramses_data(&ptl,np,infile,type); np += npadd;
		sprintf(type,"SINK"); sprintf(infile,"%sHR5.%.5d.%s.%.5d.dat",infolder,nstep,type,nowfile);
		npadd = read_ramses_data(&ptl,np,infile,type); np += npadd;
		sprintf(type,"GAS"); sprintf(infile,"%sHR5.%.5d.%s.%.5d.dat",infolder,nstep,type,nowfile);
		npadd = read_ramses_data(&ptl,np,infile,type); np += npadd;
		printf("P%d has np= %ld : %s \n", myid, np, infile);



		if(nowfile == finalfile-1){
			int src,dest;
			src = (myid+1+mid)%mid;
			dest = (myid-1+mid)%mid;
			MPI_Sendrecv(&npwrite,sizeof(size_t),MPI_BYTE,dest,0,
					&npread,sizeof(size_t),MPI_BYTE,src,0,MPI_COMM_WORLD,&status);
			{
				ptl = (FoFTPtlStruct*)realloc(ptl,
						sizeof(FoFTPtlStruct)*(np+npread));
				ReadBottomFaceContact(ptl+np,npread,linked,src,nstep,nz);
				np += npread;
			}
		}
		zmin = 2.E23;
		zmax = -2.E23;
		for(i=0;i<np;i++){
			zmin = MIN(zmin,ptl[i].z);
			zmax = MAX(zmax,ptl[i].z);
		}

		if(nowfile==initfile) zminlocal = zmin;

		ptl = (FoFTPtlStruct *)realloc(ptl,sizeof(FoFTPtlStruct)*np);
		printf("P%d: Now we have zmin=%g zmax=%g with np = %d\n",myid,zmin,zmax,np);
		if(np /NODE_HAVE_PARTICLE> ntreemax){
			ntreemax = np/NODE_HAVE_PARTICLE;
			TREE = (FoFTStruct *)realloc(TREE,sizeof(FoFTStruct)*ntreemax);
		}
//		MPI_Barrier(MPI_COMM_WORLD); if(myid==0) printf("passed 1\n");
		for(i=0;i<np;i++){
//			ptl[i].type = TYPE_PTL;
			ptl[i].sibling = &ptl[i+1];
			ptl[i].haloindx = -1;
		}
		ptl[np-1].sibling = NULL;
		FoF_Make_Tree(TREE,ptl,np,box);
		/* Local FoF and tag with a new halo number nhalo */
		nhalo = 0;
		for(i=0;i<np;i++){
			if(ptl[i].included == NO){
				p.x = ptl[i].x;
				p.y = ptl[i].y;
				p.z = ptl[i].z;
				p.link02 = ptl[i].link02;
				num=pnew_fof_link(&p,fof_link,TREE,ptl,linked,nhalo,lx,ly,lz);
				nhalo ++;
			}
		}
		printf("P%d has %ld halos from %ld particles\n",myid,nhalo,np);
		{
			HaloBound *halobound;
			int flag;
			halobound = (HaloBound *)malloc(sizeof(HaloBound)*nhalo);
			CheckHaloBound(nhalo,halobound,ptl,np,fof_link,
				zmin,zmax, zminlocal);
			printf("P%d has checked halo boundary\n",myid);fflush(stdout);
			if(nowfile != finalfile-1) {
				WriteIsolatedHalo(nhalo,halobound,ptl,linked,halofile,
						memparticlefile);
			}
			else if(nowfile == finalfile-1){
				WriteFinalHalo(nhalo,halobound,ptl,linked,halofile,
						memparticlefile);
				MPI_Finalize();
				return 0;
			}
			printf("P%d has saved isolated halos\n",myid);fflush(stdout);
			if(nowfile == initfile) flag = 0;
			else flag = 1;
			npwrite += WriteBottomFaceContact(nhalo,halobound,ptl,linked,
					flag,nstep);
			printf("P%d has saved Bottom Faced Halo %ld\n",myid,npwrite);fflush(stdout);
			np = StackUpContactParticleLeftWard(nhalo,halobound,ptl,np);
			free(halobound);
		}
	}
	MPI_Finalize();
}
