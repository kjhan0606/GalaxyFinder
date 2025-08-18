/*
 
icc -g -o mfrac mfrac.c -DINDEX -DVarPM   -DXYZDBL  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL  -lm
 
*/
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
/*
#include "header.h"
#include "tree.h"
#include "defs.h"
*/
#include "ramses.h"
#include "fof.h"





void FREAD(HaloQ *haloq, FILE *fp, DmType *dm, GasType *gas,
		SinkType *sink, StarType *star){
    fread(dm, sizeof(DmType), haloq->npdm,fp);
    fread(gas, sizeof(GasType), haloq->npgas,fp);
    fread(sink, sizeof(SinkType), haloq->npsink,fp);
    fread(star, sizeof(StarType), haloq->npstar,fp);
}

int main(int argc, char *argv[]) {

	int i,j,k;

	int nstep = atoi(argv[1]);
	char rheader[190];
	char dir[190];

	sprintf(dir,"./FoF_Data/FoF.%.5d/",nstep);
	sprintf(rheader,"%s/FoF_halo_cat.%.5d",dir,nstep);
	FILE *rhfp = fopen(rheader,"r");
	float bias,astep,anow,npower;
	char rfilename[190];
	sprintf(rfilename,"%sFoF_member_particle.%.5d",dir,nstep); 
	FILE *fp = fopen(rfilename,"r"); 
	{ 
		float bias; 
		float astep,anow; 
		int npower; 
		float size, hubble, omep,omepb,omeplam,amax;
		fread(&size,sizeof(float),1,rhfp); 
		fread(&hubble,sizeof(float),1,rhfp); 
		fread(&omep,sizeof(float),1,rhfp); 
		fread(&omepb,sizeof(float),1,rhfp); 
		fread(&omeplam,sizeof(float),1,rhfp); 
		fread(&amax,sizeof(float),1,rhfp); 
		fread(&anow,sizeof(float),1,rhfp); 
	}

	int np;
	int nx, ny;

	nx = ny = 256;

	float img[nx*ny];
	for(i=0;i<nx*ny;i++) img[i] = 0;


	float xmin,ymin,xmax,ymax;
	xmin = 10.;
	xmax = 14.3;
//	ymin = 0.005;
//	ymax = 0.05;
	ymin = 0.000;
	ymax = 0.5;
	float xstep = (xmax-xmin)/nx;
	float ystep = (ymax-ymin)/ny;




	HaloQ haloq[1000000];
	while((np=fread(haloq,sizeof(HaloQ), 1000000, rhfp)) > 0){
		int i;
		for(i=0;i<np;i++){
			dptype mfrac;
//			mfrac = (haloq[i].mstar+haloq[i].mgas+haloq[i].msink)/haloq[i].mass;
//			mfrac = (haloq[i].mstar)/haloq[i].mass;
			mfrac = (haloq[i].mgas)/haloq[i].mass;

			float x = (log10(haloq[i].mass)-xmin)/xstep;
			float y = (mfrac-ymin)/ystep;
			int ix = x;
			int iy = y;
			if(ix>=0&& ix < nx && iy >=0 && iy < ny) img[ix+nx*iy] += 1.;
		}
	}
	{
		int j;
		for(i=0;i<nx;i++){
			double sum=0;
			for(j=0;j<ny;j++){
				sum += img[i+nx*j];
			}
			if(sum>0){
				for(j=0;j<ny;j++){
					img[i+nx*j] = img[i+nx*j]/sum;
				}
			}
		}
	}
	FILE *wp = fopen(argv[2],"w");
	fwrite(&nx, sizeof(int), 1, wp);
	fwrite(&ny, sizeof(int), 1, wp);
	fwrite(img, sizeof(float), nx*ny, wp);
	fwrite(&xmin, sizeof(float), 1, wp);
	fwrite(&xmax, sizeof(float), 1, wp);
	fwrite(&ymin, sizeof(float), 1, wp);
	fwrite(&ymax, sizeof(float), 1, wp);
	fclose(wp);




	float sigma[nstep],prob[nstep];
	float minsig = -4;
	float maxsig = 4;
	float sigstep = (maxsig-minsig)/nstep;
	for(i=0;i<nstep;i++){
		float sig = minsig + sigstep*i;
		sigma[i] = sig;
		prob[i] = 1-erfc(sig)/2.;
	}



	float dist[nx*ny];

	for(i=0;i<nx;i++){
		double runningsum = 0;
		for(j=0;j<ny;j++){
			runningsum += img[i+nx*j];
		}
		if(runningsum >0){
			runningsum = 0;
			for(j=0;j<ny;j++){
				runningsum += img[i+nx*j];
				if(runningsum >1.L) runningsum = 1;
				else if(runningsum ==0) {
					dist[i+nx*j] = minsig;
					goto pos1;
				}
				for(k=0;k<nstep-1;k++){
					if(prob[k+1]>=runningsum && prob[k]<= runningsum){
						dist[i+nx*j] = sigma[k];
						goto pos1;
					}
				}
				dist[i+nx*j] = dist[i+nx*(j-1)];
pos1:           continue;
			}
		}
		else {
			for(j=0;j<ny;j++) dist[i+nx*j] = 0;
		}
	}
	wp = fopen(argv[3],"w");
	fwrite(&nx,sizeof(int), 1, wp);
	fwrite(&ny,sizeof(int), 1, wp);
	fwrite(dist,sizeof(float), nx*ny, wp);
	fclose(wp);



	wp = fopen(argv[4],"w");
	for(i=0;i<nx;i++){
		float mass = xmin + (i+0.5)*xstep;
		mass = pow(10.L,mass);
		float y1,y2,y0;
		for(j=0;j<ny;j++){
			float ytmp = ymin + ystep*(j+0.5);
			if(dist[i+nx*j]<0 && dist[i+nx*(j+1)] >= 0) y0 = ytmp;
			else if(dist[i+nx*j]<1 && dist[i+nx*(j+1)] >= 1) y1 = ytmp;
			else if(dist[i+nx*j]<2 && dist[i+nx*(j+1)] >= 2) y2 = ytmp;
		}
		fprintf(wp, "%g %g %g %g\n", mass,y0, y1, y2);
	}

	for(i=nx-1;i>=0;i--){
		float mass = xmin + (i+0.5)*xstep;
		mass = pow(10.L,mass);
		float y1,y2,y0;
		y0 = y1 = y2 = ymin;
		for(j=0;j<ny;j++){
			float ytmp = ymin + ystep*(j+0.5);
			if(dist[i+nx*j]<0 && dist[i+nx*(j+1)] >= 0) y0 = ytmp;
			else if(dist[i+nx*j]<-1 && dist[i+nx*(j+1)] >= -1) y1 = ytmp;
			else if(dist[i+nx*j]<-2 && dist[i+nx*(j+1)] >= -2) y2 = ytmp;
		}
		fprintf(wp, "%g %g %g %g\n", mass,y0, y1, y2);
	}
	fclose(wp);
}

