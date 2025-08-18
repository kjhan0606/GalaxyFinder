/*
 
icc -g -o dist dist.c -DINDEX -DVarPM   -DXYZDBL  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL  -lm
 
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
#include "ramses.h"
#include "fof.h"
*/


#define nstep 4096



int main(int argc, char *argv[]) {

	float xmin,ymin,xmax,ymax;


	int i,j,k;
	float sigma[nstep],prob[nstep];
	float minsig = -4;
	float maxsig = 4;
	float sigstep = (maxsig-minsig)/nstep;
	for(i=0;i<nstep;i++){
		float sig = minsig + sigstep*i;
		sigma[i] = sig;
		prob[i] = 1-erfc(sig)/2.;
//		printf("%g %g\n", prob[i],sigma[i]);
	}



	FILE *fp = fopen(argv[1],"r");
	int nx,ny;
	fread(&nx,sizeof(int), 1, fp);
	fread(&ny,sizeof(int), 1, fp);
	float img[nx*ny];
	fread(img, sizeof(float), nx*ny, fp);
	fread(&xmin, sizeof(float), 1, fp);
	fread(&xmax, sizeof(float), 1, fp);
	fread(&ymin, sizeof(float), 1, fp);
	fread(&ymax, sizeof(float), 1, fp);
	fclose(fp);


	float xstep, ystep;
	xstep = (xmax-xmin)/nx;
	ystep = (ymax-ymin)/ny;
//	printf("xmin/max = %g %g : ymin/max = %g %g\n",xmin,xmax,ymin,ymax);


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
	FILE *wp = fopen(argv[2],"w");
	fwrite(&nx,sizeof(int), 1, wp);
	fwrite(&ny,sizeof(int), 1, wp);
	fwrite(dist,sizeof(float), nx*ny, wp);
	fclose(wp);
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
		printf("%g %g %g %g\n", mass,y0, y1, y2);
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
		printf("%g %g %g %g\n", mass,y0, y1, y2);
	}
}


