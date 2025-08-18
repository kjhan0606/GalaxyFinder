/* 
 icc -o fofmassfunc fofmassfunc.c -lm  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL -DXYZDBL -g 
 * */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "../ramses.h"
#include "fof.h"







HaloQ *halo;


#define NH 500
float unitmass;
#define pow3(a) ((a)*(a)*(a))

int main(int argc, char *argv[]){
    int i,j,k;
    int mp;

    FILE *fp, *fp2;
    char infile[180],infile2[190];
    float boxsize,hubble,npower,omep,omepl,bias,amax,astep,anow,omepb;
    int nx,nspace;


	dptype rminmass, rmaxmass;

	rminmass = log10(5.E6);
	rmaxmass = log10(5.E14);
	int nbin= atoi(argv[2]);
	dptype binstep=(rmaxmass-rminmass)/nbin;

	dptype *massbin = (dptype*)malloc(sizeof(dptype)*nbin);
	dptype *mtotbin = (dptype*)malloc(sizeof(dptype)*nbin);
	dptype *mgasbin = (dptype*)malloc(sizeof(dptype)*nbin);
	dptype *mstarbin = (dptype*)malloc(sizeof(dptype)*nbin);
	dptype *mdmbin = (dptype*)malloc(sizeof(dptype)*nbin);
	dptype *msinkbin = (dptype*)malloc(sizeof(dptype)*nbin);
	dptype *err = (dptype*)malloc(sizeof(dptype)*nbin);

	for(i=0;i<nbin;i++){
		massbin[i] = rminmass + (rmaxmass-rminmass)/nbin * i;
		mstarbin[i] = mgasbin[i] = mdmbin[i] = msinkbin[i] = mtotbin[i] = 0;
	}



    halo=(HaloQ*) malloc(sizeof(HaloQ)*NH);
	sprintf(infile,"./FoF_Data/FoF.%.5d/FoF_halo_cat.%.5d", atoi(argv[1]), atoi(argv[1]));
    fp=fopen(infile,"r");
    fread(&boxsize,sizeof(float),1,fp); 
	fread(&hubble,sizeof(float),1,fp); 
	fread(&omep,sizeof(float),1,fp); 
	fread(&omepb,sizeof(float),1,fp); 
	fread(&omepl,sizeof(float),1,fp); 
	fread(&amax,sizeof(float),1,fp); 
	fread(&anow,sizeof(float),1,fp); 
	fprintf(stderr,"redshift = %g\n", (amax/anow)-1);

	while((mp=fread(halo,sizeof(HaloQ),NH,fp))){
		int ntmp;
		for(i=0;i<mp;i++){
//			halo[i].mass /= (hubble/100);
			if(halo[i].mass >0) {ntmp = (log10(halo[i].mass)-rminmass)/binstep; if(ntmp>=0 && ntmp < nbin) mtotbin[ntmp] += 1;}
			if(halo[i].mgas>0) {ntmp = (log10(halo[i].mgas)-rminmass)/binstep; if(ntmp>=0 && ntmp < nbin) mgasbin[ntmp] += 1;}
			if(halo[i].mdm>0) {ntmp = (log10(halo[i].mdm)-rminmass)/binstep; if(ntmp>=0 && ntmp < nbin) mdmbin[ntmp] += 1;}
			if(halo[i].mstar>0) {ntmp = (log10(halo[i].mstar)-rminmass)/binstep; if(ntmp>=0 && ntmp < nbin) mstarbin[ntmp] += 1;}
			if(halo[i].msink>0) {ntmp = (log10(halo[i].msink)-rminmass)/binstep; if(ntmp>=0 && ntmp < nbin) msinkbin[ntmp] += 1;}
		}
    }
	fclose(fp);
	dptype yzmin,yzmax;
	yzmin = 0.445312; 
	yzmax = 0.554688; 
//	dptype volume = boxsize*boxsize*boxsize*(yzmax-yzmin)*(yzmax-yzmin) *1.E6/( hubble*hubble *hubble);
	dptype volume = boxsize*boxsize*boxsize*(yzmax-yzmin)*(yzmax-yzmin);
	dptype cmass;
	for(i=0;i<nbin;i++){
		err[i] = sqrt(mtotbin[i]) / volume;
		mtotbin[i] /= volume;
		mgasbin[i] /= volume;
		mstarbin[i] /= volume;
		mdmbin[i] /= volume;
		msinkbin[i] /= volume;
	}
	FILE *wp = fopen("fofmassfunc.dat","w");
	for(i=0;i<nbin;i++){ 
		dptype binsize; 
		binsize = pow(10,massbin[i+1]) - pow(10,massbin[i]); 
		massbin[i] = (pow(10,massbin[i]) + pow(10,massbin[i+1]))*0.5; 
		for(j=i;j<nbin;j++){
			cmass += massbin[i];
		}
		printf("%lg %lg  %lg %lg %lg\n",massbin[i], mtotbin[i]/binsize*massbin[i]*log(10.),
				(mtotbin[i]+err[i])/binsize*massbin[i]*log(10.),
				(mtotbin[i]-err[i])/binsize*massbin[i]*log(10.),
				cmass);
		fprintf(wp,"%lg %lg  %lg %lg %lg\n",massbin[i], mtotbin[i]/binsize*massbin[i]*log(10.),
				(mtotbin[i]+err[i])/binsize*massbin[i]*log(10.),
				(mtotbin[i]-err[i])/binsize*massbin[i]*log(10.),
				cmass);
	}
	
}
