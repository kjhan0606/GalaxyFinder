/* 
 icc -o fofread fofread.c -lm  -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL -DXYZDBL -g 
 * */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "ramses.h"
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
    float size,hubble,npower,omep,omepl,bias,amax,astep,anow,omepb;
    int nx,nspace;


    halo=(HaloQ*) malloc(sizeof(HaloQ)*NH);
    fp=fopen(argv[1],"r");
    fp2=fopen(argv[2],"r");
        fread(&size,sizeof(float),1,fp);
        fread(&hubble,sizeof(float),1,fp);
        fread(&omep,sizeof(float),1,fp);
        fread(&omepb,sizeof(float),1,fp);
        fread(&omepl,sizeof(float),1,fp);
        fread(&amax,sizeof(float),1,fp);
        fread(&anow,sizeof(float),1,fp);

	j = 0;
    while((mp=fread(halo,sizeof(HaloQ),NH,fp))){
        for(i=0;i<mp;i++) { 
			void *p = (void *)malloc(sizeof(SinkType)*halo[i].np);
			DmType *dm = (DmType*)p;
			fread(dm, sizeof(DmType), halo[i].npdm, fp2);
			GasType *gas = (GasType*)p;
			fread(gas, sizeof(GasType), halo[i].npgas, fp2);
			SinkType *sink = (SinkType*)p;
			fread(sink, sizeof(SinkType), halo[i].npsink, fp2);
			StarType *star = (StarType*)p;
			fread(star, sizeof(StarType), halo[i].npstar, fp2);


			free(p);
			j ++;
		}
    }
	fclose(fp);
	fclose(fp2);

}
