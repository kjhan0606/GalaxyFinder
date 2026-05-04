/* NAME: rd_sink_csv
 * PURPOSE: Read RAMSES sink CSV (output_NNNNN/sink_NNNNN.csv) and populate
 *          ram->sink[] using the same SinkType layout/units as rd_sink.c.
 *          The RAMSES CSV columns (in order) are:
 *            id, msink, x, y, z, vx, vy, vz,
 *            lx, ly, lz, tform, acc_rate, del_mass,
 *            rho_gas, cs**2, etherm, vx_gas, vy_gas, vz_gas,
 *            mbh, dmfsink, level
 *          Mapping to SinkType: id->id, msink->mass, x/y/z->x/y/z,
 *          vx/vy/vz->vx/vy/vz, lx/ly/lz->Jx/Jy/Jz, tform->tbirth,
 *          del_mass->dMsmbh. Remaining SinkType fields (dMBH_coarse,
 *          dMEd_coarse, Esave, Sx,Sy,Sz, Smag, eps) are not present in
 *          the CSV and are zeroed.
 *          If the file does not exist, returns 0 with ram->nsink=0.
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "ramses.h"
#include "Memory.h"

#define SINK_CSV_LINEMAX 4096

int rd_sink_csv(RamsesType *ram, char *infile){
	FILE *fp = fopen(infile,"r");
	if(fp == NULL){
		ram->nsink = 0;
		ram->nindsink = 0;
		ram->sink = NULL;
		fprintf(stdout,"rd_sink_csv: %s not found; setting nsink=0\n", infile);
		fflush(stdout);
		return 0;
	}

	char line[SINK_CSV_LINEMAX];
	int nsink = 0;
	if(fgets(line, sizeof(line), fp) == NULL){
		fclose(fp);
		ram->nsink = 0;
		ram->nindsink = 0;
		ram->sink = NULL;
		return 0;
	}
	while(fgets(line, sizeof(line), fp)){
		char *p = line;
		while(*p==' '||*p=='\t') p++;
		if(*p=='\0'||*p=='\n'||*p=='#') continue;
		nsink++;
	}

	ram->nsink = nsink;
	ram->nindsink = nsink;
	if(nsink == 0){
		fclose(fp);
		ram->sink = NULL;
		return 0;
	}

	ram->sink = (SinkType*)Malloc(sizeof(SinkType)*nsink, PPTR(ram->sink));
	SinkType *sink = ram->sink;

	rewind(fp);
	if(fgets(line, sizeof(line), fp) == NULL){
		fclose(fp);
		ram->nsink = 0;
		Free(ram->sink);
		ram->sink = NULL;
		return 0;
	}

	int i = 0;
	while(fgets(line, sizeof(line), fp) && i < nsink){
		char *p = line;
		while(*p==' '||*p=='\t') p++;
		if(*p=='\0'||*p=='\n'||*p=='#') continue;

		int id, level;
		double msink, x, y, z, vx, vy, vz;
		double lx, ly, lz, tform, acc_rate, del_mass;
		double rho_gas, cs2, etherm, vx_gas, vy_gas, vz_gas;
		double mbh, dmfsink;
		int n = sscanf(line,
			" %d , %lf , %lf , %lf , %lf , %lf , %lf , %lf ,"
			" %lf , %lf , %lf , %lf , %lf , %lf ,"
			" %lf , %lf , %lf , %lf , %lf , %lf ,"
			" %lf , %lf , %d",
			&id, &msink, &x, &y, &z, &vx, &vy, &vz,
			&lx, &ly, &lz, &tform, &acc_rate, &del_mass,
			&rho_gas, &cs2, &etherm, &vx_gas, &vy_gas, &vz_gas,
			&mbh, &dmfsink, &level);
		if(n < 23){
			fprintf(stderr,"rd_sink_csv: short line %d (got %d fields): %s",
				i, n, line);
			continue;
		}
		sink[i].id = id;
		sink[i].mass = msink;
		sink[i].x = x; sink[i].y = y; sink[i].z = z;
		sink[i].vx = vx; sink[i].vy = vy; sink[i].vz = vz;
		sink[i].Jx = lx; sink[i].Jy = ly; sink[i].Jz = lz;
		sink[i].tbirth = tform;
		sink[i].dMsmbh = del_mass;
		sink[i].dMBH_coarse = 0.0;
		sink[i].dMEd_coarse = 0.0;
		sink[i].Esave = 0.0;
		sink[i].Sx = 0.0; sink[i].Sy = 0.0; sink[i].Sz = 0.0;
		sink[i].Smag = 0.0;
		sink[i].eps = 0.0;
		i++;
	}
	fclose(fp);

	int npartp = i;
	ram->nsink = npartp;
	ram->nindsink = npartp;

	for(i=0;i<npartp;i++){
		sink[i].mass *= ram->scale_m;
		sink[i].x *= ram->mpcscale_l;
		sink[i].y *= ram->mpcscale_l;
		sink[i].z *= ram->mpcscale_l;
		sink[i].vx *= ram->kmscale_v;
		sink[i].vy *= ram->kmscale_v;
		sink[i].vz *= ram->kmscale_v;
		sink[i].dMsmbh *= ram->scale_m;
		sink[i].dMBH_coarse *= ram->scale_m;
		sink[i].dMEd_coarse *= ram->scale_m;
		sink[i].Jx *= ram->scale_m * ram->scale_l/kpc * ram->kmscale_v;
		sink[i].Jy *= ram->scale_m * ram->scale_l/kpc * ram->kmscale_v;
		sink[i].Jz *= ram->scale_m * ram->scale_l/kpc * ram->kmscale_v;
		sink[i].Smag *= ram->scale_m * ram->scale_l/kpc * ram->kmscale_v;
	}

	fprintf(stdout,"rd_sink_csv: read %d sink particles from %s\n",
		npartp, infile);
	fflush(stdout);
	return npartp;
}
