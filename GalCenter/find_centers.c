#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<sys/types.h>
#include "Memory.h"
#include "header.h"
#include "ramses.h"
#include "tree.h"
#include "defs.h"
#include "hfind.h"
#include "galcenter.h"

typedef struct VecPos{
	dptype x,y,z,mass;
}VecPos;


extern float size;



Center Find_Center(VecPos *r, lint mp, dptype targetmass){
	lint i,j,k;
	dptype Rmax2= -9999;
	dptype cx,cy,cz,tmass;
	cx = cy = cz = tmass = 0;
	for(i=0;i<mp;i++){
		cx += r[i].x*r[i].mass;
		cy += r[i].y*r[i].mass;
		cz += r[i].z*r[i].mass;
		tmass += r[i].mass;
	}
	cx = cx/tmass;
	cy = cy/tmass;
	cz = cz/tmass;
	for(i=0;i<mp;i++){
		dptype tmpx = r[i].x - cx;
		dptype tmpy = r[i].y - cy;
		dptype tmpz = r[i].z - cz;
		dptype dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		Rmax2 = MAX(Rmax2, dist2);
	}
	dptype cox,coy,coz,omass;
	dptype cnx,cny,cnz,nmass;
	dptype Rox,Rnx,Rox2,Rnx2;
	dptype shiftr,shiftR;

	cox = cx; coy = cy; coz = cz; omass = tmass;
	Rox = sqrt(Rmax2)*pow(tmass/targetmass,0.33333333333L);

	int iter = 0;
	int criteria=1;
	do{
		Rox2 = Rox*Rox;
		cnx = cny = cnz = nmass = 0;
		for(i=0;i<mp;i++){
			dptype tmpx = r[i].x-cx;
			dptype tmpy = r[i].y-cy;
			dptype tmpz = r[i].z-cz;
			dptype dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
			if(dist2 < Rox2){
				cnx += r[i].x*r[i].mass;
				cny += r[i].y*r[i].mass;
				cnz += r[i].z*r[i].mass;
				nmass += r[i].mass;
			}
		}
		if(nmass >0) {
			cnx = cnx/nmass;
			cny = cny/nmass;
			cnz = cnz/nmass;
			shiftr = sqrt((cnx-cox)*(cnx-cox)+(cny-coy)*(cny-coy)+(cnz-coz)*(cnz-coz));
	
			Rnx = Rox*pow(targetmass/nmass, 0.3333333L);
			shiftR = Rnx/Rox-1.;
	
			if(shiftr < MIN_SHIFT && fabs(shiftR) < MIN_SHIFT) criteria = 0;
			else criteria = 1;
	
			cox = cox + (cnx-cox)*rSTEP;
			coy = coy + (cny-coy)*rSTEP;
			coz = coz + (cnz-coz)*rSTEP;
			Rox = Rox + (Rnx-Rox)*RSTEP;
		}
		else {
			Rox = Rox*2.;
		}
		iter ++;

	} while(criteria && iter < 200);


	Center center;
	center.x = cox;
	center.y = coy;
	center.z = coz;
	center.R = Rox;

	return center;

}

void find_centers(GalCenter *gal,FoFTPtlStruct *bp, lint np){
	lint mp;
	lint i,j,k;
	VecPos *r;
	dptype targetmass,nowmass;
	dptype cx,cy,cz;
	{
		dptype xmin,xmax;
		xmin = 1e20;
		xmax =-1e20;
		for(i=0;i<np;i++){
			xmin = MIN(xmin, bp[i].x);
			xmax = MAX(xmin, bp[i].x);
		}
		if(xmin < 20 && xmax > 700){
			dptype cx, tmass;
			cx = tmass = 0;
			for(i=0;i<np;i++){
				cx += bp[i].x *bp[i].mass;
				tmass += bp[i].mass;
			}
			cx = cx/tmass;
			if(cx > size*0.5) {
				for(i=0;i<np;i++){
					if(bp[i].x < size*0.5) bp[i].x += size;
				}
			}
			else {
				for(i=0;i<np;i++){
					if(bp[i].x > size*0.5) bp[i].x -= size;
				}
			}
		}

	}

	{
		mp = 0;
		for(i=0;i<np;i++){
			if(bp[i].type == TYPE_STAR){
				mp ++;
			}
		}
		if(mp>0){
			r = (VecPos*)Malloc(sizeof(VecPos)*mp,PPTR(r));
			nowmass = mp = 0;
			cx = cy = cz = 0;
			for(i=0;i<np;i++){
				if(bp[i].type == TYPE_STAR){
					r[mp].x = bp[i].x;
					r[mp].y = bp[i].y;
					r[mp].z = bp[i].z;
					r[mp].mass = bp[i].mass;
					nowmass += r[mp].mass;
					cx += r[mp].x *r[mp].mass;
					cy += r[mp].y *r[mp].mass;
					cz += r[mp].z *r[mp].mass;
					mp++;
				}
			}
			targetmass = MIN(nowmass*0.5, 1e9);
		
		}
		if(mp>30) gal->gal = Find_Center(r,mp, targetmass);
		else if(mp >0) {
			gal->gal.x = cx/nowmass; gal->gal.y = cy/nowmass; gal->gal.z = cz/nowmass;
		}
		else {
			gal->gal.x = gal->gal.y = gal->gal.z = NO_CENTER_POS;
		}
		if(mp>0) Free(r);
	}
	{
		mp = 0;
		for(i=0;i<np;i++){
			if(bp[i].type == TYPE_DM){
				mp ++;
			}
		}
		if(mp>0){
			r = (VecPos*)Malloc(sizeof(VecPos)*mp,PPTR(r));
			nowmass = mp = 0;
			cx = cy = cz = 0;
			for(i=0;i<np;i++){
				if(bp[i].type == TYPE_DM){
					r[mp].x = bp[i].x;
					r[mp].y = bp[i].y;
					r[mp].z = bp[i].z;
					r[mp].mass = bp[i].mass;
					nowmass += r[mp].mass;
					cx += r[mp].x *r[mp].mass;
					cy += r[mp].y *r[mp].mass;
					cz += r[mp].z *r[mp].mass;
					mp++;
				}
			}
	
			targetmass = MIN(nowmass*0.5, 1e9);
	
		}
		if(mp>30) gal->dmhalo = Find_Center(r,mp, targetmass);
		else if(mp >0) {
			gal->dmhalo.x = cx/nowmass; gal->dmhalo.y = cy/nowmass; gal->dmhalo.z = cz/nowmass;
		}
		else {
			gal->dmhalo.x = gal->dmhalo.y = gal->dmhalo.z = NO_CENTER_POS;
		}
		if(mp>0) Free(r);
	}
	{
		mp = 0;
		for(i=0;i<np;i++){
			if(bp[i].type == TYPE_GAS){
				mp ++;
			}
		}
		if(mp>0){
			r = (VecPos*)Malloc(sizeof(VecPos)*mp,PPTR(r));
			nowmass = mp = 0;
			cx = cy = cz = 0;
			for(i=0;i<np;i++){
				if(bp[i].type == TYPE_GAS){
					r[mp].x = bp[i].x;
					r[mp].y = bp[i].y;
					r[mp].z = bp[i].z;
					r[mp].mass = bp[i].mass;
					nowmass += r[mp].mass;
					cx += r[mp].x *r[mp].mass;
					cy += r[mp].y *r[mp].mass;
					cz += r[mp].z *r[mp].mass;
					mp++;
				}
			}
			targetmass = MIN(nowmass*0.5, 1e9);
	
		}
		if(mp>30) gal->gas = Find_Center(r,mp, targetmass);
		else if(mp >0) {
			gal->gas.x = cx/nowmass; gal->gas.y = cy/nowmass; gal->gas.z = cz/nowmass;
		}
		else {
			gal->gas.x = gal->gas.y = gal->gas.z = NO_CENTER_POS;
		}
		if(mp>0) Free(r);
	}

}
