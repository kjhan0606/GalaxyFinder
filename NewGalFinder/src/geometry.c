/* Periodic unwrap of one FoF halo into a simply connected box. */
#include "finder_internal.h"

int compare_float_ascending(const void *a,const void *b){
	const float aa = *(const float*)a;
	const float bb = *(const float*)b;
	if(aa<bb) return -1;
	else if(aa>bb) return 1;
	else return 0;
}

void unwrap_periodic_box(SimpleBasicParticleType *particles,lint n_particles,float *xinit,float *yinit,
		float *zinit,float *xmax,float *ymax,float *zmax){
	float min,max;
	SimpleBasicParticleType *tmpr;
	float *tmp,*tmparr;
	lint i,j;
	int *idx;
	float maxdist,tmpdist;
	int upbound;
	*xinit = size; *xmax = -size;
	*yinit = size; *ymax = -size;
	*zinit = size; *zmax = -size;
	for(i=0;i<n_particles;i++){
		tmpr = particles+i;
		*xinit = MIN(*xinit,tmpr->x);
		*xmax  = MAX(*xmax ,tmpr->x);
		*yinit = MIN(*yinit,tmpr->y);
		*ymax  = MAX(*ymax ,tmpr->y);
		*zinit = MIN(*zinit,tmpr->z);
		*zmax  = MAX(*zmax ,tmpr->z);
	}
	if(floor(*xinit) <= 0+nblur && ceil(*xmax) >= size-nblur){
		tmparr = (float *) Malloc(sizeof(float)*n_particles,PPTR(tmparr));
		tmp = tmparr;tmpr = particles;
		for(i=0;i<n_particles;i++) *(tmp++) = (tmpr++)->x;
		qsort(tmparr,n_particles,sizeof(float),compare_float_ascending);
		maxdist = -10.;
		for(i=0;i<n_particles-1;i++){
			tmpdist = tmparr[i+1] - tmparr[i];
			if(tmpdist > maxdist) {
				upbound = i;
				maxdist = tmpdist;
			}
		}
		for(i=0;i<n_particles;i++){
			if(particles[i].x <= tmparr[upbound]) particles[i].x += size;
		}
		/*
		for(i=0;i<=upbound;i++){
			particles[i].x += nx;
		}
		*/
		*xinit = size; *xmax = -size;
		for(i=0;i<n_particles;i++){
			*xinit = MIN(*xinit,particles[i].x);
			*xmax = MAX(*xmax,particles[i].x);
		}
		Free(tmparr);
	}
	if(floor(*yinit) <= 0+nblur && ceil(*ymax) >= size-nblur){
		tmparr = (float *) Malloc(sizeof(float)*n_particles,PPTR(tmparr));
		tmp = tmparr;
		tmpr = particles;
		for(i=0;i<n_particles;i++) *(tmp++) = tmpr++->y;
		qsort(tmparr,n_particles,sizeof(float),compare_float_ascending);
		maxdist = -10.;
		for(i=0;i<n_particles-1;i++){
			tmpdist = tmparr[i+1] - tmparr[i];
			if(tmpdist > maxdist) {
				upbound = i;
				maxdist = tmpdist;
			}
		}
		for(i=0;i<n_particles;i++){
			if(particles[i].y <= tmparr[upbound]) particles[i].y += size;
		}
		/*
		for(i=0;i<=upbound;i++){
			particles[i].y += ny;
		}
		*/
		*yinit = size; *ymax = -size;
		for(i=0;i<n_particles;i++){
			*yinit = MIN(*yinit,particles[i].y);
			*ymax = MAX(*ymax,particles[i].y);
		}
		Free(tmparr);
	}
	if(floor(*zinit) <= 0+nblur && ceil(*zmax) >= size-nblur){
		tmparr = (float *) Malloc(sizeof(float)*n_particles,PPTR(tmparr));
		tmp = tmparr;
		tmpr = particles;
		for(i=0;i<n_particles;i++) *(tmp++) = tmpr++->z;
		qsort(tmparr,n_particles,sizeof(float),compare_float_ascending);
		maxdist = -10.;
		for(i=0;i<n_particles-1;i++){
			tmpdist = tmparr[i+1] - tmparr[i];
			if(tmpdist > maxdist) {
				upbound = i;
				maxdist = tmpdist;
			}
		}
		for(i=0;i<n_particles;i++){
			if(particles[i].z <= tmparr[upbound]) particles[i].z += size;
		}
		/*
		for(i=0;i<=upbound;i++){
			particles[i].z += nz;
		}
		*/
		*zinit = size; *zmax = -size;
		for(i=0;i<n_particles;i++){
			*zinit = MIN(*zinit,particles[i].z);
			*zmax = MAX(*zmax,particles[i].z);
		}
		Free(tmparr);
	}
}
