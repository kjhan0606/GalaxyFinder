/* Iso-density shells. A shell is the set of particles and cores between two density thresholds. */
#include "finder_internal.h"

int count_members_in_radius(SimpleBasicParticleType *particles,int n_particles,int *indexes,
		float cx,float cy,float cz,float radius2){
	int i,j,k;
	int ncount;
	float tmpx,tmpy,tmpz,dist2;
	ncount = 0;
	for(i=0;i<n_particles;i++){
		tmpx = particles[indexes[i]].x-cx; 
		tmpy = particles[indexes[i]].y-cy; 
		tmpz = particles[indexes[i]].z-cz;
		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		if(dist2 < radius2) ncount ++;
	}
	return ncount;
}

int shell_particles_from_remaining(SimpleBasicParticleType  *particles,int n_particles,int shell_id, int *indexes,
		Coretype *one_core){
	int i,j,k,particle_id;
	int n_members,omem;
	float cx,cy,cz,tmpx,tmpy,tmpz,tidalr2,dist2;

	n_members = 0;
	/* Include shell particles that are not bound to any halos. */
	for(i=0;i<halo.shells[shell_id].n_particles;i++)
		if(is_bound((particle_id=halo.shells[shell_id].particle_ids[i]))==NOT) indexes[n_members++] = particle_id;
	for(i=0;i<n_particles;i++) if(is_remaining(i) != NOT) indexes[n_members++] = i;

	cx = one_core->position.x; cy = one_core->position.y; cz = one_core->position.z;
	tidalr2 = one_core->tidal_radius;
	tidalr2 = tidalr2*tidalr2;
	omem = n_members;
	n_members = 0;
	for(i=0;i<omem;i++){
		tmpx = particles[indexes[i]].x -cx;
		tmpy = particles[indexes[i]].y -cy;
		tmpz = particles[indexes[i]].z -cz;
		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		if(dist2 < tidalr2) {
			indexes[n_members++] = indexes[i];
		}
	}
	return n_members;
}

int member_indexes_in_shell(int n_particles,int galaxy_id, int *indexes,int shell_id){
	int i,j,k;
	int n_members,nlcore;

	n_members  = 0;
	for(i=0;i<n_particles;i++) if(halo.particles[i].galaxy_id == galaxy_id) indexes[n_members++] = i;
	for(i=shell_id+1;i<halo.n_shells;i++){
		nlcore = halo.shells[i].n_cores;
		for(j = 0;j<nlcore;j++)
			if(halo.shells[i].core_ids[j] == galaxy_id)
				for(k=0;k<halo.shells[i].n_particles;k++)
					indexes[n_members++] = halo.shells[i].particle_ids[k];

	}
	return n_members;
}

int shell_candidates_from_pool(SimpleBasicParticleType *particles, int shell_id, int *indexes,
		Coretype *one_core, int *remaining_pool, int nremaining){
	int i,particle_id,n_members=0;
	float cx = one_core->position.x, cy = one_core->position.y, cz = one_core->position.z;
	float tidalr2 = one_core->tidal_radius * one_core->tidal_radius;
	float tmpx,tmpy,tmpz,dist2;
	for(i=0;i<halo.shells[shell_id].n_particles;i++){
		particle_id = halo.shells[shell_id].particle_ids[i];
		if(is_bound(particle_id) != NOT) continue;
		tmpx = particles[particle_id].x - cx; tmpy = particles[particle_id].y - cy; tmpz = particles[particle_id].z - cz;
		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		if(dist2 < tidalr2) indexes[n_members++] = particle_id;
	}
	for(i=0;i<nremaining;i++){
		particle_id = remaining_pool[i];
		if(is_remaining(particle_id) == NOT) continue;
		tmpx = particles[particle_id].x - cx; tmpy = particles[particle_id].y - cy; tmpz = particles[particle_id].z - cz;
		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		if(dist2 < tidalr2) indexes[n_members++] = particle_id;
	}
	return n_members;
}

int member_indexes_from_csr(int galaxy_id, int *indexes, int shell_id,
		int *halo_csr_list, int *halo_csr_offset){
	int i,j,k,n_members,nlcore;
	int *src = halo_csr_list + halo_csr_offset[galaxy_id];
	int n = halo_csr_offset[galaxy_id+1] - halo_csr_offset[galaxy_id];
	for(i=0;i<n;i++) indexes[i] = src[i];
	n_members = n;
	for(i=shell_id+1;i<halo.n_shells;i++){
		nlcore = halo.shells[i].n_cores;
		for(j=0;j<nlcore;j++)
			if(halo.shells[i].core_ids[j] == galaxy_id)
				for(k=0;k<halo.shells[i].n_particles;k++)
					indexes[n_members++] = halo.shells[i].particle_ids[k];
	}
	return n_members;
}

int next_unenclosed_core(Coretype  *cores, int n_cores,float threshold){
	int i;
	for(i=0;i<n_cores;i++) {
		if(core_is_enclosed(cores, i) == NOT && cores[i].saddle_density > threshold) {
			return i;
		}
	}
	return -1;
}

void put_unassigned_in_last_shell(int n_particles,Coretype  *cores,int n_cores){
	int i,j,k;
	int nrest,shell_id;

	nrest = 0;
	for(i=0;i<n_particles;i++) if(is_core_particle(i) == NOT && is_shell(i) ==NOT) nrest++;
	if(nrest ==0) return;

	if(halo.n_shells >0) {
		shell_id = halo.n_shells -1;
		if(halo.shells[shell_id].particle_ids) 
			halo.shells[shell_id].particle_ids = Realloc(halo.shells[shell_id].particle_ids,sizeof(int)*(halo.shells[shell_id].n_particles+nrest));
		else halo.shells[shell_id].particle_ids = Malloc(sizeof(int)*nrest,PPTR(halo.shells[shell_id].particle_ids));

		for(i=0;i<n_particles;i++)
			if(is_core_particle(i) == NOT && is_shell(i) ==NOT) {
				mark_shell(i);
				clear_remaining(i);
				halo.shells[shell_id].particle_ids[(halo.shells[shell_id].n_particles++)] = i;
				set_galaxy_id(i,NOT_HALO_MEMBER);
			}
	}
	else if(halo.n_shells ==0){
		halo.shells[halo.n_shells].particle_ids = Malloc(sizeof(int)*nrest,PPTR(halo.shells[halo.n_shells].particle_ids));
		halo.shells[halo.n_shells].n_particles = 0;
		for(i=0;i<n_particles;i++)
			if(is_core_particle(i) == NOT && is_shell(i) ==NOT) {
				mark_shell(i);
				clear_remaining(i);
				set_galaxy_id(i,NOT_HALO_MEMBER);
				halo.shells[halo.n_shells].particle_ids[(halo.shells[halo.n_shells].n_particles++)] = i;
			}
		halo.shells[halo.n_shells].n_cores = n_cores;
		halo.shells[halo.n_shells].core_ids = Malloc(sizeof(int)*n_cores,PPTR(halo.shells[halo.n_shells].core_ids));
		for(i=0;i<n_cores;i++) halo.shells[halo.n_shells].core_ids[i] = i;
		halo.n_shells++;
	}
	else {
		fprintf(stderr,"Strange value of halo.n_shells in put_unassigned_in_last_shell\n");
		fprintf(stderr,"Ok.. Now exiting with halo.n_shells =%d\n",halo.n_shells);
		exit(999);
	}
}

void build_density_shell(int n_particles,long long *neighbor,int n_neighbors,
		Coretype *cores,int n_cores, float dthreshold){
	long long i,j,k,jj,kk;
	int core_id,*contactlist,nlist=0;
	int numenclosedpeaks;
	int now,new;
	int nowcore;

	contactlist = (int *)Malloc(sizeof(int)*n_particles,PPTR(contactlist));
	for(i=0;i<n_cores;i++) reset_core_marks(cores, i);

	for(j=0;j<n_particles;j++) clear_visited(j);

	while((nowcore = next_unenclosed_core(cores,n_cores,dthreshold))>=0){
		for(j=0;j<nlist;j++) {
	 		int jj = contactlist[j];
			clear_visited(jj);// unmarking particles (initialization)
		}
		now = nlist = 0;
		mark_visited((contactlist[nlist++]= cores[nowcore].peak_particle));
		while(now<nlist){
			kk = (long long)contactlist[now]*(long long)n_neighbors;
			for(k=0;k<n_neighbors;k++){
				new = neighbor[kk];
				if(halo.particles[new].density > dthreshold && is_visited(new) == NOT){
					mark_visited(new);
					contactlist[nlist++] = new;
				}
				kk++;
			}
			now++;
		}
		if(nlist==0) continue;
		/* erase shell or cores particles from this shell-particle list
		 * and turn on the shell flag */
		halo.shells[halo.n_shells].particle_ids = Malloc(sizeof(int)*nlist,PPTR(halo.shells[halo.n_shells].particle_ids));
		now = 0;
		for(j=0;j<nlist;j++){
			if(is_shell(contactlist[j])==NOT && is_core_particle(contactlist[j])==NOT){
				halo.shells[halo.n_shells].particle_ids[now++] = contactlist[j];
				mark_shell(contactlist[j]);
				clear_remaining(contactlist[j]);
				set_galaxy_id(contactlist[j],NOT_HALO_MEMBER);
			}
		}
		halo.shells[halo.n_shells].n_particles = now;
		halo.shells[halo.n_shells].particle_ids = Realloc(halo.shells[halo.n_shells].particle_ids,sizeof(int)*halo.shells[halo.n_shells].n_particles);

		halo.shells[halo.n_shells].n_cores = 0;
		halo.shells[halo.n_shells].core_ids = (int *)Malloc(sizeof(int)*n_cores,PPTR(halo.shells[halo.n_shells].core_ids));
		for(j=0;j<nlist;j++){
			now = contactlist[j];
			if(is_peak(now)!=NOT){
				for(core_id=0;core_id<n_cores;core_id++){
					if(cores[core_id].peak_particle == now){
						mark_core_enclosed(cores, core_id);
						halo.shells[halo.n_shells].core_ids[halo.shells[halo.n_shells].n_cores++] = core_id;
						break;
					}
				}
			}
		}
		halo.shells[halo.n_shells].core_ids = Realloc(halo.shells[halo.n_shells].core_ids,sizeof(int)*halo.shells[halo.n_shells].n_cores);
		halo.n_shells ++;
	}
	Free(contactlist);
}

#ifdef SCORE_NMEM
int compare_cores_by_density(const void *a,const void *b){
	Coresortdentype *aa,*bb;
	aa = (Coresortdentype *)a;
	bb = (Coresortdentype *)b;
	if(aa->n_members < bb->n_members) return 1;
	else if(aa->n_members > bb->n_members) return -1;
	else return 0;
}
#else
int compare_cores_by_density(const void *a,const void *b){
	Coresortdentype *aa,*bb;
	aa = (Coresortdentype *)a;
	bb = (Coresortdentype *)b;
	if(aa->density < bb->density) return 1;
	else if(aa->density > bb->density) return -1;
	else return 0;
}
#endif

void return_unbound_shell_particles(int shell_id,SimpleBasicParticleType *particles){
	int i,j,k;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<halo.shells[shell_id].n_particles;i++) 
		if(is_bound(halo.shells[shell_id].particle_ids[i]) == NOT) {
			set_galaxy_id(halo.shells[shell_id].particle_ids[i],NOT_HALO_MEMBER);
			mark_remaining(halo.shells[shell_id].particle_ids[i]);
		}
}
