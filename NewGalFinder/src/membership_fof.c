/* Per-galaxy Friends-of-Friends. OpenMP and serial paths both call the octree in ost.c.
 *
 * The stars-only pass stays on. It is a stellar isodensity cut, not a
 * boundedness test: equal-mass stars percolate near 60-80 times the cosmic
 * mean matter density. Dropping it (halo05 M_inc1) moved the background
 * stars onto the BCG. The link length is each particle's link02 from
 * gfind.c, not FOFLINK4MEMBERSHIP. Paper note: PGalF_compare/docs/02_methods.md.
 */
#include "finder_internal.h"
#include "stage_timer.h"

enum { LINK_STARS_SEPARATELY = 1 };

void link_galaxy_members(SimpleBasicParticleType *particles, int n_particles, int n_cores, Coretype *cores) {
	StageTimer stage_stars, stage_all;
	stage_begin(&stage_stars);
	if (n_particles < NOMPFoF) {
		if (LINK_STARS_SEPARATELY && MINSTELLARMASS >= 0)
			link_members_stars(particles, n_particles, n_cores, cores);
	} else if (LINK_STARS_SEPARATELY && MINSTELLARMASS >= 0) {
		link_members_openmp(particles, n_particles, n_cores, cores, TYPE_STAR);
	}
	stage_end(&stage_stars, "fof_stars");
	stage_begin(&stage_all);
	if (n_particles < NOMPFoF)
		link_members_all_species(particles, n_particles, n_cores, cores);
	else
		link_members_openmp(particles, n_particles, n_cores, cores, TYPE_ALL);
	stage_end(&stage_all, "fof_all");
}

int fof_find_root(ompFoFParticleType *fof_particles,int i){
	if(fof_particles[i].imother != i) fof_particles[i].imother = fof_find_root(fof_particles,fof_particles[i].imother);
	return fof_particles[i].imother;
}

int fof_within_link(ompFoFParticleType *fof_particles, int i,int j){
	POSTYPE tmpx = fof_particles[i].x-fof_particles[j].x;
	POSTYPE tmpy = fof_particles[i].y-fof_particles[j].y;
	POSTYPE tmpz = fof_particles[i].z-fof_particles[j].z;
	POSTYPE dist = sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
	if(dist < 0.5*( fof_particles[i].link02 + fof_particles[j].link02)){
		return 1;
	}
	else {
		return 0;
	}
	return 0;
}

void fof_unite(ompFoFParticleType *fof_particles, int i,int j){
	int root_i = fof_find_root(fof_particles, i);
	int root_j = fof_find_root(fof_particles, j);
	if(root_i != root_j) {
		if(root_i < root_j) fof_particles[root_j].imother = root_i;
		else fof_particles[root_i].imother = root_j;
	}
}

void link_members_openmp(SimpleBasicParticleType *particles, int n_particles, int n_cores, 
		Coretype *cores, int ptype){
	FoFTPtlStruct *fof_particles;
	int i,j,k;
	int core_id;
	float link_length = FOFLINK4MEMBERSHIP;
	fof_particles = (FoFTPtlStruct*) Malloc(sizeof(FoFTPtlStruct)*n_particles, PPTR(fof_particles));
	/* One pass over the halo. Each core then walks only its own members. */
	int *csr_count = (int *)Malloc(sizeof(int)*(size_t)n_cores, PPTR(csr_count));
	int *csr_off = (int *)Malloc(sizeof(int)*(size_t)(n_cores+1), PPTR(csr_off));
	for(i=0;i<n_cores;i++) csr_count[i] = 0;
	for(j=0;j<n_particles;j++){
		int h = halo.particles[j].galaxy_id;
		if(h >= 0 && h < n_cores) csr_count[h]++;
	}
	csr_off[0] = 0;
	for(i=0;i<n_cores;i++) csr_off[i+1] = csr_off[i] + csr_count[i];
	int *csr_at = (int *)Malloc(sizeof(int)*(size_t)n_cores, PPTR(csr_at));
	int *csr_idx = (int *)Malloc(sizeof(int)*(size_t)(csr_off[n_cores] > 0 ? csr_off[n_cores] : 1), PPTR(csr_idx));
	for(i=0;i<n_cores;i++) csr_at[i] = csr_off[i];
	for(j=0;j<n_particles;j++){
		int h = halo.particles[j].galaxy_id;
		if(h >= 0 && h < n_cores) csr_idx[csr_at[h]++] = j;
	}
	Free(csr_at);
	for(core_id=0;core_id<n_cores;core_id++){
		int n_members = 0;
		/* Stellar FoF only drops stars. Dark cores keep any star the
		 * bound test already accepted. */
		if(ptype == TYPE_STAR && cores[core_id].is_dark) continue;
		if(ptype == TYPE_STAR){
			for(k=csr_off[core_id]; k<csr_off[core_id+1]; k++){
				j = csr_idx[k];
				if(particles[j].type == TYPE_STAR){
					fof_particles[n_members].x = particles[j].x;
					fof_particles[n_members].y = particles[j].y;
					fof_particles[n_members].z = particles[j].z;
					fof_particles[n_members].link02 = particles[j].link02;
					fof_particles[n_members].indx = j;
					fof_particles[n_members].group_id = -1;
					n_members++;
				}
			}
		}
		else {
			for(k=csr_off[core_id]; k<csr_off[core_id+1]; k++){
				j = csr_idx[k];
				if(halo.particles[j].galaxy_id == core_id){
					fof_particles[n_members].x = particles[j].x;
					fof_particles[n_members].y = particles[j].y;
					fof_particles[n_members].z = particles[j].z;
					fof_particles[n_members].link02 = particles[j].link02;
					fof_particles[n_members].indx = j;
					fof_particles[n_members].group_id = -1;
					n_members++;
				}
			}
		}
		if(n_members == 0) continue;
        else if(n_members == 1) {
            for(j=0;j<n_members;j++){
                set_galaxy_id(fof_particles[j].indx,NOT_HALO_MEMBER);
                clear_bound(fof_particles[j].indx);
                mark_remaining(fof_particles[j].indx);
			}
			continue;
		}
		/* The octree writes one internal node per split and does not
		 * stop at nnode. 8*n_members overflowed (C297, 145 particles). */
		size_t nnode = (size_t)MAX(65*10000, n_members/2);
		FoFTStruct *TREE;
		particle *linked;
		DEBUGPRINT("Making FoF Tree in the parallel mode for C%d with n_particles= %d & nnode= %zu currentMemStack= %lld\n", 
				core_id, n_members, nnode, CurMemStack());
		TREE = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*nnode,PPTR(TREE));
		linked = (particle *)Malloc(sizeof(particle)*n_members,PPTR(linked));

		int n_accepted, n_new;
		n_accepted = 0;
		int group_id = 0;
#ifdef _OPENMP
#pragma omp parallel private(i,j,k)
#endif
		{
			int thread_id = omp_get_thread_num();
			int n_threads = omp_get_num_threads();
			/* Particles are striped, not cut into contiguous blocks.
			 * Thread t owns t, t+n_threads, ... A contiguous block can be
			 * empty while that stripe still holds particles. */
			if(thread_id < n_members){
				int stride_count = 0;
				for(i=thread_id;i<n_members;i+=n_threads){
					int next = i + n_threads;
					fof_particles[i].sibling = (next < n_members) ? fof_particles + next : NULL;
					fof_particles[i].included = NO;
					stride_count++;
				}
				int nodes_per_thread = (nnode+n_threads-1)/n_threads;
				FoFTStruct *thread_tree = TREE + nodes_per_thread*thread_id;
				build_fof_tree(thread_tree, nodes_per_thread, fof_particles+thread_id, stride_count, SERIALIZED);
			}
		}

		int iloop=0;
		int largest_count = 0;
		int assigned = 0;
		int next_seed = 0;
		do {
			int first_group = (iloop == 0);
			n_new = 0;
			if(iloop ==0){
					linked[0].x = particles[cores[core_id].peak_particle].x;
					linked[0].y = particles[cores[core_id].peak_particle].y;
					linked[0].z = particles[cores[core_id].peak_particle].z;
					linked[0].link02 = particles[cores[core_id].peak_particle].link02;
					n_accepted = 0; n_new = 1;
			}
			else {
				while(next_seed < n_members && fof_particles[next_seed].included != NO) next_seed++;
				if(next_seed < n_members){
					j = next_seed;
					linked[0].x = fof_particles[j].x;
					linked[0].y = fof_particles[j].y;
					linked[0].z = fof_particles[j].z;
					linked[0].link02 = fof_particles[j].link02;
					fof_particles[j].included = YES;
					fof_particles[j].group_id = group_id;
					n_accepted = 0; n_new = 1;
					next_seed++;
				}
			}
			iloop =1;
			if(n_new == 0) break;
			while(n_new){
#ifdef _OPENMP
#pragma omp parallel  private(j)
#endif
				{
					int thread_id = omp_get_thread_num();
		            int n_threads = omp_get_num_threads();
					int nodes_per_thread = (nnode+n_threads-1)/n_threads;
					if(thread_id < n_members){
						FoFTStruct *thread_tree = TREE + nodes_per_thread*thread_id;
						for(j=n_accepted;j<n_accepted+n_new;j++){
							particle p;
							p.x = linked[j].x;
							p.y = linked[j].y;
							p.z = linked[j].z;
							p.link02 = linked[j].link02;
							visit_fof_group(&p, link_length, thread_tree, fof_particles+thread_id);
						}
					}
				}

				n_accepted += n_new;
				n_new = 0;
				for(j=0;j<n_members;j++){
					if(fof_particles[j].included == NEW){
						linked[n_accepted+n_new].x = fof_particles[j].x;
						linked[n_accepted+n_new].y = fof_particles[j].y;
						linked[n_accepted+n_new].z = fof_particles[j].z;
						linked[n_accepted+n_new].link02 = fof_particles[j].link02;
						fof_particles[j].included = YES;
						fof_particles[j].group_id = group_id;
						n_new ++;
					}
				}
				DEBUGPRINT("C%d has %d group of n_new= %d n_accepted= %d for n_members= %d pxyz= %g %g %g\n", core_id, group_id, n_new,n_accepted, n_members,
						particles[cores[core_id].peak_particle].x,
						particles[cores[core_id].peak_particle].y,
						particles[cores[core_id].peak_particle].z
						);
			}
			/* The first queue entry is the peak coordinate, not a labeled
			 * particle. Stop when the unassigned remainder cannot beat
			 * the largest group, the same test as the serial path. */
			{
				int this_group = n_accepted - (first_group ? 1 : 0);
				if(this_group > largest_count) largest_count = this_group;
				assigned += this_group;
			}
			/* The group grown from the peak is this galaxy's isodensity
			 * component. Later groups are separate clumps in the same label;
			 * they are not a reason to keep scanning every singleton. */
			if(first_group){
				group_id = 1;
				break;
			}
			group_id ++;
			if(n_members - assigned < largest_count) break;
		} while(1);



		int jmax = 0;
		int maxcount = 0;
		if(group_id > 0){
			int *counts = (int *)calloc((size_t)group_id, sizeof(int));
			for(i=0;i<n_members;i++){
				int h = fof_particles[i].group_id;
				if(h >= 0 && h < group_id) counts[h]++;
			}
			for(j=0;j<group_id;j++){
				if(counts[j] > maxcount){
					jmax = j;
					maxcount = counts[j];
					if(maxcount > n_members/2) break;
				}
			}
			free(counts);
		}
		for(j=0;j<n_members;j++){
			if(fof_particles[j].group_id != jmax) {
				set_galaxy_id(fof_particles[j].indx,NOT_HALO_MEMBER);
                clear_bound(fof_particles[j].indx);
                mark_remaining(fof_particles[j].indx);
			}
		}
		if(ptype == TYPE_STAR) 
			DEBUGPRINT("C%d 's # of star members changes from %d to %d\n",core_id,n_members,maxcount);
		else
			DEBUGPRINT("C%d 's # of AllType members changes from %d to %d\n",core_id,n_members,maxcount);
		Free(linked);
		Free(TREE);
	}
	Free(csr_idx);
	Free(csr_off);
	Free(csr_count);
	Free(fof_particles);
}

void link_members_stars(SimpleBasicParticleType *particles,int n_particles, int n_cores, Coretype *cores){
	float link_length;
	int i,j,k;
	int n_members;
	FoFTPtlStruct *fof_particles;
	FoFTStruct *TREE;
	particle *linked,p;

	link_length = FOFLINK4MEMBERSHIP;
	fof_particles = (FoFTPtlStruct *) Malloc(sizeof(FoFTPtlStruct)*n_particles,PPTR(fof_particles));
	linked = (particle *)Malloc(sizeof(particle)*n_particles,PPTR(linked));
	size_t nnode = MAX(65*10000,n_particles);
	TREE = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*nnode,PPTR(TREE));
	for(i=0;i<n_cores;i++){
		if(cores[i].is_dark) continue;
		n_members = 0;
		for(j=0;j<n_particles;j++){
			if(halo.particles[j].galaxy_id==i && particles[j].type == TYPE_STAR){
				fof_particles[n_members].type = TYPE_PTL;
				fof_particles[n_members].x = particles[j].x;
				fof_particles[n_members].y = particles[j].y;
				fof_particles[n_members].z = particles[j].z;
				fof_particles[n_members].link02 = particles[j].link02;
				fof_particles[n_members].indx = j;
				fof_particles[n_members].sibling = &fof_particles[n_members+1];
				fof_particles[n_members].included = NO;
				n_members++;
			}
		}
		if(n_members == 0) continue;
		else if(n_members == 1) {
			for(j=0;j<n_members;j++){
				set_galaxy_id(fof_particles[j].indx,NOT_HALO_MEMBER);
				clear_bound(fof_particles[j].indx);
				mark_remaining(fof_particles[j].indx);
			}
			continue;
		}
		fof_particles[n_members-1].sibling = NULL;
		int recursiveflag;
		/*
		if(nnode>65*10000) {
			recursiveflag = PTHREAD;
		}
		else {
			recursiveflag = RECURSIVE;
		}
		*/
		DEBUGPRINT("C%d is now doing the stellar" 
				" FoF with n_members= %d with nnode= %zu\n", 
				i,n_members,nnode);
		recursiveflag = SERIALIZED;
		build_fof_tree(TREE,nnode, fof_particles,n_members,recursiveflag);

		p.x = particles[cores[i].peak_particle].x;
		p.y = particles[cores[i].peak_particle].y;
		p.z = particles[cores[i].peak_particle].z;
		p.link02 = particles[cores[i].peak_particle].link02;
		int n_accepted = collect_fof_group(&p,link_length,TREE,fof_particles,linked);
		DEBUGPRINT("C%d has n_accepted %d in nstar= %d : %g %g %g link02=%g\n", 
				i, n_accepted, n_members, p.x,p.y,p.z,p.link02);
		if(n_accepted < n_members*0.5) {
			int imax=0,mlink=0,n_new;
			/* This should be checked */
			for(k=0;k<n_members;k++) {
				fof_particles[k].sibling = &fof_particles[k+1];
				fof_particles[k].included = NO;
			}
			fof_particles[n_members-1].sibling = NULL;
		    build_fof_tree(TREE,nnode, fof_particles,n_members,recursiveflag);
			/* This should be checked */
            int run_nlink = 0;
			int residual;



			for(j=0;j<n_members;j++){

				if(fof_particles[j].included == YES) continue;

				p.x = fof_particles[j].x;
				p.y = fof_particles[j].y;
				p.z = fof_particles[j].z;
				p.link02 = fof_particles[j].link02;
				/*
				for(k=0;k<n_members;k++) fof_particles[k].included = NO;
				*/

				n_new = collect_fof_group(&p,link_length,TREE,fof_particles,linked);
#ifndef OLD
                run_nlink += n_new;
				residual = n_members-run_nlink;
               
				/* update the maximum linked group */
				if(n_new >= mlink) {
					mlink = n_new;
					imax = j;
				}
				if(residual < mlink){
					break;
				}
#else
				if(n_new >= n_members*0.5) {
					mlink = n_new;
					imax = j;
					break;
				}
				else if(n_new > mlink){
					mlink = n_new;
					imax = j;
				}
#endif
			}
			if(1){
				p.x = fof_particles[imax].x;
				p.y = fof_particles[imax].y;
				p.z = fof_particles[imax].z;
				p.link02 = fof_particles[imax].link02;
				for(k=0;k<n_members;k++) {
					fof_particles[k].sibling = &fof_particles[k+1];
					fof_particles[k].included = NO;
				}
				fof_particles[n_members-1].sibling = NULL;
				build_fof_tree(TREE,nnode, fof_particles,n_members,recursiveflag);

				n_new = collect_fof_group(&p,link_length,TREE,fof_particles,linked);
			}
		}

		k = n_members;
		for(j=0;j<n_members;j++){
			if(fof_particles[j].included == NO){
				set_galaxy_id(fof_particles[j].indx,NOT_HALO_MEMBER);
				clear_bound(fof_particles[j].indx);
				mark_remaining(fof_particles[j].indx);
				k--;
			}
		}
		DEBUGPRINT("C%d 's # of star members changes from %d to %d\n",i,n_members,k);
	}
	Free(TREE);Free(linked);Free(fof_particles);
}

void link_members_all_species(SimpleBasicParticleType *particles,int n_particles, int n_cores, Coretype *cores){
	float link_length;
	int i,j,k;
	int n_members;
	FoFTPtlStruct *fof_particles;
	FoFTStruct *TREE;
	particle *linked,p;

	link_length = FOFLINK4MEMBERSHIP;
	fof_particles = (FoFTPtlStruct *) Malloc(sizeof(FoFTPtlStruct)*n_particles,PPTR(fof_particles));
	linked = (particle *)Malloc(sizeof(particle)*n_particles,PPTR(linked));
	size_t nnode = MAX(65*10000,n_particles);
	TREE = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*nnode,PPTR(TREE));
	for(i=0;i<n_cores;i++){
		n_members = 0;
		for(j=0;j<n_particles;j++){
			if(halo.particles[j].galaxy_id==i){
				fof_particles[n_members].type = TYPE_PTL;
				fof_particles[n_members].x = particles[j].x;
				fof_particles[n_members].y = particles[j].y;
				fof_particles[n_members].z = particles[j].z;
				fof_particles[n_members].link02 = particles[j].link02;
				fof_particles[n_members].indx = j;
				fof_particles[n_members].sibling = &fof_particles[n_members+1];
				fof_particles[n_members].included = NO;
				n_members++;
			}
		}
		if(n_members == 0) continue;
		else if(n_members == 1) {
			for(j=0;j<n_members;j++){
				set_galaxy_id(fof_particles[j].indx,NOT_HALO_MEMBER);
				clear_bound(fof_particles[j].indx);
				mark_remaining(fof_particles[j].indx);
			}
			continue;
		}
		fof_particles[n_members-1].sibling = NULL;

		int recursiveflag;
		/*
		if(nnode>65*10000) {
			recursiveflag = PTHREAD;
		}
		else {
			recursiveflag = RECURSIVE;
		}
		*/
		recursiveflag = SERIALIZED;
		build_fof_tree(TREE,nnode,fof_particles,n_members,recursiveflag);

		p.x = particles[cores[i].peak_particle].x;
		p.y = particles[cores[i].peak_particle].y;
		p.z = particles[cores[i].peak_particle].z;
		p.link02 = particles[cores[i].peak_particle].link02;
		int n_accepted = collect_fof_group(&p,link_length,TREE,fof_particles,linked);
		if(n_accepted < n_members*0.5) {
			int imax=0,mlink=0,n_new;
			/* This should be checked */
			for(k=0;k<n_members;k++) {
				fof_particles[k].sibling = &fof_particles[k+1];
				fof_particles[k].included = NO;
			}
			fof_particles[n_members-1].sibling = NULL;
			build_fof_tree(TREE,nnode,fof_particles,n_members,recursiveflag);
			/* This should be checked */


            int run_nlink = 0;
			int residual;

			for(j=0;j<n_members;j++){

				if(fof_particles[j].included == YES) continue;

				p.x = fof_particles[j].x;
				p.y = fof_particles[j].y;
				p.z = fof_particles[j].z;
				p.link02 = fof_particles[j].link02;
				/*
				for(k=0;k<n_members;k++) fof_particles[k].included = NO;
				*/

				n_new = collect_fof_group(&p,link_length,TREE,fof_particles,linked);
#ifdef DEBUG
#endif

#ifndef OLD
                run_nlink += n_new;
				residual = n_members-run_nlink;
               
				/* update the maximum linked group */
				if(n_new >= mlink) {
					mlink = n_new;
					imax = j;
				}
				if(residual < mlink){
					break;
				}
#else
				if(n_new >= n_members*0.5) {
					mlink = n_new;
					imax = j;
					break;
				}
				else if(n_new > mlink){
					mlink = n_new;
					imax = j;
				}
#endif
			}
			if(mlink <= 5){
				for(j=0;j<n_members;j++) fof_particles[j].included = NO;
			}
			else {
				p.x = fof_particles[imax].x;
				p.y = fof_particles[imax].y;
				p.z = fof_particles[imax].z;
				p.link02 = fof_particles[imax].link02;
				for(k=0;k<n_members;k++) {
					fof_particles[k].sibling = &fof_particles[k+1];
					fof_particles[k].included = NO;
				}
				fof_particles[n_members-1].sibling = NULL;
				build_fof_tree(TREE,nnode,fof_particles,n_members,recursiveflag);
				n_new = collect_fof_group(&p,link_length,TREE,fof_particles,linked);
			}
		}

		k = n_members;
		for(j=0;j<n_members;j++){
			if(fof_particles[j].included == NO){
				set_galaxy_id(fof_particles[j].indx,NOT_HALO_MEMBER);
				clear_bound(fof_particles[j].indx);
				mark_remaining(fof_particles[j].indx);
				k--;
			}
		}
		DEBUGPRINT("C%d 's # of total members changes from %d to %d\n",i,n_members,k);
	}
	Free(TREE);Free(linked);Free(fof_particles);
}
