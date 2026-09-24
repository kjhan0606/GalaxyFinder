/* Density peaks, watershed basins, and persistence pruning. */
#include "finder_internal.h"

int compare_cores_by_peak_density(const void *a, const void *b){
	Coretype *aa,*bb;
	aa = (Coretype *)a;
	bb = (Coretype *)b;
	if(aa->peak_density < bb->peak_density) return 1;
	else if(aa->peak_density > bb->peak_density) return -1;
	else return 0;
}

int compare_cores_by_star_count(const void *a, const void *b){
	Coretype *aa,*bb;
	aa = (Coretype *)a;
	bb = (Coretype *)b;
	if(aa->n_stars < bb->n_stars) return 1;
	else if(aa->n_stars > bb->n_stars) return -1;
	else return 0;
}

int  merge_nearby_peaks(SimpleBasicParticleType *particles,int n_particles,Coretype *cores,int n_cores, int iflag, float fof_link){
	int i,j,k;
	if(fof_link <= 0.f) fof_link = (float)MERGINGPEAKLENGTH;
	FoFTPtlStruct *ptl = (FoFTPtlStruct *) Malloc(sizeof(FoFTPtlStruct)*n_cores,PPTR(ptl));
	particle *linked = (particle *)Malloc(sizeof(particle)*n_cores,PPTR(linked));
	size_t nnode = MAX(65*10000,n_cores);
	FoFTStruct *TREE = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*nnode,PPTR(TREE));
	if(iflag ==0) { 
		qsort(cores,n_cores,sizeof(Coretype),compare_cores_by_peak_density);
	}
	else if(iflag == 1){
		qsort(cores,n_cores,sizeof(Coretype),compare_cores_by_star_count);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<n_cores;i++){
		ptl[i].type = TYPE_PTL;
		ptl[i].x = particles[cores[i].peak_particle].x;
		ptl[i].y = particles[cores[i].peak_particle].y;
		ptl[i].z = particles[cores[i].peak_particle].z;
		ptl[i].link02 = fof_link;
		ptl[i].indx = 1;
		ptl[i].sibling = &ptl[i+1];
		ptl[i].included = NO;
	}
	ptl[n_cores-1].sibling = NULL;

	int recursiveflag;
	if(nnode > 65*10000) recursiveflag= PTHREAD;
	else recursiveflag = RECURSIVE;
	recursiveflag = SERIALIZED;
	DEBUGPRINT("Before cores merging %d\n", n_cores);
	build_fof_tree(TREE,nnode, ptl,n_cores, recursiveflag);

	for(i=0;i<n_cores;i++){
		if(ptl[i].included == NO){
			particle p;
			p.x = ptl[i].x;
			p.y = ptl[i].y;
			p.z = ptl[i].z;
			p.link02 = fof_link;
			int nlink = collect_fof_group_from(&p,fof_link,TREE,ptl,linked);
			ptl[i].included = NO;
		}
	}
	j = 0;
	if(iflag==0){
		for(i=0;i<n_cores;i++){
			if(ptl[i].included == NO){
				cores[j] = cores[i];
				j++;
			}
		}
	}
	else if(iflag==1){
		for(i=0;i<n_cores;i++){
			if(ptl[i].included == NO){
				cores[j] = cores[i];
				mark_peak(cores[i].peak_particle);
				j++;
			}
			else {
				clear_peak(cores[i].peak_particle);
			}
		}
	}

	n_cores = j;

	Free(TREE); Free(linked); Free(ptl);

	return n_cores;
}

int find_sph_density_peaks(float *den,int numneigh,long long *neighbor,int n_particles,
		Coretype **Core, int jflag, SimpleBasicParticleType *particles){
	Coretype *cores = *Core;
	long long i,j,k;
	float *me,*you;
	int iflag;
	int n_cores;
	n_cores = 0;
	DEBUGPRINT("Now before find_sph_density_peaks with %d particles with flag= %d\n", n_particles,jflag);
	if(jflag == 1){
		for(i=0;i<n_particles;i++){
			if(den[i] > PEAKTHRESHOLD && particles[i].type == TYPE_STAR)
			{
				iflag = 1;
				k = i*numneigh;
				for(j=0;j<numneigh;j++){
					if(den[neighbor[k+j]] > den[i]) {
						iflag = 0;
						break;
					}
				}
				if(iflag == 1){
					cores[n_cores].peak_particle = i;
					cores[n_cores].position.x = particles[i].x;
					cores[n_cores].position.y = particles[i].y;
					cores[n_cores].position.z = particles[i].z;
					cores[n_cores].velocity.x = particles[i].vx;
					cores[n_cores].velocity.y = particles[i].vy;
					cores[n_cores].velocity.z = particles[i].vz;
					cores[n_cores].peak_density = den[i];
					if(n_cores >= halo.max_cores-10) {
						halo.max_cores += MAXNUMCORE;
						*Core = Realloc(*Core, sizeof(Coretype)*halo.max_cores);
						cores = *Core;
					}
					/*
						printf("Warning.. insufficient cores: %d with max ncores= %d for i= %d n_particles= %d\n", n_cores, MAXNUMCORE, i, n_particles);
					*/
					n_cores ++;
				}
			}
		}
	}
	else {
		for(i=0;i<n_particles;i++){
			if(den[i] > PEAKTHRESHOLD)
			{
				iflag = 1;
				k = i*numneigh;
				for(j=0;j<numneigh;j++){
					if(den[neighbor[k+j]] > den[i]) {
						iflag = 0;
						break;
					}
				}
				if(iflag == 1){
					if(n_cores >= MAXNUMCORE){
						fprintf(stderr,"Error exceeding the number of cores: %d :: %lld   %d   \n", n_cores, i, n_particles);
						exit(999);
					}
					cores[n_cores].peak_particle = i;
					cores[n_cores].position.x = particles[i].x;
					cores[n_cores].position.y = particles[i].y;
					cores[n_cores].position.z = particles[i].z;
					cores[n_cores].velocity.x = particles[i].vx;
					cores[n_cores].velocity.y = particles[i].vy;
					cores[n_cores].velocity.z = particles[i].vz;
					cores[n_cores].peak_density = den[i];

					if(n_cores >= halo.max_cores-10) {
						halo.max_cores += MAXNUMCORE;
						cores = Realloc(cores, sizeof(Coretype)*halo.max_cores);
						*Core = cores;
					}
					n_cores ++;
				}
			}
		}
	}
	DEBUGPRINT("Now before merging peak. n_cores= %d\n", n_cores);
//	if(n_cores > 10) n_cores = merge_nearby_peaks(particles,n_particles,cores,n_cores,0);
	DEBUGPRINT("Now after merging peak\n");
	return n_cores;
}

void copy_density_into_particle_state(float *density,int n_particles,Coretype *cores,int n_cores){
	int i,now;
	for(i=0;i<n_particles;i++) {
		halo.particles[i].density = density[i];
		clear_particle_marks(i);
	}
	for(i=0;i<n_cores;i++) mark_peak(cores[i].peak_particle);
}

int merge_underpopulated_dark_cores(SimpleBasicParticleType *particles,int n_particles,Coretype *cores,int n_cores,
		long long *neighbor,int n_neighbors){
	int i,j,k;
	int num;
	int *contactlist;
	Coresortdentype *core_order;
	float minden;
	num = 0;
	for(i=0;i<n_cores;i++) if(cores[i].n_particles < MINCORENMEM) {
#ifdef DEBUG
		printf("DMSmartingFinding: %d'th cores is being erased for %d < mincorenemme= %d\n", i,cores[i].n_particles,MINCORENMEM);
#endif
		num++;
	}
#ifdef DEBUG
	printf("merge_underpopulated_dark_cores: The number of erased cores is %d from %d for mincore = %d\n",num,n_cores, MINCORENMEM);
#endif
	if(num <2) {
		num = 0;
		for(i=0;i<n_cores;i++) 
			if(cores[i].n_particles >= MINCORENMEM) cores[num++]=cores[i];
		return num;
	}
	core_order = (Coresortdentype *)Malloc(sizeof(Coresortdentype)*num,PPTR(core_order));
	num = 0;
	/* Dump the cores data to the sorted cores array */
	for(i=0;i<n_cores;i++) if(cores[i].n_particles < MINCORENMEM){
		core_order[num].n_members = cores[i].n_particles;
		core_order[num].galaxy = cores + i;
		core_order[num].density = halo.particles[cores[i].peak_particle].density; /* peak density */
		num++; /* num is the total number of cores that is underpopulated
				  while n_cores is the total number of cores. */
	}
	minden = 2.e23;
	for(i=0;i<n_particles;i++) minden = MIN(minden,halo.particles[i].density);
	/* Sort core_order in reverse order of peak density */
	qsort(core_order,num,sizeof(Coresortdentype),compare_cores_by_density);
	/* Unset the peak flag for the underpopulated cores.
	   Erase underpopulated peaks from peak (cores) list and only consider the sorted cores. */
	for(i=0;i<n_cores;i++) 
		if(cores[i].n_particles < MINCORENMEM) clear_peak(cores[i].peak_particle);
	for(i=0;i<num;i++) {
		reset_order_marks(&core_order[i]);
		reset_order_marks(&core_order[i]);
	}
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
	for(j=0;j<n_particles;j++) {
		clear_visited(j);
	}

	/* Now rock'n roll */
	contactlist = (int *)Malloc(sizeof(int)*n_particles,PPTR(contactlist));
	for(i=0;i<num;i++){/* From the highest density peak that are underpopulated */
		if(order_is_enclosed(&core_order[i]) == NOT){
			float upden = (core_order[i].galaxy)->saddle_density;
			float downden = minden;
			float denthr;
			int ncontact=0, now,new;
			int mcontact;
			float corestarmass = 0;
			int iter = 0;
			while(iter==0 || fabs(upden-downden)/downden>COREDENRESOLUTION){
 				for(j=0;j<ncontact;j++) {
					int jj = contactlist[j];
					clear_visited(jj);/* Set all the particle not to be visited */
				}
				int breakflag = ncontact = now = 0;
				mcontact = 0;

				corestarmass = 0;

				denthr = 0.25*(3.*upden+downden);
				mark_visited((contactlist[ncontact++] = (core_order[i].galaxy)->peak_particle));
				while(now<ncontact && breakflag ==0){
					long long kk = (long long)contactlist[now]*(long long)n_neighbors;
					for(k=0;k<n_neighbors;k++){
						new = neighbor[kk];
						if(halo.particles[new].density>denthr){
							if(is_visited(new) ==NOT && is_peak(new) == NOT){
								mark_visited(new);
								contactlist[ncontact++] = new;
								if(particles[new].type == TYPE_STAR) mcontact ++;
								corestarmass += particles[new].mass;
							}
							else if(is_visited(new) == NOT && is_peak(new) !=NOT){
								downden = denthr;
								breakflag = 1;
								break;
							}
						}
						kk++;
					}
					now++;
				}
				if(breakflag==0) {
					upden = denthr;
					/* This is newly inserted because we don't have to
					 * find the exact value of the threashold. Rather than that
					 * it is sufficient to satisfy the number of cores particles
					 * should be larger than and equal to "MINCORENMEM"
					 * */
					if(mcontact >= MINCORENMEM && corestarmass >= MINCORESTARMASS) break;
				}
				iter ++;
			}
			(core_order[i].galaxy)->saddle_density = (denthr = upden);

			/* Now scoop up cores particles */
//			for(j=0;j<n_particles;j++) clear_visited(j);
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				clear_visited(jj);
			}
			ncontact = now = 0;
			mcontact = 0;

			corestarmass = particles[contactlist[ncontact]].mass;

			mark_visited((contactlist[ncontact++] = (core_order[i].galaxy)->peak_particle));
			while(now<ncontact){
				long long kk = (long long)contactlist[now]*(long long)n_neighbors;
				for(k=0;k<n_neighbors;k++){
					new = neighbor[kk];
					if(halo.particles[new].density > denthr && is_visited(new) == NOT){
						mark_visited(new);
						if(particles[new].type == TYPE_STAR) mcontact ++;
						contactlist[ncontact++] = new;
						corestarmass += particles[new].mass;
					}
					kk++;
				}
				now++;
			}

			if(mcontact >= MINCORENMEM && corestarmass >= MINCORESTARMASS) {/* restore this as a meaningful peak */
#ifdef DEBUG
				printf("New merging cores is detected with n_particles= %d in %d cores\n",ncontact,i );
#endif
				mark_peak((core_order[i].galaxy)->peak_particle);
				mark_order_confirmed(&core_order[i]);
			}
			else {
				clear_peak((core_order[i].galaxy)->peak_particle);
				clear_order_confirmed(&core_order[i]);
			}
			for(j=i+1;j<num;j++){
				k = (core_order[j].galaxy)->peak_particle;
				if(is_visited(k) != NOT ) {
					/* delete this peak forever since there is no possibility 
					 * for this peak to get a meaning. */
					mark_order_enclosed(&core_order[j]);
					clear_order_confirmed(&core_order[j]);
				}
			}
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				clear_visited(jj);
			}
		}
	}
	Free(contactlist);

	for(i=0;i<n_cores;i++) mark_real_core(cores, i);
	for(i=0;i<num;i++) 
		if(order_is_confirmed(&core_order[i]) == NOT) clear_real_core(cores, ((int)(core_order[i].galaxy - cores)));
	num = 0;
	for(i=0;i<n_cores;i++){
		if(core_is_real(cores, i) != NOT){/* If this is real cores */
			cores[num++] = cores[i];
		}
	}
#ifdef DEBUG
	printf("Total number of survived cores is %d from %d\n",num,n_cores);
#endif
	n_cores = num;
	return n_cores;
}

int merge_underpopulated_cores(SimpleBasicParticleType *particles,int n_particles,Coretype *cores,int n_cores,
		long long *neighbor,int n_neighbors){
	long long i,j,k;
	int num;
	int *contactlist;
	Coresortdentype *core_order;
	float minden;
	num = 0;
	for(i=0;i<n_cores;i++) if(cores[i].n_stars < MINCORENMEM || cores[i].star_mass < MINCORESTARMASS) num++;
	DEBUGPRINT("merge_underpopulated_cores: The number of erased cores is %d from %d for mincore = %d\n",num,n_cores, MINCORENMEM);
	if(num <2) {
		num = 0;
		for(i=0;i<n_cores;i++) 
			if(cores[i].n_stars >= MINCORENMEM && cores[i].star_mass >= MINCORESTARMASS) cores[num++]=cores[i];
		return num;
	}
	core_order = (Coresortdentype *)Malloc(sizeof(Coresortdentype)*num,PPTR(core_order));
	num = 0;
	/* Dump the cores data to the sorted cores array */
	for(i=0;i<n_cores;i++) if(cores[i].n_stars < MINCORENMEM || cores[i].star_mass < MINCORESTARMASS){
		core_order[num].n_members = cores[i].n_particles;
		core_order[num].galaxy = cores + i;
		core_order[num].density = halo.particles[cores[i].peak_particle].density; /* peak density */
		DEBUGPRINT("C%lld has num= %d peakden= %g coreden= %g\n", i, cores[i].n_particles,
				halo.particles[cores[i].peak_particle].density, cores[i].saddle_density);
		num++; /* num is the total number of cores that is underpopulated
				  while n_cores is the total number of cores. */
	}
	minden = 2.e23;
	for(i=0;i<n_particles;i++) minden = MIN(minden,halo.particles[i].density);
	/* Sort core_order in reverse order of peak density */
	qsort(core_order,num,sizeof(Coresortdentype),compare_cores_by_density);
	/* Unset the peak flag for the underpopulated cores.
	   Erase underpopulated peaks from peak (cores) list and only consider the sorted cores. */
	for(i=0;i<n_cores;i++) 
		if(cores[i].n_stars < MINCORENMEM || cores[i].star_mass < MINCORESTARMASS) clear_peak(cores[i].peak_particle);
	for(i=0;i<num;i++) {
		reset_order_marks(&core_order[i]);
		reset_order_marks(&core_order[i]);
	}
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
	for(j=0;j<n_particles;j++) {
		clear_visited(j);
	}

	// Now rock'n roll 
	contactlist = (int *)Malloc(sizeof(int)*n_particles,PPTR(contactlist));
	for(i=0;i<num;i++){/* From the highest density peak that are underpopulated */
		if(order_is_enclosed(&core_order[i]) == NOT){
			float upden = (core_order[i].galaxy)->saddle_density;
			float downden = minden;
			float denthr;
			int ncontact=0, now,new;
			int mcontact;
			float corestarmass = 0;
			int iter = 0;
//			while(iter==0 || fabs(upden-downden)/denthr>COREDENRESOLUTION)
			denthr = 0.5*(upden+downden);
			DEBUGPRINT("SC%lld has upden= %g downden= %g denthr= %g res= %g\n", 
						i, upden,downden, denthr, COREDENRESOLUTION);
			do {
 				for(j=0;j<ncontact;j++) {
					int jj = contactlist[j];
					clear_visited(jj);/* Set all the particle not to be visited */
				}
				int breakflag = ncontact = now = 0;
				mcontact = 0;

				corestarmass = 0;

				denthr = 0.5*(upden+downden);
				mark_visited((contactlist[ncontact++] = (core_order[i].galaxy)->peak_particle));
				while(now<ncontact && breakflag ==0){
					long long kk = (long long)contactlist[now]*(long long)n_neighbors;
					for(k=0;k<n_neighbors;k++){
						new = neighbor[kk];
						if(halo.particles[new].density>denthr){
							if(is_visited(new) ==NOT && is_peak(new) == NOT){
								mark_visited(new);
								contactlist[ncontact++] = new;
								if(particles[new].type == TYPE_STAR) {
									mcontact ++; corestarmass += particles[new].mass;
								}
							}
							else if(is_visited(new) == NOT && is_peak(new) !=NOT){
								downden = denthr;
								breakflag = 1;
								break;
							}
						}
						kk++;
					}
					now++;
				}
				if(breakflag==0) {
					upden = denthr;
					/* This is newly inserted because we don't have to
					 * find the exact value of the threashold. Rather than that
					 * it is sufficient to satisfy the number of cores particles
					 * should be larger than and equal to "MINCORENMEM"
					 * */
					if(mcontact >= MINCORENMEM && corestarmass >= MINCORESTARMASS) break;
				}
				iter ++;
			} while(fabs(upden-downden)/denthr>COREDENRESOLUTION 
					&& fabs(minden-(core_order[i].galaxy)->saddle_density)>1.);
			(core_order[i].galaxy)->saddle_density = (denthr = upden);

			// Now scoop up cores particles
//			for(j=0;j<n_particles;j++) clear_visited(j);
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				clear_visited(jj);
			}
			ncontact = now = 0;
			mcontact = 0;

			corestarmass = particles[contactlist[ncontact]].mass;

			mark_visited((contactlist[ncontact++] = (core_order[i].galaxy)->peak_particle));
			while(now<ncontact){
				long long kk = (long long)contactlist[now]*(long long)n_neighbors;
				for(k=0;k<n_neighbors;k++){
					new = neighbor[kk];
					if(halo.particles[new].density > denthr && is_visited(new) == NOT){
						mark_visited(new);
						if(particles[new].type == TYPE_STAR) mcontact ++;
						contactlist[ncontact++] = new;
						corestarmass += particles[new].mass;
					}
					kk++;
				}
				now++;
			}

			if(mcontact >= MINCORENMEM && corestarmass >= MINCORESTARMASS) {
				// restore this as a meaningful peak 
				DEBUGPRINT("New merging cores is detected with n_particles= %d in %lld cores\n",ncontact,i );
				mark_peak((core_order[i].galaxy)->peak_particle);
				mark_order_confirmed(&core_order[i]);
			}
			else {
				clear_peak((core_order[i].galaxy)->peak_particle);
				clear_order_confirmed(&core_order[i]);
			}
			for(j=i+1;j<num;j++){
				k = (core_order[j].galaxy)->peak_particle;
				if(is_visited(k) != NOT ) {
					/* delete this peak forever since there is no possibility 
					 * for this peak to get a meaning. */
					mark_order_enclosed(&core_order[j]);
					clear_order_confirmed(&core_order[j]);
				}
			}
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				clear_visited(jj);
			}
		}
	}
	Free(contactlist);

	for(i=0;i<n_cores;i++) mark_real_core(cores, i);
	for(i=0;i<num;i++) 
		if(order_is_confirmed(&core_order[i]) == NOT) clear_real_core(cores, ((int)(core_order[i].galaxy - cores)));
	num = 0;
	for(i=0;i<n_cores;i++){
		if(core_is_real(cores, i) != NOT){/* If this is real cores */
			cores[num++] = cores[i];
		}
	}
	DEBUGPRINT("Total number of survived cores is %d from %d\n",num,n_cores);
	n_cores = num;
	return n_cores;
}

int  trim_cores_to_saddle(SimpleBasicParticleType *particles,int n_particles,long long *neighbor,
		int n_neighbors,Coretype *cores,int n_cores){
	int i,j,k;
//	float denthr,upden,downden;
	int *Tcontactlist;
	int newnmem,iflag;
	float minden;

	iflag = 0;
recycling:
	for(i=0;i<n_particles;i++) set_galaxy_id(i,NOT_HALO_MEMBER);
	minden = 1.e23;
	for(i=0;i<n_particles;i++) minden = MIN(minden,halo.particles[i].density);
	minden = MAX(DENFLOOR, minden);
	DEBUGPRINT("the minimum density %g \n", minden);

	int nthreads=1;
#ifdef _OPENMP
#pragma omp parallel
	{
		int it = omp_get_thread_num();
		if(it ==0) nthreads = omp_get_num_threads();
	}
	nthreads = MIN(nthreads, MAXTHREADS);
#endif



	size_t numlinkingwatershedding= MIN(n_particles, MAXNUMWATERSHEDDING);

	Tcontactlist = (int *)Malloc(sizeof(int)*numlinkingwatershedding*nthreads,
			PPTR(Tcontactlist));

	/* peak_to_core[p] = cores index whose peak particle is p, else -1.
	 * Built once on the master thread before the parallel watershed loop;
	 * read-only inside the loop so it is race-free. Used by each thread to
	 * record which peer peak it first touched (-> cores[i].merge_into) for
	 * the later H-maxima / persistence-based prune. */
	int *peak_to_core = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(peak_to_core));
	for(j=0;j<n_particles;j++) peak_to_core[j] = -1;
	for(i=0;i<n_cores;i++){
		peak_to_core[cores[i].peak_particle] = i;
		cores[i].merge_into = -1;
	}

#ifdef _OPENMP
#pragma omp parallel for private(j,i)
#endif
	for(j=0;j<n_particles;j++) {
		for(i=0;i<nthreads;i++) clear_visited_thread(j,i);
	}

	for(i=0;i<n_particles;i++) clear_particle_marks(i);
	for(i=0;i<n_cores;i++) mark_peak(cores[i].peak_particle);


	/*
#ifdef _OPENMP
#pragma omp parallel private(i,j,k) num_threads(nthreads)
#endif
	{
		int it = 0;
#ifdef _OPENMP
		it = omp_get_thread_num();
#endif
		int *contactlist = Tcontactlist + (long)it*numlinkingwatershedding;
		for(i=it;i<n_cores;i+=nthreads){
		*/
	{
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) num_threads(nthreads) schedule(dynamic)
#endif
		for(i=0;i<n_cores;i++){
			int it = omp_get_thread_num();
			int *contactlist = Tcontactlist + (long long)it*numlinkingwatershedding;
			// Now find cores density threshold 
			float upden = halo.particles[cores[i].peak_particle].density;
			float downden = minden;
			float denthr;
			int mcontact;
			int ncontact=0, now;
			float _coreres = COREDENRESOLUTION;
			do{
				// initialization before a search for the cores density 
 				for(j=0;j<ncontact;j++) {
					int jj = contactlist[j];
					clear_visited_thread(jj,it);// Set all the particle not to be visited 
				}
				int breakflag = 0; 
				ncontact = now = 0;
				// Trial value 
				denthr = 0.5*(upden+downden); 
				// Now peak particle is included. 
				contactlist[ncontact++] = cores[i].peak_particle;
				mark_visited_thread(cores[i].peak_particle,it); 
				while(now < ncontact && breakflag ==0){
					long long kk = (long long)contactlist[now]*(long long)n_neighbors;
					for(k=0;k<n_neighbors;k++){
						int new = neighbor[kk];
						if(halo.particles[new].density>denthr) {
							if(is_visited_thread(new,it) == NOT && is_peak(new) == NOT){
								mark_visited_thread(new,it);
								contactlist[ncontact++] = new;
							}
							else if(is_visited_thread(new,it) == NOT && is_peak(new) != NOT){
								/* Record the peer peak we first touched at
								 * this (sub-saddle) density. As bisection
								 * tightens denthr toward the true saddle,
								 * later iterations overwrite this with the
								 * peer that survives the highest threshold,
								 * which is the immediate saddle neighbour. */
								cores[i].merge_into = peak_to_core[new];
								breakflag = 1;
								break;
							}
						}
						kk ++;
					}
					now++;
				}
				if(breakflag==0) upden = denthr;
				else if(breakflag==1) downden = denthr;
				/* Adaptive resolution: when denthr is within 1% of peak density
				 * (peak tightly squeezed by neighbors), tighten the relative tolerance
				 * 100x (1e-3 -> 1e-5) to better resolve the watershed boundary. */
				{
					float _peakden = halo.particles[cores[i].peak_particle].density;
					_coreres = (_peakden > 0 && denthr / _peakden > 0.99)
						? COREDENRESOLUTION * 0.01 : COREDENRESOLUTION;
				}
			}while(fabs((upden-downden)/denthr) > _coreres);
			cores[i].saddle_density = (denthr = upden);
			/* Now scoop up cores particles */
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				clear_visited_thread(jj,it);
			}
			{
				int breakflag = 0; 
				ncontact = now = 0;
				// Trial value 
				// Now peak particle is included. 
				contactlist[ncontact++] = cores[i].peak_particle;
				mark_visited_thread(cores[i].peak_particle,it); 
				while(now < ncontact){
					long long kk = (long long)contactlist[now]*(long long)n_neighbors;
					for(k=0;k<n_neighbors;k++){
						int new = neighbor[kk];
						if(halo.particles[new].density>denthr) {
							if(is_visited_thread(new,it) == NOT && is_peak(new) == NOT){
								mark_visited_thread(new,it);
								contactlist[ncontact++] = new;
							}
							else if(is_visited_thread(new,it) == NOT && is_peak(new) != NOT){
								/* Backup capture: if bisection converged at
								 * the first try (no prior breakflag), we still
								 * see a peer here. */
								if(cores[i].merge_into < 0)
									cores[i].merge_into = peak_to_core[new];
								breakflag = 1;
								break;
							}
						}
						kk ++;
					}
					now++;
				}
				if(breakflag ==1) { //no coredensity is found. discard this peak.
					/* Clear marks before dropping the list. Zeroing ncontact
					 * first left every flooded particle visited, so a later
					 * core on this thread could not step onto them. */
					for(j=0;j<ncontact;j++)
						clear_visited_thread(contactlist[j],it);
					ncontact = 0;
				}
			}
			/* Now scoop up cores particles */
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				clear_visited_thread(jj,it);
			}

			/*
			ncontact = now = 0;
			float corestarmass = 0;
			mcontact = 0;
			if(particles[cores[i].peak_particle].type == TYPE_STAR) {
				mcontact ++;
				corestarmass += particles[cores[i].peak_particle].mass;
			}
			mark_visited_thread((contactlist[ncontact++] = cores[i].peak_particle),it);
			while(now<ncontact){
				long long kk = (long long)contactlist[now]*(long long)n_neighbors;
				for(k=0;k<n_neighbors;k++){
					int new = neighbor[kk];
					if(halo.particles[new].density > denthr && is_visited_thread(new,it) == NOT){

						if(particles[new].type== TYPE_STAR) {
							mcontact ++;
							corestarmass += particles[new].mass;
						}

						mark_visited_thread(new,it);
						contactlist[ncontact++] = new;
					}
					kk++;
				}
				now++;
			}
			*/

			mcontact = 0;
			float corestarmass = 0;
			for(j=0;j<ncontact;j++){
				int jj = contactlist[j];
				if(particles[jj].type == TYPE_STAR){
					mcontact ++;
					corestarmass += particles[jj].mass;
				}
			}

			for(j=0;j<ncontact;j++) {// make sure these lines are thread-safe
				mark_core_particle(contactlist[j]);
				mark_bound(contactlist[j]);
				clear_remaining(contactlist[j]);
				set_galaxy_id(contactlist[j],i);
			}

			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				clear_visited_thread(jj,it);
			}

			/* Now cores[i].center is temporarily considered as the center for tidal radius calculation.*/
			/* But in case of kinetic energy measurement the center of mass is instead used. */
			if(iflag != 0)
			{
				double cx,cy,cz;
				double cvx,cvy,cvz;
				dptype tmass=0,amass=0;
				cx = cy = cz = cvx = cvy = cvz = 0;
				if(MINCORESTARMASS>0){
					for(j=0;j<ncontact;j++){
						if(particles[contactlist[j]].type == TYPE_STAR){
							amass = particles[contactlist[j]].mass;
							tmass += amass;
							cx += particles[contactlist[j]].x*amass;
							cy += particles[contactlist[j]].y*amass;
							cz += particles[contactlist[j]].z*amass;
							cvx += particles[contactlist[j]].vx*amass;
							cvy += particles[contactlist[j]].vy*amass;
							cvz += particles[contactlist[j]].vz*amass;
						}
					}
				}
				else {
					for(j=0;j<ncontact;j++){
						amass = particles[contactlist[j]].mass;
						tmass += amass;
						cx += particles[contactlist[j]].x*amass;
						cy += particles[contactlist[j]].y*amass;
						cz += particles[contactlist[j]].z*amass;
						cvx += particles[contactlist[j]].vx*amass;
						cvy += particles[contactlist[j]].vy*amass;
						cvz += particles[contactlist[j]].vz*amass;
					}
				}
				cx = cx/tmass;
				cy = cy/tmass;
				cz = cz/tmass;
				cvx = cvx/tmass;
				cvy = cvy/tmass;
				cvz = cvz/tmass;
				cores[i].position.x = cx; cores[i].position.y = cy; cores[i].position.z = cz;
				cores[i].velocity.x = cvx; cores[i].velocity.y = cvy; cores[i].velocity.z = cvz;
			}
			/*
			if(i==1){ 
				FILE *wwp = fopen("C1.out","w");
				for(j=0;j<ncontact;j++) fprintf(wwp,"%g %g %g %g %d\n",
						particles[contactlist[j]].x/TSC_CELL_SIZE+NCELLBUFF/2.,
						particles[contactlist[j]].y/TSC_CELL_SIZE+NCELLBUFF/2.,
						particles[contactlist[j]].z/TSC_CELL_SIZE+NCELLBUFF/2.,
						halo.particles[contactlist[j]].density,
						contactlist[j]
						);
				fclose(wwp);
			}
			if(i==697){ 
				FILE *wwp = fopen("C697.out","w");
				for(j=0;j<ncontact;j++) fprintf(wwp,"%g %g %g %g %d\n",
						particles[contactlist[j]].x/TSC_CELL_SIZE+NCELLBUFF/2.,
						particles[contactlist[j]].y/TSC_CELL_SIZE+NCELLBUFF/2.,
						particles[contactlist[j]].z/TSC_CELL_SIZE+NCELLBUFF/2.,
						halo.particles[contactlist[j]].density,
						contactlist[j]
						);
				fclose(wwp);
			}
			*/
			cores[i].n_particles = ncontact;
			cores[i].n_stars = mcontact;
			cores[i].star_mass = corestarmass;
			float xpix,ypix,zpix;
			xpix = cores[i].position.x/TSC_CELL_SIZE+NCELLBUFF/2.;
			ypix = cores[i].position.y/TSC_CELL_SIZE+NCELLBUFF/2.;
			zpix = cores[i].position.z/TSC_CELL_SIZE+NCELLBUFF/2.;
			DEBUGPRINT("C%d has number %d in %d with mstar= %g den= %g/%g at %g %g %g :: %g %g %g\n",i,
					cores[i].n_stars, cores[i].n_particles, cores[i].star_mass,
					cores[i].peak_density, cores[i].saddle_density,
					cores[i].position.x, cores[i].position.y, cores[i].position.z, xpix,ypix,zpix);
		}
	}
	Free(Tcontactlist);
	Free(peak_to_core);
	DEBUGPRINT("Now Found the cores densities for %d cores\n",n_cores);
	for(i=0;i<n_cores;i++){
		DEBUGPRINT("C%d has nstar= %d npall= %d xyz= %g %g %g\n", 
				i, cores[i].n_stars,cores[i].n_particles, cores[i].position.x, cores[i].position.y, cores[i].position.z);
	}
	if(iflag ==0){
		// Deleting peaks of cores particles less than minimum number 
		iflag = 1;
#ifdef OLD
		newnmem  =0;
		for(i=0;i<n_cores;i++){
			if(cores[i].n_stars >= MINCORENMEM && cores[i].star_mass >= MINCORESTARMASS) cores[newnmem++] = cores[i];
		}
		n_cores = newnmem;
		for(i=0;i<n_particles;i++) clear_particle_marks(i);
		for(i=0;i<n_cores;i++) mark_peak(cores[i].peak_particle);
#else
//		if(MINSTELLARMASS >0) n_cores = merge_underpopulated_cores(particles,n_particles,cores,n_cores,neighbor,n_neighbors);
//		else  n_cores = merge_underpopulated_dark_cores(particles,n_particles,cores,n_cores,neighbor,n_neighbors);
//		if(MINSTELLARMASS<=0) n_cores = merge_underpopulated_dark_cores(particles,n_particles,cores,n_cores,neighbor,n_neighbors);

//		DEBUGPRINT("Now before merging peak\n");
//		if(n_cores > 10) n_cores = merge_nearby_peaks(particles,n_particles,cores,n_cores,1);
//		DEBUGPRINT("Now after merging peak\n");

#endif

//		goto recycling; // go and restart again .
	}

	/* Re-cull cores whose refined ncontact (overwritten cores[i].n_particles above)
	 * fell below MINCORENMEM. Without this, breakflag=1 cases (peak squeezed
	 * by neighbor in bisected watershed) leave cores with nummem=0..few in
	 * the array, wasting shell-loop work and breaking MINCORENMEM invariant. */
	{
		int *_remap = (int *)Malloc(sizeof(int)*n_cores, PPTR(_remap));
		int _newnum = 0;
		for(i=0;i<n_cores;i++) {
			if(cores[i].n_particles >= MINCORENMEM) {
				_remap[i] = _newnum;
				if(_newnum != i) cores[_newnum] = cores[i];
				_newnum++;
			} else {
				_remap[i] = -1;
				clear_peak(cores[i].peak_particle);
			}
		}
		if(_newnum < n_cores) {
			for(j=0;j<n_particles;j++) {
				if(halo.particles[j].galaxy_id >= 0) {
					int _newid = _remap[halo.particles[j].galaxy_id];
					if(_newid < 0) {
						set_galaxy_id(j, NOT_HALO_MEMBER);
						clear_core_particle(j);
						clear_bound(j);
					} else if(_newid != halo.particles[j].galaxy_id) {
						set_galaxy_id(j, _newid);
					}
				}
			}
			/* Re-map merge_into through the same _remap so PrunePersistence
			 * sees valid post-cull peer indices. A peer that got culled
			 * (mapped to -1) loses its merge_into pointer for the survivor. */
			for(j=0;j<_newnum;j++) {
				int mi = cores[j].merge_into;
				cores[j].merge_into = (mi >= 0) ? _remap[mi] : -1;
			}
			DEBUGPRINT("trim_cores_to_saddle re-cull: %d -> %d cores (MINCORENMEM=%d)\n",
					n_cores, _newnum, MINCORENMEM);
			n_cores = _newnum;
		}
		Free(_remap);
	}

	return n_cores;
}

int _pp_cmp_ascpeak(const void *a, const void *b) {
    int ia = *(const int *)a, ib = *(const int *)b;
    float da = halo.particles[halo.persistence_order[ia].peak_particle].density;
    float db = halo.particles[halo.persistence_order[ib].peak_particle].density;
    if(da < db) return -1;
    if(da > db) return  1;
    return 0;
}

int _pp_find(int *p, int x) {
    while(p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

int prune_shallow_peaks(int n_particles, Coretype *cores, int n_cores,
                                  SimpleBasicParticleType *particles, float tau)
{
    int i, j;
    if(tau <= 0.f || n_cores <= 1) return n_cores;

    int *parent = (int *)Malloc(sizeof(int)*(size_t)n_cores, PPTR(parent));
    int *order  = (int *)Malloc(sizeof(int)*(size_t)n_cores, PPTR(order));
    for(i=0;i<n_cores;i++){ parent[i] = i; order[i] = i; }

    halo.persistence_order = cores;
    qsort(order, (size_t)n_cores, sizeof(int), _pp_cmp_ascpeak);

    int n_merge_events = 0;
    for(j=0;j<n_cores;j++){
        int c = order[j];
        int peer = cores[c].merge_into;
        if(peer < 0 || peer >= n_cores) continue;
        int rc = _pp_find(parent, c);
        int rp = _pp_find(parent, peer);
        if(rc == rp) continue;
        /* the lower-peak root is the merge candidate (child). */
        float dc = halo.particles[cores[rc].peak_particle].density;
        float dp = halo.particles[cores[rp].peak_particle].density;
        int child, parnt;
        if(dc <  dp) { child = rc; parnt = rp; }
        else         { child = rp; parnt = rc; }
        float peak_d   = halo.particles[cores[child].peak_particle].density;
        float saddle_d = cores[child].saddle_density;
        if(peak_d <= 0.f) continue;
        float prom = (peak_d - saddle_d) / peak_d;
        if(prom < tau){
            parent[child] = parnt;
            /* The merged object touches its neighbour at the child's saddle.
             * Leaving the parent's higher saddle drops those particles. */
            if(cores[child].saddle_density < cores[parnt].saddle_density)
                cores[parnt].saddle_density = cores[child].saddle_density;
            n_merge_events++;
        }
    }

    /* Compact: surviving roots get new indices 0..new_n-1, in original order. */
    int *new_id = (int *)Malloc(sizeof(int)*(size_t)n_cores, PPTR(new_id));
    for(i=0;i<n_cores;i++) new_id[i] = -1;
    int new_n = 0;
    for(i=0;i<n_cores;i++){
        if(_pp_find(parent, i) == i) new_id[i] = new_n++;
    }
    for(i=0;i<n_cores;i++){
        if(new_id[i] < 0) new_id[i] = new_id[_pp_find(parent, i)];
    }

    /* Clear WP_PEAK on absorbed peaks. */
    for(i=0;i<n_cores;i++){
        if(_pp_find(parent, i) != i) clear_peak(cores[i].peak_particle);
    }

    /* Remap particle membership. */
    for(i=0;i<n_particles;i++){
        int h = halo.particles[i].galaxy_id;
        if(h >= 0 && h < n_cores){
            halo.particles[i].galaxy_id = new_id[h];
        }
    }

    if(new_n < n_cores){
        Coretype *tmp = (Coretype *)Malloc(sizeof(Coretype)*(size_t)new_n, PPTR(tmp));
        for(i=0;i<n_cores;i++){
            if(_pp_find(parent, i) == i) tmp[new_id[i]] = cores[i];
        }
        for(i=0;i<new_n;i++){
            cores[i] = tmp[i];
            cores[i].n_particles   = 0;
            cores[i].n_stars  = 0;
            cores[i].star_mass = 0.f;
            cores[i].merge_into = -1; /* stale post-merge; not used downstream */
        }
        Free(tmp);

        /* Rebuild per-cores particle stats from halo.particles[].galaxy_id. */
        for(i=0;i<n_particles;i++){
            int h = halo.particles[i].galaxy_id;
            if(h >= 0 && h < new_n && is_core_particle(i)){
                cores[h].n_particles++;
                if(particles[i].type == TYPE_STAR){
                    cores[h].n_stars++;
                    cores[h].star_mass += particles[i].mass;
                }
            }
        }
    }

    DEBUGPRINT("prune_shallow_peaks: %d -> %d cores (tau=%.3f, %d merge events)\n",
               n_cores, new_n, (double)tau, n_merge_events);

    Free(new_id); Free(order); Free(parent);
    return new_n;
}
