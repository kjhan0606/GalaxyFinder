/* One FoF halo in, galaxy ids out. This is the only entry the MPI driver calls. */
#include "finder_internal.h"
#include "stage_timer.h"

int find_galaxies(FoFTPtlStruct *raw, lint n_particles,lint *galaxy_ids){
	int i,j;
	if(halo.max_cores < MAXNUMCORE) halo.max_cores = MAXNUMCORE;
	float xinit,yinit,zinit;
	float xmax,ymax,zmax;
	long long *neighbor = NULL;
	int n_neighbors,n_cores;
	float *density = NULL;
	Coretype *cores;
	SimpleBasicParticleType *particles;

	particles = (SimpleBasicParticleType*)Malloc(sizeof(SimpleBasicParticleType)*n_particles,PPTR(particles));
	double Mvir = copy_raw_particles(raw,particles,n_particles);
	double halobgpotent = background_potential(Mvir); 

	for(i=0;i<n_particles;i++){
		galaxy_ids[i] = NOT_HALO_MEMBER;
	}
//	n_neighbors = MIN(MINCORENMEM,n_particles);
	n_neighbors = MIN(NUMNEIGHBOR,n_particles);
	if(1){
		unwrap_periodic_box(particles,n_particles,&xinit,&yinit, &zinit,&xmax,&ymax,&zmax);
	}
	{
		int    nstar_gate  = count_stars(particles,n_particles);
		float  mstar_gate  = total_star_mass(particles,n_particles);
		float  mdm_gate    = total_dm_mass(particles,n_particles);
		int    do_grid     = ((nstar_gate > NUMNEIGHBOR) && (mstar_gate >= MINSTELLARMASS))
		                     || (mdm_gate >= (float)MINDMMASS);
		/* A DMO halo cannot produce a stellar peak when DM_DENSITY_WEIGHT is
		 * zero. Skip the empty 3-kpc stellar FFT and enter the dedicated,
		 * adaptively capped DMO peak finder below. */
		if(nstar_gate == 0 && mdm_gate > 0.f){
			cores = (Coretype*)Malloc(sizeof(Coretype)*halo.max_cores,PPTR(cores));
			n_cores = 0;
		}
		else if(!do_grid){
			neighbor = (long long*)Malloc(sizeof(long long)*n_particles*n_neighbors,PPTR(neighbor));
			density = (float*)Malloc(sizeof(float)*n_particles,PPTR(density));
			void findsphdensity(SimpleBasicParticleType *,int ,long long *, int , float *);
			findsphdensity(particles,n_particles,neighbor,n_neighbors,density);
			cores = (Coretype*)Malloc(sizeof(Coretype)*halo.max_cores,PPTR(cores));
			n_cores = find_sph_density_peaks(density,n_neighbors,neighbor,n_particles,&cores,0,particles);
		}
		else {
			neighbor = (long long*)Malloc(sizeof(long long)*n_particles*n_neighbors,PPTR(neighbor));
			density = (float*)Malloc(sizeof(float)*n_particles,PPTR(density));
			cores = (Coretype*)Malloc(sizeof(Coretype)*halo.max_cores,PPTR(cores));
#ifdef ADV
			/* Unified star+DM density. DM_DENSITY_WEIGHT == 0 short-circuits
			 * the DM grid in find_peaks_with_thresholds, making this byte-equivalent
			 * to the old stellar-only call. */
			void find_density_peaks(SimpleBasicParticleType *, int, int, float *,
					Coretype **, int *, int, long long *);
			find_density_peaks(particles,n_particles,n_neighbors,density, &cores, &n_cores, halo.max_cores,
					neighbor);
#else
			neighbor = (long long*)Malloc(sizeof(long long)*n_particles*(long)n_neighbors,PPTR(neighbor));
			void starfindsphdensity(SimpleBasicParticleType *,int ,long long *, int , float *);
			starfindsphdensity(particles,n_particles,neighbor,n_neighbors,density);
			n_cores = find_sph_density_peaks(density,n_neighbors,neighbor,n_particles,&cores,1,particles);
#endif
		}
		DEBUGPRINT("density calculates\n");
	}
	{
		DEBUGPRINT("%d n_cores detected\n",n_cores);
		for(i=0;i<n_cores;i++){
			DEBUGPRINT("C%d has ipeak= %d\n", i, cores[i].peak_particle);
			cores[i].is_dark = 0; /* default; dark cores get this set explicitly later */
		}
		halo.particles = (ParticleState *)Malloc(sizeof(ParticleState)*n_particles,PPTR(halo.particles));
		memset(halo.particles, 0, sizeof(ParticleState)*(size_t)n_particles);
	}
	if(n_cores == 0) {
		int nstar_now = count_stars(particles, n_particles);
		if(nstar_now == 0 && total_dm_mass(particles, n_particles) > 0.f){
			DEBUGPRINT("DMO halo: no stellar peaks, seeding on dark matter\n");
			for(i=0;i<n_particles;i++) set_galaxy_id(i, NOT_HALO_MEMBER);
			n_cores = assign_dmo_halos(particles, (int)n_particles, cores);
			if(density != NULL) Free(density);
			if(n_cores <= 0){
				for(i=0;i<n_particles;i++) galaxy_ids[i] = 0;
				Free(halo.particles);
				if(neighbor != NULL) Free(neighbor);
				Free(cores);
				Free(particles);
				return 1;
			}
			goto gogo;
		}
		for(i=0;i<n_particles;i++) galaxy_ids[i] = 0;
		return 1;
	}
renumcore :
	if(n_cores ==1) {
		int n_members;
		/*
		Free(neighbor);

		samples = (SampledParticle*)Malloc(sizeof(SampledParticle)*n_particles,PPTR(samples));
		for(i=0;i<n_particles;i++){
			samples[i].position.x = particles[i].x; samples[i].position.y = particles[i].y; samples[i].position.z = particles[i].z;
			samples[i].velocity.x = particles[i].vx; samples[i].velocity.y = particles[i].vy; samples[i].velocity.z = particles[i].vz;
			samples[i].mass = particles[i].mass;
			samples[i].particle = particles+i;
		}
		for(i=0;i<n_particles;i++) set_galaxy_id(i,NOT_HALO_MEMBER);
		n_members = unbind_isolated_halo(samples,n_particles,particles,0);
		Free(samples);

		link_members_all_species(particles,n_particles,n_cores,cores);
		for(i=0;i<n_particles;i++){
#ifdef NOBACKGROUND
			galaxy_ids[i] = 0;
#else
			if(halo.particles[i].galaxy_id==0) galaxy_ids[i] = 0;
			else galaxy_ids[i] = NOT_HALO_MEMBER;
#endif
		}

		Free(halo.particles);
		*/
		for(i=0;i<n_particles;i++) galaxy_ids[i] = 0;
		return 1;
	}
	else if(n_cores > 1) {
		copy_density_into_particle_state(density,n_particles,cores,n_cores);

		{
			StageTimer stage_trim;
			stage_begin(&stage_trim);
			n_cores = trim_cores_to_saddle(particles,n_particles,neighbor,n_neighbors,cores,n_cores);
			stage_end(&stage_trim, "watershed");
		}
		Free(density);

		/* H-maxima / persistence-based prune of weak watershed peaks.
		 * Suppresses BCG / cluster-centre over-segmentation by folding cores
		 * whose (peak-saddle)/peak < PERSISTENCE_TAU into the higher peer they
		 * first touched in trim_cores_to_saddle's bisection. PERSISTENCE_TAU<=0
		 * keeps the pre-prune behaviour. Runs BEFORE the "watershed" stage
		 * dump so downstream stages see the pruned cores list. */
		if((float)PERSISTENCE_TAU > 0.f && n_cores > 1){
			StageTimer stage_persist;
			stage_begin(&stage_persist);
			n_cores = prune_shallow_peaks(n_particles, cores, n_cores, particles,
					(float)PERSISTENCE_TAU);
			stage_end(&stage_persist, "persistence");
		}

		/* DEBUG: stage "watershed" dump — cores immediately after
		 * trim_cores_to_saddle's MINCORENMEM re-cull and the persistence prune,
		 * before any shell-loop / boundedness / membership-FoF
		 * post-processing.  Gated by NEWGAL_DUMP_STAGES env var. */
		write_stage_snapshot(raw, n_particles, n_cores, cores, "watershed");
		if(n_cores ==1) {
			goto renumcore;
		}
		else if(n_cores ==0) {
			for(i=0;i<n_particles;i++) galaxy_ids[i] = 0;
			return 1;
			/*
			Free(halo.particles);Free(density);Free(neighbor);
			Free(particles);
			return n_cores;
			*/
		}
		DEBUGPRINT("total number of cores changes to %d\n",n_cores);fflush(stdout);

		/* Post-hoc dark-galaxy classification on the unified-density cores.
		 * find_density_peaks peaks at rho_star + DM_DENSITY_WEIGHT*rho_DM, so a
		 * single cores list contains both luminous and pure-DM peaks. A cores
		 * is flagged dark if its stellar content is below the FoF-gate scale
		 * (too few stars or too little stellar mass to be a galaxy). When
		 * DM_DENSITY_WEIGHT==0 no pure-DM peak rises above PEAKTHRESHOLD,
		 * so this loop simply marks nothing. */
		{
			int n_dark = 0;
			for(i=0;i<n_cores;i++){
				if(cores[i].n_stars < MINCORENMEM
						|| cores[i].star_mass < (float)MINSTELLARMASS){
					cores[i].is_dark = 1;
					n_dark++;
				}
			}
			DEBUGPRINT("dark-galaxy classification: %d/%d cores marked dark\n",
					n_dark, n_cores);
		}

		n_cores = assign_members_from_watershed(particles, (int)n_particles, neighbor, n_neighbors, cores, n_cores);
		write_stage_snapshot(raw, n_particles, n_cores, cores, "post_shell");
	}
gogo:
		for(i=0;i<n_particles;i++){
			galaxy_ids[i] = halo.particles[i].galaxy_id;
		}

		int *hid_num = (int*)Malloc(sizeof(int)*n_cores, PPTR(hid_num));
		int *hid_new = (int*)Malloc(sizeof(int)*n_cores, PPTR(hid_new));
		for(i=0;i<n_cores;i++){
			hid_num[i] = 0;
			hid_new[i] = -1;
		}
		for(i=0;i<n_particles;i++){
			if(halo.particles[i].galaxy_id >=0) hid_num[halo.particles[i].galaxy_id]++;
		}
		int jpeak = 0;
		int numpeak = 0;
		j= 0;
		for(i=0;i<n_cores;i++){
			if(hid_num[i] > numpeak){
				numpeak = hid_num[i];
				jpeak = i;
			}
			if(hid_num[i] >0) {
				hid_new[i] = j;
				j++;
			}
		}
		n_cores = j;
		for(i=0;i<n_particles;i++){
			if(galaxy_ids[i]!= NOT_HALO_MEMBER) galaxy_ids[i] = hid_new[galaxy_ids[i]];
#ifdef NOBACKGROUND
			else galaxy_ids[i] = hid_new[jpeak];
#endif
		}


		Free(hid_new);Free(hid_num);



		Free(halo.particles);
		if(neighbor != NULL) Free(neighbor);
	/*
	*/

	Free(cores);
	Free(particles);
	return n_cores;
}
