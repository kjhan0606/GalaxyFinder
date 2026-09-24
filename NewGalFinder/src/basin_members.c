/* Membership from the watershed basins.
 *
 * Peaks and saddles stay as trim_cores_to_saddle left them. Each particle
 * keeps the peak it climbs to, so a parallel flood overwrite cannot hand a
 * satellite particle to a neighbour, and a failed climb keeps the flood id.
 * Particles below every saddle are not members yet. They join a galaxy only
 * when they sit inside that galaxy's Roche sphere and are bound to its own
 * members. The host halo mass is not a potential source, and the NFW
 * background floor is not added on this test: that floor binds every slow
 * particle out to the box edge.
 */
#define _GNU_SOURCE
#include "finder_internal.h"
#include "stage_timer.h"

static int basin_owner_cmp(const void *a, const void *b) {
	float ka = ((const float *)a)[1];
	float kb = ((const float *)b)[1];
	if(ka < kb) return -1;
	if(ka > kb) return 1;
	return 0;
}

static int climb_core(int start, long long *neighbor, int n_neighbors,
		int *peak_to_core, int *dest, int *stack, int *seen, int generation,
		int n_particles) {
	int nstack = 0;
	int owner = -1;
	int p = start;
	while(1){
		int k, best;
		float bestd;
		long long kk;
		if(p < 0 || p >= n_particles){
			owner = -1;
			break;
		}
		if(dest[p] != -2){
			owner = dest[p];
			break;
		}
		if(seen[p] == generation){
			owner = -1;
			break;
		}
		seen[p] = generation;
		stack[nstack++] = p;
		if(is_peak(p)){
			owner = peak_to_core[p];
			break;
		}
		best = -1;
		bestd = halo.particles[p].density;
		kk = (long long)p * (long long)n_neighbors;
		for(k=0;k<n_neighbors;k++){
			int nb = (int)neighbor[kk + k];
			float den;
			if(nb < 0 || nb >= n_particles) continue;
			den = halo.particles[nb].density;
			if(den > bestd){
				bestd = den;
				best = nb;
			}
		}
		if(best < 0){
			owner = -1;
			break;
		}
		p = best;
	}
	{
		int i;
		for(i=0;i<nstack;i++) dest[stack[i]] = owner;
	}
	return owner;
}

static void assign_uphill_basins(SimpleBasicParticleType *particles, int n_particles,
		long long *neighbor, int n_neighbors, Coretype *cores, int n_cores) {
	int i, generation;
	int *peak_to_core, *dest, *stack, *seen, *flood;
	int n_kept = 0, n_moved = 0, n_released = 0, n_gained = 0;
	(void)particles;
	peak_to_core = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(peak_to_core));
	dest = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(dest));
	stack = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(stack));
	seen = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(seen));
	flood = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(flood));
	for(i=0;i<n_particles;i++){
		int h = halo.particles[i].galaxy_id;
		peak_to_core[i] = -1;
		dest[i] = -2;
		seen[i] = 0;
		flood[i] = (h >= 0 && h < n_cores) ? h : -1;
	}
	for(i=0;i<n_cores;i++){
		int p = cores[i].peak_particle;
		if(p >= 0 && p < n_particles) peak_to_core[p] = i;
	}
	generation = 1;
	for(i=0;i<n_particles;i++){
		int owner, flood_id;
		if(dest[i] == -2){
			if(generation == 0x7fffffff){
				int j;
				for(j=0;j<n_particles;j++) seen[j] = 0;
				generation = 1;
			}
			climb_core(i, neighbor, n_neighbors, peak_to_core, dest, stack, seen,
					generation, n_particles);
			generation++;
		}
		owner = dest[i];
		if(owner >= 0 && !(halo.particles[i].density > cores[owner].saddle_density))
			owner = -1;
		flood_id = flood[i];
		if(owner < 0 && flood_id >= 0
				&& halo.particles[i].density > cores[flood_id].saddle_density)
			owner = flood_id;
		if(owner == flood_id && owner >= 0) n_kept++;
		else if(owner >= 0 && flood_id >= 0) n_moved++;
		else if(owner >= 0) n_gained++;
		else if(flood_id >= 0) n_released++;
		if(owner >= 0){
			set_galaxy_id(i, owner);
			mark_core_particle(i);
			mark_bound(i);
			clear_remaining(i);
		}
		else {
			set_galaxy_id(i, NOT_HALO_MEMBER);
			clear_core_particle(i);
			clear_bound(i);
			mark_remaining(i);
		}
	}
	for(i=0;i<n_cores;i++){
		int p = cores[i].peak_particle;
		if(p < 0 || p >= n_particles) continue;
		set_galaxy_id(p, i);
		mark_peak(p);
		mark_core_particle(p);
		mark_bound(p);
		clear_remaining(p);
	}
	DEBUGPRINT("basin ownership: kept %d moved %d gained %d released %d\n",
			n_kept, n_moved, n_gained, n_released);
	Free(flood);
	Free(seen);
	Free(stack);
	Free(dest);
	Free(peak_to_core);
}

static void refresh_core_counts(SimpleBasicParticleType *particles, int n_particles,
		Coretype *cores, int n_cores) {
	int i;
	for(i=0;i<n_cores;i++){
		cores[i].n_particles = 0;
		cores[i].n_stars = 0;
		cores[i].star_mass = 0.f;
	}
	for(i=0;i<n_particles;i++){
		int h = halo.particles[i].galaxy_id;
		if(h < 0 || h >= n_cores) continue;
		cores[h].n_particles++;
		if(particles[i].type == TYPE_STAR){
			cores[h].n_stars++;
			cores[h].star_mass += particles[i].mass;
		}
	}
}

static void member_frame(SimpleBasicParticleType *particles, int *members, int n_members,
		Coretype *core, int satellite, double *cx, double *cy, double *cz,
		double *cvx, double *cvy, double *cvz) {
	int i, n_nucleus = 0;
	double tmass = 0;
	*cx = *cy = *cz = *cvx = *cvy = *cvz = 0;
	if(satellite){
		float nucleus_r2 = (float)NUCLEUS_RADIUS * (float)NUCLEUS_RADIUS;
		for(i=0;i<n_members;i++){
			int id = members[i];
			float dx, dy, dz, dist2, m;
			if(particles[id].type != TYPE_STAR) continue;
			dx = particles[id].x - core->position.x;
			dy = particles[id].y - core->position.y;
			dz = particles[id].z - core->position.z;
			dist2 = dx*dx + dy*dy + dz*dz;
			if(dist2 > nucleus_r2) continue;
			m = particles[id].mass;
			n_nucleus++;
			tmass += m;
			*cx += particles[id].x * m;
			*cy += particles[id].y * m;
			*cz += particles[id].z * m;
			*cvx += particles[id].vx * m;
			*cvy += particles[id].vy * m;
			*cvz += particles[id].vz * m;
		}
		if(n_nucleus >= NUCLEUS_MIN_STARS && tmass > 0){
			*cx /= tmass; *cy /= tmass; *cz /= tmass;
			*cvx /= tmass; *cvy /= tmass; *cvz /= tmass;
			return;
		}
		*cx = core->position.x; *cy = core->position.y; *cz = core->position.z;
		*cvx = core->velocity.x; *cvy = core->velocity.y; *cvz = core->velocity.z;
		return;
	}
	for(i=0;i<n_members;i++){
		int id = members[i];
		float m;
		if(particles[id].type != TYPE_STAR) continue;
		m = particles[id].mass;
		tmass += m;
		*cx += particles[id].x * m;
		*cy += particles[id].y * m;
		*cz += particles[id].z * m;
		*cvx += particles[id].vx * m;
		*cvy += particles[id].vy * m;
		*cvz += particles[id].vz * m;
	}
	if(tmass > 0){
		*cx /= tmass; *cy /= tmass; *cz /= tmass;
		*cvx /= tmass; *cvy /= tmass; *cvz /= tmass;
	}
	else {
		*cx = core->position.x; *cy = core->position.y; *cz = core->position.z;
		*cvx = core->velocity.x; *cvy = core->velocity.y; *cvz = core->velocity.z;
	}
}

static int retain_self_bound(SimpleBasicParticleType *particles, int n_particles,
		Coretype *cores, int core_id, int central) {
	int i, n_members, n_stars, onmem, iter;
	int *members;
	SampledParticle *samples;
	members = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(members));
	n_members = 0;
	n_stars = 0;
	for(i=0;i<n_particles;i++){
		if(halo.particles[i].galaxy_id != core_id) continue;
		members[n_members++] = i;
		if(particles[i].type == TYPE_STAR) n_stars++;
	}
	if(n_members == 0){
		Free(members);
		return 0;
	}
	samples = (SampledParticle *)Malloc(sizeof(SampledParticle)*(size_t)n_members, PPTR(samples));
	pack_samples(samples, members, n_members, particles);
	onmem = n_members;
	for(iter=0; iter<BOUNDITER; iter++){
		unsigned char *bndflag;
		int nbound, w;
		bndflag = (unsigned char *)Malloc(sizeof(unsigned char)*(size_t)onmem, PPTR(bndflag));
		nbound = count_self_bound(samples, onmem, bndflag, cores+core_id,
				(core_id != central) || n_stars < 1);
		if(nbound == onmem || nbound == 0){
			if(nbound == 0) onmem = 0;
			Free(bndflag);
			break;
		}
		w = 0;
		for(i=0;i<onmem;i++){
			if(bndflag[i] == BOUND) samples[w++] = samples[i];
		}
		onmem = w;
		Free(bndflag);
	}
	for(i=0;i<n_members;i++){
		int id = members[i];
		set_galaxy_id(id, NOT_HALO_MEMBER);
		clear_core_particle(id);
		clear_bound(id);
		mark_remaining(id);
	}
	for(i=0;i<onmem;i++){
		int id = (int)(samples[i].particle - particles);
		set_galaxy_id(id, core_id);
		mark_core_particle(id);
		mark_bound(id);
		clear_remaining(id);
	}
	Free(samples);
	Free(members);
	return onmem;
}

typedef struct RadiusMass {
	float radius;
	float mass;
} RadiusMass;

typedef struct RadialPoint {
	float r;
	float mass;
	float in_host;
	int id;
} RadialPoint;

static int radial_point_cmp(const void *a, const void *b) {
	float ra = ((const RadialPoint *)a)->r;
	float rb = ((const RadialPoint *)b)->r;
	if(ra < rb) return -1;
	if(ra > rb) return 1;
	return 0;
}

static int radius_mass_cmp(const void *a, const void *b) {
	float ra = ((const RadiusMass *)a)->radius;
	float rb = ((const RadiusMass *)b)->radius;
	if(ra < rb) return -1;
	if(ra > rb) return 1;
	return 0;
}

static void set_basin_tidal_radii(SimpleBasicParticleType *particles, int n_particles,
		Coretype *cores, int n_cores, int central, double *central_grow) {
	int i, ncent;
	int *members;
	RadiusMass *shell;
	double *prefix;
	double cx, cy, cz, mcent, rfar, rvir, rgrow;
	members = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(members));
	ncent = indexes_of_galaxy(n_particles, central, members);
	shell = (RadiusMass *)Malloc(sizeof(RadiusMass)*(size_t)(ncent > 0 ? ncent : 1), PPTR(shell));
	prefix = (double *)Malloc(sizeof(double)*(size_t)(ncent+1), PPTR(prefix));
	cx = cores[central].position.x;
	cy = cores[central].position.y;
	cz = cores[central].position.z;
	mcent = 0;
	rfar = 0;
	{
		double *tmass = (double *)Malloc(sizeof(double)*(size_t)n_cores, PPTR(tmass));
		for(i=0;i<n_cores;i++) tmass[i] = 0;
		for(i=0;i<n_particles;i++){
			int h = halo.particles[i].galaxy_id;
			if(h >= 0 && h < n_cores) tmass[h] += particles[i].mass;
		}
		for(i=0;i<ncent;i++){
			int id = members[i];
			double dx = particles[id].x - cx;
			double dy = particles[id].y - cy;
			double dz = particles[id].z - cz;
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			shell[i].radius = (float)r;
			shell[i].mass = particles[id].mass;
			if(r > rfar) rfar = r;
		}
		mcent = tmass[central];
		if(ncent > 1) qsort(shell, (size_t)ncent, sizeof(RadiusMass), radius_mass_cmp);
		prefix[0] = 0;
		for(i=0;i<ncent;i++) prefix[i+1] = prefix[i] + shell[i].mass;
		rvir = virial_radius(mcent > 0 ? mcent : 0);
		rgrow = 4.0 * rvir;
		if(rfar * 1.5 > rgrow) rgrow = rfar * 1.5;
		if(!(rgrow > 0)) rgrow = rfar;
		cores[central].tidal_radius = (float)rgrow;
		*central_grow = rgrow;
		for(i=0;i<n_cores;i++){
			double dx, dy, dz, dist, menc, factor;
			int lo, hi, mid;
			if(i == central) continue;
			dx = cores[i].position.x - cx;
			dy = cores[i].position.y - cy;
			dz = cores[i].position.z - cz;
			dist = sqrt(dx*dx + dy*dy + dz*dz);
			lo = 0;
			hi = ncent;
			while(lo < hi){
				mid = lo + (hi-lo)/2;
				if((double)shell[mid].radius < dist) lo = mid + 1;
				else hi = mid;
			}
			menc = (lo > 0) ? prefix[lo] : 0;
			if(!(menc > 0)) menc = mcent;
			if(!(dist > 0) || !(menc > 0) || !(tmass[i] > 0)){
				cores[i].tidal_radius = 0.f;
				continue;
			}
			factor = pow(tmass[i]/menc, 1.0/3.0);
			if(factor > 0.5) factor = 0.5;
			cores[i].tidal_radius = (float)(dist * factor);
		}
		Free(tmass);
	}
	DEBUGPRINT("basin tidal: central C%d grow=%g M=%g members=%d\n",
			central, rgrow, mcent, ncent);
	Free(prefix);
	Free(shell);
	Free(members);
}

static void potential_from_members(int n_target, Vector3d *target, int n_source,
		Vector3d *source, float *source_mass, float *penergy) {
	int i;
	for(i=0;i<n_target;i++) penergy[i] = 0.f;
	if(n_target <= 0 || n_source <= 0) return;
	if(n_source < 800){
		float epsilon2 = (float)EPSILON * (float)EPSILON;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(static)
#endif
		for(i=0;i<n_target;i++){
			int j;
			float potent = 0.f;
			for(j=0;j<n_source;j++){
				float dx = target[i].x - source[j].x;
				float dy = target[i].y - source[j].y;
				float dz = target[i].z - source[j].z;
				float dist2 = dx*dx + dy*dy + dz*dz;
				if(dist2 > 0.f) potent += -source_mass[j]/sqrtf(dist2+epsilon2);
			}
			penergy[i] = potent * (float)potentfact;
		}
		return;
	}
	{
		size_t nnode = (n_source >= 200000)
			? (size_t)MAX(65*10000, n_source/2)
			: (size_t)MAX(64, n_source*8);
		TStruct *TREE = (TStruct *)Malloc(sizeof(TStruct)*nnode, PPTR(TREE));
		TPtlStruct *ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*(size_t)n_source, PPTR(ptl));
		float theta2 = 0.5f;
		for(i=0;i<n_source;i++){
			ptl[i].type = TYPE_PTL;
			ptl[i].x = source[i].x;
			ptl[i].y = source[i].y;
			ptl[i].z = source[i].z;
			ptl[i].mass = source_mass[i];
			ptl[i].sibling = &ptl[i+1];
		}
		ptl[n_source-1].sibling = NULL;
		build_force_tree(TREE, nnode, ptl, (size_t)n_source, theta2, SERIALIZED);
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(guided)
#endif
		for(i=0;i<n_target;i++){
			particle p;
			p.x = target[i].x;
			p.y = target[i].y;
			p.z = target[i].z;
			p.link02 = 0;
			p.indx = i;
			penergy[i] = tree_potential(&p, theta2, TREE, ptl) * (float)potentfact;
		}
		Free(ptl);
		Free(TREE);
	}
}

static int claim_self_bound_envelope(SimpleBasicParticleType *particles, int n_particles,
		Coretype *cores, int core_id, int central, double central_grow) {
	int i, n_source, n_target, n_claim, satellite;
	int *members, *targets;
	double cx, cy, cz, cvx, cvy, cvz;
	double r2;
	Vector3d *source, *target;
	float *source_mass, *penergy;
	members = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(members));
	targets = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(targets));
	n_source = indexes_of_galaxy(n_particles, core_id, members);
	if(n_source < 2){
		Free(targets);
		Free(members);
		return 0;
	}
	satellite = (core_id != central);
	member_frame(particles, members, n_source, cores+core_id, satellite,
			&cx, &cy, &cz, &cvx, &cvy, &cvz);
	if(satellite) r2 = (double)cores[core_id].tidal_radius * (double)cores[core_id].tidal_radius;
	else r2 = central_grow * central_grow;
	n_target = 0;
	if(r2 > 0){
		for(i=0;i<n_particles;i++){
			double dx, dy, dz;
			if(halo.particles[i].galaxy_id >= 0) continue;
			dx = particles[i].x - cx;
			dy = particles[i].y - cy;
			dz = particles[i].z - cz;
			if(dx*dx + dy*dy + dz*dz > r2) continue;
			targets[n_target++] = i;
		}
	}
	if(n_target == 0){
		Free(targets);
		Free(members);
		return 0;
	}
	source = (Vector3d *)Malloc(sizeof(Vector3d)*(size_t)n_source, PPTR(source));
	source_mass = (float *)Malloc(sizeof(float)*(size_t)n_source, PPTR(source_mass));
	target = (Vector3d *)Malloc(sizeof(Vector3d)*(size_t)n_target, PPTR(target));
	penergy = (float *)Malloc(sizeof(float)*(size_t)n_target, PPTR(penergy));
	for(i=0;i<n_source;i++){
		int id = members[i];
		source[i].x = particles[id].x;
		source[i].y = particles[id].y;
		source[i].z = particles[id].z;
		source_mass[i] = particles[id].mass;
	}
	for(i=0;i<n_target;i++){
		int id = targets[i];
		target[i].x = particles[id].x;
		target[i].y = particles[id].y;
		target[i].z = particles[id].z;
	}
	potential_from_members(n_target, target, n_source, source, source_mass, penergy);
	n_claim = 0;
	for(i=0;i<n_target;i++){
		int id = targets[i];
		float dx = (particles[id].x - cx) * r1kineticfact;
		float dy = (particles[id].y - cy) * r1kineticfact;
		float dz = (particles[id].z - cz) * r1kineticfact;
		float vx = dx + (particles[id].vx - cvx) * r2kineticfact;
		float vy = dy + (particles[id].vy - cvy) * r2kineticfact;
		float vz = dz + (particles[id].vz - cvz) * r2kineticfact;
		float ke = 0.5f*(vx*vx + vy*vy + vz*vz);
		if(ke > -penergy[i]) continue;
		set_galaxy_id(id, core_id);
		mark_bound(id);
		clear_remaining(id);
		clear_core_particle(id);
		n_claim++;
	}
	DEBUGPRINT("basin grow: C%d claimed %d of %d inside r=%g (sources %d)\n",
			core_id, n_claim, n_target,
			satellite ? cores[core_id].tidal_radius : (float)central_grow, n_source);
	Free(penergy);
	Free(target);
	Free(source_mass);
	Free(source);
	Free(targets);
	Free(members);
	return n_claim;
}

static int host_by_peak(Coretype *cores, int n_cores) {
	int i, best = 0;
	float den = -1.f;
	for(i=0;i<n_cores;i++){
		if(cores[i].peak_density > den){
			den = cores[i].peak_density;
			best = i;
		}
	}
	return best;
}

static double enclosed_mass(const RadiusMass *shell, const double *prefix, int n, double dist) {
	int lo = 0, hi = n, mid;
	while(lo < hi){
		mid = lo + (hi-lo)/2;
		if((double)shell[mid].radius < dist) lo = mid + 1;
		else hi = mid;
	}
	return (lo > 0) ? prefix[lo] : 0;
}

static int collect_ball(const int *head, const int *next, int nx, int ny, int nz,
		float xmin, float ymin, float zmin, float cell,
		SimpleBasicParticleType *particles, float cx, float cy, float cz, float radius,
		int *out) {
	int ix0, ix1, iy0, iy1, iz0, iz1, ix, iy, iz, nout;
	float r2;
	if(!(radius > 0.f)) return 0;
	r2 = radius * radius;
	ix0 = (int)((cx - radius - xmin) / cell);
	ix1 = (int)((cx + radius - xmin) / cell);
	iy0 = (int)((cy - radius - ymin) / cell);
	iy1 = (int)((cy + radius - ymin) / cell);
	iz0 = (int)((cz - radius - zmin) / cell);
	iz1 = (int)((cz + radius - zmin) / cell);
	if(ix0 < 0) ix0 = 0;
	if(iy0 < 0) iy0 = 0;
	if(iz0 < 0) iz0 = 0;
	if(ix1 >= nx) ix1 = nx - 1;
	if(iy1 >= ny) iy1 = ny - 1;
	if(iz1 >= nz) iz1 = nz - 1;
	nout = 0;
	for(iz=iz0; iz<=iz1; iz++){
		for(iy=iy0; iy<=iy1; iy++){
			for(ix=ix0; ix<=ix1; ix++){
				int p = head[(iz * ny + iy) * nx + ix];
				for(; p >= 0; p = next[p]){
					float dx = particles[p].x - cx;
					float dy = particles[p].y - cy;
					float dz = particles[p].z - cz;
					if(dx*dx + dy*dy + dz*dz <= r2) out[nout++] = p;
				}
			}
		}
	}
	return nout;
}

static int keep_particle(SimpleBasicParticleType *particles, int id, Coretype *core) {
	float dx, dy, dz, r2;
	if(id == core->peak_particle) return 1;
	if(particles[id].type != TYPE_STAR) return 0;
	dx = particles[id].x - core->position.x;
	dy = particles[id].y - core->position.y;
	dz = particles[id].z - core->position.z;
	r2 = (float)NUCLEUS_RADIUS * (float)NUCLEUS_RADIUS;
	return dx*dx + dy*dy + dz*dz <= r2;
}

/* Drop particles that are not bound to this set. The peak and the nucleus
 * stars stay, so a later iteration still has the satellite's own mass. */
static int retain_bound_set(SimpleBasicParticleType *particles, int *ids, int n,
		Coretype *core, int satellite) {
	int iter, ncur;
	ncur = n;
	for(iter=0; iter<BOUNDITER && ncur >= 2; iter++){
		int i, nkeep;
		double cx, cy, cz, cvx, cvy, cvz;
		Vector3d *pos;
		float *mass, *pe;
		member_frame(particles, ids, ncur, core, satellite, &cx, &cy, &cz, &cvx, &cvy, &cvz);
		pos = (Vector3d *)Malloc(sizeof(Vector3d)*(size_t)ncur, PPTR(pos));
		mass = (float *)Malloc(sizeof(float)*(size_t)ncur, PPTR(mass));
		pe = (float *)Malloc(sizeof(float)*(size_t)ncur, PPTR(pe));
		for(i=0;i<ncur;i++){
			int id = ids[i];
			pos[i].x = particles[id].x;
			pos[i].y = particles[id].y;
			pos[i].z = particles[id].z;
			mass[i] = particles[id].mass;
		}
		potential_from_members(ncur, pos, ncur, pos, mass, pe);
		nkeep = 0;
		for(i=0;i<ncur;i++){
			int id = ids[i];
			float dx, dy, dz, vx, vy, vz, ke;
			if(keep_particle(particles, id, core)){
				ids[nkeep++] = id;
				continue;
			}
			dx = (particles[id].x - cx) * r1kineticfact;
			dy = (particles[id].y - cy) * r1kineticfact;
			dz = (particles[id].z - cz) * r1kineticfact;
			vx = dx + (particles[id].vx - cvx) * r2kineticfact;
			vy = dy + (particles[id].vy - cvy) * r2kineticfact;
			vz = dz + (particles[id].vz - cvz) * r2kineticfact;
			ke = 0.5f * (vx*vx + vy*vy + vz*vz);
			if(!(ke > -pe[i])) ids[nkeep++] = id;
		}
		Free(pe);
		Free(mass);
		Free(pos);
		if(nkeep == ncur) break;
		ncur = nkeep;
	}
	return ncur;
}

static void clear_initial_members(int core_id, const int *csr, const int *csr_off) {
	int j;
	for(j=csr_off[core_id]; j<csr_off[core_id+1]; j++){
		int id = csr[j];
		set_galaxy_id(id, NOT_HALO_MEMBER);
		clear_core_particle(id);
		clear_bound(id);
		mark_remaining(id);
	}
}

static void commit_members(SimpleBasicParticleType *particles, int *ids, int n, int core_id) {
	int i;
	(void)particles;
	for(i=0;i<n;i++){
		set_galaxy_id(ids[i], core_id);
		mark_bound(ids[i]);
		clear_remaining(ids[i]);
	}
}

static void decide_satellite(SimpleBasicParticleType *particles, int n_particles,
		Coretype *cores, int core_id, int host,
		const int *head, const int *next, int nx, int ny, int nz,
		float xmin, float ymin, float zmin, float cell,
		const RadiusMass *hshell, const double *hprefix, int nhost,
		const int *csr, const int *csr_off, int *work) {
	int i, nball, ninside, nbound;
	double hx, hy, hz, sx, sy, sz, dx, dy, dz, dist, rt, m_self, m_host;
	float search;
	hx = cores[host].position.x;
	hy = cores[host].position.y;
	hz = cores[host].position.z;
	sx = cores[core_id].position.x;
	sy = cores[core_id].position.y;
	sz = cores[core_id].position.z;
	dx = sx - hx; dy = sy - hy; dz = sz - hz;
	dist = sqrt(dx*dx + dy*dy + dz*dz);
	search = (float)(0.5 * dist);
	nball = collect_ball(head, next, nx, ny, nz, xmin, ymin, zmin, cell,
			particles, (float)sx, (float)sy, (float)sz, search, work);
	{
		int w = 0;
		for(i=0;i<nball;i++){
			int id = work[i];
			int h = halo.particles[id].galaxy_id;
			if(h < 0 || h == core_id) work[w++] = id;
		}
		nball = w;
	}
	rt = 0.5 * dist;
	m_self = 0;
	m_host = 0;
	/* One radial sort. The eight Jacobi updates then read prefix sums. */
	if(nball > 0){
		RadialPoint *radial = (RadialPoint *)Malloc(sizeof(RadialPoint)*(size_t)nball, PPTR(radial));
		double *prefix = (double *)Malloc(sizeof(double)*(size_t)(nball+1), PPTR(prefix));
		double *prefix_in = (double *)Malloc(sizeof(double)*(size_t)(nball+1), PPTR(prefix_in));
		int k;
		for(k=0;k<nball;k++){
			int id = work[k];
			double px = particles[id].x - sx;
			double py = particles[id].y - sy;
			double pz = particles[id].z - sz;
			double qx = particles[id].x - hx;
			double qy = particles[id].y - hy;
			double qz = particles[id].z - hz;
			radial[k].r = (float)sqrt(px*px + py*py + pz*pz);
			radial[k].mass = particles[id].mass;
			radial[k].in_host = (sqrt(qx*qx + qy*qy + qz*qz) < dist) ? radial[k].mass : 0.f;
			radial[k].id = id;
		}
		qsort(radial, (size_t)nball, sizeof(RadialPoint), radial_point_cmp);
		prefix[0] = 0;
		prefix_in[0] = 0;
		for(k=0;k<nball;k++){
			prefix[k+1] = prefix[k] + radial[k].mass;
			prefix_in[k+1] = prefix_in[k] + radial[k].in_host;
			work[k] = radial[k].id;
		}
		for(i=0;i<8;i++){
			int lo = 0, hi = nball, mid;
			double m, m_in, menc, factor;
			while(lo < hi){
				mid = lo + (hi-lo)/2;
				if((double)radial[mid].r < rt) lo = mid + 1;
				else hi = mid;
			}
			m = prefix[lo];
			m_in = prefix_in[lo];
			menc = enclosed_mass(hshell, hprefix, nhost, dist) - m_in;
			if(!(menc > 0)) menc = enclosed_mass(hshell, hprefix, nhost, dist);
			if(!(m > 0) || !(menc > 0) || !(dist > 0)){
				rt = 0;
				m_self = m;
				m_host = menc;
				break;
			}
			factor = pow(m / (3.0 * menc), 1.0/3.0);
			if(factor > 0.5) factor = 0.5;
			rt = dist * factor;
			m_self = m;
			m_host = menc;
		}
		Free(prefix_in);
		Free(prefix);
		Free(radial);
	}
	ninside = 0;
	for(i=0;i<nball;i++){
		int id = work[i];
		double px = particles[id].x - sx;
		double py = particles[id].y - sy;
		double pz = particles[id].z - sz;
		if(sqrt(px*px + py*py + pz*pz) < rt || keep_particle(particles, id, cores+core_id))
			work[ninside++] = id;
	}
	if(ninside == 0 && cores[core_id].peak_particle >= 0){
		work[ninside++] = cores[core_id].peak_particle;
	}
	nbound = retain_bound_set(particles, work, ninside, cores+core_id, 1);
	clear_initial_members(core_id, csr, csr_off);
	commit_members(particles, work, nbound, core_id);
	cores[core_id].tidal_radius = (float)rt;
	cores[core_id].n_particles = nbound;
	DEBUGPRINT("C%d satellite rt=%g m=%g M=%g d=%g bound=%d from ball=%d c %g %g %g\n",
			core_id, rt, m_self, m_host, dist, nbound, nball,
			sx+halo.origin.x, sy+halo.origin.y, sz+halo.origin.z);
	(void)n_particles;
}

static int append_bound_candidates(SimpleBasicParticleType *particles,
		int *source, int nsrc, int *cand, int ncand, Coretype *core) {
	int i, nadd;
	double cx, cy, cz, cvx, cvy, cvz;
	Vector3d *spos, *tpos;
	float *smass, *pe;
	if(ncand <= 0 || nsrc <= 0) return nsrc;
	member_frame(particles, source, nsrc, core, 0, &cx, &cy, &cz, &cvx, &cvy, &cvz);
	spos = (Vector3d *)Malloc(sizeof(Vector3d)*(size_t)nsrc, PPTR(spos));
	smass = (float *)Malloc(sizeof(float)*(size_t)nsrc, PPTR(smass));
	tpos = (Vector3d *)Malloc(sizeof(Vector3d)*(size_t)ncand, PPTR(tpos));
	pe = (float *)Malloc(sizeof(float)*(size_t)ncand, PPTR(pe));
	for(i=0;i<nsrc;i++){
		int id = source[i];
		spos[i].x = particles[id].x;
		spos[i].y = particles[id].y;
		spos[i].z = particles[id].z;
		smass[i] = particles[id].mass;
	}
	for(i=0;i<ncand;i++){
		int id = cand[i];
		tpos[i].x = particles[id].x;
		tpos[i].y = particles[id].y;
		tpos[i].z = particles[id].z;
	}
	potential_from_members(ncand, tpos, nsrc, spos, smass, pe);
	nadd = 0;
	{
		int w = 0;
		for(i=0;i<ncand;i++){
			int id = cand[i];
			float dx = (particles[id].x - cx) * r1kineticfact;
			float dy = (particles[id].y - cy) * r1kineticfact;
			float dz = (particles[id].z - cz) * r1kineticfact;
			float vx = dx + (particles[id].vx - cvx) * r2kineticfact;
			float vy = dy + (particles[id].vy - cvy) * r2kineticfact;
			float vz = dz + (particles[id].vz - cvz) * r2kineticfact;
			float ke = 0.5f * (vx*vx + vy*vy + vz*vz);
			if(ke > -pe[i]){
				cand[w++] = id;
				continue;
			}
			source[nsrc + nadd] = id;
			nadd++;
		}
		ncand = w;
	}
	Free(pe);
	Free(tpos);
	Free(smass);
	Free(spos);
	(void)ncand;
	return nsrc + nadd;
}

static void decide_host(SimpleBasicParticleType *particles, int n_particles,
		Coretype *cores, int host,
		const int *head, const int *next, int nx, int ny, int nz,
		float xmin, float ymin, float zmin, float cell,
		const int *csr, const int *csr_off, int *work) {
	int i, nsrc, ncand, nbound, round;
	int *cand;
	double hx, hy, hz, m_basin, rfar, rsearch, mall;
	hx = cores[host].position.x;
	hy = cores[host].position.y;
	hz = cores[host].position.z;
	m_basin = 0;
	rfar = 0;
	mall = 0;
	nsrc = 0;
	for(i=0;i<n_particles;i++){
		mall += particles[i].mass;
		if(halo.particles[i].galaxy_id != host) continue;
		{
			double dx = particles[i].x - hx;
			double dy = particles[i].y - hy;
			double dz = particles[i].z - hz;
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			m_basin += particles[i].mass;
			if(r > rfar) rfar = r;
			work[nsrc++] = i;
		}
	}
	nsrc = retain_bound_set(particles, work, nsrc, cores+host, 0);
	rsearch = 4.0 * virial_radius(m_basin > 0 ? m_basin : 0);
	if(1.5 * rfar > rsearch) rsearch = 1.5 * rfar;
	{
		double rhalo = virial_radius(mall > 0 ? mall : 0);
		if(rhalo > 0 && rsearch > rhalo) rsearch = rhalo;
	}
	if(!(rsearch > 0)) rsearch = rfar;
	cand = work + n_particles;
	ncand = collect_ball(head, next, nx, ny, nz, xmin, ymin, zmin, cell,
			particles, (float)hx, (float)hy, (float)hz, (float)rsearch, cand);
	{
		int w = 0;
		for(i=0;i<ncand;i++){
			if(halo.particles[cand[i]].galaxy_id < 0) cand[w++] = cand[i];
		}
		ncand = w;
	}
	/* Source stays the host's own mass. Unassigned particles are tests,
	 * not part of the potential, until they pass. */
	for(round=0; round<3 && ncand > 0 && nsrc > 0; round++){
		int nprev = nsrc;
		nsrc = append_bound_candidates(particles, work, nsrc, cand, ncand, cores+host);
		/* Failures were packed at the front of cand. */
		ncand -= nsrc - nprev;
		if(nsrc == nprev) break;
	}
	nbound = nsrc;
	clear_initial_members(host, csr, csr_off);
	commit_members(particles, work, nbound, host);
	cores[host].tidal_radius = (float)rsearch;
	cores[host].n_particles = nbound;
	DEBUGPRINT("C%d host search=%g basin_mass=%g bound=%d c %g %g %g\n",
			host, rsearch, m_basin, nbound,
			hx+halo.origin.x, hy+halo.origin.y, hz+halo.origin.z);
}

struct DarkSort {
	SimpleBasicParticleType *particles;
	float x, y, z;
};

static int dark_nearer(const void *a, const void *b, void *arg) {
	const struct DarkSort *s = (const struct DarkSort *)arg;
	int ia = *(const int *)a;
	int ib = *(const int *)b;
	double ax = s->particles[ia].x - s->x;
	double ay = s->particles[ia].y - s->y;
	double az = s->particles[ia].z - s->z;
	double bx = s->particles[ib].x - s->x;
	double by = s->particles[ib].y - s->y;
	double bz = s->particles[ib].z - s->z;
	double da = ax*ax + ay*ay + az*az;
	double db = bx*bx + by*by + bz*bz;
	if(da < db) return -1;
	if(da > db) return 1;
	return 0;
}

static void gauss_axis(const float *in, float *out, int nx, int ny, int nz,
		int axis, const float *ker, int hw) {
	int x, y, z, t;
	for(z=0;z<nz;z++){
		for(y=0;y<ny;y++){
			for(x=0;x<nx;x++){
				double s = 0;
				for(t=-hw;t<=hw;t++){
					int xx=x, yy=y, zz=z;
					if(axis==0) xx = x+t;
					else if(axis==1) yy = y+t;
					else zz = z+t;
					if(xx<0||yy<0||zz<0||xx>=nx||yy>=ny||zz>=nz) continue;
					s += ker[t+hw] * in[((long)zz*ny+yy)*nx+xx];
				}
				out[((long)z*ny+y)*nx+x] = (float)s;
			}
		}
	}
}

/* Dark cores are sought after stellar satellites have claimed their
 * members and before the host adopts the leftovers. Seeds are maxima of
 * the DM density. m is the mass above the host shell density, and the
 * potential uses only that excess. */
static int add_dark_galaxies(SimpleBasicParticleType *particles, int n_particles,
		Coretype *cores, int n_cores, int host,
		const int *head, const int *next, int nx, int ny, int nz,
		float xmin, float ymin, float zmin, float qcell,
		const RadiusMass *hshell, const double *hprefix, int nhost,
		int *work) {
	int i, n_dm, n_peak, n_keep, ix, iy, iz;
	int gnx, gny, gnz;
	long ncell;
	float cell, span, sx, sy, sz;
	float *den, *tmp, *ker;
	int hw;
	double hx, hy, hz;
	RadiusMass *dmshell;
	double *dmprefix;
	typedef struct { float x, y, z, den; } DarkPeak;
	DarkPeak *peaks;
	const float smooth = (float)DM_GAUSSIAN_SMOOTHING_LENGTH;
	const float dedup = (float)STAR_DM_DEDUP_LENGTH;
	const int min_dm = DM_MINCORENMEM;
	const int min_nucleus = 20;
	float qxmin = xmin, qymin = ymin, qzmin = zmin;
	int planted = 0;

	if(n_cores > 0 && host >= 0 && host < n_cores){
		hx = cores[host].position.x;
		hy = cores[host].position.y;
		hz = cores[host].position.z;
	}
	else {
		hx = hy = hz = 0;
	}
	{
		int seen = 0;
		float x0=0,y0=0,z0=0,x1=0,y1=0,z1=0;
		for(i=0;i<n_particles;i++){
			if(particles[i].type != TYPE_DM) continue;
			if(!seen){
				x0=x1=particles[i].x; y0=y1=particles[i].y; z0=z1=particles[i].z;
				seen=1;
			}
			if(particles[i].x<x0) x0=particles[i].x;
			if(particles[i].y<y0) y0=particles[i].y;
			if(particles[i].z<z0) z0=particles[i].z;
			if(particles[i].x>x1) x1=particles[i].x;
			if(particles[i].y>y1) y1=particles[i].y;
			if(particles[i].z>z1) z1=particles[i].z;
		}
		if(!seen) return n_cores;
		span = x1-x0;
		if(y1-y0>span) span=y1-y0;
		if(z1-z0>span) span=z1-z0;
		cell = 0.006f;
		gnx = (int)((x1-x0)/cell)+1;
		gny = (int)((y1-y0)/cell)+1;
		gnz = (int)((z1-z0)/cell)+1;
		while((long)gnx*(long)gny*(long)gnz > 80000000L){
			cell *= 1.25f;
			gnx = (int)((x1-x0)/cell)+1;
			gny = (int)((y1-y0)/cell)+1;
			gnz = (int)((z1-z0)/cell)+1;
		}
		if(gnx<1) gnx=1;
		if(gny<1) gny=1;
		if(gnz<1) gnz=1;
		ncell = (long)gnx*(long)gny*(long)gnz;
		(void)span;
		/* Query-grid origin stays in the function arguments. Density-grid
		 * origin is x0,y0,z0 and is stored back into xmin only after the
		 * query origin has been copied. */
		xmin = x0; ymin = y0; zmin = z0;
	}
	DEBUGPRINT("dark grid %d %d %d cell=%g\n", gnx, gny, gnz, cell);
	den = (float *)Malloc(sizeof(float)*(size_t)ncell, PPTR(den));
	tmp = (float *)Malloc(sizeof(float)*(size_t)ncell, PPTR(tmp));
	{
		long c;
		for(c=0;c<ncell;c++) den[c]=0.f;
	}
	for(i=0;i<n_particles;i++){
		int gx, gy, gz;
		long c;
		if(particles[i].type != TYPE_DM) continue;
		gx = (int)((particles[i].x - xmin)/cell);
		gy = (int)((particles[i].y - ymin)/cell);
		gz = (int)((particles[i].z - zmin)/cell);
		if(gx<0) gx=0;
		if(gy<0) gy=0;
		if(gz<0) gz=0;
		if(gx>=gnx) gx=gnx-1;
		if(gy>=gny) gy=gny-1;
		if(gz>=gnz) gz=gnz-1;
		c = ((long)gz*gny + gy)*gnx + gx;
		den[c] += particles[i].mass;
	}
	{
		double sigma = (double)smooth / (double)cell;
		double norm = 0;
		hw = (int)ceil(3.0*sigma);
		if(hw<1) hw=1;
		if(hw>24) hw=24;
		ker = (float *)Malloc(sizeof(float)*(size_t)(2*hw+1), PPTR(ker));
		for(i=-hw;i<=hw;i++){
			double w = exp(-0.5*(i*i)/(sigma*sigma));
			ker[i+hw] = (float)w;
			norm += w;
		}
		for(i=0;i<2*hw+1;i++) ker[i] = (float)(ker[i]/norm);
		gauss_axis(den, tmp, gnx, gny, gnz, 0, ker, hw);
		gauss_axis(tmp, den, gnx, gny, gnz, 1, ker, hw);
		gauss_axis(den, tmp, gnx, gny, gnz, 2, ker, hw);
		Free(ker);
		{
			long c;
			for(c=0;c<ncell;c++) den[c]=tmp[c];
		}
	}
	Free(tmp);

	n_peak = 0;
	peaks = (DarkPeak *)Malloc(sizeof(DarkPeak)*4096, PPTR(peaks));
	for(iz=1;iz<gnz-1;iz++){
		for(iy=1;iy<gny-1;iy++){
			for(ix=1;ix<gnx-1;ix++){
				long c = ((long)iz*gny+iy)*gnx+ix;
				float v = den[c];
				int dx, dy, dz, maxed;
				if(!(v>0.f)) continue;
				maxed = 1;
				for(dz=-1;dz<=1 && maxed;dz++){
					for(dy=-1;dy<=1 && maxed;dy++){
						for(dx=-1;dx<=1;dx++){
							if(dx==0&&dy==0&&dz==0) continue;
							if(den[c + ((long)dz*gny+dy)*gnx+dx] >= v){
								maxed = 0;
								break;
							}
						}
					}
				}
				if(!maxed) continue;
				if(n_peak < 4096){
					peaks[n_peak].x = xmin + (ix+0.5f)*cell;
					peaks[n_peak].y = ymin + (iy+0.5f)*cell;
					peaks[n_peak].z = zmin + (iz+0.5f)*cell;
					peaks[n_peak].den = v;
					n_peak++;
				} else {
					int w, worst = 0;
					for(w=1;w<n_peak;w++){
						if(peaks[w].den < peaks[worst].den) worst = w;
					}
					if(v > peaks[worst].den){
						peaks[worst].x = xmin + (ix+0.5f)*cell;
						peaks[worst].y = ymin + (iy+0.5f)*cell;
						peaks[worst].z = zmin + (iz+0.5f)*cell;
						peaks[worst].den = v;
					}
				}
			}
		}
	}
	Free(den);
	DEBUGPRINT("dark maxima %d\n", n_peak);
	if(n_cores <= 0){
		int ibest = 0, pnear = -1, nlab = 0;
		double bestd = 1.e300, tmass = 0;
		double cx = 0, cy = 0, cz = 0, cvx = 0, cvy = 0, cvz = 0;
		float r2;
		if(n_peak <= 0){
			Free(peaks);
			return 0;
		}
		for(i=1;i<n_peak;i++){
			if(peaks[i].den > peaks[ibest].den) ibest = i;
		}
		r2 = smooth * smooth;
		for(i=0;i<n_particles;i++){
			double dx, dy, dz, d2;
			if(particles[i].type != TYPE_DM) continue;
			dx = particles[i].x - peaks[ibest].x;
			dy = particles[i].y - peaks[ibest].y;
			dz = particles[i].z - peaks[ibest].z;
			d2 = dx*dx + dy*dy + dz*dz;
			if(d2 < bestd){ bestd = d2; pnear = i; }
			if(d2 > r2) continue;
			set_galaxy_id(i, 0);
			mark_bound(i);
			clear_remaining(i);
			tmass += particles[i].mass;
			cx += particles[i].x * particles[i].mass;
			cy += particles[i].y * particles[i].mass;
			cz += particles[i].z * particles[i].mass;
			cvx += particles[i].vx * particles[i].mass;
			cvy += particles[i].vy * particles[i].mass;
			cvz += particles[i].vz * particles[i].mass;
			nlab++;
		}
		memset(cores, 0, sizeof(Coretype));
		cores[0].is_dark = 1;
		cores[0].peak_particle = pnear;
		cores[0].position.x = peaks[ibest].x;
		cores[0].position.y = peaks[ibest].y;
		cores[0].position.z = peaks[ibest].z;
		cores[0].peak_density = peaks[ibest].den;
		cores[0].n_particles = nlab;
		if(tmass > 0){
			cores[0].velocity.x = (float)(cvx / tmass);
			cores[0].velocity.y = (float)(cvy / tmass);
			cores[0].velocity.z = (float)(cvz / tmass);
		}
		else if(pnear >= 0){
			cores[0].velocity.x = particles[pnear].vx;
			cores[0].velocity.y = particles[pnear].vy;
			cores[0].velocity.z = particles[pnear].vz;
		}
		n_cores = 1;
		host = 0;
		planted = 1;
		hx = peaks[ibest].x;
		hy = peaks[ibest].y;
		hz = peaks[ibest].z;
		DEBUGPRINT("DMO host peak den=%g labeled=%d c %g %g %g\n",
				peaks[ibest].den, nlab,
				hx+halo.origin.x, hy+halo.origin.y, hz+halo.origin.z);
		(void)cx; (void)cy; (void)cz;
	}

	n_dm = 0;
	dmshell = (RadiusMass *)Malloc(sizeof(RadiusMass)*(size_t)n_particles, PPTR(dmshell));
	for(i=0;i<n_particles;i++){
		double dx, dy, dz;
		if(particles[i].type != TYPE_DM) continue;
		dx = particles[i].x - hx;
		dy = particles[i].y - hy;
		dz = particles[i].z - hz;
		dmshell[n_dm].radius = (float)sqrt(dx*dx+dy*dy+dz*dz);
		dmshell[n_dm].mass = particles[i].mass;
		n_dm++;
	}
	if(n_dm>1) qsort(dmshell, (size_t)n_dm, sizeof(RadiusMass), radius_mass_cmp);
	dmprefix = (double *)Malloc(sizeof(double)*(size_t)(n_dm+1), PPTR(dmprefix));
	dmprefix[0] = 0;
	for(i=0;i<n_dm;i++) dmprefix[i+1] = dmprefix[i] + dmshell[i].mass;

	n_keep = 0;
	{
		/* Jacobi mass must be centred on the host. A DMO call arrives with
		 * a profile centred elsewhere, so rebuild it after the host is known. */
		RadiusMass *local_shell = NULL;
		double *local_prefix = NULL;
		const RadiusMass *use_shell = hshell;
		const double *use_prefix = hprefix;
		int use_n = nhost;
		if(planted){
			int q;
			local_shell = (RadiusMass *)Malloc(sizeof(RadiusMass)*(size_t)n_particles, PPTR(local_shell));
			local_prefix = (double *)Malloc(sizeof(double)*(size_t)(n_particles+1), PPTR(local_prefix));
			for(q=0;q<n_particles;q++){
				double dx = particles[q].x - hx;
				double dy = particles[q].y - hy;
				double dz = particles[q].z - hz;
				local_shell[q].radius = (float)sqrt(dx*dx + dy*dy + dz*dz);
				local_shell[q].mass = particles[q].mass;
			}
			qsort(local_shell, (size_t)n_particles, sizeof(RadiusMass), radius_mass_cmp);
			local_prefix[0] = 0;
			for(q=0;q<n_particles;q++) local_prefix[q+1] = local_prefix[q] + local_shell[q].mass;
			use_shell = local_shell;
			use_prefix = local_prefix;
			use_n = n_particles;
		}
		(void)use_shell; (void)use_prefix; (void)use_n;
		/* The satellite loop below still calls enclosed_mass on the caller's
		 * profile. For a stellar host that profile is already centred. For
		 * DMO the caller's profile is rebuilt into use_*. Point the names
		 * the loop already uses at that buffer by writing through a copy
		 * only when we allocated one. */
		if(local_shell){
			hshell = local_shell;
			hprefix = local_prefix;
			nhost = use_n;
		}
	for(i=0;i<n_peak;i++){
		int k, nball, nsrc, nbound, ndm_bound, n_nucleus;
		int stellar_hit;
		double dx, dy, dz, dist, rho, m_smooth, rt, m_ex;
		double cx, cy, cz, cvx, cvy, cvz, tmass;
		float search;
		dx = peaks[i].x - hx;
		dy = peaks[i].y - hy;
		dz = peaks[i].z - hz;
		dist = sqrt(dx*dx+dy*dy+dz*dz);
		if(!(dist > 2.0*(double)smooth)) continue;
		stellar_hit = 0;
		for(k=0;k<n_cores;k++){
			double cut, px, py, pz, pd;
			cut = cores[k].tidal_radius > dedup ? cores[k].tidal_radius : dedup;
			px = peaks[i].x - cores[k].position.x;
			py = peaks[i].y - cores[k].position.y;
			pz = peaks[i].z - cores[k].position.z;
			pd = sqrt(px*px+py*py+pz*pz);
			if(pd < cut){ stellar_hit = 1; break; }
		}
		if(stellar_hit) continue;
		{
			double rin = 0.8*dist, rout = 1.2*dist;
			double minn = enclosed_mass(dmshell, dmprefix, n_dm, rin);
			double mout = enclosed_mass(dmshell, dmprefix, n_dm, rout);
			double vol = (4.0/3.0)*3.141592653589793*(rout*rout*rout - rin*rin*rin);
			rho = (vol>0) ? (mout-minn)/vol : 0;
		}
		if(!(rho>0)) continue;
		search = smooth;
		nball = collect_ball(head, next, nx, ny, nz, qxmin, qymin, qzmin, qcell,
				particles, peaks[i].x, peaks[i].y, peaks[i].z, search, work);
		m_smooth = 0;
		for(k=0;k<nball;k++){
			int id = work[k];
			if(particles[id].type==TYPE_DM) m_smooth += particles[id].mass;
		}
		{
			double vol = (4.0/3.0)*3.141592653589793*(double)smooth*(double)smooth*(double)smooth;
			if(!(m_smooth/vol > 3.0*rho)) continue;
		}
		search = (float)(0.5*dist);
		nball = collect_ball(head, next, nx, ny, nz, qxmin, qymin, qzmin, qcell,
				particles, peaks[i].x, peaks[i].y, peaks[i].z, search, work);
		{
			int w = 0;
			for(k=0;k<nball;k++){
				if(halo.particles[work[k]].galaxy_id < 0) work[w++] = work[k];
			}
			nball = w;
		}
		if(nball < min_dm) continue;
		rt = 0.5*dist;
		m_ex = 0;
		{
			int iter, capped = 0;
			for(iter=0; iter<8; iter++){
				double mraw = 0, menc, factor, mbg;
				for(k=0;k<nball;k++){
					int id = work[k];
					double px = particles[id].x - peaks[i].x;
					double py = particles[id].y - peaks[i].y;
					double pz = particles[id].z - peaks[i].z;
					if(sqrt(px*px+py*py+pz*pz) < rt) mraw += particles[id].mass;
				}
				mbg = (4.0/3.0)*3.141592653589793*rt*rt*rt*rho;
				m_ex = mraw - mbg;
				if(!(m_ex>0)){ m_ex = 0; break; }
				menc = enclosed_mass(hshell, hprefix, nhost, dist);
				if(menc > mraw) menc -= mraw;
				factor = pow(m_ex/(3.0*menc), 1.0/3.0);
				if(factor >= 0.5){ capped = 1; break; }
				rt = dist * factor;
			}
			if(capped || !(m_ex>0) || rt < (double)smooth) continue;
		}
		stellar_hit = 0;
		for(k=0;k<n_cores;k++){
			double px = peaks[i].x - cores[k].position.x;
			double py = peaks[i].y - cores[k].position.y;
			double pz = peaks[i].z - cores[k].position.z;
			if(sqrt(px*px+py*py+pz*pz) < rt){ stellar_hit = 1; break; }
		}
		if(stellar_hit) continue;
		nsrc = 0;
		for(k=0;k<nball;k++){
			int id = work[k];
			double px = particles[id].x - peaks[i].x;
			double py = particles[id].y - peaks[i].y;
			double pz = particles[id].z - peaks[i].z;
			if(sqrt(px*px+py*py+pz*pz) < rt) work[nsrc++] = id;
		}
		if(nsrc < min_dm || nsrc > 300000) continue;
		{
			struct DarkSort ctx;
			ctx.particles = particles;
			ctx.x = peaks[i].x;
			ctx.y = peaks[i].y;
			ctx.z = peaks[i].z;
			qsort_r(work, (size_t)nsrc, sizeof(int), dark_nearer, &ctx);
		}
		{
			double acc = 0;
			int n_ex = 0;
			for(k=0;k<nsrc;k++){
				acc += particles[work[k]].mass;
				n_ex++;
				if(acc >= m_ex) break;
			}
			nsrc = n_ex;
		}
		n_nucleus = 0;
		tmass = 0;
		cx = cy = cz = cvx = cvy = cvz = 0;
		for(k=0;k<nsrc;k++){
			int id = work[k];
			double px, py, pz, r2;
			if(particles[id].type != TYPE_DM) continue;
			px = particles[id].x - peaks[i].x;
			py = particles[id].y - peaks[i].y;
			pz = particles[id].z - peaks[i].z;
			r2 = px*px+py*py+pz*pz;
			if(r2 > (double)NUCLEUS_RADIUS*(double)NUCLEUS_RADIUS) continue;
			tmass += particles[id].mass;
			cx += particles[id].x * particles[id].mass;
			cy += particles[id].y * particles[id].mass;
			cz += particles[id].z * particles[id].mass;
			cvx += particles[id].vx * particles[id].mass;
			cvy += particles[id].vy * particles[id].mass;
			cvz += particles[id].vz * particles[id].mass;
			n_nucleus++;
		}
		if(n_nucleus < min_nucleus){
			n_nucleus = 0;
			tmass = 0;
			cx = cy = cz = cvx = cvy = cvz = 0;
			for(k=0;k<nsrc;k++){
				int id = work[k];
				if(particles[id].type != TYPE_DM) continue;
				tmass += particles[id].mass;
				cx += particles[id].x * particles[id].mass;
				cy += particles[id].y * particles[id].mass;
				cz += particles[id].z * particles[id].mass;
				cvx += particles[id].vx * particles[id].mass;
				cvy += particles[id].vy * particles[id].mass;
				cvz += particles[id].vz * particles[id].mass;
				n_nucleus++;
			}
		}
		if(!(tmass>0) || n_nucleus < min_nucleus) continue;
		cx/=tmass; cy/=tmass; cz/=tmass;
		cvx/=tmass; cvy/=tmass; cvz/=tmass;
		{
			Vector3d *pos;
			float *pmass, *pe;
			int nkeep = 0;
			pos = (Vector3d *)Malloc(sizeof(Vector3d)*(size_t)nsrc, PPTR(pos));
			pmass = (float *)Malloc(sizeof(float)*(size_t)nsrc, PPTR(pmass));
			pe = (float *)Malloc(sizeof(float)*(size_t)nsrc, PPTR(pe));
			for(k=0;k<nsrc;k++){
				int id = work[k];
				pos[k].x = particles[id].x;
				pos[k].y = particles[id].y;
				pos[k].z = particles[id].z;
				pmass[k] = particles[id].mass;
			}
			potential_from_members(nsrc, pos, nsrc, pos, pmass, pe);
			ndm_bound = 0;
			for(k=0;k<nsrc;k++){
				int id = work[k];
				float ex, ey, ez, vx, vy, vz, ke;
				double px = particles[id].x - peaks[i].x;
				double py = particles[id].y - peaks[i].y;
				double pz = particles[id].z - peaks[i].z;
				int in_nucleus = particles[id].type==TYPE_DM
					&& px*px+py*py+pz*pz <= (double)NUCLEUS_RADIUS*(double)NUCLEUS_RADIUS;
				ex = (particles[id].x - cx) * r1kineticfact;
				ey = (particles[id].y - cy) * r1kineticfact;
				ez = (particles[id].z - cz) * r1kineticfact;
				vx = ex + (particles[id].vx - cvx) * r2kineticfact;
				vy = ey + (particles[id].vy - cvy) * r2kineticfact;
				vz = ez + (particles[id].vz - cvz) * r2kineticfact;
				ke = 0.5f*(vx*vx+vy*vy+vz*vz);
				if(in_nucleus || !(ke > -pe[k])){
					work[nkeep++] = id;
					if(particles[id].type==TYPE_DM) ndm_bound++;
				}
			}
			nbound = nkeep;
			Free(pe); Free(pmass); Free(pos);
		}
		if(ndm_bound < min_dm) continue;
		if(nsrc>0 && (double)nbound / (double)nsrc < 0.5) continue;
		if(n_cores >= halo.max_cores-1) break;
		{
			Coretype *core = cores + n_cores;
			int p = -1;
			memset(core, 0, sizeof(*core));
			for(k=0;k<nbound;k++){
				if(particles[work[k]].type==TYPE_DM){ p = work[k]; break; }
			}
			core->is_dark = 1;
			core->peak_particle = p;
			core->position.x = (float)cx;
			core->position.y = (float)cy;
			core->position.z = (float)cz;
			core->velocity.x = (float)cvx;
			core->velocity.y = (float)cvy;
			core->velocity.z = (float)cvz;
			core->peak_density = peaks[i].den;
			core->saddle_density = (float)rho;
			core->tidal_radius = (float)rt;
			core->n_particles = nbound;
			commit_members(particles, work, nbound, n_cores);
			DEBUGPRINT("dark C%d bound %d dm %d rt=%g d=%g m=%g c %g %g %g\n",
					n_cores, nbound, ndm_bound, rt, dist, m_ex,
					cx+halo.origin.x, cy+halo.origin.y, cz+halo.origin.z);
			n_cores++;
			n_keep++;
		}
	}
		if(local_shell){
			Free(local_prefix);
			Free(local_shell);
		}
	}
	DEBUGPRINT("dark galaxies kept %d\n", n_keep);
	Free(dmprefix);
	Free(dmshell);
	Free(peaks);
	(void)qcell;
	return n_cores;
}

int assign_members_from_watershed(SimpleBasicParticleType *particles, int n_particles,
		long long *neighbor, int n_neighbors, Coretype *cores, int n_cores) {
	int i, host, nx, ny, nz, ncell;
	int *head, *next, *csr, *csr_off, *filled, *work;
	float *order;
	float xmin, ymin, zmin, xmax, ymax, zmax, cell, span;
	RadiusMass *hshell;
	double *hprefix;
	assign_uphill_basins(particles, n_particles, neighbor, n_neighbors, cores, n_cores);
	refresh_core_counts(particles, n_particles, cores, n_cores);
	host = host_by_peak(cores, n_cores);
	DEBUGPRINT("membership host C%d peak_density=%g (labels are not m or M)\n",
			host, cores[host].peak_density);

	csr_off = (int *)Malloc(sizeof(int)*(size_t)(n_cores+1), PPTR(csr_off));
	filled = (int *)Malloc(sizeof(int)*(size_t)n_cores, PPTR(filled));
	for(i=0;i<n_cores;i++) filled[i] = 0;
	for(i=0;i<n_particles;i++){
		int h = halo.particles[i].galaxy_id;
		if(h >= 0 && h < n_cores) filled[h]++;
	}
	csr_off[0] = 0;
	for(i=0;i<n_cores;i++) csr_off[i+1] = csr_off[i] + filled[i];
	csr = (int *)Malloc(sizeof(int)*(size_t)(csr_off[n_cores] > 0 ? csr_off[n_cores] : 1), PPTR(csr));
	for(i=0;i<n_cores;i++) filled[i] = 0;
	for(i=0;i<n_particles;i++){
		int h = halo.particles[i].galaxy_id;
		if(h >= 0 && h < n_cores) csr[csr_off[h] + filled[h]++] = i;
	}
	Free(filled);

	xmin = xmax = particles[0].x;
	ymin = ymax = particles[0].y;
	zmin = zmax = particles[0].z;
	for(i=1;i<n_particles;i++){
		if(particles[i].x < xmin) xmin = particles[i].x;
		if(particles[i].y < ymin) ymin = particles[i].y;
		if(particles[i].z < zmin) zmin = particles[i].z;
		if(particles[i].x > xmax) xmax = particles[i].x;
		if(particles[i].y > ymax) ymax = particles[i].y;
		if(particles[i].z > zmax) zmax = particles[i].z;
	}
	span = xmax - xmin;
	if(ymax - ymin > span) span = ymax - ymin;
	if(zmax - zmin > span) span = zmax - zmin;
	cell = 0.01f;
	if(span > 0.f && span / cell > 250.f) cell = span / 250.f;
	nx = (int)((xmax - xmin) / cell) + 1;
	ny = (int)((ymax - ymin) / cell) + 1;
	nz = (int)((zmax - zmin) / cell) + 1;
	if(nx < 1) nx = 1;
	if(ny < 1) ny = 1;
	if(nz < 1) nz = 1;
	ncell = nx * ny * nz;
	head = (int *)Malloc(sizeof(int)*(size_t)ncell, PPTR(head));
	next = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(next));
	for(i=0;i<ncell;i++) head[i] = -1;
	for(i=0;i<n_particles;i++){
		int ix = (int)((particles[i].x - xmin) / cell);
		int iy = (int)((particles[i].y - ymin) / cell);
		int iz = (int)((particles[i].z - zmin) / cell);
		int c;
		if(ix < 0) ix = 0;
		if(iy < 0) iy = 0;
		if(iz < 0) iz = 0;
		if(ix >= nx) ix = nx - 1;
		if(iy >= ny) iy = ny - 1;
		if(iz >= nz) iz = nz - 1;
		c = (iz * ny + iy) * nx + ix;
		next[i] = head[c];
		head[c] = i;
	}

	hshell = (RadiusMass *)Malloc(sizeof(RadiusMass)*(size_t)n_particles, PPTR(hshell));
	hprefix = (double *)Malloc(sizeof(double)*(size_t)(n_particles+1), PPTR(hprefix));
	{
		double hx = cores[host].position.x;
		double hy = cores[host].position.y;
		double hz = cores[host].position.z;
		for(i=0;i<n_particles;i++){
			double dx = particles[i].x - hx;
			double dy = particles[i].y - hy;
			double dz = particles[i].z - hz;
			hshell[i].radius = (float)sqrt(dx*dx + dy*dy + dz*dz);
			hshell[i].mass = particles[i].mass;
		}
	}
	qsort(hshell, (size_t)n_particles, sizeof(RadiusMass), radius_mass_cmp);
	hprefix[0] = 0;
	for(i=0;i<n_particles;i++) hprefix[i+1] = hprefix[i] + hshell[i].mass;

	work = (int *)Malloc(sizeof(int)*(size_t)n_particles*2, PPTR(work));
	order = (float *)Malloc(sizeof(float)*(size_t)n_cores*2, PPTR(order));
	for(i=0;i<n_cores;i++){
		order[2*i] = (float)i;
		order[2*i+1] = (i == host) ? 1.e30f : -cores[i].peak_density;
	}
	qsort(order, (size_t)n_cores, sizeof(float)*2, basin_owner_cmp);
	{
		StageTimer stage_sats;
		stage_begin(&stage_sats);
	for(i=0;i<n_cores;i++){
		int core_id = (int)order[2*i];
		if(core_id == host) continue;
		decide_satellite(particles, n_particles, cores, core_id, host,
				head, next, nx, ny, nz, xmin, ymin, zmin, cell,
				hshell, hprefix, n_particles, csr, csr_off, work);
	}
		stage_end(&stage_sats, "stellar_satellites");
	}
	{
		StageTimer stage_dark;
		stage_begin(&stage_dark);
	n_cores = add_dark_galaxies(particles, n_particles, cores, n_cores, host,
			head, next, nx, ny, nz, xmin, ymin, zmin, cell,
			hshell, hprefix, n_particles, work);
		stage_end(&stage_dark, "dark_galaxies");
	}
	{
		StageTimer stage_host;
		stage_begin(&stage_host);
	decide_host(particles, n_particles, cores, host,
			head, next, nx, ny, nz, xmin, ymin, zmin, cell,
			csr, csr_off, work);
		stage_end(&stage_host, "host_membership");
	}
	Free(order);
	Free(work);
	Free(hprefix);
	Free(hshell);
	Free(next);
	Free(head);
	Free(csr);
	Free(csr_off);

	refresh_core_counts(particles, n_particles, cores, n_cores);
	/* Stellar linking is the isodensity cut. Failures stay in the pool.
	 * The host pass has already finished, so it does not adopt them. */
	link_galaxy_members(particles, n_particles, n_cores, cores);
	refresh_core_counts(particles, n_particles, cores, n_cores);
	return n_cores;
}

/* A friends-of-friends halo with no stars. The brightest dark-matter
 * density peak is the host. The other peaks are satellites, with the
 * same excess-mass and Jacobi test used for star-poor cores. */
int assign_dmo_halos(SimpleBasicParticleType *particles, int n_particles,
		Coretype *cores) {
	int i, host, nx, ny, nz, ncell, n_cores;
	int *head, *next, *csr, *csr_off, *filled, *work;
	float xmin, ymin, zmin, xmax, ymax, zmax, cell, span;
	RadiusMass dummy;
	double pref[2];
	StageTimer stage_dark, stage_host;
	xmin = xmax = particles[0].x;
	ymin = ymax = particles[0].y;
	zmin = zmax = particles[0].z;
	for(i=1;i<n_particles;i++){
		if(particles[i].x < xmin) xmin = particles[i].x;
		if(particles[i].y < ymin) ymin = particles[i].y;
		if(particles[i].z < zmin) zmin = particles[i].z;
		if(particles[i].x > xmax) xmax = particles[i].x;
		if(particles[i].y > ymax) ymax = particles[i].y;
		if(particles[i].z > zmax) zmax = particles[i].z;
	}
	span = xmax - xmin;
	if(ymax - ymin > span) span = ymax - ymin;
	if(zmax - zmin > span) span = zmax - zmin;
	cell = 0.01f;
	if(span > 0.f && span / cell > 250.f) cell = span / 250.f;
	nx = (int)((xmax - xmin) / cell) + 1;
	ny = (int)((ymax - ymin) / cell) + 1;
	nz = (int)((zmax - zmin) / cell) + 1;
	if(nx < 1) nx = 1;
	if(ny < 1) ny = 1;
	if(nz < 1) nz = 1;
	ncell = nx * ny * nz;
	head = (int *)Malloc(sizeof(int)*(size_t)ncell, PPTR(head));
	next = (int *)Malloc(sizeof(int)*(size_t)n_particles, PPTR(next));
	for(i=0;i<ncell;i++) head[i] = -1;
	for(i=0;i<n_particles;i++){
		int ix = (int)((particles[i].x - xmin) / cell);
		int iy = (int)((particles[i].y - ymin) / cell);
		int iz = (int)((particles[i].z - zmin) / cell);
		int c;
		if(ix < 0) ix = 0;
		if(iy < 0) iy = 0;
		if(iz < 0) iz = 0;
		if(ix >= nx) ix = nx - 1;
		if(iy >= ny) iy = ny - 1;
		if(iz >= nz) iz = nz - 1;
		c = (iz * ny + iy) * nx + ix;
		next[i] = head[c];
		head[c] = i;
	}
	work = (int *)Malloc(sizeof(int)*(size_t)n_particles*2, PPTR(work));
	dummy.radius = 0.f;
	dummy.mass = 0.f;
	pref[0] = pref[1] = 0;
	stage_begin(&stage_dark);
	n_cores = add_dark_galaxies(particles, n_particles, cores, 0, 0,
			head, next, nx, ny, nz, xmin, ymin, zmin, cell,
			&dummy, pref, 1, work);
	stage_end(&stage_dark, "dmo_peaks");
	if(n_cores <= 0){
		Free(work); Free(next); Free(head);
		return 0;
	}
	host = 0;
	csr_off = (int *)Malloc(sizeof(int)*(size_t)(n_cores+1), PPTR(csr_off));
	filled = (int *)Malloc(sizeof(int)*(size_t)n_cores, PPTR(filled));
	for(i=0;i<n_cores;i++) filled[i] = 0;
	for(i=0;i<n_particles;i++){
		int h = halo.particles[i].galaxy_id;
		if(h >= 0 && h < n_cores) filled[h]++;
	}
	csr_off[0] = 0;
	for(i=0;i<n_cores;i++) csr_off[i+1] = csr_off[i] + filled[i];
	csr = (int *)Malloc(sizeof(int)*(size_t)(csr_off[n_cores] > 0 ? csr_off[n_cores] : 1), PPTR(csr));
	for(i=0;i<n_cores;i++) filled[i] = 0;
	for(i=0;i<n_particles;i++){
		int h = halo.particles[i].galaxy_id;
		if(h >= 0 && h < n_cores) csr[csr_off[h] + filled[h]++] = i;
	}
	stage_begin(&stage_host);
	decide_host(particles, n_particles, cores, host,
			head, next, nx, ny, nz, xmin, ymin, zmin, cell,
			csr, csr_off, work);
	stage_end(&stage_host, "dmo_host");
	Free(filled); Free(csr); Free(csr_off);
	Free(work); Free(next); Free(head);
	refresh_core_counts(particles, n_particles, cores, n_cores);
	link_galaxy_members(particles, n_particles, n_cores, cores);
	refresh_core_counts(particles, n_particles, cores, n_cores);
	DEBUGPRINT("DMO catalogue cores %d\n", n_cores);
	return n_cores;
}
