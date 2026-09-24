/* Tidal radii and energy tests. Shell growth uses the host potential; the final pass uses self-energy. */
#include "finder_internal.h"

double nfw_concentration(double mass){
	double mstar= 1.3E13L;
	double redshift;
	double res;
	redshift = amax/a - 1.L;
	res = 14.5 *pow(mass/mstar,-0.15)/(1.L+redshift);
	return res;
}

double nfw_alpha(double c,double s){
	double res;
	{
		res = 2.L-(c*s/(1.L+c*s))*(c*s/(1.L+c*s));
		res = res/(log(1.L+c*s)-c*s/(1.L+c*s));
		return res;
	}
}

double circular_velocity_squared(double mass, double Rv){
	double res;
	res = G*mass/Rv;
	return res;
}

double nfw_shape_factor(double c){
	double res;
	res =  1./(log(1.L+c)-c/(1.L+c));
	return res;
}

double virial_radius(double Mvir){
	double res;
	res = pow(3.*Mvir/(4.L*M_PI*vv*RHOC),0.3333333333333L);
	return res;
}

double background_potential(double Mvir){
	double massincgs = Mvir*Msun;
	double c = nfw_concentration(Mvir);
	double gc = nfw_shape_factor(c);
	double Rvir = virial_radius(Mvir);
	double Vcir2 = circular_velocity_squared(Mvir, Rvir)*potentfact;
	double halobgpotent = Vcir2*(1-gc*log(1+c));
	return halobgpotent;
}

double tidal_radius_factor(int n_particles, int myid, int mp, int uid, int mpall, Coretype *cores, float Mvir,
		double mratio, double bmass){
	double alpha,ellipse_contribution;
	double w2,wx,wy,wz,dx,dy,dz,dist2,dvx,dvy,dvz;
	double dist,pot,res,s,Rvir,v2,nfwc;
	double kpdratio;

	dx = cores[myid].position.x - cores[uid].position.x;
	dy = cores[myid].position.y - cores[uid].position.y;
	dz = cores[myid].position.z - cores[uid].position.z;
	dist = sqrt(dx*dx+dy*dy+dz*dz);
	dvx = cores[myid].velocity.x - cores[uid].velocity.x;
	dvy = cores[myid].velocity.y - cores[uid].velocity.y;
	dvz = cores[myid].velocity.z - cores[uid].velocity.z;
	v2 = dvx*dvx+dvy*dvy+dvz*dvz;
	wx = dy*dvz-dz*dvy;
	wy = dz*dvx-dx*dvz;
	wz = dx*dvy-dy*dvx;
	w2 = (wx*wx+wy*wy+wz*wz)/(dist*dist);
	pot = Mvir/dist;
	kpdratio = r2kineticfact*r2kineticfact/potentfact;
	ellipse_contribution = w2/pot * kpdratio;
	/* This is for the case when the minor body is not bound to the major body. */
	/* This is the case when the pair do not bound to each other. 
	 * In this case the contribution of the orbital motion to the
	 * tidal radius is assumed to be negligible. */
/*
#error This should be modified further.
*/

	Rvir = pow(3.*Mvir/(4.L*M_PI*vv*RHOC),0.3333333333333L);
	/* Change to the simulation length scale */
	/*
	Rvir = Rvir/size*rng;
	*/
	s = dist/Rvir;

	nfwc = nfw_concentration(bmass);
#ifdef OLD_BETA
	if(ellipse_contribution > 2.L){
		ellipse_contribution = 0.L;
	}
#else
	if(s>1) {
		double fact;
		fact = v2/(bmass/Rvir)*kpdratio;
		if(fact > 2.L) ellipse_contribution = 0.L;
	}
	else {
		double fact;
		fact = v2*s/(nfw_shape_factor(nfwc)*log(1.L+nfwc*s)*bmass/Rvir )*kpdratio;
		if(fact > 2.L) ellipse_contribution = 0.L;
	}
#endif
//	alpha = nfw_alpha(nfw_concentration(bmass),s);
	alpha = nfw_alpha(nfwc,s);
	res = pow(mratio/(ellipse_contribution+alpha),0.33333333333333L);
#ifdef DEBUG
	/*
	printf(" ellipse_c= %g alpha= %g c= %g s= %g\n",ellipse_contribution,alpha,nfw_concentration(mpall),s);
	*/
#endif
	return res;
}

void self_halo_potential(int n_particles,Vector3d *r,float *mass,float *penergy){
    TStruct *TREE;
    TPtlStruct *ptl;
    particle p;
    float theta2=0.5;
	dptype Mvir = 0;
    int i,j,k;
	if(n_particles < 100){
		float epsilon2 = EPSILON * EPSILON;
		for(i=0;i<n_particles;i++){
			Mvir += mass[i];
		}
		double halobgpotent= background_potential(Mvir);
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
		for(i=0;i<n_particles;i++){
			float potent = 0;
			for(j=0;j<n_particles;j++){
				float ptlmass = mass[j]; 
				float tmpx = r[i].x -  r[j].x;
				float tmpy = r[i].y -  r[j].y;
				float tmpz = r[i].z -  r[j].z;
				float dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  
				if(dist2 >0) potent += -ptlmass/sqrt(dist2+epsilon2);
			}
			penergy[i] = potent*potentfact + halobgpotent;
		}
	}
	else { 
		size_t nnode = MAX(65*10000,n_particles*0.5);
		TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE)); 
		ptl = (TPtlStruct *) Malloc(sizeof(TPtlStruct)*n_particles,PPTR(ptl)); 
		for(i=0;i<n_particles;i++){ 
			ptl[i].type = TYPE_PTL; 
			ptl[i].x = r[i].x; ptl[i].y = r[i].y; ptl[i].z = r[i].z; 
			ptl[i].mass = mass[i]; 
			ptl[i].sibling = &ptl[i+1]; 
			Mvir += mass[i]; 
		}
		double halobgpotent= background_potential(Mvir);
		ptl[n_particles-1].sibling = NULL;
	
		int recursiveflag;
		if(nnode > 65*10000) recursiveflag = PTHREAD;
		else recursiveflag = RECURSIVE;

		recursiveflag = SERIALIZED;
		build_force_tree(TREE,nnode,ptl,n_particles,theta2,recursiveflag);
#ifdef _OPENMP
#pragma omp parallel for private(i,p) schedule(guided)
#endif
		for(i=0;i<n_particles;i++){ 
			p.x = ptl[i].x; p.y = ptl[i].y; p.z = ptl[i].z; 
			penergy[i] = tree_potential(&p,theta2,TREE,ptl)*potentfact + halobgpotent; 
		} 
		Free(ptl); 
		Free(TREE);
	}
}

void external_halo_potential(int nend,Vector3d *r,float *mass, float *penergy, int snp, Vector3d *sr, float *smass){
    TStruct *TREE;
    TPtlStruct *ptl;
    particle p;
    float theta2=0.5;
    int i,j,k;
    int n_particles;

	if(nend < 1000 ){
		float epsilon2 = EPSILON * EPSILON;
		dptype Mvir = 0;
		for(i=0;i<snp;i++){
			Mvir += smass[i];
		}
		double halobgpotent= background_potential(Mvir);
		if(nend >= MAXTHREADS){
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
			for(i=0;i<nend;i++){
				float potent = 0;
				for(j=0;j<snp;j++){
					float ptlmass = smass[j]; 
					float tmpx = r[i].x -  sr[j].x;
					float tmpy = r[i].y -  sr[j].y;
					float tmpz = r[i].z -  sr[j].z;
					float dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  
					if(dist2 >0) potent += -ptlmass/sqrt(dist2+epsilon2);
				}
				penergy[i] = potent*potentfact + halobgpotent;
			}
		}
		else {
			for(i=0;i<nend;i++) penergy[i] = 0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j) reduction(+:penergy[:nend])
#endif
			for(j=0;j<snp;j++){
				float ptlmass = smass[j];
				for(i=0;i<nend;i++){
					float potent;
					float tmpx = r[i].x -  sr[j].x;
					float tmpy = r[i].y -  sr[j].y;
					float tmpz = r[i].z -  sr[j].z;
					float dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  
					if(dist2 >0) potent = -ptlmass/sqrt(dist2+epsilon2);
					penergy[i] += potent*potentfact + halobgpotent;
				}
			}
		}
	}
	else {
		size_t nnode = MAX(65*10000,snp*0.5);
	    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE)); 
		ptl = (TPtlStruct *) Malloc(sizeof(TPtlStruct)*snp,PPTR(ptl)); 
		dptype Mvir = 0; 
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(+:Mvir)
#endif
		for(i=0;i<snp;i++){ 
			ptl[i].type = TYPE_PTL; 
			ptl[i].x = sr[i].x; ptl[i].y = sr[i].y; ptl[i].z = sr[i].z; 
			ptl[i].mass = smass[i]; 
			ptl[i].sibling = &ptl[i+1]; 
			Mvir += smass[i]; 
		} 
		ptl[snp-1].sibling = NULL; 
		double halobgpotent =  background_potential(Mvir); 

		int recursiveflag;
		if(nnode > 65*10000) recursiveflag = PTHREAD;
		else recursiveflag = RECURSIVE;
		recursiveflag = SERIALIZED;
		build_force_tree(TREE,nnode,ptl,snp,theta2,recursiveflag);
#ifdef _OPENMP
#pragma omp parallel for private(i,p) schedule(dynamic)
#endif
		for(i=0;i<nend;i++){ 
			p.x = r[i].x; p.y = r[i].y; p.z = r[i].z; 
			penergy[i] = tree_potential(&p,theta2,TREE,ptl)*potentfact + halobgpotent; 
		} 
		Free(ptl); 
		Free(TREE);
	}
}

int count_bound_against_members(SampledParticle *candidate_samples,int tnp,unsigned char *bndflag,
		SampledParticle *member_samples,int snp, Coretype *one_core){
	int i,j,k,nbound;
	float *kenergy,*penergy;
	double cx,cy,cz,cvx,cvy,cvz;
	float distx,disty,distz;
	float tmpvx,tmpvy,tmpvz;
	Vector3d *tr,*tvr;
	Vector3d *sr;
	float *mass;

	tr = (Vector3d*)Malloc(sizeof(Vector3d)*tnp,PPTR(tr));
	tvr = (Vector3d*)Malloc(sizeof(Vector3d)*tnp,PPTR(tvr));
	mass = (float*)Malloc(sizeof(float)*tnp,PPTR(mass));
	kenergy = (float*)Malloc(sizeof(float)*tnp,PPTR(kenergy));
	penergy = (float*)Malloc(sizeof(float)*tnp,PPTR(penergy));

	/*
	cx = cy = cz = cvx = cvy = cvz = 0.L;
	for(i=0;i<tnp;i++){
		cx += (tr[i].x = candidate_samples[i].position.x);
		cy += (tr[i].y = candidate_samples[i].position.y);
		cz += (tr[i].z = candidate_samples[i].position.z);
		cvx += (tvr[i].x = candidate_samples[i].velocity.x);
		cvy += (tvr[i].y = candidate_samples[i].velocity.y);
		cvz += (tvr[i].z = candidate_samples[i].velocity.z);
	}
	cx = cx/(double)tnp; cy = cy/(double)tnp; cz = cz/(double)tnp;
	cvx = cvx/(double)tnp; cvy = cvy/(double)tnp; cvz = cvz/(double)tnp;
	*/
	for(i=0;i<tnp;i++){
		tr[i].x = candidate_samples[i].position.x;
		tr[i].y = candidate_samples[i].position.y;
		tr[i].z = candidate_samples[i].position.z;
		tvr[i].x = candidate_samples[i].velocity.x;
		tvr[i].y = candidate_samples[i].velocity.y;
		tvr[i].z = candidate_samples[i].velocity.z;
		mass[i] = candidate_samples[i].mass;
	}
	cx = one_core->position.x; cy = one_core->position.y; cz = one_core->position.z;
	cvx = one_core->velocity.x; cvy = one_core->velocity.y; cvz = one_core->velocity.z;

#ifdef _OPENMP
#pragma omp parallel for private(i,distx,disty, distz,tmpvx,tmpvy,tmpvz)
#endif
	for(i=0;i<tnp;i++){
		distx = (tr[i].x-cx)*r1kineticfact;
		disty = (tr[i].y-cy)*r1kineticfact;
		distz = (tr[i].z-cz)*r1kineticfact;
		tmpvx = distx + (tvr[i].x-cvx)*r2kineticfact;
		tmpvy = disty + (tvr[i].y-cvy)*r2kineticfact;
		tmpvz = distz + (tvr[i].z-cvz)*r2kineticfact;
		kenergy[i] = 0.5*(tmpvx*tmpvx+tmpvy*tmpvy+tmpvz*tmpvz);
	}

	sr = (Vector3d*)Malloc(sizeof(Vector3d)*snp,PPTR(sr));
	float *smass = (float*)Malloc(sizeof(float)*snp,PPTR(smass));
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	for(i=0;i<snp;i++){
		sr[i].x = member_samples[i].position.x; sr[i].y = member_samples[i].position.y; sr[i].z = member_samples[i].position.z;
		smass[i] = member_samples[i].mass;
	}
	if(tnp>0) external_halo_potential(tnp,tr,mass,penergy,snp,sr, smass);
	Free(smass);
	Free(sr);

	nbound = 0;
	for(i=0;i<tnp;i++){
		if(kenergy[i] > -penergy[i]) {
			bndflag[i] = UNBOUND;
		}
		else {
			bndflag[i] = BOUND;
			nbound ++;
		}
	}
	Free(penergy); Free(kenergy);
	Free(mass); 
	Free(tvr);Free(tr);
	return nbound;
}

int claim_shell_particles(SampledParticle *candidate_samples,int tnp, SampledParticle *member_samples, int snp,
		SimpleBasicParticleType *particles,int n_particles,Coretype *one_core){
	int i,j,k;
	unsigned char *bndflag;
	int n_members,omem;
	bndflag = (unsigned char *)Malloc(sizeof(unsigned char)*tnp,PPTR(bndflag));
	for(i=0;i<tnp;i++) bndflag[i] = BOUND;
	omem = count_bound_against_members(candidate_samples,tnp,bndflag,member_samples,snp,one_core);
	n_members = 0;
	for(i=0;i<tnp;i++) {
		if(bndflag[i] == BOUND) {
			mark_bound(((int)((candidate_samples)[i].particle - particles)));
			candidate_samples[n_members++] = candidate_samples[i];
		}
		else {
			mark_remaining(((int)((candidate_samples)[i].particle - particles)));
			clear_bound(((int)((candidate_samples)[i].particle - particles)));
			/* */
			set_galaxy_id(((int)((candidate_samples)[i].particle - particles)),NOT_HALO_MEMBER);
		}
	}
	Free(bndflag);
	return n_members;
}

int count_self_bound(SampledParticle *samples,int n_particles,unsigned char *bndflag, Coretype *one_core, int satellite_frame){
	int i,j,k,nbound;
	float *kenergy,*penergy;
	double cx,cy,cz,cvx,cvy,cvz;
	float distx,disty,distz;
	float tmpvx,tmpvy,tmpvz;
	Vector3d *r,*vr;
	float *mass;

	r = (Vector3d*)Malloc(sizeof(Vector3d)*n_particles,PPTR(r));
	vr = (Vector3d*)Malloc(sizeof(Vector3d)*n_particles,PPTR(vr));
	mass = (float*)Malloc(sizeof(float)*n_particles,PPTR(mass));
	kenergy = (float*)Malloc(sizeof(float)*n_particles,PPTR(kenergy));
	penergy = (float*)Malloc(sizeof(float)*n_particles,PPTR(penergy));
	for(i=0;i<n_particles;i++){
		r[i].x = samples[i].position.x;
		r[i].y = samples[i].position.y;
		r[i].z = samples[i].position.z;
		vr[i].x = samples[i].velocity.x;
		vr[i].y = samples[i].velocity.y;
		vr[i].z = samples[i].velocity.z;
		mass[i] = samples[i].mass;
	}
	{
		dptype tmass = 0;
		cx = cy = cz = cvx = cvy = cvz = 0.L;
		if(satellite_frame && one_core){
			int n_nucleus = 0;
			float nucleus_r2 = (float)NUCLEUS_RADIUS * (float)NUCLEUS_RADIUS;
			for(i=0;i<n_particles;i++){
				float dx,dy,dz,dist2;
				if(samples[i].particle->type != TYPE_STAR) continue;
				dx = r[i].x - one_core->position.x;
				dy = r[i].y - one_core->position.y;
				dz = r[i].z - one_core->position.z;
				dist2 = dx*dx + dy*dy + dz*dz;
				if(dist2 > nucleus_r2) continue;
				n_nucleus++;
				tmass += mass[i];
				cx += r[i].x*mass[i];
				cy += r[i].y*mass[i];
				cz += r[i].z*mass[i];
				cvx += vr[i].x*mass[i];
				cvy += vr[i].y*mass[i];
				cvz += vr[i].z*mass[i];
			}
			/* Own-member stars only. Too few nucleus stars must not fall
			 * back onto a host-contaminated candidate COM. */
			if(n_nucleus >= NUCLEUS_MIN_STARS && tmass > 0){
				cx /= tmass; cy /= tmass; cz /= tmass;
				cvx /= tmass; cvy /= tmass; cvz /= tmass;
			}
			else {
				cx = one_core->position.x; cy = one_core->position.y; cz = one_core->position.z;
				cvx = one_core->velocity.x; cvy = one_core->velocity.y; cvz = one_core->velocity.z;
			}
		}
		else {
			for(i=0;i<n_particles;i++){
				if(samples[i].particle->type == TYPE_STAR){
					cx +=  r[i].x*samples[i].mass;
					cy +=  r[i].y*samples[i].mass;
					cz +=  r[i].z*samples[i].mass;
					cvx += vr[i].x*samples[i].mass;
					cvy += vr[i].y*samples[i].mass;
					cvz += vr[i].z*samples[i].mass;
					tmass += samples[i].mass;
				}
			}
			cx = cx/tmass; cy = cy/tmass; cz = cz/tmass;
			cvx = cvx/tmass; cvy = cvy/tmass; cvz = cvz/tmass;
		}
	}

	for(i=0;i<n_particles;i++){
		distx = (r[i].x-cx)*r1kineticfact;
		disty = (r[i].y-cy)*r1kineticfact;
		distz = (r[i].z-cz)*r1kineticfact;
		tmpvx = distx + (vr[i].x-cvx)*r2kineticfact;
		tmpvy = disty + (vr[i].y-cvy)*r2kineticfact;
		tmpvz = distz + (vr[i].z-cvz)*r2kineticfact;
		kenergy[i] = 0.5*(tmpvx*tmpvx+tmpvy*tmpvy+tmpvz*tmpvz);
	}
	self_halo_potential(n_particles,r,mass,penergy);
	nbound = 0;
	for(i=0;i<n_particles;i++){
		if(kenergy[i] > -penergy[i]) {
			bndflag[i] = UNBOUND;
		}
		else {
			bndflag[i] = BOUND;
			nbound ++;
		}
	}
	/*
#ifdef DEBUG
	for(j=0;j<n_particles;j++){
		printf("P%d of n_particles = %d is ke %g and pe %g bnd? = %d\n",j,n_particles,
				kenergy[j],penergy[j],bndflag[j]);
	}
#endif
*/
	Free(penergy); Free(kenergy);
	Free(mass);
	Free(vr);Free(r);
	DEBUGPRINT(" Total bound : %d  among %d\n", nbound,n_particles);
	return nbound;
}

int unbind_isolated_halo(SampledParticle *samples,int n_particles,SimpleBasicParticleType *particles,int galaxy_id){
	int i,j,k;
	unsigned char *bndflag;
	int n_members,onmem;
	bndflag = (unsigned char *)Malloc(sizeof(unsigned char)*n_particles,PPTR(bndflag));
	for(i=0;i<n_particles;i++) bndflag[i] = BOUND;
	onmem = n_particles;
	for(i=0;i<BOUNDITER;i++){
		n_members = 0;
		for(j=0;j<onmem;j++){
			if(bndflag[j] == BOUND){
				samples[n_members++] = samples[j];
			}
		}
		if(n_members==0) {
			Free(bndflag);
			return n_members;
		}
		n_members = count_self_bound(samples,(k=n_members),bndflag,NULL,0);
#ifdef DEBUG
		printf("iterating to find bound particles %d:    %d from %d\n",i,n_members,onmem);
#endif
		/*
		n_members = count_self_bound(samples,(k=n_members),bndflag);
		*/
		if(n_members == onmem  || n_members ==0) {
			break;
		}
		onmem = n_members;
	}
	Free(bndflag);
	for(i=0;i<n_members;i++) set_galaxy_id(((int)((samples)[i].particle - particles)),galaxy_id);
//	for(i=0;i<n_particles;i++) set_galaxy_id(i,galaxy_id);
#ifdef DEBUG
	printf("alone halo loses member particles from %d to %d\n",n_particles,n_members);
#endif
	return n_members;
}

int unbind_one_galaxy(int nkp,SampledParticle *samples,int n_particles,SimpleBasicParticleType *particles,int galaxy_id, Coretype *cores, int satellite_frame){
	int i,j,k;
	unsigned char *bndflag;
	int n_members;
	bndflag = (unsigned char*)Malloc(sizeof(unsigned char)*n_particles,PPTR(bndflag));
	for(i=0;i<nkp;i++) bndflag[i] = BOUND;

	n_members = count_self_bound(samples,nkp,bndflag,cores,satellite_frame);

	for(i=0;i<nkp;i++) {
		if(bndflag[i]==BOUND) {
			set_galaxy_id(((int)((samples)[i].particle - particles)),galaxy_id);
			clear_remaining(((int)((samples)[i].particle - particles)));
		}
		else if(halo.particles[((int)((samples)[i].particle - particles))].galaxy_id == galaxy_id){
			set_galaxy_id(((int)((samples)[i].particle - particles)),NOT_HALO_MEMBER);
			mark_remaining(((int)((samples)[i].particle - particles)));
		}
	}
	Free(bndflag);
	return n_members;
}

int compare_cores_by_mass(const void *a,const void *b){
	/* Sort cores ascending by current cumulative member mass (tmass).
	 * Previously sorted by member count (n_members); switched so the shell
	 * loop processes the lightest cores first regardless of how
	 * non-uniform per-particle stellar masses are. n_members and tmass
	 * agree only when star-particle masses are uniform. */
	Coresorttype *aa,*bb;
	aa = (Coresorttype *)a;
	bb = (Coresorttype *)b;
	if(aa->tmass < bb->tmass) return -1;
	else if(aa->tmass > bb->tmass) return +1;
	else return 0;
}

void init_nfw_tidal_table(void){
}

void set_tidal_radii(Coretype *cores,int n_cores,
		Coresorttype *core_order,int ncore,SimpleBasicParticleType *particles,int n_particles){
	int i,j,k;
#ifdef OLD_TIDAL
	if(iflag ==1) {
		void mkRtidal(void);
		(void)mkRtidal();
		iflag = 0;
	}
#endif

	MemberCloud *core2member;
	core2member = (MemberCloud *)Malloc(sizeof(MemberCloud)*n_cores,PPTR(core2member));
	for(i=0;i<n_cores;i++) {
		core2member[i].n = 0;
		core2member[i].total_mass = 0;
	}
	for(i=0;i<n_particles;i++) if(halo.particles[i].galaxy_id>=0) (core2member[halo.particles[i].galaxy_id].n)++;
//	DEBUGPRINT(" C2 cores has %d members\n", core2member[2].n);
	for(i=0;i<n_cores;i++){
		if(core2member[i].n ==0) {
			core2member[i].positions = NULL;
			core2member[i].masses = NULL;
		}
		else {
			core2member[i].positions = (Vector3d*) Malloc(sizeof(Vector3d)*core2member[i].n,PPTR(core2member[i].positions));
			core2member[i].masses = (float*) Malloc(sizeof(float)*core2member[i].n,PPTR(core2member[i].masses));
			core2member[i].n = 0;
			core2member[i].total_mass = 0;
		}
	}
	for(i=0;i<n_particles;i++)
		if(halo.particles[i].galaxy_id>=0) {
			k = halo.particles[i].galaxy_id;
			((core2member[k].positions)+core2member[k].n)->x = particles[i].x;
			((core2member[k].positions)+core2member[k].n)->y = particles[i].y;
			((core2member[k].positions)+core2member[k].n)->z = particles[i].z;
			*((core2member[k].masses)+core2member[k].n) = particles[i].mass;
			(core2member[k].n)++;
			core2member[k].total_mass += particles[i].mass;
		}
#ifdef _OPENMP
#pragma omp parallel private(i,k)
#endif
	{
#ifdef _OPENMP
		int it = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int central_id = central_stellar_core(cores, n_cores);
		for(i=it;i<ncore;i+=nthreads)
#else
		int central_id = central_stellar_core(cores, n_cores);
		for(i=0;i<ncore;i++)
#endif
		{
			int coreID = ((int)(core_order[i].galaxy - cores));
			int npCore = core2member[coreID].n;
			cores[coreID].n_particles = npCore;
			cores[coreID].tidal_radius = MAX_TIDAL_R;
			if(coreID == central_id){
				DEBUGPRINT("C%d keeps the halo tidal radius, star_mass=%g\n",
						coreID, cores[coreID].star_mass);
				continue;
			}
	
			double coreX = cores[coreID].position.x; double coreY = cores[coreID].position.y; double coreZ = cores[coreID].position.z;
			if(npCore>0){
				int jdcore;
				/* Bidirectional pair scan with explicit mass gate. The previous
				 * "jdcore=i+1" form relied on compare_cores_by_mass's n_members-ordering to act
				 * as the host/satellite ordering, but n_members and tmass disagree
				 * when star-particle masses are non-uniform — leaving the
				 * largest-n_members (not necessarily largest-mass) cores with
				 * Rtidal=MAX_TIDAL_R and applying the tidal_radius_factor formula
				 * (which assumes mratio < 1) to backwards pairs. Now every
				 * pair (i,j) is examined and the tidal bound is applied only
				 * when self is the lighter member, decoupling the loop order
				 * from the physics. */
				for(jdcore=0;jdcore<ncore;jdcore++){
					if(jdcore == i) continue;
					int countCoreID = ((int)(core_order[jdcore].galaxy - cores));
					int bnpCore = core2member[countCoreID].n;
					double bmass = core2member[countCoreID].total_mass;
					/* The stellar central is a host even if this core has more total mass. */
					if(countCoreID != central_id && core2member[coreID].total_mass >= bmass) continue;
					double bcoreX = cores[countCoreID].position.x; double bcoreY = cores[countCoreID].position.y; double bcoreZ = cores[countCoreID].position.z;
					double r2 = (bcoreX-coreX)*(bcoreX-coreX)+
						(bcoreY-coreY)*(bcoreY-coreY)+
						(bcoreZ-coreZ)*(bcoreZ-coreZ);
					if(bnpCore>0){
#ifndef OLD_TIDAL
						float Mvir=0;
						int inum = 0;
						for(k=0;k<bnpCore;k++){
							double tmpx = ((core2member[countCoreID].positions+k)->x) - bcoreX;
							double tmpy = ((core2member[countCoreID].positions+k)->y) - bcoreY;
							double tmpz = ((core2member[countCoreID].positions+k)->z) - bcoreZ;
							double dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
							if(dist2< r2) {
								Mvir += core2member[countCoreID].masses[k];
								inum++;
							}
						}
						if(inum >0){
							int mtmp = core2member[coreID].total_mass /Mvir * NUM_MASS;
							if(mtmp >= NUM_MASS) {
								cores[coreID].tidal_radius = sqrt(r2)*0.5; //  middle 
							}
							else {
								double tidal_radius;
								double mratio = core2member[coreID].total_mass/Mvir;
								tidal_radius  = tidal_radius_factor(npCore,coreID,inum,countCoreID,bnpCore,cores, Mvir, mratio,bmass);
								tidal_radius = MIN(tidal_radius,0.5);
								cores[coreID].tidal_radius = MIN(cores[coreID].tidal_radius, sqrt(r2)*tidal_radius);
							}
						}
						else{
							cores[coreID].tidal_radius = MIN(cores[coreID].tidal_radius, MAX_TIDAL_R);
						}
#else
						float mratio = core2member[coreID].total_mass/bmass;
						float Rvir = pow(3.*bmass/(4.L*M_PI*vv*RHOC),0.3333333333333L);
						float dist_over_Rv = sqrt(r2)/Rvir;
						float nfw_c = nfw_concentration(core2member[countCoreID].total_mass);
						float nfw_rtidal(float, float, float);
						float tidal_radius = nfw_rtidal(mratio, dist_over_Rv, nfw_c);
						tidal_radius = MIN(tidal_radius,0.5);
						cores[coreID].tidal_radius = MIN(cores[coreID].tidal_radius, sqrt(r2)*tidal_radius);
	
#endif
					}
				}
			}
			else {
				cores[coreID].tidal_radius = 0.;
			}
			DEBUGPRINT("C%d has tidal radius %g for mass = %g\n",coreID,cores[coreID].tidal_radius,core2member[coreID].total_mass);
		}
	}

	for(i=n_cores-1;i>=0;i--) {
		Free(core2member[i].positions);
		Free(core2member[i].masses);
	}
	Free(core2member);
}

int indexes_of_galaxy(int n_particles,int galaxy_id, int *indexes){
	int i,j,k;
	int n_members,nlcore;
	n_members  = 0;
	for(i=0;i<n_particles;i++) if(halo.particles[i].galaxy_id == galaxy_id) indexes[n_members++] = i;
	return n_members;
}

void pack_samples_parallel(SampledParticle *samples, int *list, int nlist,SimpleBasicParticleType *particles){
	int ii;
#ifdef _OPENMP
#pragma omp parallel for private(ii) schedule(guided)
#endif
	for(ii=0;ii<nlist;ii++){
		samples[ii].position.x = particles[list[ii]].x;
		samples[ii].position.y = particles[list[ii]].y;
		samples[ii].position.z = particles[list[ii]].z;
		samples[ii].velocity.x = particles[list[ii]].vx;
		samples[ii].velocity.y = particles[list[ii]].vy;
		samples[ii].velocity.z = particles[list[ii]].vz;
		samples[ii].mass = particles[list[ii]].mass;
		samples[ii].particle = particles+list[ii];
	}
}

int	list_tidal_candidates_adaptive(SimpleBasicParticleType *particles,int n_particles,int hid,int *list,
		Coretype *cores, Coresorttype *core_order,int now,int n_cores){
	int i,j,k,num,onum,inum;
	int id;
	float cx,cy,cz;
	float rt2;
	float tmpx,tmpy,tmpz,dist2;
	cx = cores[hid].position.x; cy = cores[hid].position.y; cz = cores[hid].position.z;
	rt2 = (cores[hid].tidal_radius)*(cores[hid].tidal_radius);
	
	num = 0;
	for(i=0;i<n_particles;i++){
		if(halo.particles[i].galaxy_id==hid){
			list[num++] = i;
		}
	}
	onum = num;
	for(i=0;i<n_particles;i++){
		if(halo.particles[i].galaxy_id != hid){
			tmpx =  particles[i].x-cx; tmpy =  particles[i].y-cy; tmpz =  particles[i].z-cz;
			dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
			if(dist2 < rt2) list[num++] = i;
		}
	}
    inum = onum;
	for(i=onum;i<num;i++){
		id = halo.particles[list[i]].galaxy_id;
		if(id==NOT_HALO_MEMBER) {
			list[inum++] = list[i];
		}
		else {
			for(j=now+1;j<n_cores;j++){
				if(id == ((int)(core_order[j].galaxy - cores))){
					list[inum++] = list[i];
					break;
				}
			}
		}
	}
	num = inum;
	return num;
}

int central_stellar_core(Coretype *cores, int n_cores){
	int i, best;
	float mstar;
	best = 0;
	mstar = -1.f;
	for(i=0;i<n_cores;i++){
		if(cores[i].star_mass > mstar){
			mstar = cores[i].star_mass;
			best = i;
		}
	}
	return best;
}

int buried_in_heavier_core(Coretype *cores, Coresorttype *order, int n_order, int self_index){
	Coretype *self;
	float self_mass, r_self, sx, sy, sz;
	int h;
	if(self_index < 0 || self_index >= n_order) return 0;
	self = order[self_index].galaxy;
	self_mass = order[self_index].tmass;
	r_self = self->tidal_radius;
	if(!(r_self > 0.f)) return 0;
	sx = self->position.x; sy = self->position.y; sz = self->position.z;
	for(h=0;h<n_order;h++){
		Coretype *host;
		float dx, dy, dz, r_host, cap;
		if(h == self_index) continue;
		if(!(order[h].tmass > self_mass)) continue;
		host = order[h].galaxy;
		dx = host->position.x - sx;
		dy = host->position.y - sy;
		dz = host->position.z - sz;
		r_host = host->tidal_radius;
		/* The heaviest core keeps MAX_TIDAL_R, which would mark every
		 * lighter core as buried. Cap the host sphere at 10 satellite radii. */
		cap = 10.f * r_self;
		if(!(r_host < cap)) r_host = cap;
		if(dx*dx + dy*dy + dz*dz < r_host*r_host) return 1;
	}
	return 0;
}

int list_own_members(int n_particles, int galaxy_id, int *list){
	int i, n;
	n = 0;
	for(i=0;i<n_particles;i++){
		if(halo.particles[i].galaxy_id == galaxy_id) list[n++] = i;
	}
	return n;
}

int	list_tidal_candidates(SimpleBasicParticleType *particles,int n_particles,int hid,int *list,
		Coretype *cores, Coresorttype *core_order,int now,int n_cores){
	int i,j,k,num,onum,inum;
	int id;
	float cx,cy,cz;
	float rt2;
	float tmpx,tmpy,tmpz,dist2;
	cx = cores[hid].position.x; cy = cores[hid].position.y; cz = cores[hid].position.z;
	rt2 = (cores[hid].tidal_radius)*(cores[hid].tidal_radius);
	
	num = 0;
	for(i=0;i<n_particles;i++){
		tmpx =  particles[i].x-cx; tmpy =  particles[i].y-cy; tmpz =  particles[i].z-cz;
		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		/* Include all the particles within the tidal radius irrepective of
		 * whether they are cores particles to other halo candidates. */
		if(dist2 < rt2) {
			list[num++] = i;
		}
		/* Erase from member list halo particles that are beyond the tidal radius. */
		else if(halo.particles[i].galaxy_id == hid){
			set_galaxy_id(i,NOT_HALO_MEMBER);
			mark_remaining(i);
		}
	}
	inum = 0;
	for(i=0;i<num;i++){
		id = halo.particles[list[i]].galaxy_id;
		/* If it is the member particle or free, simply include it to the list. */
		if(id==NOT_HALO_MEMBER || id == hid) {
			list[inum++] = list[i];
		}
		/* If it is not a cores particle of another halo and it is 
		 * a member particle of another larger system, then include it. */
		else if(is_core_particle(list[i]) == NOT) {
			for(j=now+1;j<n_cores;j++){
				if(id == ((int)(core_order[j].galaxy - cores))){
					list[inum++] = list[i];
					break;
				}
			}
		}
	}
	num = inum;
	return num;
}
