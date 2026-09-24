/* Copy the raw FoF particle record and, when asked, write an intermediate catalogue. */
#include "finder_internal.h"

double copy_raw_particles(FoFTPtlStruct *raw,SimpleBasicParticleType *particles,lint n_particles){
	lint i;
	double xmin,ymin,zmin;
	double tmass = 0;
	xmin = ymin = zmin =  1.E30L;
	for(i=0;i<n_particles;i++){
		xmin = MIN(xmin,raw[i].x);
		ymin = MIN(ymin,raw[i].y);
		zmin = MIN(zmin,raw[i].z);
	}
	halo.origin.x = xmin;
	halo.origin.y = ymin;
	halo.origin.z = zmin;
	for(i=0;i<n_particles;i++){
		particles[i].type = raw[i].type;
		particles[i].mass = raw[i].mass;
		particles[i].x = raw[i].x - xmin;
		particles[i].y = raw[i].y - ymin;
		particles[i].z = raw[i].z - zmin;
		particles[i].vx = raw[i].vx;
		particles[i].vy = raw[i].vy;
		particles[i].vz = raw[i].vz;
		particles[i].indx = raw[i].indx;
		particles[i].link02 = raw[i].link02;
		tmass += raw[i].mass;
	}
	return tmass;
}

float total_star_mass(SimpleBasicParticleType *particles, int n_particles){
    float starmass=0;
    int i,j,k;
    for(i=0;i<n_particles;i++){
        if(particles[i].type == TYPE_STAR)
            starmass += particles[i].mass;
    }
    return starmass;

}

int count_stars(SimpleBasicParticleType *particles, int n_particles){
    int starmass=0;
    int i,j,k;
    for(i=0;i<n_particles;i++){
        if(particles[i].type == TYPE_STAR)
            starmass += 1;
    }
    return starmass;

}

float total_dm_mass(SimpleBasicParticleType *particles, int n_particles){
    float dmmass=0;
    int i;
    for(i=0;i<n_particles;i++){
        if(particles[i].type == TYPE_DM)
            dmmass += particles[i].mass;
    }
    return dmmass;
}

int _dump_stage_enabled(const char *tag){
    const char *list = getenv("NEWGAL_DUMP_STAGES");
    if(!list || !*list) return 0;
    size_t tag_len = strlen(tag);
    const char *p = list;
    while(*p){
        const char *e = strchr(p, ',');
        size_t n = e ? (size_t)(e - p) : strlen(p);
        if(n == tag_len && strncmp(p, tag, n) == 0) return 1;
        if(!e) break;
        p = e + 1;
    }
    return 0;
}

void write_stage_snapshot(FoFTPtlStruct *raw, lint n_particles,
                                int n_cores, Coretype *cores,
                                const char *tag){
    const char *outdir = getenv("NEWGAL_DUMP_CORES_DIR");
    if(outdir == NULL || *outdir == '\0') return;
    if(!_dump_stage_enabled(tag)) return;

    long min_np = 1000000L;
    const char *env_minnp = getenv("NEWGAL_DUMP_CORES_MIN_NP");
    if(env_minnp && *env_minnp){
        long v = atol(env_minnp);
        if(v > 0) min_np = v;
    }
    if((long)n_particles < min_np) return;

    char fwrp[512], fwlist[512], fwbp[512], fwflag[512];
    snprintf(fwrp,   sizeof(fwrp),   "%s/GALFIND.DATA.cores.%s.%d",    outdir, tag, myid);
    snprintf(fwlist, sizeof(fwlist), "%s/GALCATALOG.LIST.cores.%s.%d", outdir, tag, myid);
    snprintf(fwbp,   sizeof(fwbp),   "%s/background_ptl.cores.%s.%d",  outdir, tag, myid);
    snprintf(fwflag, sizeof(fwflag), "%s/iscore_flag.cores.%s.%d",     outdir, tag, myid);

    FILE *wrp   = fopen(fwrp,   "ab");
    FILE *wlist = fopen(fwlist, "ab");
    FILE *wbp   = fopen(fwbp,   "ab");
    FILE *wflag = fopen(fwflag, "ab");
    if(!wrp || !wlist || !wbp || !wflag){
        fprintf(stderr,"write_stage_snapshot: fopen failed in %s\n", outdir);
        if(wrp)   fclose(wrp);
        if(wlist) fclose(wlist);
        if(wbp)   fclose(wbp);
        if(wflag) fclose(wflag);
        return;
    }

    /* HaloInfo: aggregate over all particles (matches GetHaloInfo). */
    HaloInfo hinfo;
    hinfo.nsub = n_cores;
    hinfo.ndm = hinfo.ngas = hinfo.nsink = hinfo.nstar = 0;
    hinfo.totm = hinfo.mdm = hinfo.mgas = hinfo.msink = hinfo.mstar = 0;
    hinfo.x = hinfo.y = hinfo.z = 0;
    hinfo.vx = hinfo.vy = hinfo.vz = 0;
    {
        lint i;
        for(i=0;i<n_particles;i++){
            dptype mass=0, vx=0, vy=0, vz=0;
            if(raw[i].type == TYPE_DM){
                hinfo.ndm++;
                mass = raw[i].p.dm.mass; hinfo.mdm += mass;
                vx = raw[i].p.dm.vx; vy = raw[i].p.dm.vy; vz = raw[i].p.dm.vz;
            } else if(raw[i].type == TYPE_SINK){
                hinfo.nsink++;
                mass = raw[i].p.sink.mass; hinfo.msink += mass;
                vx = raw[i].p.sink.vx; vy = raw[i].p.sink.vy; vz = raw[i].p.sink.vz;
            } else if(raw[i].type == TYPE_STAR){
                hinfo.nstar++;
                mass = raw[i].p.star.mass; hinfo.mstar += mass;
                vx = raw[i].p.star.vx; vy = raw[i].p.star.vy; vz = raw[i].p.star.vz;
            } else if(raw[i].type == TYPE_GAS){
                hinfo.ngas++;
                mass = raw[i].p.gas.mass; hinfo.mgas += mass;
                vx = raw[i].p.gas.vx; vy = raw[i].p.gas.vy; vz = raw[i].p.gas.vz;
            }
            hinfo.totm += mass;
            hinfo.x  += mass*raw[i].x;
            hinfo.y  += mass*raw[i].y;
            hinfo.z  += mass*raw[i].z;
            hinfo.vx += mass*vx;
            hinfo.vy += mass*vy;
            hinfo.vz += mass*vz;
        }
    }
    hinfo.npall = (int)n_particles;
    if(hinfo.totm > 0){
        hinfo.x/=hinfo.totm; hinfo.y/=hinfo.totm; hinfo.z/=hinfo.totm;
        hinfo.vx/=hinfo.totm; hinfo.vy/=hinfo.totm; hinfo.vz/=hinfo.totm;
    }
    fwrite(&hinfo, sizeof(HaloInfo), 1, wrp);
    fwrite(&hinfo, sizeof(HaloInfo), 1, wlist);

    size_t maxbytes = sizeof(FoFTPtlStruct) * (size_t)n_particles;
    void *subdata = malloc(maxbytes);
    unsigned char *flagbuf = (unsigned char *)malloc((size_t)n_particles);
    if(!subdata || !flagbuf){
        fprintf(stderr,"write_stage_snapshot: malloc failed\n");
        if(subdata) free(subdata);
        if(flagbuf) free(flagbuf);
        fclose(wrp); fclose(wlist); fclose(wbp); fclose(wflag);
        return;
    }

    /* Per-cores SubInfo + particle payload. */
    {
        int k; lint i;
        for(k=0;k<n_cores;k++){
            SubInfo s;
            size_t ndm=0, ngas=0, nsink=0, nstar=0;
            for(i=0;i<n_particles;i++){
                if(halo.particles[i].galaxy_id != k) continue;
                if(raw[i].type == TYPE_DM)        ndm++;
                else if(raw[i].type == TYPE_GAS)  ngas++;
                else if(raw[i].type == TYPE_SINK) nsink++;
                else if(raw[i].type == TYPE_STAR) nstar++;
            }
            s.npdm = (int)ndm; s.npgas = (int)ngas;
            s.npsink = (int)nsink; s.npstar = (int)nstar;
            s.npall = (int)(ndm + ngas + nsink + nstar);
            s.totm = s.mdm = s.mgas = s.msink = s.mstar = 0;
            s.x = s.y = s.z = s.vx = s.vy = s.vz = 0;

            char *p_dm   = (char*)subdata;
            char *p_gas  = p_dm   + ndm  * sizeof(DmType);
            char *p_sink = p_gas  + ngas * sizeof(GasType);
            char *p_star = p_sink + nsink* sizeof(SinkType);
            unsigned char *f_dm   = flagbuf;
            unsigned char *f_gas  = f_dm   + ndm;
            unsigned char *f_sink = f_gas  + ngas;
            unsigned char *f_star = f_sink + nsink;

            for(i=0;i<n_particles;i++){
                if(halo.particles[i].galaxy_id != k) continue;
                dptype mass=0;
                unsigned char b = (is_core_particle(i) != NOT) ? 1 : 0;
                if(raw[i].type == TYPE_DM){
                    memcpy(p_dm, &raw[i].p.dm, sizeof(DmType));
                    p_dm += sizeof(DmType);
                    *f_dm++ = b;
                    mass = raw[i].p.dm.mass; s.mdm += mass;
                } else if(raw[i].type == TYPE_GAS){
                    memcpy(p_gas, &raw[i].p.gas, sizeof(GasType));
                    p_gas += sizeof(GasType);
                    *f_gas++ = b;
                    mass = raw[i].p.gas.mass; s.mgas += mass;
                } else if(raw[i].type == TYPE_SINK){
                    memcpy(p_sink, &raw[i].p.sink, sizeof(SinkType));
                    p_sink += sizeof(SinkType);
                    *f_sink++ = b;
                    mass = raw[i].p.sink.mass; s.msink += mass;
                } else if(raw[i].type == TYPE_STAR){
                    memcpy(p_star, &raw[i].p.star, sizeof(StarType));
                    p_star += sizeof(StarType);
                    *f_star++ = b;
                    mass = raw[i].p.star.mass; s.mstar += mass;
                }
                s.totm += mass;
                s.x  += mass*raw[i].x;
                s.y  += mass*raw[i].y;
                s.z  += mass*raw[i].z;
                s.vx += mass*raw[i].vx;
                s.vy += mass*raw[i].vy;
                s.vz += mass*raw[i].vz;
            }
            if(s.totm > 0){
                s.x/=s.totm; s.y/=s.totm; s.z/=s.totm;
                s.vx/=s.totm; s.vy/=s.totm; s.vz/=s.totm;
            }
            size_t bytes = ndm  * sizeof(DmType)   + ngas  * sizeof(GasType)
                         + nsink* sizeof(SinkType) + nstar * sizeof(StarType);
            fwrite(&s, sizeof(SubInfo), 1, wrp);
            fwrite(&s, sizeof(SubInfo), 1, wlist);
            fwrite(subdata, 1, bytes, wrp);
            fwrite(flagbuf, 1, s.npall, wflag);
        }
    }

    /* Background: halo.particles[i].galaxy_id not in [0, n_cores). */
    {
        SubInfo s;
        size_t ndm=0, ngas=0, nsink=0, nstar=0;
        lint i;
        for(i=0;i<n_particles;i++){
            int h = halo.particles[i].galaxy_id;
            if(h >= 0 && h < n_cores) continue;
            if(raw[i].type == TYPE_DM)        ndm++;
            else if(raw[i].type == TYPE_GAS)  ngas++;
            else if(raw[i].type == TYPE_SINK) nsink++;
            else if(raw[i].type == TYPE_STAR) nstar++;
        }
        s.npdm = (int)ndm; s.npgas = (int)ngas;
        s.npsink = (int)nsink; s.npstar = (int)nstar;
        s.npall = (int)(ndm + ngas + nsink + nstar);
        s.totm = s.mdm = s.mgas = s.msink = s.mstar = 0;
        s.x = s.y = s.z = s.vx = s.vy = s.vz = 0;

        char *p_dm   = (char*)subdata;
        char *p_gas  = p_dm   + ndm  * sizeof(DmType);
        char *p_sink = p_gas  + ngas * sizeof(GasType);
        char *p_star = p_sink + nsink* sizeof(SinkType);
        unsigned char *f_dm   = flagbuf;
        unsigned char *f_gas  = f_dm   + ndm;
        unsigned char *f_sink = f_gas  + ngas;
        unsigned char *f_star = f_sink + nsink;

        for(i=0;i<n_particles;i++){
            int h = halo.particles[i].galaxy_id;
            if(h >= 0 && h < n_cores) continue;
            dptype mass=0;
            unsigned char b = (is_core_particle(i) != NOT) ? 1 : 0;
            if(raw[i].type == TYPE_DM){
                memcpy(p_dm, &raw[i].p.dm, sizeof(DmType));
                p_dm += sizeof(DmType);
                *f_dm++ = b;
                mass = raw[i].p.dm.mass; s.mdm += mass;
            } else if(raw[i].type == TYPE_GAS){
                memcpy(p_gas, &raw[i].p.gas, sizeof(GasType));
                p_gas += sizeof(GasType);
                *f_gas++ = b;
                mass = raw[i].p.gas.mass; s.mgas += mass;
            } else if(raw[i].type == TYPE_SINK){
                memcpy(p_sink, &raw[i].p.sink, sizeof(SinkType));
                p_sink += sizeof(SinkType);
                *f_sink++ = b;
                mass = raw[i].p.sink.mass; s.msink += mass;
            } else if(raw[i].type == TYPE_STAR){
                memcpy(p_star, &raw[i].p.star, sizeof(StarType));
                p_star += sizeof(StarType);
                *f_star++ = b;
                mass = raw[i].p.star.mass; s.mstar += mass;
            }
            s.totm += mass;
            s.x  += mass*raw[i].x;
            s.y  += mass*raw[i].y;
            s.z  += mass*raw[i].z;
            s.vx += mass*raw[i].vx;
            s.vy += mass*raw[i].vy;
            s.vz += mass*raw[i].vz;
        }
        if(s.totm > 0){
            s.x/=s.totm; s.y/=s.totm; s.z/=s.totm;
            s.vx/=s.totm; s.vy/=s.totm; s.vz/=s.totm;
        }
        size_t bytes = ndm  * sizeof(DmType)   + ngas  * sizeof(GasType)
                     + nsink* sizeof(SinkType) + nstar * sizeof(StarType);
        fwrite(&s, sizeof(SubInfo), 1, wbp);
        fwrite(subdata, 1, bytes, wbp);
        fwrite(flagbuf, 1, s.npall, wflag);
    }

    free(subdata);
    free(flagbuf);
    fflush(wrp); fflush(wlist); fflush(wbp); fflush(wflag);
    fclose(wrp); fclose(wlist); fclose(wbp); fclose(wflag);

    fprintf(stderr,
            "[write_stage_snapshot rank=%d stage=%s] n_particles=%ld n_cores=%d -> %s\n",
            myid, tag, (long)n_particles, n_cores, outdir);
    fflush(stderr);

    const char *exit_after = getenv("NEWGAL_DUMP_CORES_EXIT_AFTER");
    if(exit_after && *exit_after && atoi(exit_after) != 0){
        fprintf(stderr,
                "[write_stage_snapshot rank=%d stage=%s] NEWGAL_DUMP_CORES_EXIT_AFTER "
                "set, calling MPI_Abort to stop the run.\n", myid, tag);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
}
