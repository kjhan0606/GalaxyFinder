/* Shared types for one FoF halo.
 *
 * Species and tree-node kinds stay in tree.h. They still share the value 1
 * (TYPE_PTL and TYPE_STAR). Do not compare a tree-node type with a species.
 */
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include "Memory.h"
#include "ramses.h"
#include "tree.h"
#include "defs.h"
#include "params.h"
#include "header.h"
#include "generated_macros.h"

typedef struct Vec3d {
    double x, y, z;
} Vec3d;
typedef Vec3d Vector3d;

typedef struct SampledParticle {
    Vec3d position;
    Vec3 velocity;
    float mass;
    SimpleBasicParticleType *particle;
} SampledParticle;

typedef struct Coresorttype {
    int n_members;
    float tmass;
    Coretype *galaxy;
} Coresorttype;

typedef struct Coresortdentype {
    int n_members;
    float density;
    Coretype *galaxy;
    unsigned char flag;
} Coresortdentype;

typedef struct Shelltype {
    int n_cores;
    int *core_ids;
    int n_particles;
    int *particle_ids;
} Shelltype;

typedef struct MemberCloud {
    Vec3d *positions;
    float *masses;
    int n;
    float total_mass;
} MemberCloud;

typedef struct ompLinkedList {
    int bpid;
    particle p;
} ompLinkedList;

#include "particle_state.h"

static inline void pack_samples(SampledParticle *samples, const int *indexes,
                                int n_indexes, SimpleBasicParticleType *particles) {
    int i;
    for (i = 0; i < n_indexes; i++) {
        samples[i].position.x = particles[indexes[i]].x;
        samples[i].position.y = particles[indexes[i]].y;
        samples[i].position.z = particles[indexes[i]].z;
        samples[i].velocity.x = particles[indexes[i]].vx;
        samples[i].velocity.y = particles[indexes[i]].vy;
        samples[i].velocity.z = particles[indexes[i]].vz;
        samples[i].mass = particles[indexes[i]].mass;
        samples[i].particle = particles + indexes[i];
    }
}

/* Scratch for the one FoF halo this rank is identifying. */
typedef struct HaloScratch {
    ParticleState *particles;
    Shelltype shells[1000000];
    int n_shells;
    int max_cores;
    Vec3d origin;
    Coretype *persistence_order;
} HaloScratch;

extern HaloScratch halo;

static inline void clear_particle_marks(int i) { halo.particles[i].marks = 0; }
static inline void clear_thread_marks(int i) {
    int t;
    for (t = 0; t < MAXTHREADS; t++) halo.particles[i].thread_marks[t] = 0;
}
static inline void set_galaxy_id(int i, int galaxy) { halo.particles[i].galaxy_id = galaxy; }
static inline void mark_peak(int i)      { halo.particles[i].marks |= PARTICLE_IS_PEAK; }
static inline void clear_peak(int i)     { halo.particles[i].marks &= (unsigned char)~PARTICLE_IS_PEAK; }
static inline int  is_peak(int i)        { return halo.particles[i].marks & PARTICLE_IS_PEAK; }
static inline void toggle_peak(int i)    { halo.particles[i].marks ^= PARTICLE_IS_PEAK; }
static inline void mark_visited(int i)   { halo.particles[i].marks |= PARTICLE_IS_VISITED; }
static inline void clear_visited(int i)  { halo.particles[i].marks &= (unsigned char)~PARTICLE_IS_VISITED; }
static inline int  is_visited(int i)     { return halo.particles[i].marks & PARTICLE_IS_VISITED; }
static inline void toggle_visited(int i) { halo.particles[i].marks ^= PARTICLE_IS_VISITED; }
static inline void mark_core_particle(int i)   { halo.particles[i].marks |= PARTICLE_IS_CORE; }
static inline void clear_core_particle(int i)  { halo.particles[i].marks &= (unsigned char)~PARTICLE_IS_CORE; }
static inline int  is_core_particle(int i)     { return halo.particles[i].marks & PARTICLE_IS_CORE; }
static inline void toggle_core(int i)          { halo.particles[i].marks ^= PARTICLE_IS_CORE; }
static inline void mark_shell(int i)    { halo.particles[i].marks |= PARTICLE_IS_SHELL; }
static inline void clear_shell(int i)   { halo.particles[i].marks &= (unsigned char)~PARTICLE_IS_SHELL; }
static inline int  is_shell(int i)      { return halo.particles[i].marks & PARTICLE_IS_SHELL; }
static inline void toggle_shell(int i)  { halo.particles[i].marks ^= PARTICLE_IS_SHELL; }
static inline void mark_bound(int i)    { halo.particles[i].marks |= PARTICLE_IS_BOUND; }
static inline void clear_bound(int i)   { halo.particles[i].marks &= (unsigned char)~PARTICLE_IS_BOUND; }
static inline int  is_bound(int i)      { return halo.particles[i].marks & PARTICLE_IS_BOUND; }
static inline void toggle_bound(int i)  { halo.particles[i].marks ^= PARTICLE_IS_BOUND; }
static inline void mark_remaining(int i)   { halo.particles[i].marks |= PARTICLE_IS_REMAINING; }
static inline void clear_remaining(int i)  { halo.particles[i].marks &= (unsigned char)~PARTICLE_IS_REMAINING; }
static inline int  is_remaining(int i)     { return halo.particles[i].marks & PARTICLE_IS_REMAINING; }
static inline void toggle_remaining(int i) { halo.particles[i].marks ^= PARTICLE_IS_REMAINING; }
static inline void mark_visited_thread(int i, int t)  { halo.particles[i].thread_marks[t] |= THREAD_IS_VISITED; }
static inline void clear_visited_thread(int i, int t) { halo.particles[i].thread_marks[t] &= (unsigned char)~THREAD_IS_VISITED; }
static inline int  is_visited_thread(int i, int t)    { return halo.particles[i].thread_marks[t] & THREAD_IS_VISITED; }
static inline void toggle_visited_thread(int i, int t){ halo.particles[i].thread_marks[t] ^= THREAD_IS_VISITED; }

extern int myid, nid;
extern double onesolarmass;
extern double com2real, real2com, potentfact;
extern double pntmass;
extern double r1kineticfact, r2kineticfact;
extern float amax, a, rng, size, hubble;
extern int ng, nspace;
extern int nx, ny, nz;
extern float omep, omeplam, acoeff[8];
extern float epsilon;
extern float m_tidal[NUM_MASS], r_tidal[NUM_MASS];

#include "generated_protos.h"
