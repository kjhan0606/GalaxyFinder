# PGalF Data Structures Reference

This document describes all major data structures used in PGalF.

## Table of Contents

1. [Particle Types](#1-particle-types)
2. [Tree Structures](#2-tree-structures)
3. [Working Structures](#3-working-structures)
4. [Output Structures](#4-output-structures)
5. [Memory Layout](#5-memory-layout)

---

## 1. Particle Types

All particle types are defined in `ramses.h`.

### 1.1 Base Particle Type (PmType / DmType / StarType)

Used for dark matter and star particles:

```c
typedef struct PmType {
    dptype x, y, z;          // Position [cMpc/h]
    dptype vx, vy, vz;       // Velocity [km/s physical]
    dptype mass;             // Mass [Msun/h]
    idtype id;               // Unique particle ID
    int levelp;              // AMR refinement level
    familytype family, tag;  // Particle type flags

#ifdef OUTPUT_PARTICLE_POTENTIAL
    dptype potent;           // Gravitational potential
#endif

#ifndef NBODY
    dptype tp, zp;           // Birth time, metallicity
#ifdef NCHEM
    dptype chem[NCHEM];      // Chemical abundances
#endif
    dptype mass0;            // Initial mass [Msun/h]
#ifdef FBK
    dptype fbk;              // Feedback energy
#endif
#endif
} PmType;

typedef PmType DmType;       // Dark matter alias
typedef PmType StarType;     // Star particle alias
```

**Size**: 72-128 bytes depending on compile flags

**Family codes**:
| Code | Type |
|------|------|
| 0 | Dark Matter |
| 1 | Star |
| 2 | Sink/AGN |
| 4 | Gas |

### 1.2 Gas Cell Type (GasType)

```c
typedef struct GasType {
    dptype x, y, z, cellsize;  // Position and cell size [cMpc/h]
    float vx, vy, vz;          // Velocity [km/s physical]
    dptype den;                // Gas density [simulation units]
    float temp;                // Temperature [K/μ/density]
    float metallicity;         // Metallicity Z
#ifdef NCHEM
    float chem[NCHEM];         // Chemical abundances
#endif
#ifdef NION
    float ion[NION];           // Ionization fractions
#endif
    float mass;                // Cell mass [Msun/h]
    dptype potent;             // Gravitational potential
    dptype fx, fy, fz;         // Gravitational force
} GasType;
```

**Size**: 88-120 bytes depending on compile flags

### 1.3 Sink/AGN Type (SinkType)

```c
typedef struct SinkType {
    dptype x, y, z;            // Position [cMpc/h]
    dptype vx, vy, vz;         // Velocity [km/s]
    dptype mass;               // Mass [Msun/h]
    dptype tbirth;             // Birth time
    dptype Jx, Jy, Jz;         // Angular momentum [Msun/h km/s kpc]
    dptype Sx, Sy, Sz;         // Spin direction (unit vector)
    dptype dMsmbh;             // SMBH accretion rate
    dptype dMBH_coarse;        // Coarse accretion rate
    dptype dMEd_coarse;        // Eddington rate
    dptype Esave;              // Stored energy
    dptype Smag;               // Spin magnitude
    dptype eps;                // Softening length
    int id;                    // Unique ID
} SinkType;
```

**Size**: 136 bytes

---

## 2. Tree Structures

### 2.1 FoF Tree Node (FoFTStruct)

Used for Friends-of-Friends grouping:

```c
typedef struct FoFTStruct {
    unsigned int type;         // TYPE_TREE (0) for nodes
    void *sibling;             // Next node at same level
    POSTYPE L;                 // Node side length
    POSTYPE dist2;             // Distance² to query point
    POSTYPE dist;              // Distance to query point
    POSTYPE x0, y0, z0;        // Node lower corner
    POSTYPE monox, monoy, monoz;  // Center of mass
    int Nparticle;             // Number of particles in node
    void *daughter;            // First child (node or particle)
    float maxlink02, minlink02;  // Linking length range
} FoFTStruct;
```

**Size**: 72-88 bytes

### 2.2 FoF Tree Particle (FoFTPtlStruct)

Unified particle structure for tree operations:

```c
typedef struct FoFTPtlStruct {
    unsigned int type;         // TYPE_STAR=1, TYPE_DM=2, TYPE_GAS=3, TYPE_SINK=4
    void *sibling;             // Next particle in linked list
    enum boolean included;     // YES/NO/ENCLOSE for FoF status
    dptype x, y, z;            // Position
    dptype vx, vy, vz;         // Velocity
    dptype mass;               // Mass
    float link02;              // FoF linking length squared
    union {
        DmType dm;
        StarType star;
        SinkType sink;
        GasType gas;
    } p;                       // Original particle data
    size_t haloindx;           // Parent halo index
    int indx;                  // Particle index
} FoFTPtlStruct;
```

**Size**: 216-280 bytes (largest of union members)

### 2.3 k-d Tree Node (for Nearest Neighbor)

```c
typedef struct TStruct {
    GENERALTPtlPOINTER;        // type, sibling
    void *daughter;            // Child pointer
    dptype L;                  // Node size
    dptype dist, dist2;        // Distance metrics
    dptype x0, y0, z0;         // Node position
    dptype mass;               // Total mass
    int Nparticle;             // Particle count
    dptype monox, monoy, monoz;  // Center of mass
    dptype trQ;                // Trace of quadrupole
    dptype quad[6];            // Quadrupole tensor
    dptype dist_over_thetasq;  // Opening angle criterion
} TStruct;
```

---

## 3. Working Structures

### 3.1 Working Particle (WorkingParticle)

Per-particle flags and density during processing:

```c
typedef struct WorkingParticle {
    float den;                 // Computed density
    int haloid;                // Assigned subhalo ID (-1 = unassigned)
    unsigned char flag;        // Status flags
    unsigned char flag2[MAXTHREADS];  // Thread-local flags
} WorkingParticle;
```

**Size**: 8 + MAXTHREADS bytes

**Flag bit definitions**:
```c
#define WP_RESET     0x00  // No flags set
#define WP_PEAK      0x01  // Is a density peak
#define WP_VISITED   0x02  // Visited during water-shedding
#define WP_CORE      0x04  // Belongs to core
#define WP_SHELL     0x08  // Belongs to shell
#define WP_BOUND     0x10  // Gravitationally bound
#define WP_REMAINING 0x20  // Remaining after processing
```

**Flag manipulation macros**:
```c
SET_PEAK(i)       // Set peak flag
UNSET_PEAK(i)     // Clear peak flag
IS_PEAK(i)        // Test peak flag
TOGGLE_PEAK(i)    // Toggle peak flag
```

### 3.2 Core Type (Coretype)

Represents an identified density peak/galaxy core:

```c
typedef struct Coretype {
    int peak;              // Index of peak particle
    int nummem;            // Number of member particles
    int numstar;           // Number of star particles
    float starmass;        // Total stellar mass
    float coredensity;     // Density threshold for core
    float cx, cy, cz;      // Center of mass position
    float cvx, cvy, cvz;   // Center of mass velocity
    float Rtidal;          // Tidal radius
    float density;         // Peak density
    unsigned char flag;    // Status flags
} Coretype;
```

**Size**: 52 bytes

**Core flag definitions**:
```c
#define CORE_RESET     0x00  // No flags
#define CORE_ENCLOSED  0x01  // Enclosed within another core
#define CORE_CONFIRMED 0x02  // Confirmed as real core
#define REAL_CORE      0x04  // Definitely a real core
```

### 3.3 Core Sort Type (Coresorttype)

For sorting cores by membership:

```c
typedef struct Coresorttype {
    int nmem;              // Number of members
    Coretype *core;        // Pointer to core
} Coresorttype;
```

### 3.4 Simple Particle Type (SimpleBasicParticleType)

Lightweight particle structure for tree operations:

```c
typedef struct SimpleBasicParticleType {
    unsigned int type;
    dptype mass, x, y, z;
    dptype vx, vy, vz;
    dptype link02;
    idtype indx;
    struct SimpleBasicParticleType *bp;  // Linked list pointer
} SimpleBasicParticleType;
```

---

## 4. Output Structures

### 4.1 Halo Information (HaloInfo)

Summary of a FoF halo with identified substructures:

```c
typedef struct HaloInfo {
    int nsub;              // Number of subhalos
    int ndm, nstar, nsink, ngas, npall;  // Particle counts
    dptype totm;           // Total mass [Msun/h]
    dptype mdm, mgas, msink, mstar;      // Component masses
    dptype x, y, z;        // Center of mass [cMpc/h]
    dptype vx, vy, vz;     // Bulk velocity [km/s]
} HaloInfo;
```

**Size**: 96 bytes

### 4.2 Subhalo Information (SubInfo)

Properties of an individual subhalo:

```c
typedef struct SubInfo {
    int npdm, npgas, npsink, npstar, npall;  // Particle counts
    dptype totm;           // Total mass [Msun/h]
    dptype mdm, mgas, msink, mstar;          // Component masses
    dptype x, y, z;        // Center of mass [cMpc/h]
    dptype vx, vy, vz;     // Bulk velocity [km/s]
} SubInfo;
```

**Size**: 88 bytes

### 4.3 Halo Catalog Entry (HaloQ)

Input catalog entry for FoF halos:

```c
typedef struct HaloQ {
    size_t np, npstar, npgas, npdm, npsink;  // Particle counts
    POSTYPE x, y, z;       // Halo center [cMpc/h]
    float mass, mstar, mgas, mdm, msink;     // Masses
    float vx, vy, vz;      // Bulk velocity [km/s]
} HaloQ;
```

**Size**: 76-92 bytes

---

## 5. Memory Layout

### 5.1 Custom Memory Allocator

The code uses a stack-based memory pool to avoid fragmentation:

```c
// Memory.h
typedef long long INT8;

INT8 Make_Total_Memory(INT8 size);    // Initialize pool
void *Malloc(INT8 size, void **ptr);  // Allocate from stack
void *Realloc(void *ptr, INT8 size);  // Resize allocation
void Free(void *ptr);                  // Release to stack
INT8 CheckAvailableMemory(void);      // Query free space
long CurMemStack(void);               // Current stack position
void InitialOldMemStack(long pos);    // Reset stack to position
```

### 5.2 Memory Layout Example

For processing a halo with N particles:

```
+------------------+ ← Stack base (Make_Total_Memory)
|  Global arrays   |
+------------------+
|  FoFTPtlStruct   | N × 280 bytes
|  (bp array)      |
+------------------+
|  WorkingParticle | N × (8 + MAXTHREADS) bytes
|  (wp array)      |
+------------------+
|  Coretype array  | MAXNUMCORE × 52 bytes
+------------------+
|  k-d tree nodes  | ~2N × 80 bytes
+------------------+
|  FoF tree nodes  | ~2N × 80 bytes
+------------------+
|  TSC density     | nx × ny × nz × 4 bytes
|  grid            |
+------------------+
|  FFTW buffers    | 2 × nx × ny × nz × 4 bytes
+------------------+
|  Temporary       | Variable
|  arrays          |
+------------------+ ← Stack top
```

### 5.3 Memory Estimation

For a halo with N particles and grid size G³:

| Component | Memory (bytes) |
|-----------|----------------|
| Particle data | N × 280 |
| Working flags | N × 72 |
| Core array | 52 MB (1M cores) |
| k-d tree | N × 160 |
| FoF tree | N × 160 |
| TSC grid | G³ × 4 |
| FFT buffers | G³ × 8 |
| **Total** | **~680N + 12G³ + 52M** |

Example: 10M particles, 512³ grid:
- Particles: 6.8 GB
- Grid: 1.6 GB
- Trees: 3.2 GB
- **Total: ~12 GB**

---

## Type Definitions Summary

| Type | Defined In | Description |
|------|------------|-------------|
| `dptype` | ramses.h | double (position/mass) |
| `idtype` | ramses.h | int or long long (particle ID) |
| `familytype` | ramses.h | char (particle family) |
| `POSTYPE` | tree.h | double or float (tree positions) |
| `lint` | header.h | int (local index) |
| `indxtype` | defs.h | long (global index) |
| `INT8` | Memory.h | long long (memory sizes) |

---

## Preprocessor Configurations

| Flag | Effect |
|------|--------|
| `XYZDBL` | Use double for positions |
| `DOUBLE_PRECISION` | Use long long for IDs |
| `READ_SINK` | Enable sink particles |
| `NCHEM=N` | Track N chemical elements |
| `NENER=N` | Track N radiation groups |
| `NBODY` | Pure N-body (no stars) |
| `OUTPUT_PARTICLE_POTENTIAL` | Store potentials |
