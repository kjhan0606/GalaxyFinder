# opFoF Data Structures Reference

This document provides a comprehensive reference for all data structures used in opFoF, including their memory layouts, field descriptions, and usage patterns.

---

## Table of Contents

1. [Particle Structures](#particle-structures)
2. [Halo Structures](#halo-structures)
3. [Tree Structures](#tree-structures)
4. [Simulation Parameters](#simulation-parameters)
5. [Memory Management](#memory-management)
6. [MPI Communication Structures](#mpi-communication-structures)
7. [Type Definitions](#type-definitions)
8. [Memory Layout Diagrams](#memory-layout-diagrams)

---

## Particle Structures

### FoFTPtlStruct - Main Particle Structure

The primary structure representing a particle in the FoF algorithm.

```c
typedef struct FoFTPtlStruct {
    unsigned int type;           // Particle type identifier
    POSTYPE x, y, z;             // Position coordinates
    POSTYPE link02;              // Linking length for this particle
    union {
        DmType dm;               // Dark matter properties
        StarType star;           // Stellar properties
        GasType gas;             // Gas properties
        SinkType sink;           // Sink (black hole) properties
    };
    size_t haloindx;             // Assigned halo index (0 = unlinked)
} FoFTPtlStruct;
```

**Field Details:**

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `type` | unsigned int | 4 bytes | Particle type: 1=STAR, 2=DM, 3=GAS, 4=SINK |
| `x, y, z` | POSTYPE | 8 bytes each (if XYZDBL) | Position in Mpc/h |
| `link02` | POSTYPE | 8 bytes | Particle-specific linking length |
| `union` | varies | ~40 bytes | Type-specific properties |
| `haloindx` | size_t | 8 bytes | Halo membership index |

**Total Size:** ~80-96 bytes per particle (depending on alignment)

### DmType - Dark Matter Particle

```c
typedef struct DmType {
    float vx, vy, vz;            // Velocity components (km/s)
    float mass;                  // Particle mass (Msun/h)
    long long id;                // Unique particle identifier
    int level;                   // AMR refinement level
} DmType;
```

**Field Details:**

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `vx, vy, vz` | float | 4 bytes each | Velocity in km/s |
| `mass` | float | 4 bytes | Mass in Msun/h |
| `id` | long long | 8 bytes | Unique ID (supports >4 billion particles) |
| `level` | int | 4 bytes | AMR level (0 = coarsest) |

### StarType - Stellar Particle

```c
typedef struct StarType {
    float vx, vy, vz;            // Velocity components (km/s)
    float mass;                  // Stellar mass (Msun/h)
    long long id;                // Unique particle identifier
    float age;                   // Stellar age (Gyr)
    float metallicity;           // Metal mass fraction (Z)
} StarType;
```

**Additional Fields:**

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `age` | float | 4 bytes | Time since formation |
| `metallicity` | float | 4 bytes | Z = M_metals / M_total |

### GasType - Gas Particle

```c
typedef struct GasType {
    float vx, vy, vz;            // Velocity components (km/s)
    float mass;                  // Gas mass (Msun/h)
    float density;               // Gas density (Msun/h / (Mpc/h)³)
    float temperature;           // Gas temperature (K)
    float metallicity;           // Metal mass fraction
    float sfr;                   // Star formation rate (Msun/yr)
} GasType;
```

**Additional Fields:**

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `density` | float | 4 bytes | Local gas density |
| `temperature` | float | 4 bytes | Gas temperature |
| `sfr` | float | 4 bytes | Star formation rate |

### SinkType - Sink (Black Hole) Particle

```c
typedef struct SinkType {
    float vx, vy, vz;            // Velocity components (km/s)
    float mass;                  // Sink mass (Msun/h)
    float birth_time;            // Formation time (scale factor)
    float lx, ly, lz;            // Angular momentum vector
    float mdot;                  // Accretion rate (Msun/yr)
} SinkType;
```

**Additional Fields:**

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `birth_time` | float | 4 bytes | Scale factor at formation |
| `lx, ly, lz` | float | 4 bytes each | Specific angular momentum |
| `mdot` | float | 4 bytes | Mass accretion rate |

---

## Halo Structures

### HaloQ - Halo Properties Structure

The main output structure containing halo properties.

```c
typedef struct HaloQ {
    size_t np;                   // Total particle count
    size_t npstar;               // Star particle count
    size_t npgas;                // Gas particle count
    size_t npdm;                 // Dark matter particle count
    size_t npsink;               // Sink particle count
    POSTYPE x, y, z;             // Center of mass position (Mpc/h)
    double mass;                 // Total halo mass (Msun/h)
    double mstar;                // Total stellar mass
    double mgas;                 // Total gas mass
    double mdm;                  // Total dark matter mass
    double msink;                // Total sink mass
    float vx, vy, vz;            // Center of mass velocity (km/s)
} HaloQ;
```

**Field Details:**

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `np` | size_t | 8 bytes | Total number of member particles |
| `npstar` | size_t | 8 bytes | Number of star particles |
| `npgas` | size_t | 8 bytes | Number of gas particles |
| `npdm` | size_t | 8 bytes | Number of DM particles |
| `npsink` | size_t | 8 bytes | Number of sink particles |
| `x, y, z` | POSTYPE | 8 bytes each | COM position |
| `mass` | double | 8 bytes | Total mass |
| `mstar` | double | 8 bytes | Stellar mass |
| `mgas` | double | 8 bytes | Gas mass |
| `mdm` | double | 8 bytes | DM mass |
| `msink` | double | 8 bytes | Sink mass |
| `vx, vy, vz` | float | 4 bytes each | COM velocity |

**Total Size:** 128 bytes per halo

### HaloBound - Boundary Halo Tracking

Used to track halos that intersect domain boundaries.

```c
typedef struct HaloBound {
    size_t nmem;                 // Number of member particles
    POSTYPE zmin;                // Minimum z coordinate
    POSTYPE zmax;                // Maximum z coordinate
    int boundflag;               // Boundary contact flag
} HaloBound;
```

**boundflag Values:**

| Value | Meaning |
|-------|---------|
| 0 | Internal halo (no boundary contact) |
| 1 | Touches lower Z boundary |
| 2 | Touches upper Z boundary |
| 3 | Touches both boundaries (spans domain) |

---

## Tree Structures

### FoFTStruct - Tree Node Structure

Represents a node in the oct-tree spatial index.

```c
typedef struct FoFTStruct {
    POSTYPE dist;                // Maximum extent from center
    POSTYPE monox;               // Center x coordinate
    POSTYPE monoy;               // Center y coordinate
    POSTYPE monoz;               // Center z coordinate
    void *daughter;              // Pointer to first child/particle list
    void *sibling;               // Pointer to next sibling node
    int Nparticle;               // Number of particles in subtree
} FoFTStruct;
```

**Field Details:**

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `dist` | POSTYPE | 8 bytes | Max distance from center to any particle |
| `monox,y,z` | POSTYPE | 8 bytes each | Geometric center |
| `daughter` | void* | 8 bytes | Child pointer or particle list |
| `sibling` | void* | 8 bytes | Next sibling in tree |
| `Nparticle` | int | 4 bytes | Particle count in subtree |

**Total Size:** 52 bytes per node (with padding: 56 bytes)

### Tree Node Types

```
Internal Node (Nparticle > NODE_HAVE_PARTICLE):
    daughter -> first child node
    sibling -> next sibling node

Leaf Node (Nparticle <= NODE_HAVE_PARTICLE):
    daughter -> linked list of particles
    sibling -> next sibling node
```

### Particle Linking in Leaves

In leaf nodes, particles are linked via sibling pointers:

```c
// Traversing particles in a leaf node
FoFTPtlStruct *p = (FoFTPtlStruct *)leaf->daughter;
while (p != NULL) {
    // Process particle p
    p = (FoFTPtlStruct *)p->sibling_ptr;  // pseudo-code
}
```

---

## Simulation Parameters

### PMHeader - Simulation Header

Contains global simulation parameters read from RAMSES info files.

```c
typedef struct PMHeader {
    double omega_m;              // Matter density parameter
    double omega_l;              // Dark energy density parameter
    double omega_b;              // Baryon density parameter
    double H0;                   // Hubble constant (km/s/Mpc)
    double boxlen_ini;           // Initial box size (Mpc/h)
    double aexp;                 // Current scale factor
    double amax;                 // Maximum scale factor
    double unit_l;               // Length unit conversion
    double unit_d;               // Density unit conversion
    double unit_t;               // Time unit conversion
    int ncpu;                    // Number of RAMSES CPUs used
    int nstep;                   // Snapshot number
} PMHeader;
```

### Cosmological Parameters

```c
// Critical density at z=0
#define RHO_CRIT 2.7755e11  // Msun/h / (Mpc/h)³

// Mean matter density
rho_m = RHO_CRIT * omega_m;

// Mean inter-particle separation
d_mean = pow(mass_particle / rho_m, 1.0/3.0);

// Linking length (0.2 of mean separation)
link_length = 0.2 * d_mean;
```

---

## Memory Management

### Memory Pool Structure

opFoF uses a custom memory allocator for efficiency:

```c
// Global memory pool (Memory.c)
static char *MemoryPool;         // Base pointer
static size_t PoolSize;          // Total pool size
static size_t CurrentOffset;     // Current allocation position

// Allocation stack
static void *SMalloc[MAXSTACK];  // Allocation pointers
static size_t SizeMalloc[MAXSTACK]; // Allocation sizes
static int StackTop;             // Stack pointer
```

### Memory Allocation Pattern

```
Pool Layout:
+------------------+------------------+------------------+----
| Allocation 1     | Allocation 2     | Allocation 3     | ...
+------------------+------------------+------------------+----
^                  ^                  ^
SMalloc[0]         SMalloc[1]         SMalloc[2]
```

### Memory Functions

```c
// Initialize memory pool
void InitMemory(size_t size_mb);

// Allocate from pool
void *MyMalloc(size_t size);

// Free (stack-based, must free in reverse order)
void MyFree(void *ptr);

// Get current usage
size_t GetMemoryUsage(void);
```

---

## MPI Communication Structures

### Boundary Particle Exchange

```c
// Structure for boundary particle data
typedef struct BoundaryParticle {
    FoFTPtlStruct particle;      // Full particle data
    size_t original_halo;        // Halo index on source rank
} BoundaryParticle;

// Buffer for MPI communication
typedef struct MPIBuffer {
    void *data;                  // Data buffer
    size_t count;                // Number of elements
    MPI_Datatype type;           // MPI data type
} MPIBuffer;
```

### Communication Pattern

```
Rank Layout (Z-slab decomposition):

Rank 0: [zmin=0, zmax=L/P)
        ↓ send lower boundary
Rank 1: [zmin=L/P, zmax=2L/P)
        ↓ send lower boundary  ↑ receive upper boundary
Rank 2: [zmin=2L/P, zmax=3L/P)
        ...
```

---

## Type Definitions

### Position Type

```c
#ifdef XYZDBL
    typedef double POSTYPE;      // Double precision positions
#else
    typedef float POSTYPE;       // Single precision positions
#endif
```

### Particle Type Constants

```c
#define TYPE_STAR   1
#define TYPE_DM     2
#define TYPE_GAS    3
#define TYPE_SINK   4

// Type checking macros
#define IS_DM(p)    ((p)->type == TYPE_DM)
#define IS_STAR(p)  ((p)->type == TYPE_STAR)
#define IS_GAS(p)   ((p)->type == TYPE_GAS)
#define IS_SINK(p)  ((p)->type == TYPE_SINK)
```

### Index Types

```c
typedef size_t halo_index_t;     // Halo index type (64-bit)
typedef long long particle_id_t; // Particle ID type (64-bit)

#define UNLINKED 0               // Unassigned halo index
```

### Status Flags

```c
// Halo boundary flags
#define BOUND_NONE   0           // No boundary contact
#define BOUND_LOWER  1           // Lower Z boundary
#define BOUND_UPPER  2           // Upper Z boundary
#define BOUND_BOTH   3           // Both boundaries

// Processing status
#define STATUS_PENDING   0
#define STATUS_PROCESSED 1
#define STATUS_OUTPUT    2
```

---

## Memory Layout Diagrams

### Particle Array Layout

```
Memory Layout for N particles:

Address:  0x0000    0x0060    0x00C0    0x0120    ...
          +--------+--------+--------+--------+
          |Ptl[0]  |Ptl[1]  |Ptl[2]  |Ptl[3]  | ...
          +--------+--------+--------+--------+
          96 bytes each (with alignment)

Single Particle (FoFTPtlStruct):
Offset    Field           Size
0x00      type            4 bytes
0x04      (padding)       4 bytes
0x08      x               8 bytes
0x10      y               8 bytes
0x18      z               8 bytes
0x20      link02          8 bytes
0x28      union (dm/star/..) ~40 bytes
0x50      haloindx        8 bytes
0x58      (padding)       8 bytes
0x60      --- next particle ---
```

### Tree Node Layout

```
Oct-tree Structure:

                    [Root Node]
                    Nparticle=1000
                   /    |    \    \
            [Oct0]  [Oct1] ... [Oct7]
            Np=100  Np=200     Np=50
           /   \
      [Sub0] [Sub1]  ... (continue until Np <= 8)
      Np=8   Np=12

Leaf Node:
+------------------+
| FoFTStruct       |
|   dist, mono...  |
|   daughter ------+---> [Particle 1] -> [Particle 2] -> ... -> NULL
|   sibling        |
|   Nparticle=8    |
+------------------+
```

### Halo Catalog File Layout

```
File: FoF_halo_cat.NNNNN

Offset      Content                     Size
0x0000      Header (7 floats)           28 bytes
0x001C      HaloQ[0]                    128 bytes
0x009C      HaloQ[1]                    128 bytes
0x011C      HaloQ[2]                    128 bytes
...
EOF

Total size = 28 + 128 * num_halos bytes
```

### Member Particle File Layout

```
File: FoF_member_particle.NNNNN

Offset      Content                     Size
0x0000      Header (7 floats)           28 bytes
0x001C      Halo 0 particles
            - DM particles (npdm)
            - Gas particles (npgas)
            - Sink particles (npsink)
            - Star particles (npstar)
...         Halo 1 particles
...

Particle order within each halo:
[DM][DM]...[GAS][GAS]...[SINK][SINK]...[STAR][STAR]...
```

---

## Best Practices

### Memory Efficiency

1. **Use appropriate types**: Use `float` for velocities (sufficient precision), `double` for positions and masses

2. **Minimize padding**: Arrange structure fields by size (largest first)

3. **Pool allocation**: Use the memory pool for large allocations to avoid fragmentation

### Type Safety

```c
// Always check particle type before accessing union
void process_particle(FoFTPtlStruct *p) {
    switch(p->type) {
        case TYPE_DM:
            process_dm(p->dm);
            break;
        case TYPE_STAR:
            process_star(p->star);
            break;
        // ...
    }
}
```

### Endianness Considerations

```c
// opFoF uses native byte order
// When reading on different architecture, byte-swap may be needed

#ifdef SWAP_ENDIAN
    void swap_haloq(HaloQ *h) {
        h->np = bswap_64(h->np);
        // ... swap all fields
    }
#endif
```

---

## Version Compatibility

| Version | Changes |
|---------|---------|
| 1.0 | Initial structure definitions |
| 1.1 | Added sink particle support |
| 1.2 | Changed haloindx to size_t (64-bit) |
| 1.3 | Added metallicity to GasType |
