# opFoF API Reference

This document provides a complete reference for all public functions, macros, and interfaces in opFoF.

---

## Table of Contents

1. [Core FoF Functions](#core-fof-functions)
2. [Tree Functions](#tree-functions)
3. [I/O Functions](#io-functions)
4. [Memory Management](#memory-management)
5. [MPI Communication](#mpi-communication)
6. [Utility Functions](#utility-functions)
7. [Timing Functions](#timing-functions)
8. [Macros and Constants](#macros-and-constants)

---

## Core FoF Functions

### FoF_Make_Tree

Constructs the oct-tree spatial index from particle data.

```c
void FoF_Make_Tree(FoFTPtlStruct *particles, int nparticles);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `particles` | `FoFTPtlStruct*` | Array of particle data |
| `nparticles` | `int` | Number of particles |

**Returns:** None (modifies global tree structure)

**Description:**
Builds a hierarchical oct-tree for efficient spatial queries. The tree is used by the FoF linking functions to find neighbors within the linking length.

**Example:**
```c
FoFTPtlStruct *ptl = load_particles("snapshot.dat", &nptl);
FoF_Make_Tree(ptl, nptl);
```

**Notes:**
- Must be called before `new_fof_link()` or `pnew_fof_link()`
- Tree is stored in global variables
- Call `Free_Tree()` when done to release memory

---

### new_fof_link

Links particles using FoF algorithm (non-periodic boundaries).

```c
int new_fof_link(FoFTPtlStruct *particles, int nparticles,
                 POSTYPE link_length);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `particles` | `FoFTPtlStruct*` | Array of particles with tree built |
| `nparticles` | `int` | Number of particles |
| `link_length` | `POSTYPE` | FoF linking length in Mpc/h |

**Returns:** `int` - Number of halos identified

**Description:**
Performs Friends-of-Friends linking on particles without periodic boundary conditions. Each particle's `haloindx` field is set to its assigned halo ID (1-indexed, 0 means unlinked).

**Example:**
```c
FoF_Make_Tree(particles, nptl);
int nhalos = new_fof_link(particles, nptl, 0.2 * mean_sep);
printf("Found %d halos\n", nhalos);
```

---

### pnew_fof_link

Links particles using FoF algorithm (periodic boundaries).

```c
int pnew_fof_link(FoFTPtlStruct *particles, int nparticles,
                  POSTYPE link_length, POSTYPE BoxSize);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `particles` | `FoFTPtlStruct*` | Array of particles |
| `nparticles` | `int` | Number of particles |
| `link_length` | `POSTYPE` | FoF linking length |
| `BoxSize` | `POSTYPE` | Simulation box size |

**Returns:** `int` - Number of halos identified

**Description:**
Same as `new_fof_link()` but handles periodic boundary conditions. Particles near box edges are checked for links across the periodic boundary.

**Example:**
```c
int nhalos = pnew_fof_link(particles, nptl, link_len, 100.0);
```

---

### haloproperty

Calculates properties for a single halo.

```c
void haloproperty(FoFTPtlStruct *linked, size_t nlinked,
                  HaloQ *halo, POSTYPE BoxSize, int pflag);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `linked` | `FoFTPtlStruct*` | Array of particles in halo |
| `nlinked` | `size_t` | Number of particles |
| `halo` | `HaloQ*` | Output structure for halo properties |
| `BoxSize` | `POSTYPE` | Box size for periodic handling |
| `pflag` | `int` | Periodic flag (0=non-periodic, 1=periodic) |

**Returns:** None (fills `halo` structure)

**Description:**
Computes center-of-mass position, velocity, total mass, and component masses for a halo. Handles periodic boundary wrapping if `pflag=1`.

**Example:**
```c
HaloQ halo;
haloproperty(linked_particles, nlinked, &halo, 100.0, 1);
printf("Halo mass: %.3e Msun/h\n", halo.mass);
```

---

### WriteIsolatedHalo

Writes a complete halo (not on domain boundary) to output.

```c
void WriteIsolatedHalo(FoFTPtlStruct *linked, size_t nlinked,
                       FILE *fcat, FILE *fmem, int pflag, POSTYPE BoxSize);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `linked` | `FoFTPtlStruct*` | Halo particles |
| `nlinked` | `size_t` | Particle count |
| `fcat` | `FILE*` | Catalog file pointer |
| `fmem` | `FILE*` | Member particle file pointer |
| `pflag` | `int` | Periodic flag |
| `BoxSize` | `POSTYPE` | Box size |

**Returns:** None

---

### WriteFinalHalo

Writes final merged boundary halo to output.

```c
void WriteFinalHalo(FoFTPtlStruct *linked, size_t nlinked,
                    FILE *fcat, FILE *fmem);
```

---

### CheckHaloBound

Checks if halo intersects domain boundaries.

```c
int CheckHaloBound(FoFTPtlStruct *linked, size_t nlinked,
                   POSTYPE zmin, POSTYPE zmax, POSTYPE link_len);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `linked` | `FoFTPtlStruct*` | Halo particles |
| `nlinked` | `size_t` | Particle count |
| `zmin` | `POSTYPE` | Domain lower Z boundary |
| `zmax` | `POSTYPE` | Domain upper Z boundary |
| `link_len` | `POSTYPE` | Linking length |

**Returns:** `int` - Boundary flag (0=internal, 1=lower, 2=upper, 3=both)

---

## Tree Functions

### TreeWalk

Walks the tree to find neighbors within radius.

```c
void TreeWalk(FoFTStruct *node, FoFTPtlStruct *particle,
              POSTYPE radius, FoFTPtlStruct **neighbors, int *nneighbors);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `node` | `FoFTStruct*` | Tree node to start from |
| `particle` | `FoFTPtlStruct*` | Query particle |
| `radius` | `POSTYPE` | Search radius |
| `neighbors` | `FoFTPtlStruct**` | Output neighbor array |
| `nneighbors` | `int*` | Output neighbor count |

---

### Free_Tree

Releases tree memory.

```c
void Free_Tree(void);
```

---

## I/O Functions

### ReadRamsesParticles

Reads particle data from RAMSES format files.

```c
int ReadRamsesParticles(int nstep, int file_idx, FoFTPtlStruct **particles,
                        int *nparticles, PMHeader *header);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `nstep` | `int` | Snapshot number |
| `file_idx` | `int` | File index within snapshot |
| `particles` | `FoFTPtlStruct**` | Output particle array |
| `nparticles` | `int*` | Output particle count |
| `header` | `PMHeader*` | Output header info |

**Returns:** `int` - 0 on success, -1 on error

**Example:**
```c
FoFTPtlStruct *ptl;
int nptl;
PMHeader header;

if (ReadRamsesParticles(100, 0, &ptl, &nptl, &header) < 0) {
    fprintf(stderr, "Failed to read particles\n");
    exit(1);
}
```

---

### ReadHeader

Reads simulation header/info file.

```c
int ReadHeader(int nstep, int file_idx, PMHeader *header);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `nstep` | `int` | Snapshot number |
| `file_idx` | `int` | File index |
| `header` | `PMHeader*` | Output header structure |

**Returns:** `int` - 0 on success, -1 on error

---

### WriteHaloCatalog

Writes halo catalog header.

```c
void WriteHaloCatalog(FILE *fp, PMHeader *header);
```

---

### WriteMemberHeader

Writes member particle file header.

```c
void WriteMemberHeader(FILE *fp, PMHeader *header);
```

---

## Memory Management

### InitMemory

Initializes the memory pool.

```c
void InitMemory(size_t size_mb);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `size_mb` | `size_t` | Pool size in megabytes |

**Example:**
```c
InitMemory(17000);  // 17 GB pool
```

---

### MyMalloc

Allocates memory from pool.

```c
void *MyMalloc(size_t size);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `size` | `size_t` | Number of bytes to allocate |

**Returns:** `void*` - Pointer to allocated memory, NULL on failure

**Example:**
```c
FoFTPtlStruct *particles = MyMalloc(nptl * sizeof(FoFTPtlStruct));
if (!particles) {
    fprintf(stderr, "Allocation failed\n");
    exit(1);
}
```

---

### MyFree

Frees memory back to pool.

```c
void MyFree(void *ptr);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `ptr` | `void*` | Pointer to free |

**Notes:**
- Must free in reverse order of allocation (stack-based)
- Freeing out of order will cause errors

---

### GetMemoryUsage

Returns current memory pool usage.

```c
size_t GetMemoryUsage(void);
```

**Returns:** `size_t` - Bytes currently allocated

---

### FreeMemory

Releases entire memory pool.

```c
void FreeMemory(void);
```

---

## MPI Communication

### BIG_MPI_Send

Sends large messages in chunks.

```c
int BIG_MPI_Send(void *buf, size_t count, MPI_Datatype datatype,
                 int dest, int tag, MPI_Comm comm);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `buf` | `void*` | Send buffer |
| `count` | `size_t` | Number of elements |
| `datatype` | `MPI_Datatype` | Element type |
| `dest` | `int` | Destination rank |
| `tag` | `int` | Message tag |
| `comm` | `MPI_Comm` | Communicator |

**Returns:** `int` - MPI error code

**Description:**
Standard MPI_Send may fail for very large messages. This function breaks the message into chunks of ~1GB and sends them sequentially.

---

### BIG_MPI_Recv

Receives large messages in chunks.

```c
int BIG_MPI_Recv(void *buf, size_t count, MPI_Datatype datatype,
                 int source, int tag, MPI_Comm comm, MPI_Status *status);
```

---

### WriteBottomFaceContact

Sends boundary particles to neighboring rank.

```c
void WriteBottomFaceContact(FoFTPtlStruct *particles, size_t nparticles,
                            size_t *halo_ids, int dest_rank);
```

---

### ReadBottomFaceContact

Receives boundary particles from neighboring rank.

```c
int ReadBottomFaceContact(FoFTPtlStruct **particles, size_t *nparticles,
                          size_t **halo_ids, int source_rank);
```

---

### StackUpContactParticleLeftWard

Accumulates boundary particles for leftward transfer.

```c
void StackUpContactParticleLeftWard(FoFTPtlStruct *particles, size_t nparticles,
                                    size_t halo_id);
```

---

## Utility Functions

### ComputeLinkLength

Calculates FoF linking length from particle mass.

```c
POSTYPE ComputeLinkLength(double mass, double omega_m, double b_factor);
```

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `mass` | `double` | Particle mass in Msun/h |
| `omega_m` | `double` | Matter density parameter |
| `b_factor` | `double` | Linking factor (typically 0.2) |

**Returns:** `POSTYPE` - Linking length in Mpc/h

**Example:**
```c
POSTYPE link_len = ComputeLinkLength(1e10, 0.3089, 0.2);
```

---

### PeriodicDistance

Calculates minimum image distance.

```c
POSTYPE PeriodicDistance(POSTYPE x1, POSTYPE y1, POSTYPE z1,
                         POSTYPE x2, POSTYPE y2, POSTYPE z2,
                         POSTYPE BoxSize);
```

**Returns:** `POSTYPE` - Distance with periodic wrapping

---

### SortParticlesByZ

Sorts particles by Z coordinate.

```c
void SortParticlesByZ(FoFTPtlStruct *particles, int nparticles);
```

---

### zhpsort / xhpsort / yhpsort

Heap sort by specific coordinate.

```c
void zhpsort(FoFTPtlStruct *particles, size_t *indices, int n);
void xhpsort(FoFTPtlStruct *particles, size_t *indices, int n);
void yhpsort(FoFTPtlStruct *particles, size_t *indices, int n);
```

---

## Timing Functions

### InitTimer

Initializes timing system.

```c
void InitTimer(void);
```

---

### StartTimer

Starts a named timer.

```c
void StartTimer(const char *name);
```

**Example:**
```c
StartTimer("TreeBuild");
FoF_Make_Tree(ptl, nptl);
StopTimer("TreeBuild");
```

---

### StopTimer

Stops a named timer.

```c
void StopTimer(const char *name);
```

---

### GetTimer

Gets elapsed time for a timer.

```c
double GetTimer(const char *name);
```

**Returns:** `double` - Elapsed time in seconds

---

### PrintTimers

Prints all timer values.

```c
void PrintTimers(void);
```

---

### WallTime

Returns current wall clock time.

```c
double WallTime(void);
```

**Returns:** `double` - Time in seconds since epoch

---

## Macros and Constants

### Type Identifiers

```c
#define TYPE_STAR   1    // Stellar particle
#define TYPE_DM     2    // Dark matter particle
#define TYPE_GAS    3    // Gas particle
#define TYPE_SINK   4    // Sink (black hole) particle
```

### Tree Constants

```c
#define NODE_HAVE_PARTICLE 8  // Max particles per leaf
```

### Physical Constants

```c
#define RHO_CRIT 2.7755e11    // Critical density (Msun/h / (Mpc/h)³)
```

### Boundary Flags

```c
#define BOUND_NONE   0    // Internal halo
#define BOUND_LOWER  1    // Touches lower boundary
#define BOUND_UPPER  2    // Touches upper boundary
#define BOUND_BOTH   3    // Touches both boundaries
```

### Memory Macros

```c
#define NMEG 17000L           // Default pool size (MB)
#define MAXSTACK 100          // Maximum allocation stack depth
```

### Utility Macros

```c
// Type checking
#define IS_DM(p)    ((p)->type == TYPE_DM)
#define IS_STAR(p)  ((p)->type == TYPE_STAR)
#define IS_GAS(p)   ((p)->type == TYPE_GAS)
#define IS_SINK(p)  ((p)->type == TYPE_SINK)

// Min/Max
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// Square
#define SQ(x) ((x)*(x))
```

---

## Error Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| -1 | Generic error |
| -2 | File not found |
| -3 | Memory allocation failed |
| -4 | MPI error |
| -5 | Invalid parameter |

---

## Thread Safety

opFoF is **not thread-safe** by default. The following are not thread-safe:

- Global tree structure
- Memory pool
- File I/O operations

For hybrid MPI+OpenMP parallelization, each thread must:
1. Have its own tree
2. Use thread-local memory
3. Synchronize file access

---

## Version History

| Version | API Changes |
|---------|-------------|
| 1.0 | Initial API |
| 1.1 | Added `BIG_MPI_Send/Recv` |
| 1.2 | Added sink particle support |
| 1.3 | Changed `haloindx` to `size_t` |
