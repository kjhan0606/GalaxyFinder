# PGalF - Parallel Galaxy Finder

**PGalF** (Parallel Galaxy Finder) is a high-performance MPI+OpenMP code for identifying galaxies and subhalos within Friends-of-Friends (FoF) halos from cosmological N-body/hydrodynamic simulations. It uses density peak finding with water-shedding algorithms and physical tidal radius calculations to robustly identify bound substructures.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Algorithm Description](#algorithm-description)
- [Installation](#installation)
- [Usage](#usage)
- [Input/Output Formats](#inputoutput-formats)
- [Configuration Parameters](#configuration-parameters)
- [Code Architecture](#code-architecture)
- [Performance Considerations](#performance-considerations)
- [File Descriptions](#file-descriptions)
- [References](#references)

---

## Overview

PGalF processes FoF halo catalogs and identifies substructures (galaxies/subhalos) within each halo using the following approach:

1. **Density Estimation**: Compute particle densities using TSC (Triangular-Shaped Cloud) interpolation with Gaussian smoothing
2. **Peak Finding**: Identify local density maxima as potential galaxy cores
3. **Water-shedding**: Group particles to their nearest density peaks following density gradients
4. **Boundedness Test**: Iteratively remove unbound particles using energy criteria
5. **Tidal Radius**: Calculate physical tidal radii using NFW profile models
6. **Final Assignment**: Assign remaining particles to subhalos via FoF linking

### Parallelization Strategy

- **MPI Master-Slave Pattern**: Rank 0 handles I/O and work distribution; worker ranks process individual halos
- **OpenMP Threading**: Parallel density calculations and tree operations within each halo
- Scales to thousands of MPI processes with 8-64 OpenMP threads per process

---

## Features

- Multi-component support: Dark matter, gas, stars, and sink/AGN particles
- Hybrid MPI+OpenMP parallelization for large-scale simulations
- Physical tidal radius calculation using NFW profiles
- Iterative boundedness determination
- Handles periodic boundary conditions
- Stack-based memory management for billion-particle simulations
- Supports RAMSES simulation format (adaptable to GADGET)

---

## Algorithm Description

### Step 1: Density Calculation

Particles are mapped onto a 3D grid using **Triangular-Shaped Cloud (TSC)** interpolation:

```
Cell size: TSC_CELL_SIZE = 0.004 cMpc/h
Smoothing: Gaussian_Smoothing_Length = 0.008 cMpc/h
```

The density field is smoothed using FFT-based Gaussian convolution (FFTW3 library).

### Step 2: Peak Identification

Local density maxima are identified where:
- Density > `PEAKTHRESHOLD` (2000 h^2 Msun/ckpc^3)
- Particle has highest density among its `NUMNEIGHBOR` (32) nearest neighbors

### Step 3: Water-Shedding

Each particle is assigned to the nearest peak by following the density gradient:
1. Build k-d tree for nearest neighbor queries
2. For each particle, find the neighbor with highest density
3. Recursively trace to the density peak
4. Particles converging to the same peak form a proto-core

### Step 4: Peak Merging

Underpopulated peaks are merged with nearby dominant peaks:
- Merge if separation < `MERGINGPEAKLENGTH` (4e-3 cMpc/h)
- Merge if member count < `MINCORENMEM` (30 particles)

### Step 5: Boundedness Iteration

For each core, test particle binding over `BOUNDITER` (4) iterations:

```
For each iteration:
    For each particle in core:
        Calculate kinetic energy (KE) in CoM frame
        Calculate potential energy (PE) via tree-based force
        If KE + PE > 0: mark as unbound, move to shell
    Recalculate core center-of-mass
    Break if no particles removed
```

### Step 6: Tidal Radius Calculation

The tidal radius is computed using an NFW profile lookup table:

```c
r_tidal = f(m/M, d/R_vir, c)
```

Where:
- `m/M`: subhalo-to-host mass ratio
- `d/R_vir`: distance from host center in virial radii
- `c`: host halo concentration parameter

### Step 7: Final Membership Assignment

Remaining unassigned particles are linked to cores using FoF with linking length `FOFLINK4MEMBERSHIP` (0.005 cMpc/h).

---

## Installation

### Prerequisites

- MPI implementation (OpenMPI, MPICH, or Intel MPI)
- Intel compilers (icc/icx) with OpenMP support
- FFTW3 library (single precision, with OpenMP)
- GNU Make

### Build Instructions

```bash
# Edit Makefile to set FFTW path
vim Makefile
# Set: FFTW = /path/to/fftw3

# Build the code
make clean
make all

# Verify executable
ls -la gfind.exe
```

### Compiler Flags

The default compilation uses:

```makefile
CC = mpiicx
FC = mpiifx
OPT = -O3 -qopenmp -DINDEX -DVarPM -DXYZDBL -DNENER=0 -DNPRE=8 -DREAD_SINK -DNCHEM=3 -DADV
LIBS = -lfftw3f_omp -lfftw3f -lm
```

Key preprocessor definitions:
| Flag | Description |
|------|-------------|
| `-DXYZDBL` | Use double precision for positions |
| `-DREAD_SINK` | Enable sink/AGN particle support |
| `-DNCHEM=3` | Track 3 chemical elements |
| `-DNMEG=90000L` | Allocate 90 GB memory per worker |

---

## Usage

### Basic Execution

```bash
# Single process (debugging)
./gfind.exe <snapshot_number>

# Example: process snapshot 5
./gfind.exe 5

# Parallel execution with 64 MPI processes
mpirun -np 64 ./gfind.exe 5

# With OpenMP threads
export OMP_NUM_THREADS=8
mpirun -np 64 ./gfind.exe 5
```

### Restart from Offset

For restarting or processing specific portions:

```bash
./gfind.exe <snapshot> <header_offset> <data_offset>

# Example: restart from specific file positions
./gfind.exe 5 12345678 98765432
```

### Expected Directory Structure

```
./FoF_Data/
└── FoF.XXXXX/
    ├── FoF_halo_cat.XXXXX       # Halo catalog (header)
    ├── FoF_member_particle.XXXXX # Particle data
    ├── GALCATALOG.LIST.XXXXX    # Output: ASCII catalog
    ├── GALFIND.DATA.XXXXX       # Output: Binary data
    └── background_ptl.XXXXX     # Output: Unassigned particles
```

---

## Input/Output Formats

### Input Files

#### FoF_halo_cat.XXXXX (Binary Header)

```c
float size;      // Box size [cMpc/h]
float hubble;    // Hubble parameter (H0/100)
float omep;      // Omega_matter
float omepb;     // Omega_baryon
float omeplam;   // Omega_Lambda
float amax;      // Maximum scale factor
float anow;      // Current scale factor

// Followed by array of HaloQ structures:
struct HaloQ {
    size_t np, npstar, npgas, npdm, npsink;  // Particle counts
    float x, y, z;           // Halo center [cMpc/h]
    float mass, mstar, mgas, mdm, msink;  // Masses [Msun/h]
    float vx, vy, vz;        // Bulk velocity [km/s]
};
```

#### FoF_member_particle.XXXXX (Binary Data)

For each halo, particles are stored sequentially:
1. `npdm` DmType structures
2. `npgas` GasType structures
3. `npsink` SinkType structures
4. `npstar` StarType structures

### Output Files

#### GALCATALOG.LIST.XXXXX (Binary Catalog)

```c
// Per FoF halo:
struct HaloInfo {
    int nsub;                    // Number of subhalos
    int ndm, nstar, nsink, ngas, npall;  // Particle counts
    double totm, mdm, mgas, msink, mstar; // Masses [Msun/h]
    double x, y, z, vx, vy, vz;  // Center-of-mass position/velocity
};

// Per subhalo (repeated nsub times):
struct SubInfo {
    int npdm, npgas, npsink, npstar, npall;
    double totm, mdm, mgas, msink, mstar;
    double x, y, z, vx, vy, vz;
};
```

#### GALFIND.DATA.XXXXX (Binary with Particles)

Same structure as GALCATALOG.LIST but includes full particle data after each SubInfo.

---

## Configuration Parameters

All physical parameters are in `params.h`:

### Core Algorithm Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `BOUNDITER` | 4 | Iterations for boundedness testing |
| `PEAKTHRESHOLD` | 2000 | Minimum density for peaks [h^2 Msun/ckpc^3] |
| `MINCORENMEM` | 30 | Minimum particles per confirmed core |
| `NUMNEIGHBOR` | 32 | Neighbors for density kernel |
| `NSHELLDIVIDE` | 10 | Number of iso-density shells |
| `MERGINGPEAKLENGTH` | 4e-3 | Peak merge separation [cMpc/h] |

### Density Grid Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TSC_CELL_SIZE` | 0.004 | TSC grid cell size [cMpc/h] |
| `Gaussian_Smoothing_Length` | 0.008 | Gaussian smoothing scale [cMpc/h] |
| `NCELLBUFF` | 10 | Boundary buffer cells (must be even) |
| `DENFLOOR` | 1.0 | Minimum density floor |

### Performance Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MAXTHREADS` | 64 | Maximum OpenMP threads |
| `DEEPSIZE` | 1024 | OMP parallelization threshold |
| `NOMPFoF` | 500000 | Particle count for OMP FoF |
| `FOFLINK4MEMBERSHIP` | 0.005 | Final FoF linking length [cMpc/h] |

### Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MAXNUMCORE` | 1,000,000 | Maximum cores per halo |
| `MAXNUMWATERSHEDDING` | 100,000,000 | Maximum water-shed iterations |

---

## Code Architecture

### MPI Communication Pattern

```
Rank 0 (Master)                    Ranks 1-N (Workers)
     |                                   |
     |-- Read halo from file             |
     |-- Wait for READY signal  <--------|-- Send READY
     |-- Send particle data    --------->|-- Receive particles
     |                                   |-- subhalo_den()
     |-- Wait for WRITING      <---------|-- Send WRITING + results
     |-- Write output                    |
     |-- Loop...                         |-- Loop...
```

### Key MPI Tags

| Tag | Value | Purpose |
|-----|-------|---------|
| `READY` | 1 | Worker ready signal |
| `WRITING` | 3 | Worker sending results |
| `NP_TAG` | 17 | Particle count |
| `R_TAG` | 32 | Position/particle data |
| `PT2H_TAG` | 55 | Particle-to-halo mapping |
| `MPEAK_TAG` | 96 | Number of identified peaks |

### Memory Management

Custom stack-based allocator (`Memory2.c`) to handle large datasets:

```c
Make_Total_Memory(size);  // Initialize pool
Malloc(size, &ptr);       // Allocate
Realloc(ptr, new_size);   // Resize
Free(ptr);                // Release
```

Worker ranks allocate `NMEG` (90 GB) by default; master allocates 80 GB.

---

## Performance Considerations

### Scalability

- **Strong scaling**: Up to ~1000 MPI ranks for 10^8 particle halos
- **Weak scaling**: Linear with number of halos

### Memory Requirements

Per worker rank:
- Base: ~1 GB for code and structures
- Per halo: ~200 bytes/particle during processing
- Peak: 2-3x particle data for tree structures

### Optimization Tips

1. **Balance OMP threads vs MPI ranks**: Typically 8-16 threads optimal
2. **Set NMEG appropriately**: Match to available memory per node
3. **Large halos dominate**: Consider sorting halos by size for load balance
4. **TSC_CELL_SIZE**: Smaller values increase resolution but cost O(N^3) memory

---

## File Descriptions

### Core Source Files

| File | Lines | Description |
|------|-------|-------------|
| `gfind.c` | 899 | Main program, MPI master-slave loop, I/O |
| `subhaloden.mod6.c` | ~2600 | Core algorithms: density, peaks, binding |
| `ost.c` | ~1000 | Barnes-Hut tree for FoF linking |
| `nnost.c` | ~1300 | k-d tree for nearest neighbor queries |
| `tsc_omp2.c` | 507 | TSC density interpolation (OpenMP) |
| `gsmooth.c` | 116 | FFT-based Gaussian smoothing |
| `mkRtidal.c` | 129 | NFW tidal radius lookup |
| `force_spline.mod3.c` | 189 | Gravitational force via spline |

### Header Files

| File | Description |
|------|-------------|
| `params.h` | Physical and algorithm parameters |
| `ramses.h` | RAMSES data structures (DmType, GasType, etc.) |
| `tree.h` | Tree structures (FoFTStruct, Coretype) |
| `header.h` | Physical constants, MPI tags |
| `defs.h` | Global variable definitions |
| `Memory.h` | Memory allocator interface |

### Support Files

| File | Description |
|------|-------------|
| `Memory2.c` | Stack-based memory allocator |
| `utils.c` | MPI utilities (BIG_MPI_Send/Recv) |
| `nrutil.c` | Numerical Recipes utilities |
| `b2l.c` | Big-to-Little Endian conversions |
| `spline.mod2.f` | Fortran spline interpolation |
| `Treewalk.near.c` | Tree traversal for neighbors |

---

## Particle Types

The code handles four particle families:

| Type | ID | Structure | Key Fields |
|------|----|-----------|------------|
| Dark Matter | 0 | `DmType` | x,y,z,vx,vy,vz,mass,id |
| Stars | 1 | `StarType` | + tp (birth time), zp (metallicity), chem[] |
| Sinks/AGN | 2 | `SinkType` | + Jx,Jy,Jz (angular momentum), dMsmbh |
| Gas | 4 | `GasType` | + den, temp, metallicity, cellsize |

---

## Physical Units

| Quantity | Unit |
|----------|------|
| Position | comoving Mpc/h |
| Velocity | km/s (peculiar) |
| Mass | Msun/h |
| Density | h^2 Msun/ckpc^3 |
| Time | Scale factor or Gyr |

---

## References

- **FoF Algorithm**: Davis et al. (1985), ApJ, 292, 371
- **NFW Profile**: Navarro, Frenk & White (1997), ApJ, 490, 493
- **Subhalo Finding**: Springel et al. (2001), MNRAS, 328, 726
- **RAMSES Code**: Teyssier (2002), A&A, 385, 337

---

## Additional Documentation

Detailed documentation is available in the `docs/` directory:

| Document | Description |
|----------|-------------|
| **[ALGORITHM.md](docs/ALGORITHM.md)** | In-depth algorithm descriptions: density estimation, peak finding, water-shedding, boundedness testing, tidal radius |
| **[DATA_STRUCTURES.md](docs/DATA_STRUCTURES.md)** | Complete reference for particle types, tree structures, output formats, and memory layout |
| **[CONFIGURATION.md](docs/CONFIGURATION.md)** | Parameter tuning guide for different resolutions and science cases |
| **[TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md)** | Common issues, error messages, and solutions |
| **[EXAMPLES.md](docs/EXAMPLES.md)** | Practical usage examples, job scripts, and Python code for reading outputs |

### Related Components

This package is part of the PGalF (Parallel Galaxy Finder) suite:

- **NewDD**: Domain decomposition preprocessor for RAMSES data (see `../NewDD/README.md`)
- **GADGET Adaptation**: Guide for adapting to GADGET HDF5 format (see `../NewDD/README4GADGET.md`)

---

## License

This code is provided for scientific research purposes. Please cite the appropriate papers when using this software for publications.

---

## Contact

For questions or bug reports, please contact the development team or open an issue in the repository.
