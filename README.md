# PGalF Suite: Galaxy Finding Pipeline for RAMSES Simulations

A comprehensive suite for identifying galaxies and halos from RAMSES cosmological simulation outputs. The pipeline processes raw simulation data through domain decomposition, halo finding, and galaxy identification.

## Pipeline Overview

```
RAMSES Output → NewDD → opFoF → NewGalFinder → Galaxy Catalogs
               (Slab)   (Halo)   (Galaxy)
```

## Components

| Directory | Description |
|-----------|-------------|
| **NewDD** | RAMSES domain decomposition into Z-directional slabs for parallel processing |
| **opFoF** | MPI-parallel Friends-of-Friends (FoF) halo finder |
| **NewGalFinder** | Parallel galaxy/subhalo finder using density peaks and water-shedding |
| **GalCenter** | Galaxy center identification utilities |

## Quick Start

For a dark-matter-only RAMSES snapshot written with `LONGINT` particle IDs
and particle potentials, build NewDD and opFoF with one matching ABI:

```bash
MPI_CC=mpicc ./build_dmo_pipeline.sh
```

Do not mix a NewDD slab catalogue and an opFoF executable built with
different particle-layout flags. opFoF first links particles within each MPI
slab domain, exchanges only the face bands with adjacent ranks, and propagates
component labels until they converge. Resolved halo pieces then move through
an adjacent-rank ring to their unique owner; no `MPI_Alltoall[v]` or rank-0
particle gather is used. This permits a halo to span multiple MPI domains
without duplicate ownership or fragmentation.

The distributed path was regression-tested on matched one- and two-fluid DMO
catalogues with 8, 16, 32, and 64 MPI ranks. For each fluid representation,
all rank counts produced identical particle membership, zero duplicate member
IDs, and the same halo and member totals.

### 0. Configure Build System
```bash
# Default configuration (Intel OneAPI mpiicx/mpiifx compilers, optimized build)
./configure

# Debug build
./configure --debug

# Custom compilers and FFTW path
./configure --cc=mpicc --fc=mpifort --fftw=/path/to/fftw

# See all options
./configure --help
```

### 1. Build All Components
```bash
make all          # Build everything
# Or build individually:
make galcenter    # Build GalCenter only
make galfinder    # Build NewGalFinder only
make newdd        # Build NewDD only
```

### 2. Domain Decomposition
```bash
mpirun -np 8 ./NewDD/newdd.exe <snapshot> <nsplit>
```

### 3. Halo Finding
```bash
cd opFoF
make this
mpirun -np 8 ./opfof.exe <snapshot> <nfiles>
```

### 4. Galaxy Finding
```bash
mpirun -np 64 ./NewGalFinder/gfind.exe <snapshot>
```

## Build System

The project uses a `configure` script to generate Makefiles for all subdirectories from `Makefile.in` templates.

| File | Description |
|------|-------------|
| `configure` | Build configuration script |
| `Makefile.in` | Top-level Makefile template |
| `GalCenter/Makefile.in` | GalCenter Makefile template |
| `NewGalFinder/Makefile.in` | NewGalFinder Makefile template |
| `NewDD/Makefile.in` | NewDD Makefile template |

### Configure Options

| Option | Default | Description |
|--------|---------|-------------|
| `--cc=CC` | `mpiicx` | C compiler for GalCenter/NewDD |
| `--fc=FC` | `mpiifx` | Fortran compiler for GalCenter/NewDD |
| `--galfinder-cc=CC` | `mpiicx` | C compiler for NewGalFinder |
| `--galfinder-fc=FC` | `mpiifx` | Fortran compiler for NewGalFinder |
| `--fftw=PATH` | `/home/kjhan/local` | FFTW installation path |
| `--opt=FLAGS` | `-O3` | Optimization flags |
| `--debug` | - | Use `-g` debug flags |
| `--openmp=FLAGS` | `-qopenmp` | OpenMP flags |
| `--no-openmp` | - | Disable OpenMP |

Per-component options (e.g., `--galcenter-nmeg=N`, `--galfinder-nchem=N`, `--newdd-ndust=N`) are also available. Run `./configure --help` for the full list.

## Shared Files

- `ramses.h` - RAMSES data structure definitions
- `params.h` - Algorithm parameters
- `libmyram.a` - Custom memory management library

## Documentation

Each subdirectory contains detailed documentation:
- [NewDD/README.md](NewDD/README.md) - Domain decomposition details
- [opFoF/README.md](opFoF/README.md) - FoF algorithm and usage
- [NewGalFinder/README.md](NewGalFinder/README.md) - Galaxy finder algorithm

## Requirements

- MPI (OpenMPI, MPICH, or Intel MPI)
- Intel or GCC compilers with OpenMP support
- FFTW3 library (for NewGalFinder)

## Output Units

| Quantity | Unit |
|----------|------|
| Position | comoving Mpc/h |
| Velocity | km/s |
| Mass | Msun/h |

## License

For scientific research purposes.
