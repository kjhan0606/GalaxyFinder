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

### 1. Domain Decomposition
```bash
cd NewDD
make all
mpirun -np 8 ./newdd.exe <snapshot> <nsplit>
```

### 2. Halo Finding
```bash
cd opFoF
make this
mpirun -np 8 ./opfof.exe <snapshot> <nfiles>
```

### 3. Galaxy Finding
```bash
cd NewGalFinder
make all
mpirun -np 64 ./gfind.exe <snapshot>
```

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
