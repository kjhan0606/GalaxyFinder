# GalCenter Module

A post-processing utility for computing galaxy center positions from the PGalF (Parallel Galaxy Finder) pipeline.

## Overview

GalCenter calculates precise center-of-mass positions for galaxies identified by the NewGalFinder halo-finding stage. It processes binary galaxy data files and computes three distinct center positions for each galaxy:

- **Stellar Center** - Center of mass of star particles
- **Dark Matter Center** - Center of mass of dark matter particles
- **Gas Center** - Center of mass of gas particles

## Pipeline Context

```
RAMSES Output → NewDD → opFoF → NewGalFinder → GalCenter → Galaxy Catalogs
                                                    ↑
                                              (this module)
```

## Files

| File | Description |
|------|-------------|
| `galcenter.c` | Main MPI master program - distributes galaxy data to workers |
| `find_centers.c` | Core algorithm for computing center positions |
| `galcenter.h` | Local header defining Center and GalCenter structs |
| `tree.h` | Tree data structure definitions |
| `check.c` | Binary file validation utility |
| `checkcen.c` | Output center sanity checker |
| `Makefile` | Build configuration |
| `Rules.make` | Compiler settings and flags |

### Symlinked Headers (from ../NewGalFinder/)

- `Memory.h` - Memory management interface
- `header.h` - Common constants and macros
- `defs.h` - Data structure definitions (HaloInfo, SubInfo, etc.)
- `hfind.h` - Galaxy finder function declarations
- `utils.c` - Shared utility functions

## Building

### Using configure (Recommended)

```bash
# From the top-level GalaxyFinder directory:
./configure                               # Default settings
./configure --debug                       # Debug build
./configure --galcenter-nmeg=10000        # Custom NMEG

# Build GalCenter
make galcenter
```

#### Key configure options for GalCenter

| Option | Default | Description |
|--------|---------|-------------|
| `--cc=CC` | `mpicc` | C compiler |
| `--fc=FC` | `mpifort` | Fortran compiler |
| `--opt=FLAGS` | `-O3` | Optimization flags |
| `--debug` | - | Use `-g` debug flags |
| `--openmp=FLAGS` | `-qopenmp` | OpenMP flags |
| `--galcenter-nmeg=N` | `5000` | Memory pool size (MB) |
| `--galcenter-nchem=N` | `9` | Number of chemical species |
| `--galcenter-ndust=N` | `4` | Number of dust species |

Run `./configure --help` for the full list of options.

### Direct Build (Alternative)

```bash
cd GalCenter
make all      # Build galcenter.exe
make new      # Clean and rebuild
make clean    # Remove object files and executables
```

**Requirements:**
- MPI compiler (mpicc or mpiicx)
- OpenMP support
- Math library (-lm)

## Usage

```bash
mpirun -np <N> ./galcenter.exe <snapshot_number>
```

**Example:**
```bash
mpirun -np 8 ./galcenter.exe 50    # Process snapshot 50
```

## Input/Output

**Input:**
- `FoF_Data/FoF.NNNNN/GALFIND.DATA.NNNNN` - Binary galaxy catalog from NewGalFinder

**Output:**
- `FoF_Data/FoF.NNNNN/GALFIND.CENTER.NNNNN` - Binary file containing GalCenter structs

## Algorithm

The center-finding algorithm uses iterative refinement:

1. Calculate initial mass-weighted center-of-mass
2. Estimate radius based on mass ratio
3. Iteratively refine (up to 200 iterations):
   - Select particles within current sphere
   - Recalculate center using enclosed particles
   - Update radius based on enclosed mass
   - Converge when position shift < 0.005 and radius shift < 0.5%

**Special handling:**
- Particles < 30: Simple mass-weighted center
- Particles = 0: Marked as NO_CENTER_POS (-999999)
- Periodic boundaries: Automatic wrapping for edge-crossing halos

## Data Structures

```c
typedef struct GalCenter {
    idtype haloid;      // Parent halo ID
    idtype subgalid;    // Subgal/subhalo ID within halo
    Center gal;         // Star center (x, y, z, R)
    Center dmhalo;      // Dark matter center
    Center gas;         // Gas center
} GalCenter;
```

## MPI Architecture

Master-slave design:
- **Master (rank 0):** Reads catalog, distributes work, collects results, writes output
- **Workers (rank 1+):** Receive particle data, compute centers, return results

## Units

- Position: comoving Mpc/h
- Velocity: km/s
- Mass: Msun/h

## Configuration

Key parameters are set by the top-level `./configure` script:
- `NMEG=5000L` - Memory per processor (change with `--galcenter-nmeg=N`)
- OpenMP enabled via `-qopenmp` (change with `--openmp=FLAGS` or `--no-openmp`)
- Optimization via `-O3` (change with `--opt=FLAGS` or `--debug`)

## Related Modules

- `../NewGalFinder/` - Upstream galaxy finder
- `../opFoF/` - Friends-of-Friends halo finder
- `../NewDD/` - Domain decomposition
