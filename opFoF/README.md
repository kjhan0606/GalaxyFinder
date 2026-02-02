# opFoF: OpenMP/MPI-Parallel Friends-of-Friends Halo Finder

**opFoF** (Open Friends-of-Friends) is a high-performance, MPI-parallel halo finder designed to identify dark matter halos and their constituent particles (dark matter, stars, gas, and sinks) from cosmological N-body and hydrodynamic simulations. It is specifically optimized for reading **NewDD** (New Domain Decomposition) output data in RAMSES format.

---

## Quick Links

| Document | Description |
|----------|-------------|
| [Quick Start Guide](docs/QUICKSTART.md) | Get running in 5 minutes |
| [Algorithm Details](docs/ALGORITHM.md) | In-depth FoF algorithm description |
| [Data Structures](docs/DATA_STRUCTURES.md) | Complete structure reference |
| [Configuration Guide](docs/CONFIGURATION.md) | All configuration options |
| [Examples](docs/EXAMPLES.md) | Usage examples and workflows |
| [Output Format](docs/OUTPUT_FORMAT.md) | Binary format specification |
| [API Reference](docs/API.md) | Function documentation |
| [Troubleshooting](docs/TROUBLESHOOTING.md) | Common issues and solutions |

---

## Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Directory Structure](#directory-structure)
4. [Dependencies & Requirements](#dependencies--requirements)
5. [Compilation](#compilation)
6. [Usage](#usage)
7. [Input Format](#input-format)
8. [Output Format](#output-format)
9. [Algorithm Description](#algorithm-description)
10. [Configuration Options](#configuration-options)
11. [Utility Programs](#utility-programs)
12. [Performance Considerations](#performance-considerations)
13. [Troubleshooting](#troubleshooting)
14. [References](#references)

---

## Overview

The Friends-of-Friends (FoF) algorithm is a standard method in computational cosmology for identifying gravitationally bound structures (halos) in N-body simulations. The algorithm links particles that are within a specified linking length `b` (typically 0.2 times the mean inter-particle separation) into groups.

**opFoF** implements this algorithm with:
- MPI parallelization for distributed memory systems
- Tree-based spatial indexing for O(N log N) performance
- Support for periodic and non-periodic boundary conditions
- Multi-component particle handling (DM, stars, gas, sinks)
- Efficient memory management for large datasets

---

## Features

### Core Capabilities
- **MPI-Parallel Processing**: Distributes particle data across multiple MPI ranks for scalable halo finding
- **Tree-Based Spatial Indexing**: Uses oct-tree structures for efficient neighbor searches
- **Multi-Component Support**: Handles dark matter, stellar, gas, and sink (black hole) particles
- **Periodic Boundary Conditions**: Correctly identifies halos that wrap around simulation box edges
- **Domain Decomposition**: Z-slab based partitioning with boundary particle exchange
- **Large-Scale Support**: Handles billion-particle simulations with optimized memory management

### Particle Types Supported
| Type ID | Component | Properties |
|---------|-----------|------------|
| 2 | Dark Matter | position, velocity, mass, id, level |
| 1 | Stars | position, velocity, mass, id, age, metallicity |
| 3 | Gas | position, velocity, mass, density, temperature |
| 4 | Sinks | position, velocity, mass, birth time, angular momentum |

### Output Properties
For each identified halo:
- Total particle count and counts by component
- Total mass and mass by component
- Center-of-mass position (x, y, z)
- Center-of-mass velocity (vx, vy, vz)
- Individual member particle data

---

## Directory Structure

```
opFoF/
├── opfof.c                    # Main program entry point
├── Treewalk.fof.ordered.c     # Active tree walking implementation
├── Treewalk.fof.ordered.org.c # Original tree walking (backup)
├── Treewalk.fof.c             # Legacy tree walking variant
├── fof.c                      # Test/benchmark program
├── fof.h                      # Core data structure definitions
├── pmheader.h                 # Particle & simulation structures
├── Memory.c                   # Memory allocator implementation
├── Memory.h                   # Memory management interface
├── Memory2.c                  # Alternative memory manager
├── Time.c                     # Timing utilities
├── Time.h                     # Timing interface
├── ramses2read.c              # RAMSES format data reader
├── read.c                     # Generic particle reading routines
├── header.c                   # Header/metadata file I/O
├── migrate.c                  # Particle migration between MPI ranks
├── kjhmigrate.c               # Additional migration code
├── Makefile                   # Build configuration
├── Rules.make                 # Compilation flags and rules
├── params.h -> ../params.h    # Parameter definitions (symlink)
├── ramses.h -> ../ramses.h    # RAMSES structures (symlink)
│
├── # Analysis & Post-Processing Tools
├── fofmassfunc.c              # Mass function calculator
├── checkfof.c                 # Output validation utility
├── fofread.c                  # FoF catalog reader
├── how2read.c                 # Example reader for output format
├── how2readfof.c              # Detailed output format documentation
│
├── # Compiled Executables
├── opfof.exe                  # Main executable
├── fofread                    # Catalog reader
├── fofmassfunc                # Mass function tool
├── checkfof                   # Validation tool
├── how2readfof                # Output parser
│
└── mfrac/                     # Mass fraction analysis
    ├── mfrac.c                # Mass fraction calculator
    ├── dist.c                 # Distribution calculator
    └── check.c                # Validation checker
```

---

## Dependencies & Requirements

### Software Requirements
| Requirement | Version | Notes |
|-------------|---------|-------|
| MPI | OpenMPI/MPICH/Intel MPI | Intel MPI recommended |
| C Compiler | Intel ICC or GCC | Intel ICC preferred for optimization |
| Fortran Compiler | Intel IFORT or gfortran | Required for some components |
| Make | GNU Make 3.8+ | Build system |

### Hardware Requirements
- **Architecture**: 64-bit (uses `long long` indices for >2 billion particles)
- **Memory**: ~17 GB RAM per MPI rank (configurable via `-DNMEG`)
- **Storage**: Sufficient for input snapshots and output catalogs
- **Network**: High-bandwidth interconnect recommended for large MPI jobs

### Library Dependencies
- **libm**: Standard math library
- **libmyram.a**: Custom memory manager (included, pre-compiled)

---

## Compilation

### Quick Start

```bash
cd opFoF

# Clean previous build
make clean

# Build the main executable
make this

# Or rebuild everything from scratch
make all
```

### Compilation Flags

The `Rules.make` file contains all compilation flags:

```makefile
# Compiler selection (Intel MPI)
CC = mpiicc
F77 = mpiifort

# Optimization flags
OPT = -DINTEL -g

# Common flags
COMFLAGS = -DINDEX -DVarPM -DXYZDBL

# Feature flags
CDFLAGS = -DWGROUPSIZE=8 -DNMEG=17000L -DINCLUDE_TREE_FORCE \
          -D_LARGE_FILES -DSAVESLICE -DPMSEEDFORCE -DQUADHILBERT \
          -DNENER=0 -DNPRE=8 -DREAD_SINK -DNCHEM=9 -DNDUST=4
```

### Flag Descriptions

| Flag | Description |
|------|-------------|
| `-DINTEL` | Enable Intel compiler-specific optimizations |
| `-DINDEX` | Enable particle ID indexing |
| `-DVarPM` | Support variable particle mass |
| `-DXYZDBL` | Use double precision for coordinates |
| `-DWGROUPSIZE=8` | MPI work group size (8 processes per group) |
| `-DNMEG=17000L` | Memory pool size in MB (17 GB default) |
| `-DINCLUDE_TREE_FORCE` | Include tree force calculations |
| `-D_LARGE_FILES` | Support files larger than 2 GB |
| `-DREAD_SINK` | Enable sink particle reading |
| `-DQUADHILBERT` | Use quad precision for Hilbert curves |
| `-DNCHEM=9` | Number of chemical elements tracked |
| `-DNDUST=4` | Number of dust species |

### Customizing Memory Allocation

To change the memory pool size, modify `-DNMEG` in `Rules.make`:

```makefile
# For 32 GB per rank
CDFLAGS = ... -DNMEG=32000L ...

# For 8 GB per rank
CDFLAGS = ... -DNMEG=8000L ...
```

### Using GCC Instead of Intel

Modify `Rules.make`:

```makefile
CC = mpicc
F77 = mpif77
OPT = -g -O2  # Remove -DINTEL
```

---

## Usage

### Command-Line Syntax

opFoF supports two execution modes:

#### Mode 1: Single Snapshot Analysis

```bash
mpirun -np <nprocs> ./opfof.exe <snapshot_number> <num_files>
```

**Parameters:**
- `<snapshot_number>`: Snapshot index (e.g., 100 for snapshot 00100)
- `<num_files>`: Number of data files per snapshot

**Example:**
```bash
# Analyze snapshot 150 with 64 data files using 8 MPI ranks
mpirun -np 8 ./opfof.exe 150 64
```

#### Mode 2: Multiple Snapshot Range

```bash
mpirun -np <nprocs> ./opfof.exe <start_snap> <end_snap> <num_files>
```

**Parameters:**
- `<start_snap>`: First snapshot number
- `<end_snap>`: Last snapshot number
- `<num_files>`: Number of data files per snapshot

**Example:**
```bash
# Analyze snapshots 100-200 with 128 files each using 32 MPI ranks
mpirun -np 32 ./opfof.exe 100 200 128
```

### Input Data Location

By default, opFoF looks for input data in the current working directory with the following naming convention:

```
SN.<NNNNN>.<TYPE>.<MMMMM>.dat    # Particle data files
SN.<NNNNN>.<MMMMM>.info          # Info/header files
```

Where:
- `<NNNNN>` = 5-digit snapshot number (zero-padded)
- `<TYPE>` = DM, STAR, GAS, or SINK
- `<MMMMM>` = 5-digit file index (zero-padded)

### SLURM Batch Script Example

```bash
#!/bin/bash
#SBATCH --job-name=opfof
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=140G

module load intel-mpi

cd /path/to/data
mpirun -np 32 /path/to/opfof.exe 150 512
```

---

## Input Format

### RAMSES Binary Format

opFoF reads particle data in RAMSES binary format, produced by NewDD domain decomposition:

#### File Naming Convention

```
SN.00150.DM.00000.dat      # Dark matter file 0
SN.00150.DM.00001.dat      # Dark matter file 1
...
SN.00150.STAR.00000.dat    # Star particles file 0
SN.00150.GAS.00000.dat     # Gas particles file 0
SN.00150.SINK.00000.dat    # Sink particles file 0
SN.00150.00000.info        # Metadata/header file
```

#### Info File Format

The `.info` file contains cosmological parameters in text format:

```
define BoxSize  100.0       # Box size in Mpc/h
define Omega_m  0.3089      # Matter density parameter
define Omega_l  0.6911      # Dark energy density parameter
define Omega_b  0.0486      # Baryon density parameter
define H0       67.74       # Hubble parameter in km/s/Mpc
define aexp     0.5         # Current scale factor
define amax     1.0         # Maximum scale factor
```

#### Particle Data Structure

Each particle file contains binary arrays of particle properties:

**Dark Matter:**
```c
struct DmType {
    float vx, vy, vz;      // Velocity components
    float mass;            // Particle mass
    long long id;          // Unique particle ID
    int level;             // Refinement level
};
```

**Stars:**
```c
struct StarType {
    float vx, vy, vz;      // Velocity components
    float mass;            // Stellar mass
    long long id;          // Unique particle ID
    float age;             // Stellar age
    float metallicity;     // Metal content
};
```

**Gas:**
```c
struct GasType {
    float vx, vy, vz;      // Velocity components
    float mass;            // Gas mass
    float density;         // Gas density
    float temperature;     // Gas temperature
};
```

**Sinks:**
```c
struct SinkType {
    float vx, vy, vz;      // Velocity components
    float mass;            // Sink mass
    float birth_time;      // Formation time
    float lx, ly, lz;      // Angular momentum
};
```

---

## Output Format

### Output Files

opFoF produces two binary output files per snapshot:

```
FoF_halo_cat.<NNNNN>          # Halo catalog
FoF_member_particle.<NNNNN>   # Member particle data
```

### Halo Catalog Format (FoF_halo_cat)

#### Header Section (7 floats)

```c
struct CatalogHeader {
    float size;       // Box size in Mpc/h
    float hubble;     // Hubble parameter (H0/100)
    float omep;       // Omega matter
    float omepb;      // Omega baryon
    float omeplam;    // Omega lambda
    float amax;       // Maximum scale factor
    float anow;       // Current scale factor
};
```

#### Halo Records (HaloQ structure)

```c
struct HaloQ {
    size_t np;        // Total particle count
    size_t npstar;    // Star particle count
    size_t npgas;     // Gas particle count
    size_t npdm;      // Dark matter particle count
    size_t npsink;    // Sink particle count
    POSTYPE x, y, z;  // Center-of-mass position (Mpc/h)
    double mass;      // Total halo mass (Msun/h)
    double mstar;     // Stellar mass
    double mgas;      // Gas mass
    double mdm;       // Dark matter mass
    double msink;     // Sink mass
    float vx, vy, vz; // Center-of-mass velocity (km/s)
};
```

### Member Particle Format (FoF_member_particle)

Contains the same header followed by particle data grouped by halo:

```
[Header - 7 floats]
[Halo 1 DM particles]
[Halo 1 Gas particles]
[Halo 1 Sink particles]
[Halo 1 Star particles]
[Halo 2 DM particles]
...
```

### Reading Output in C

```c
#include <stdio.h>
#include "fof.h"

int main() {
    FILE *fp = fopen("FoF_halo_cat.00150", "rb");

    // Read header
    float header[7];
    fread(header, sizeof(float), 7, fp);

    printf("Box size: %.2f Mpc/h\n", header[0]);
    printf("Scale factor: %.4f\n", header[6]);

    // Read halos
    HaloQ halo;
    int count = 0;
    while(fread(&halo, sizeof(HaloQ), 1, fp) == 1) {
        if(halo.np > 30) {  // Filter small halos
            printf("Halo %d: Np=%zu, M=%.3e Msun/h, pos=(%.2f, %.2f, %.2f)\n",
                   count++, halo.np, halo.mass, halo.x, halo.y, halo.z);
        }
    }

    fclose(fp);
    return 0;
}
```

### Reading Output in Python

```python
import numpy as np
import struct

def read_fof_catalog(filename):
    """Read opFoF halo catalog."""
    with open(filename, 'rb') as f:
        # Read header
        header = struct.unpack('7f', f.read(28))
        box_size, hubble, omega_m, omega_b, omega_l, amax, anow = header

        # Define HaloQ structure
        # Note: Adjust dtype based on your POSTYPE (float32 or float64)
        halo_dtype = np.dtype([
            ('np', 'u8'), ('npstar', 'u8'), ('npgas', 'u8'),
            ('npdm', 'u8'), ('npsink', 'u8'),
            ('x', 'f8'), ('y', 'f8'), ('z', 'f8'),  # POSTYPE=double
            ('mass', 'f8'), ('mstar', 'f8'), ('mgas', 'f8'),
            ('mdm', 'f8'), ('msink', 'f8'),
            ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4')
        ])

        # Read all halos
        halos = np.fromfile(f, dtype=halo_dtype)

    return {
        'header': {
            'box_size': box_size, 'hubble': hubble,
            'omega_m': omega_m, 'omega_b': omega_b, 'omega_l': omega_l,
            'amax': amax, 'anow': anow
        },
        'halos': halos
    }

# Example usage
data = read_fof_catalog('FoF_halo_cat.00150')
halos = data['halos']

# Filter halos with > 100 particles
massive = halos[halos['np'] > 100]
print(f"Found {len(massive)} halos with >100 particles")
print(f"Most massive halo: {massive['mass'].max():.3e} Msun/h")
```

---

## Algorithm Description

### Friends-of-Friends Method

The FoF algorithm identifies halos by linking particles that are separated by less than a linking length `b`:

$$b = 0.2 \times \bar{d}$$

where $\bar{d}$ is the mean inter-particle separation:

$$\bar{d} = \left( \frac{m_p}{\rho_{crit} \times \Omega_m} \right)^{1/3}$$

### Implementation Details

#### 1. Tree Building (`FoF_Make_Tree`)

The algorithm constructs an oct-tree for efficient spatial queries:

```
1. Start with all particles in root node
2. Recursively subdivide nodes with > 8 particles
3. Each internal node stores:
   - Bounding extent (dist)
   - Center of mass (monox, monoy, monoz)
   - Pointers to children/siblings
4. Leaf nodes contain particle pointers
```

#### 2. Particle Linking (`new_fof_link` / `pnew_fof_link`)

```
1. For each unlinked particle (seed):
   a. Initialize linked[] array with seed
   b. Tree walk to find neighbors within link length
   c. Add unlinked neighbors to linked[]
   d. Recursively search from new members
   e. Assign halo index to all linked particles
2. Repeat until all particles are assigned
```

#### 3. Periodic Boundary Handling

For periodic simulations, the algorithm:
- Checks for halo members across box boundaries
- Uses `fof_open_periodic()` to unwrap coordinates
- Properly calculates center-of-mass for wrapped halos

#### 4. MPI Domain Decomposition

```
1. Divide simulation volume into Z-slabs
2. Each MPI rank handles subset of files
3. Initial boundary: zmin = myid * (box/nprocs)
4. Process internal halos first
5. Exchange boundary particles via MPI_Send/Recv
6. Merge boundary halos across rank boundaries
```

#### 5. Halo Property Calculation

For each identified halo:
```c
// Center of mass
x_cm = sum(m_i * x_i) / sum(m_i)
y_cm = sum(m_i * y_i) / sum(m_i)
z_cm = sum(m_i * z_i) / sum(m_i)

// Velocity (mass-weighted)
vx_cm = sum(m_i * vx_i) / sum(m_i)
vy_cm = sum(m_i * vy_i) / sum(m_i)
vz_cm = sum(m_i * vz_i) / sum(m_i)

// Component masses
M_total = sum(m_i)
M_DM = sum(m_i) for type==2
M_star = sum(m_i) for type==1
M_gas = sum(m_i) for type==3
M_sink = sum(m_i) for type==4
```

---

## Configuration Options

### Key Parameters (Hardcoded)

| Parameter | Default | Location | Description |
|-----------|---------|----------|-------------|
| `fof_link` | 0.2 | opfof.c | Linking length in units of mean separation |
| `NODE_HAVE_PARTICLE` | 8 | fof.h | Max particles per tree leaf node |
| `WGROUPSIZE` | 8 | Rules.make | MPI work group size |
| `NMEG` | 17000 | Rules.make | Memory pool size (MB) |

### Modifying Linking Length

The linking length is calculated in `ramses2read.c`:

```c
// link02 = 0.2 * (mass / (2.7755e11 * omega_m))^(1/3)
link02 = 0.2 * pow(mass / (2.7755e11 * omega_m), 1.0/3.0);
```

To change the linking length, modify the `0.2` factor.

### Minimum Halo Size

The default minimum halo size for output is 30 particles. To modify:

```c
// In output routines
if(haloq.np > 30) {  // Change 30 to desired minimum
    // Output halo
}
```

---

## Utility Programs

### fofread

Read and display halo catalog:

```bash
./fofread FoF_halo_cat.00150
```

### fofmassfunc

Calculate halo mass function:

```bash
./fofmassfunc FoF_halo_cat.00150 > mass_function.dat
```

Output columns: `log(M)`, `dn/dlogM`, `N(>M)`

### checkfof

Validate output file integrity:

```bash
./checkfof FoF_halo_cat.00150
```

### how2readfof

Parse and display detailed output format:

```bash
./how2readfof FoF_halo_cat.00150
```

### mfrac (Mass Fraction Tools)

```bash
cd mfrac

# Calculate mass fractions
./mfrac ../FoF_halo_cat.00150

# Calculate distributions
./dist ../FoF_halo_cat.00150

# Validate results
./check ../FoF_halo_cat.00150
```

---

## Performance Considerations

### Scaling Guidelines

| Simulation Size | Recommended MPI Ranks | Memory per Rank |
|-----------------|----------------------|-----------------|
| 256³ particles | 8 | 4 GB |
| 512³ particles | 32 | 8 GB |
| 1024³ particles | 128 | 16 GB |
| 2048³ particles | 512 | 32 GB |

### Optimization Tips

1. **Match ranks to files**: Best performance when `nfiles` is divisible by `nprocs`

2. **Memory allocation**: Set `-DNMEG` slightly above expected peak usage

3. **I/O bandwidth**: Use parallel file system (Lustre, GPFS) for large runs

4. **Load balancing**: Ensure uniform particle distribution across Z-slabs

5. **MPI tuning**: For Intel MPI:
   ```bash
   export I_MPI_SHM_LMT=shm
   export I_MPI_DAPL_UD=enable
   ```

### Common Performance Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Memory overflow | Insufficient `-DNMEG` | Increase memory pool |
| Slow I/O | Serial file access | Use parallel filesystem |
| Load imbalance | Non-uniform data | Adjust domain decomposition |
| MPI timeout | Large boundary exchange | Increase timeout limits |

---

## Troubleshooting

### Common Errors

**"Cannot allocate memory"**
```
Solution: Increase -DNMEG in Rules.make and recompile
```

**"File not found: SN.00100..."**
```
Solution: Verify input file naming convention matches expected format
Check: ls SN.00100.*.dat
```

**"MPI_Send: Message truncated"**
```
Solution: The code uses BIG_MPI_Send for large messages.
If still failing, check MPI buffer limits.
```

**"Segmentation fault in tree walk"**
```
Solution: Usually indicates corrupted input data.
Verify input files with checkfof utility.
```

### Debug Mode

Compile with debug symbols:

```makefile
OPT = -DINTEL -g -O0 -DDEBUG
```

Enable verbose output:
```bash
export OPFOF_VERBOSE=1
mpirun -np 8 ./opfof.exe 150 64
```

---

## References

### Algorithm

- Davis, M. et al. (1985). "The evolution of large-scale structure in a universe dominated by cold dark matter." ApJ, 292, 371-394.
- [FoF algorithm description](https://en.wikipedia.org/wiki/Friends-of-friends_algorithm)

### Related Documentation

- [NewDD README](../NewDD/README.md) - Domain decomposition tool
- [PGalF README](../NewGalFinder/README.md) - Galaxy finder using FoF output
- [RAMSES format](https://www.ics.uzh.ch/~teyssier/ramses/RAMSES.html)

### Citation

If you use opFoF in your research, please cite:

```bibtex
@software{opfof,
  author = {Han, K. J.},
  title = {opFoF: Open Friends-of-Friends Halo Finder},
  year = {2024},
  url = {https://github.com/your-repo/opFoF}
}
```

---

## Contact

For questions, bug reports, or feature requests:
- Create an issue in the repository
- Contact the maintainers directly

---

## License

[Specify your license here]

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial release |
| 1.1 | 2025 | Multi-component support |
| 1.2 | 2025 | Improved MPI scalability |
