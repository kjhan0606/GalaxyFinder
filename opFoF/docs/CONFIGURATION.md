# opFoF Configuration Guide

This document provides comprehensive information about configuring opFoF for different simulation setups, hardware configurations, and scientific requirements.

---

## Table of Contents

1. [Compile-Time Configuration](#compile-time-configuration)
2. [Runtime Configuration](#runtime-configuration)
3. [Memory Configuration](#memory-configuration)
4. [MPI Configuration](#mpi-configuration)
5. [Input Data Configuration](#input-data-configuration)
6. [Output Configuration](#output-configuration)
7. [Algorithm Parameters](#algorithm-parameters)
8. [Performance Tuning](#performance-tuning)
9. [Platform-Specific Settings](#platform-specific-settings)
10. [Configuration Examples](#configuration-examples)

---

## Compile-Time Configuration

### Rules.make Configuration

All compile-time options are controlled through `Rules.make`:

```makefile
# Compiler selection
CC = mpiicc                  # Intel MPI C compiler
F77 = mpiifort               # Intel MPI Fortran compiler

# Optimization level
OPT = -DINTEL -g             # Debug mode with Intel optimizations

# Common flags (always enabled)
COMFLAGS = -DINDEX -DVarPM -DXYZDBL

# Feature flags
CDFLAGS = -DWGROUPSIZE=8 -DNMEG=17000L -DINCLUDE_TREE_FORCE \
          -D_LARGE_FILES -DSAVESLICE -DPMSEEDFORCE -DQUADHILBERT \
          -DNENER=0 -DNPRE=8 -DREAD_SINK -DNCHEM=9 -DNDUST=4
```

### Flag Reference Table

#### Essential Flags

| Flag | Default | Description |
|------|---------|-------------|
| `-DINDEX` | ON | Enable particle ID indexing |
| `-DVarPM` | ON | Variable particle mass support |
| `-DXYZDBL` | ON | Double precision coordinates |

#### Memory Flags

| Flag | Default | Description |
|------|---------|-------------|
| `-DNMEG=N` | 17000L | Memory pool size in MB |
| `-D_LARGE_FILES` | ON | Support files > 2GB |

#### Particle Type Flags

| Flag | Default | Description |
|------|---------|-------------|
| `-DREAD_SINK` | ON | Enable sink particle reading |
| `-DNCHEM=N` | 9 | Number of chemical species |
| `-DNDUST=N` | 4 | Number of dust species |

#### Parallelization Flags

| Flag | Default | Description |
|------|---------|-------------|
| `-DWGROUPSIZE=N` | 8 | MPI work group size |
| `-DQUADHILBERT` | ON | Quad-precision Hilbert curves |

#### Algorithm Flags

| Flag | Default | Description |
|------|---------|-------------|
| `-DINCLUDE_TREE_FORCE` | ON | Include tree force calculations |
| `-DPMSEEDFORCE` | ON | Enable PM seed force |
| `-DSAVESLICE` | ON | Save vertical slices |

### Enabling/Disabling Features

To disable a feature, remove or comment out the flag:

```makefile
# Disable sink particle reading
# CDFLAGS = ... -DREAD_SINK ...  # Commented out
CDFLAGS = ... # -DREAD_SINK removed
```

To enable conditional compilation:

```makefile
# Enable debug output
CDFLAGS = ... -DDEBUG -DVERBOSE ...
```

### Precision Configuration

#### Double Precision Positions (Recommended)

```makefile
COMFLAGS = -DINDEX -DVarPM -DXYZDBL
```

This sets `POSTYPE` to `double`, providing:
- 15-16 significant digits
- Essential for large boxes (>1 Gpc)
- Slightly higher memory usage

#### Single Precision Positions

```makefile
COMFLAGS = -DINDEX -DVarPM
# Remove -DXYZDBL
```

This sets `POSTYPE` to `float`:
- 6-7 significant digits
- Suitable for small boxes (<500 Mpc)
- 33% memory reduction for positions

---

## Runtime Configuration

### Command-Line Arguments

```bash
# Single snapshot mode
./opfof.exe <snapshot_number> <num_files>

# Multi-snapshot mode
./opfof.exe <start_snap> <end_snap> <num_files>
```

### Environment Variables

```bash
# OpenMP threads (if OpenMP enabled)
export OMP_NUM_THREADS=4

# Memory limits
export MALLOC_MMAP_THRESHOLD_=1048576

# MPI settings
export I_MPI_SHM_LMT=shm
export I_MPI_FABRICS=shm:ofi

# Debug output
export OPFOF_DEBUG=1
export OPFOF_VERBOSE=1
```

---

## Memory Configuration

### Memory Pool Size

The `-DNMEG` flag sets the memory pool size in megabytes:

```makefile
# Examples:
-DNMEG=8000L     # 8 GB per MPI rank
-DNMEG=17000L    # 17 GB per MPI rank (default)
-DNMEG=32000L    # 32 GB per MPI rank
-DNMEG=64000L    # 64 GB per MPI rank
```

### Memory Estimation

Calculate required memory per rank:

```
Memory = (N_particles / N_ranks) × bytes_per_particle
       + tree_overhead
       + communication_buffers
       + safety_margin

Where:
- bytes_per_particle ≈ 100 bytes
- tree_overhead ≈ 60 bytes × N_particles
- communication_buffers ≈ 10% of particle memory
- safety_margin ≈ 20%
```

**Example calculation:**
```
Simulation: 1024³ particles = 1.07 billion particles
MPI ranks: 128
Particles per rank: 8.4 million

Memory per rank:
  Particles: 8.4M × 100 bytes = 840 MB
  Tree nodes: 8.4M × 60 bytes = 504 MB
  Comm buffers: 134 MB
  Safety: 296 MB
  Total: ~1.8 GB

Recommended -DNMEG=4000L (with overhead)
```

### Memory Monitoring

Enable memory tracking:

```makefile
CDFLAGS = ... -DMEMORY_DEBUG ...
```

This prints memory usage statistics during execution.

---

## MPI Configuration

### Work Group Size

```makefile
-DWGROUPSIZE=N
```

Recommendations:
- **Small clusters (< 64 cores)**: `-DWGROUPSIZE=4`
- **Medium clusters (64-256 cores)**: `-DWGROUPSIZE=8`
- **Large clusters (> 256 cores)**: `-DWGROUPSIZE=16`

### MPI Implementation Settings

#### Intel MPI

```bash
# Environment setup
source /opt/intel/oneapi/setvars.sh

# Recommended settings
export I_MPI_SHM_LMT=shm
export I_MPI_DAPL_UD=enable
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=omp
```

#### OpenMPI

```bash
# Use with GCC
mpirun --bind-to core --map-by socket -np 32 ./opfof.exe 100 128
```

#### MPICH

```bash
# Use hydra launcher
mpiexec -n 32 -bind-to core ./opfof.exe 100 128
```

### Process Placement

For optimal performance:

```bash
# One MPI rank per NUMA node
mpirun -np 4 -ppn 1 ./opfof.exe 100 64

# Multiple ranks per node
mpirun -np 32 -ppn 8 ./opfof.exe 100 256
```

---

## Input Data Configuration

### File Naming Convention

opFoF expects input files in this format:

```
SN.<NNNNN>.<TYPE>.<MMMMM>.dat
SN.<NNNNN>.<MMMMM>.info
```

Where:
- `NNNNN` = 5-digit zero-padded snapshot number
- `TYPE` = DM, STAR, GAS, or SINK
- `MMMMM` = 5-digit zero-padded file index

### Modifying File Patterns

To change the naming convention, edit `ramses2read.c`:

```c
// Current pattern
sprintf(filename, "SN.%05d.%s.%05d.dat", nstep, type_str, file_idx);

// Alternative pattern (example)
sprintf(filename, "snapshot_%03d/%s_particles_%03d.bin", nstep, type_str, file_idx);
```

### Data Directory

By default, opFoF reads from the current working directory. To specify a different path:

```c
// In opfof.c or ramses2read.c
#define DATA_DIR "/path/to/simulation/data/"
sprintf(filename, DATA_DIR "SN.%05d.%s.%05d.dat", ...);
```

### Supported Particle Types

Configure which particle types to read:

```c
// In opfof.c
#define READ_DM    1    // Always read dark matter
#define READ_STAR  1    // Read stellar particles
#define READ_GAS   1    // Read gas particles
#define READ_SINK  1    // Read sink particles (requires -DREAD_SINK)
```

---

## Output Configuration

### Output Files

Default output file names:

```
FoF_halo_cat.<NNNNN>          # Halo catalog
FoF_member_particle.<NNNNN>   # Member particles
```

### Modifying Output Path

```c
// In Treewalk.fof.ordered.c
#define OUTPUT_DIR "/path/to/output/"
sprintf(catfile, OUTPUT_DIR "FoF_halo_cat.%05d", nstep);
sprintf(memfile, OUTPUT_DIR "FoF_member_particle.%05d", nstep);
```

### Minimum Halo Size

Only halos with more than N particles are written:

```c
// In output routines
#define MIN_HALO_PARTICLES 30  // Default threshold

if (halo.np > MIN_HALO_PARTICLES) {
    WriteHalo(&halo);
}
```

To change:
```c
#define MIN_HALO_PARTICLES 50   // More restrictive
#define MIN_HALO_PARTICLES 10   // Less restrictive
```

### Output Verbosity

Control output detail level:

```c
#define OUTPUT_ALL_PROPERTIES    // Include all halo properties
#define OUTPUT_MEMBER_PARTICLES  // Include member particle data
// #define OUTPUT_PARTICLE_IDS   // Include particle IDs (large files)
```

---

## Algorithm Parameters

### Linking Length

The FoF linking length is set in `ramses2read.c`:

```c
// Default: b = 0.2 × mean inter-particle separation
double b_factor = 0.2;
link02 = b_factor * pow(mass / (2.7755e11 * omega_m), 1.0/3.0);
```

To modify:
```c
// More aggressive linking (larger halos)
double b_factor = 0.25;

// Less aggressive linking (smaller halos)
double b_factor = 0.15;
```

### Tree Parameters

```c
// Maximum particles per leaf node
#define NODE_HAVE_PARTICLE 8

// Tree opening angle (for future force calculations)
#define TREE_OPENING_ANGLE 0.5
```

### Boundary Buffer Zone

```c
// Extra zone near boundaries for halo detection
#define BOUNDARY_BUFFER (2.0 * link_length)
```

---

## Performance Tuning

### Optimization Levels

```makefile
# Debug mode (slow, with symbols)
OPT = -DINTEL -g -O0

# Balanced (debugging possible, good performance)
OPT = -DINTEL -g -O2

# Maximum performance
OPT = -DINTEL -O3 -xHost -ipo

# Aggressive (may affect numerical stability)
OPT = -DINTEL -O3 -xHost -ipo -fast -no-prec-div
```

### Vectorization

Enable auto-vectorization:
```makefile
OPT = -DINTEL -O3 -xHost -qopt-report=5
```

### Loop Unrolling

```makefile
OPT = -DINTEL -O3 -funroll-loops
```

### Profile-Guided Optimization

```bash
# Step 1: Build with profiling
make clean
OPT="-DINTEL -O2 -prof-gen" make this

# Step 2: Run with representative data
mpirun -np 8 ./opfof.exe 100 64

# Step 3: Rebuild with profile data
make clean
OPT="-DINTEL -O3 -prof-use" make this
```

---

## Platform-Specific Settings

### Intel Clusters

```makefile
CC = mpiicc
F77 = mpiifort
OPT = -DINTEL -O3 -xHost -ipo
```

### AMD EPYC Systems

```makefile
CC = mpicc
F77 = mpif77
OPT = -O3 -march=znver2 -mtune=znver2
COMFLAGS = -DINDEX -DVarPM -DXYZDBL  # Remove -DINTEL
```

### IBM Power Systems

```makefile
CC = mpixlc
F77 = mpixlf
OPT = -O3 -qarch=pwr9 -qtune=pwr9
COMFLAGS = -DINDEX -DVarPM -DXYZDBL
```

### NVIDIA GPU Clusters (Future)

```makefile
# Placeholder for GPU support
CC = nvcc
NVFLAGS = -arch=sm_80 -O3
```

### Cray Systems

```makefile
CC = cc  # Cray compiler wrapper
F77 = ftn
OPT = -O3 -hfp3
```

---

## Configuration Examples

### Example 1: Small Test Run

For quick testing on a workstation:

```makefile
# Rules.make
CC = mpicc
F77 = mpif77
OPT = -g -O0 -DDEBUG

COMFLAGS = -DINDEX -DVarPM -DXYZDBL
CDFLAGS = -DWGROUPSIZE=2 -DNMEG=4000L -D_LARGE_FILES -DREAD_SINK
```

```bash
# Run command
mpirun -np 2 ./opfof.exe 0 4
```

### Example 2: Production Run on HPC Cluster

For large-scale production:

```makefile
# Rules.make
CC = mpiicc
F77 = mpiifort
OPT = -DINTEL -O3 -xHost -ipo

COMFLAGS = -DINDEX -DVarPM -DXYZDBL
CDFLAGS = -DWGROUPSIZE=16 -DNMEG=30000L -DINCLUDE_TREE_FORCE \
          -D_LARGE_FILES -DPMSEEDFORCE -DQUADHILBERT \
          -DREAD_SINK -DNCHEM=9 -DNDUST=4
```

```bash
# SLURM job script
#!/bin/bash
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --mem=240G

module load intel-mpi
srun ./opfof.exe 100 200 1024
```

### Example 3: Memory-Constrained System

For systems with limited RAM:

```makefile
# Rules.make
CC = mpicc
F77 = mpif77
OPT = -O2

COMFLAGS = -DINDEX -DVarPM  # Remove -DXYZDBL for smaller memory
CDFLAGS = -DWGROUPSIZE=4 -DNMEG=8000L -D_LARGE_FILES
```

### Example 4: Dark Matter Only

For DM-only simulations:

```makefile
# Rules.make
CDFLAGS = -DWGROUPSIZE=8 -DNMEG=16000L -D_LARGE_FILES \
          -DQUADHILBERT
# Remove -DREAD_SINK, -DNCHEM, -DNDUST
```

### Example 5: High-Resolution Cosmological Box

For 4096³ simulations:

```makefile
# Rules.make
CC = mpiicc
F77 = mpiifort
OPT = -DINTEL -O3 -xHost -ipo -parallel

COMFLAGS = -DINDEX -DVarPM -DXYZDBL
CDFLAGS = -DWGROUPSIZE=32 -DNMEG=60000L -DINCLUDE_TREE_FORCE \
          -D_LARGE_FILES -DQUADHILBERT -DREAD_SINK
```

```bash
# Run on 512 MPI ranks
mpirun -np 512 ./opfof.exe 100 2048
```

---

## Troubleshooting Configuration

### Common Configuration Errors

**Error: "undefined reference to `MPI_Init`"**
```
Solution: Ensure MPI compiler wrappers are used
  CC = mpiicc    (Intel MPI)
  CC = mpicc     (OpenMPI/MPICH)
```

**Error: "NMEG value too large"**
```
Solution: Reduce memory pool or use long suffix
  -DNMEG=17000L  (correct)
  -DNMEG=17000   (may overflow on 32-bit systems)
```

**Warning: "implicit declaration of function"**
```
Solution: Ensure all headers are included
  Check #include statements in source files
```

**Error: "cannot allocate memory"**
```
Solution:
  1. Reduce -DNMEG value
  2. Use more MPI ranks
  3. Check system memory limits (ulimit -v)
```

### Validation

After changing configuration, validate with:

```bash
# Check compilation
make clean && make this 2>&1 | tee build.log
grep -i error build.log
grep -i warning build.log

# Quick test run
mpirun -np 2 ./opfof.exe 0 2
```
