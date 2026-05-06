# PGalF Configuration Guide

This guide provides detailed instructions for configuring PGalF parameters for different simulation types and science goals.

## Table of Contents

1. [Parameter Overview](#1-parameter-overview)
2. [Resolution-Dependent Settings](#2-resolution-dependent-settings)
3. [Science-Case Configurations](#3-science-case-configurations)
4. [Memory Optimization](#4-memory-optimization)
5. [Parallel Execution Tuning](#5-parallel-execution-tuning)
6. [Compile-Time Options](#6-compile-time-options)
7. [Dark-Galaxy Detection](#9-dark-galaxy-detection)

---

## 1. Parameter Overview

All runtime parameters are defined in `params.h`. Changing these requires recompilation.

### Parameter Categories

| Category | Parameters | Impact |
|----------|------------|--------|
| Core Algorithm | BOUNDITER, PEAKTHRESHOLD, MINCORENMEM | Galaxy detection sensitivity |
| Density Grid | TSC_CELL_SIZE, Gaussian_Smoothing_Length | Resolution vs memory tradeoff |
| Performance | MAXTHREADS, DEEPSIZE, NOMPFoF | Parallelization efficiency |
| Limits | MAXNUMCORE, MAXNUMWATERSHEDDING | Memory allocation bounds |

---

## 2. Resolution-Dependent Settings

### 2.1 High-Resolution Simulations (dx < 100 pc)

For zoom-in simulations or high-resolution cosmological boxes:

```c
// params.h settings for high resolution
#define TSC_CELL_SIZE 0.001              // 1 kpc/h cells
#define Gaussian_Smoothing_Length 0.002  // 2 kpc/h smoothing
#define PEAKTHRESHOLD 5000.L             // Higher density threshold
#define MINCORENMEM 50                   // More particles required
#define NUMNEIGHBOR 64                   // More neighbors for density
#define MERGINGPEAKLENGTH 1.e-3          // Smaller merge distance
```

**Rationale**:
- Smaller cells resolve fine structure
- Higher threshold avoids noise peaks
- More neighbors improve density estimates

### 2.2 Standard Resolution (100 pc < dx < 1 kpc)

Default settings work well for typical cosmological simulations:

```c
// params.h default settings
#define TSC_CELL_SIZE 0.004              // 4 kpc/h cells
#define Gaussian_Smoothing_Length 0.008  // 8 kpc/h smoothing
#define PEAKTHRESHOLD 2000.L             // Standard threshold
#define MINCORENMEM 30                   // 30 particles minimum
#define NUMNEIGHBOR 32                   // Standard neighbors
#define MERGINGPEAKLENGTH 4.e-3          // 4 kpc/h merge
```

### 2.3 Low-Resolution Simulations (dx > 1 kpc)

For large-volume, low-resolution simulations:

```c
// params.h settings for low resolution
#define TSC_CELL_SIZE 0.010              // 10 kpc/h cells
#define Gaussian_Smoothing_Length 0.020  // 20 kpc/h smoothing
#define PEAKTHRESHOLD 1000.L             // Lower threshold
#define MINCORENMEM 20                   // Fewer particles
#define NUMNEIGHBOR 32                   // Keep standard
#define MERGINGPEAKLENGTH 1.e-2          // 10 kpc/h merge
```

### 2.4 Resolution Quick Reference

| Simulation Type | TSC_CELL_SIZE | Smoothing | PEAKTHRESHOLD | MINCORENMEM |
|----------------|---------------|-----------|---------------|-------------|
| Ultra-high res (<50 pc) | 0.0005 | 0.001 | 10000 | 100 |
| High res (50-200 pc) | 0.001 | 0.002 | 5000 | 50 |
| Standard (200 pc-1 kpc) | 0.004 | 0.008 | 2000 | 30 |
| Low res (1-5 kpc) | 0.010 | 0.020 | 1000 | 20 |
| Very low res (>5 kpc) | 0.020 | 0.040 | 500 | 15 |

---

## 3. Science-Case Configurations

### 3.1 Dwarf Galaxy Detection

To detect low-mass dwarf galaxies:

```c
#define PEAKTHRESHOLD 500.L       // Lower threshold for faint systems
#define MINCORENMEM 15            // Accept smaller galaxies
#define MINSTELLARMASS 1.e5       // Include very small stellar systems
#define BOUNDITER 6               // More iterations for accuracy
```

**Note**: Lower thresholds may increase false positives. Verify with visual inspection.

### 3.2 Cluster Galaxy Finding

For galaxy clusters with many subhalos:

```c
#define MAXNUMCORE 5000000        // Allow more cores per halo
#define NSHELLDIVIDE 20           // Finer shell division
#define MERGINGPEAKLENGTH 2.e-3   // Smaller merge to preserve substructure
#define FOFLINK4MEMBERSHIP 0.003  // Tighter linking
```

### 3.3 High-z Galaxy Finding (z > 3)

For high-redshift simulations where galaxies are smaller:

```c
#define TSC_CELL_SIZE 0.002       // Finer grid (smaller galaxies)
#define Gaussian_Smoothing_Length 0.004
#define PEAKTHRESHOLD 3000.L      // Higher threshold (denser universe)
#define MINCORENMEM 20            // Accept smaller systems
```

### 3.4 AGN Host Identification

When focusing on AGN/sink particle hosts:

```c
// Ensure sink particles are included
// Compile with -DREAD_SINK flag

#define MINCORENMEM 10            // Small cores around AGN
#define BOUNDITER 4               // Standard iterations
```

### 3.5 Satellite Galaxy Studies

For studying satellite populations:

```c
#define NSHELLDIVIDE 15           // More shells for better membership
#define BOUNDITER 5               // More unbinding iterations
#define FOFLINK4MEMBERSHIP 0.004  // Moderate linking
```

---

## 4. Memory Optimization

### 4.1 Memory Budget Calculation

Estimate memory requirements before running:

```
Memory (GB) ≈ (N_particles × 0.7 + N_grid³ × 12) / 10⁹
```

Where:
- N_particles = maximum particles per halo
- N_grid = max(halo_size / TSC_CELL_SIZE) for each dimension

### 4.2 Reducing Memory Usage

For memory-constrained systems:

```c
// Increase cell size to reduce grid memory
#define TSC_CELL_SIZE 0.008       // Larger cells

// Reduce maximum threads
#define MAXTHREADS 16             // Fewer thread-local arrays

// Reduce core limit
#define MAXNUMCORE 100000         // Limit core allocation
```

### 4.3 Memory Allocation (Makefile)

Adjust per-worker memory in Makefile:

```makefile
# For systems with 64 GB RAM per node, 8 workers per node
DFLAGS = -DNMEG=8000L    # 8 GB per worker

# For systems with 128 GB RAM per node, 4 workers per node
DFLAGS = -DNMEG=30000L   # 30 GB per worker

# For large memory nodes (256+ GB)
DFLAGS = -DNMEG=90000L   # 90 GB per worker (default)
```

### 4.4 Grid Size Limits

The TSC grid size is limited by available memory:

| TSC_CELL_SIZE | Max Halo Size | Grid Memory |
|---------------|---------------|-------------|
| 0.002 | 0.5 Mpc/h | ~8 GB |
| 0.004 | 1.0 Mpc/h | ~8 GB |
| 0.008 | 2.0 Mpc/h | ~8 GB |
| 0.010 | 2.5 Mpc/h | ~8 GB |

For very large halos (>2 Mpc/h), increase TSC_CELL_SIZE accordingly.

---

## 5. Parallel Execution Tuning

### 5.1 MPI vs OpenMP Balance

Optimal balance depends on halo size distribution:

| Scenario | MPI Ranks | OMP Threads | Rationale |
|----------|-----------|-------------|-----------|
| Many small halos | High (64+) | Low (4) | Distribute halos widely |
| Few large halos | Low (8-16) | High (16-32) | Parallel within halos |
| Mixed distribution | Medium (32) | Medium (8) | Balanced approach |

### 5.2 OpenMP Parameters

```c
// For large halos with many particles
#define DEEPSIZE 512              // Start OMP earlier in tree
#define NOMPFoF 100000            // Lower threshold for OMP FoF

// For small halos (reduce overhead)
#define DEEPSIZE 2048             // Delay OMP parallelization
#define NOMPFoF 1000000           // Higher threshold
```

### 5.3 Job Script Example

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

export OMP_NUM_THREADS=8
export OMP_STACKSIZE=512M
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun ./gfind.exe 100
```

### 5.4 Load Balancing Considerations

Halos are processed in file order. For better load balance:

1. **Pre-sort halos by size** (in FoF output)
2. Use dynamic MPI scheduling (already implemented)
3. Monitor with timing output:

```bash
# Enable timing in code
# Compile with -DDEBUG=1 -DLOG=1
mpirun -np 32 ./gfind.exe 100 2>&1 | tee timing.log
```

---

## 6. Compile-Time Options

### 6.1 Using configure

The Makefile is generated by the top-level `configure` script. Use configure options instead of editing the Makefile directly:

```bash
# From the GalaxyFinder top-level directory:
./configure [options]
make galfinder
```

Key configure options for NewGalFinder:

| Option | Default | Description |
|--------|---------|-------------|
| `--galfinder-cc=CC` | `mpiicx` | C compiler |
| `--galfinder-fc=FC` | `mpiifx` | Fortran compiler |
| `--fftw=PATH` | `/home/kjhan/local` | FFTW path |
| `--opt=FLAGS` | `-O3` | Optimization flags |
| `--debug` | - | Use `-g` debug flags |
| `--openmp=FLAGS` | `-qopenmp` | OpenMP flags |
| `--galfinder-nmeg=N` | `90000` | Memory pool size (MB) |
| `--galfinder-nchem=N` | `3` | Chemical species count |
| `--galfinder-fftw-libs=L` | `-lfftw3f_omp -lfftw3f` | FFTW link flags |

### 6.2 Makefile Flags Reference

| Flag | Effect |
|------|--------|
| `-DXYZDBL` | Double precision positions (recommended) |
| `-DREAD_SINK` | Include sink/AGN particles |
| `-DNCHEM=N` | Track N chemical elements |
| `-DNENER=N` | Track N radiation energy groups |
| `-DNPRE=8` | RAMSES precision indicator |
| `-DDEBUG=1` | Enable debug output |
| `-DLOG=1` | Enable logging |
| `-DVarPM` | Variable particle mass |
| `-DINDEX` | Include particle indices |

### 6.3 Common Configurations

**Standard RAMSES with stars and AGN**:
```bash
./configure --galfinder-nchem=3
```

**Debug build**:
```bash
./configure --debug
```

### 6.4 Recompilation After Changes

After modifying `params.h`:

```bash
make clean
make galfinder
```

After changing configure options:

```bash
./configure [new options]
make clean
make galfinder
```

---

## 7. Parameter Validation

### 7.1 Consistency Checks

Ensure parameters are consistent:

```c
// Smoothing must be ≥ 2 × cell size
assert(Gaussian_Smoothing_Length >= 2 * TSC_CELL_SIZE);

// NUMNEIGHBOR must be ≤ MAX_NUM_NEAR (in tree.h)
assert(NUMNEIGHBOR <= MAX_NUM_NEAR);  // MAX_NUM_NEAR = 32 by default

// NCELLBUFF must be even
assert(NCELLBUFF % 2 == 0);
```

### 7.2 Diagnostic Output

Monitor key metrics during runs:

```
P0 has omep=0.3 rng=1 r1kineticfact=XXX potentfact=XXX
Now passing through 1000: 50000 particles
```

Check that:
- `potentfact` is reasonable (~10^-10 to 10^-8)
- Particle counts match expectations
- Processing time per halo is consistent

---

## 8. Quick Configuration Templates

### Template 1: Default (Copy-Paste Ready)

```c
// params.h - Default configuration
#define BOUNDITER 4
#define COREDENRESOLUTION (1.e-3)
#define PEAKTHRESHOLD 2000.L
#define MINCORENMEM 30
#define NUMNEIGHBOR 32
#define NSHELLDIVIDE 10
#define MERGINGPEAKLENGTH 4.e-3
#define MINCORESTARMASS -1
#define MINSTELLARMASS 2.e6
#define MAXNUMCORE 1000000
#define TSC_CELL_SIZE 0.004
#define Gaussian_Smoothing_Length 0.008
#define DENFLOOR 1.
#define NCELLBUFF 10
#define DEEPSIZE 1024
#define NOMPFoF 500000
#define FOFLINK4MEMBERSHIP 0.005
#define MAXNUMWATERSHEDDING 100000000L
#define MAXTHREADS 64
```

### Template 2: High-Resolution Zoom

```c
// params.h - High-resolution zoom simulation
#define BOUNDITER 5
#define COREDENRESOLUTION (1.e-4)
#define PEAKTHRESHOLD 5000.L
#define MINCORENMEM 50
#define NUMNEIGHBOR 64
#define NSHELLDIVIDE 15
#define MERGINGPEAKLENGTH 1.e-3
#define MINCORESTARMASS -1
#define MINSTELLARMASS 1.e5
#define MAXNUMCORE 500000
#define TSC_CELL_SIZE 0.001
#define Gaussian_Smoothing_Length 0.002
#define DENFLOOR 1.
#define NCELLBUFF 10
#define DEEPSIZE 512
#define NOMPFoF 200000
#define FOFLINK4MEMBERSHIP 0.002
#define MAXNUMWATERSHEDDING 100000000L
#define MAXTHREADS 64
```

### Template 3: Large Volume Survey

```c
// params.h - Large volume, low resolution
#define BOUNDITER 3
#define COREDENRESOLUTION (1.e-2)
#define PEAKTHRESHOLD 1000.L
#define MINCORENMEM 20
#define NUMNEIGHBOR 32
#define NSHELLDIVIDE 8
#define MERGINGPEAKLENGTH 1.e-2
#define MINCORESTARMASS -1
#define MINSTELLARMASS 1.e7
#define MAXNUMCORE 2000000
#define TSC_CELL_SIZE 0.010
#define Gaussian_Smoothing_Length 0.020
#define DENFLOOR 1.
#define NCELLBUFF 10
#define DEEPSIZE 2048
#define NOMPFoF 1000000
#define FOFLINK4MEMBERSHIP 0.010
#define MAXNUMWATERSHEDDING 100000000L
#define MAXTHREADS 64
```

---

## 9. Dark-Galaxy Detection

Detection of DM-only subhalos (dark galaxies) is built into the standard
pipeline. Peak finding always runs on the unified density field

```
ρ_total = ρ_star + DM_DENSITY_WEIGHT × ρ_DM
```

with each component smoothed at its own scale. Cores are classified as dark
post-hoc when their stellar content falls below the FoF-gate scale. See
[ALGORITHM.md §8](ALGORITHM.md#8-dark-galaxy-detection) for the
algorithmic details.

The feature is controlled at runtime by `DM_DENSITY_WEIGHT` in `params.h`:

- `DM_DENSITY_WEIGHT = 0.1f` (default) — DM grid contributes; pure-DM peaks
  are detected and classified post-hoc as dark galaxies.
- `DM_DENSITY_WEIGHT = 0.f` — the DM TSC, DM Gaussian smoothing, and combine
  loop are all short-circuited; the run is byte-equivalent to the
  stellar-only finder. No compile flag toggle is required.

### 9.1 Active parameters (unified path)

| Parameter | Default | Effect |
|-----------|---------|--------|
| `DM_DENSITY_WEIGHT` | `0.1f` | Weight of the smoothed DM grid in `ρ_total`. Set to `0.f` to disable the DM contribution entirely (DM TSC + smoothing are short-circuited). |
| `DM_GAUSSIAN_SMOOTHING_LENGTH` | `0.012` (cMpc/h) | Smoothing scale applied to the DM grid. Should track the inter-particle separation of a 50–100 DM particle clump at the lightest DM species. |
| `DM_TSC_CELL_SIZE` | `TSC_CELL_SIZE` | TSC cell size for the DM grid. Increase to lower memory cost on DM-rich halos. |
| `MINDMMASS` | `10 × MINSTELLARMASS` | DM mass at which a stellar-empty halo is allowed to enter the grid path. Acts as a FoF-level gate only — does not select per-core. |

### 9.2 Legacy parameters (stand-alone DM-only path)

These are still defined in `params.h` for callers of `lagFindDarkCore` (a
stand-alone DM-only peak finder kept for back-compat), but they are **not**
referenced by `subhalo_den()` under the unified path.

| Parameter | Default | Note |
|-----------|---------|------|
| `DM_PEAKTHRESHOLD` | `1.e3` | DM-only grid threshold (legacy). |
| `DM_MERGINGPEAKLENGTH` | `5.e-3` (cMpc/h) | DM-only peak merge distance (legacy). |
| `DM_MINCORENMEM` | `50` | Minimum DM particles per dark core (legacy). |
| `STAR_DM_DEDUP_LENGTH` | `5.e-3` (cMpc/h) | Star/DM peak dedup distance from the removed dual-pass approach. |

### 9.3 Tuning recipes

**Want exactly the stellar-only result, but keep the unified call path:**
set `DM_DENSITY_WEIGHT 0.f`. The DM grid is allocated only when nonzero,
so this is the cheapest configuration and is byte-equivalent to
`lagFindStellarCore` for stellar peaks.

**Cluster cores where DM swamps galaxies:** lower
`DM_DENSITY_WEIGHT` to ~0.05 or raise `PEAKTHRESHOLD`. The
`ρ_total = ρ_star + 0.1 ρ_DM` design intentionally puts a pure-DM
ρ_DM = 1000 peak at the same combined density as a luminous core with
ρ_star = 100, ρ_DM = 0; tighten this if your cluster has many DM
concentrations between bright galaxies.

**Known dark subhalos are missed:** raise `DM_DENSITY_WEIGHT` toward 0.2–0.3,
or shrink `DM_GAUSSIAN_SMOOTHING_LENGTH` so individual subhalos are not
washed out into a smooth halo background. Watch for spurious peaks.

**High-resolution simulations:** scale `DM_GAUSSIAN_SMOOTHING_LENGTH` with
the DM particle mass — roughly halve it when the DM particle mass drops by
8×, the same factor you would use for `Gaussian_Smoothing_Length` against
the stellar mass resolution.
