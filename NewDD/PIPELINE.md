# PGalF Pipeline Integration Guide

This document describes how NewDD integrates with the PGalF (Parallel Galaxy Finder) pipeline.

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        PGalF Pipeline Workflow                               │
└─────────────────────────────────────────────────────────────────────────────┘

  RAMSES Simulation          GADGET Simulation
        │                          │
        └──────────┬───────────────┘
                   │
                   ▼
         ┌─────────────────┐
         │     NewDD       │  ← Domain Decomposition
         │  (This Package) │     (Z-axis slab partitioning)
         └────────┬────────┘
                  │
                  │  Output: FoF_Data/NewDD.XXXXX/
                  │  - DM particles per slab
                  │  - Star particles per slab
                  │  - Gas cells per slab
                  │  - Sink particles per slab
                  │
                  ▼
         ┌─────────────────┐
         │      opFoF      │  ← Halo Finding
         │   (OpenMP FoF)  │     (Friends-of-Friends algorithm)
         └────────┬────────┘
                  │
                  │  Output: Halo catalogs
                  │  - Halo positions, masses
                  │  - Particle membership
                  │
                  ▼
         ┌─────────────────┐
         │  NewGalFinder   │  ← Galaxy Identification
         │                 │     (Subhalo/galaxy detection)
         └────────┬────────┘
                  │
                  │  Output: Galaxy catalogs
                  │  - Galaxy positions
                  │  - Stellar masses
                  │  - Gas content
                  │
                  ▼
         ┌─────────────────┐
         │   GalCenter     │  ← Galaxy Center Finding
         │                 │     (Precise center calculation)
         └────────┬────────┘
                  │
                  │  Output: Final catalogs
                  │  - Refined galaxy centers
                  │  - Galaxy properties
                  │
                  ▼
          Science Analysis
```

---

## Pipeline Components

### 1. NewDD (Domain Decomposition)

**Purpose**: Reorganize simulation data into spatially-sorted slabs for efficient parallel processing.

**Input**:
- RAMSES: `output_XXXXX/` directories
- GADGET: HDF5 snapshot files (with modification)

**Output**: `FoF_Data/NewDD.XXXXX/`
- `SN.XXXXX.[TYPE].[SLAB].dat` - Binary data files
- `SN.XXXXX.[SLAB].info` - Metadata files

**Key Features**:
- Z-axis slab decomposition
- Unit conversion to physical units
- MPI parallel I/O

---

### 2. opFoF (Halo Finder)

**Purpose**: Identify dark matter halos using the Friends-of-Friends algorithm.

**Input**: NewDD output (`SN.XXXXX.DM.*.dat`)

**Key Parameters**:
- Linking length (typically 0.2 * mean particle separation)
- Minimum halo mass (particle count threshold)

**Output**:
- Halo catalogs with positions, masses, velocities
- Particle-to-halo membership

---

### 3. NewGalFinder (Galaxy Identification)

**Purpose**: Identify galaxies within halos using stellar and gas components.

**Input**:
- opFoF halo catalogs
- NewDD star and gas data

**Key Features**:
- Subhalo detection
- Galaxy property calculation
- Merger tree construction

---

### 4. GalCenter (Center Finding)

**Purpose**: Calculate precise galaxy centers using advanced algorithms.

**Input**: NewGalFinder galaxy catalogs

**Methods**:
- Shrinking sphere
- Potential minimum
- Density peak

---

## Data Flow Specifications

### NewDD → opFoF Interface

**File Format**: Binary arrays

```
SN.XXXXX.DM.[SLAB].dat
└── Array of DmType structures
    ├── x, y, z      (double) - Position [cMpc/h]
    ├── vx, vy, vz   (double) - Velocity [km/s]
    ├── mass         (double) - Mass [Msun/h]
    └── id           (long)   - Particle ID
```

**Coordinate Convention**:
- Positions: Comoving Mpc/h
- Velocities: Physical km/s
- Masses: Solar masses/h

### NewDD → NewGalFinder Interface

**Star Data**: `SN.XXXXX.STAR.[SLAB].dat`
```
StarType structure:
├── Position, velocity, mass (same as DM)
├── tp       - Birth time (scale factor)
├── zp       - Metallicity (mass fraction)
├── chem     - Chemical abundances (if NCHEM > 0)
├── mass0    - Initial mass [Msun/h]
├── birth_d  - Dust birth parameter
└── partp    - Particle pointer/index
```

**Gas Data**: `SN.XXXXX.GAS.[SLAB].dat`
```
GasType structure:
├── Position, velocity
├── cellsize    - Cell size [cMpc/h]
├── den         - Density [code units]
├── temp        - Temperature [K/mu]
├── metallicity - Metal fraction
├── chem        - Chemical composition (if NCHEM > 0)
├── dust        - Dust species mass fractions (if NDUST > 0)
├── mass        - Cell mass [Msun/h]
├── potent      - Gravitational potential
└── fx, fy, fz  - Gravitational force
```

---

## Running the Full Pipeline

### Step 1: Domain Decomposition (NewDD)

```bash
cd NewDD
mpirun -np 32 ./newdd.exe 100 256
```

### Step 2: Halo Finding (opFoF)

```bash
cd ../opFoF
mpirun -np 32 ./opfof.exe 100 256
```

### Step 3: Galaxy Finding (NewGalFinder)

```bash
cd ../NewGalFinder
mpirun -np 32 ./newgalfinder.exe 100
```

### Step 4: Center Finding (GalCenter)

```bash
cd ../GalCenter
./galcenter.exe 100
```

---

## Parameter Consistency

Ensure consistent parameters across pipeline components:

| Parameter | NewDD | opFoF | NewGalFinder |
|-----------|-------|-------|--------------|
| Snapshot | `ISTEP` | `ISTEP` | `ISTEP` |
| Nsplit | `NSPLIT` | `NSPLIT` | - |
| Box size | from info | from info | from info |
| Cosmology | from info | from info | from info |

---

## Batch Processing Multiple Snapshots

### Shell Script Example

```bash
#!/bin/bash
# process_pipeline.sh

SNAPSHOTS="50 60 70 80 90 100"
NSPLIT=128
NP=32

for SNAP in $SNAPSHOTS; do
    echo "Processing snapshot $SNAP..."

    # Step 1: NewDD
    cd /path/to/NewDD
    mpirun -np $NP ./newdd.exe $SNAP $NSPLIT

    # Step 2: opFoF
    cd /path/to/opFoF
    mpirun -np $NP ./opfof.exe $SNAP $NSPLIT

    # Step 3: NewGalFinder
    cd /path/to/NewGalFinder
    mpirun -np $NP ./newgalfinder.exe $SNAP

    # Step 4: GalCenter
    cd /path/to/GalCenter
    ./galcenter.exe $SNAP

    echo "Snapshot $SNAP complete."
done
```

### SLURM Job Array

```bash
#!/bin/bash
#SBATCH --job-name=pgalf
#SBATCH --array=50,60,70,80,90,100
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --time=04:00:00

SNAP=$SLURM_ARRAY_TASK_ID
NSPLIT=128

cd /path/to/NewDD
srun ./newdd.exe $SNAP $NSPLIT

cd /path/to/opFoF
srun ./opfof.exe $SNAP $NSPLIT

# Continue with other steps...
```

---

## Directory Structure

Recommended project layout:

```
PGalF_Project/
├── output_00050/          # RAMSES snapshots
├── output_00060/
├── ...
│
├── NewDD/                 # Domain decomposition
│   ├── newdd.exe
│   └── ...
│
├── opFoF/                 # Halo finder
│   ├── opfof.exe
│   └── ...
│
├── NewGalFinder/          # Galaxy finder
│   ├── newgalfinder.exe
│   └── ...
│
├── GalCenter/             # Center finder
│   ├── galcenter.exe
│   └── ...
│
├── FoF_Data/              # Pipeline outputs
│   ├── NewDD.00050/       # NewDD output
│   ├── NewDD.00060/
│   ├── Halos.00050/       # opFoF output
│   ├── Galaxies.00050/    # NewGalFinder output
│   └── Centers.00050/     # GalCenter output
│
└── Analysis/              # Science analysis scripts
    ├── plot_halos.py
    └── ...
```

---

## Troubleshooting Pipeline Issues

### Data Consistency Check

```python
import numpy as np

# Check NewDD output
# Adjust NCHEM value to match your compilation flags
NCHEM = 9
dm_dtype = np.dtype([('x','f8'),('y','f8'),('z','f8'),
                     ('vx','f8'),('vy','f8'),('vz','f8'),
                     ('mass','f8'),('id','i8'),
                     ('levelp','i4'),('family','i1'),('tag','i1'),
                     ('tp','f8'),('zp','f8'),('chem','f8',NCHEM),
                     ('mass0','f8'),('birth_d','f8'),('partp','i4')])

dm = np.fromfile('FoF_Data/NewDD.00050/SN.00050.DM.00000.dat',
                 dtype=dm_dtype)

# Sanity checks
print(f"Particles: {len(dm)}")
print(f"Position range: [{dm['x'].min():.3f}, {dm['x'].max():.3f}]")
print(f"Mass range: [{dm['mass'].min():.3e}, {dm['mass'].max():.3e}]")
print(f"Family values: {np.unique(dm['family'])}")  # Should be [1] for DM
```

### Common Issues

| Stage | Issue | Solution |
|-------|-------|----------|
| NewDD → opFoF | Wrong particle count | Check info file matches binary |
| opFoF → NewGalFinder | Missing halos | Lower minimum mass threshold |
| NewGalFinder → GalCenter | No galaxies | Check star particle availability |

---

## Performance Optimization

### Pipeline Bottlenecks

1. **NewDD**: I/O bound (reading RAMSES files)
   - Use SSD storage
   - Optimize WGROUPSIZE

2. **opFoF**: Compute bound (FoF linking)
   - Increase OpenMP threads
   - Optimize linking length

3. **NewGalFinder**: Memory bound (particle sorting)
   - Ensure sufficient RAM
   - Process fewer slabs per job

### Recommended Resources

| Component | CPUs | Memory | Time (typical) |
|-----------|------|--------|----------------|
| NewDD | 32 | 64 GB | 30 min |
| opFoF | 32 | 128 GB | 1-2 hr |
| NewGalFinder | 16 | 64 GB | 30 min |
| GalCenter | 4 | 16 GB | 10 min |

---

## See Also

- [README.md](README.md) - Full NewDD documentation
- [CONFIGURATION.md](CONFIGURATION.md) - Advanced configuration
- [README4GADGET.md](README4GADGET.md) - GADGET support
