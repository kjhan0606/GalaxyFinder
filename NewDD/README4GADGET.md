# Adapting NewDD for GADGET HDF5 Format

This guide explains how to modify the NewDD package to read GADGET-style HDF5 simulation data instead of RAMSES format, enabling the same slab decomposition functionality for GADGET simulations.

---

## Table of Contents

1. [Overview](#overview)
2. [GADGET HDF5 File Structure](#gadget-hdf5-file-structure)
3. [Data Mapping: GADGET to NewDD](#data-mapping-gadget-to-newdd)
4. [Required Modifications](#required-modifications)
5. [New Files to Create](#new-files-to-create)
6. [Makefile Changes](#makefile-changes)
7. [Implementation Details](#implementation-details)
8. [Unit Conversion](#unit-conversion)
9. [Example Implementation](#example-implementation)
10. [Testing](#testing)

---

## Overview

The NewDD package currently reads RAMSES simulation outputs, which use:
- Fortran binary format files
- AMR (Adaptive Mesh Refinement) grid structure
- Multiple files per CPU domain

GADGET simulations use:
- HDF5 format (modern versions) or Fortran binary (legacy)
- Particle-based SPH (no AMR grid)
- Single or multiple files per snapshot

The key modifications involve:
1. Creating new HDF5 reading routines
2. Adapting the metadata parsing
3. Mapping GADGET particle types to NewDD structures
4. Adjusting unit conversions

---

## GADGET HDF5 File Structure

### File Naming Convention

```
snapshot_XXX.hdf5              # Single file
snapshot_XXX.0.hdf5            # Multiple files (file 0)
snapshot_XXX.1.hdf5            # Multiple files (file 1)
...
```

### HDF5 Group Structure

```
snapshot_XXX.hdf5
├── Header                     # Simulation metadata (attributes)
├── PartType0/                 # Gas particles
│   ├── Coordinates            # Position [N, 3]
│   ├── Velocities             # Velocity [N, 3]
│   ├── Masses                 # Mass [N] (or single value in header)
│   ├── ParticleIDs            # Unique IDs [N]
│   ├── InternalEnergy         # Thermal energy [N]
│   ├── Density                # SPH density [N]
│   ├── SmoothingLength        # SPH smoothing length [N]
│   ├── ElectronAbundance      # (optional) [N]
│   ├── NeutralHydrogenAbundance # (optional) [N]
│   ├── StarFormationRate      # (optional) [N]
│   └── Metallicity            # (optional) [N] or [N, NumMetals]
├── PartType1/                 # Dark matter (halo)
│   ├── Coordinates
│   ├── Velocities
│   ├── Masses                 # (often uniform, stored in header)
│   └── ParticleIDs
├── PartType2/                 # Disk particles (optional)
├── PartType3/                 # Bulge particles (optional)
├── PartType4/                 # Stars
│   ├── Coordinates
│   ├── Velocities
│   ├── Masses
│   ├── ParticleIDs
│   ├── StellarFormationTime   # Scale factor of formation
│   ├── Metallicity
│   └── InitialMass            # (optional)
└── PartType5/                 # Black holes/sinks
    ├── Coordinates
    ├── Velocities
    ├── Masses
    ├── ParticleIDs
    ├── BH_Mass                 # True BH mass
    └── BH_Mdot                 # Accretion rate
```

### Header Attributes

```
Header/
├── NumPart_ThisFile           # [6] Particle count per type in this file
├── NumPart_Total              # [6] Total particles per type (low bits)
├── NumPart_Total_HighWord     # [6] Total particles per type (high bits)
├── MassTable                  # [6] Particle masses (if uniform)
├── Time                       # Scale factor (cosmological) or time
├── Redshift                   # Current redshift
├── BoxSize                    # Simulation box size
├── NumFilesPerSnapshot        # Number of files per snapshot
├── Omega0                     # Matter density parameter
├── OmegaLambda                # Dark energy density parameter
├── HubbleParam                # Hubble parameter h = H0/100
├── Flag_Sfr                   # Star formation enabled
├── Flag_Cooling               # Cooling enabled
├── Flag_StellarAge            # Stellar ages recorded
├── Flag_Metals                # Metallicity recorded
├── Flag_Feedback              # Feedback enabled
└── UnitLength_in_cm           # Length unit in cm
    UnitMass_in_g              # Mass unit in grams
    UnitVelocity_in_cm_per_s   # Velocity unit in cm/s
```

---

## Data Mapping: GADGET to NewDD

### Particle Type Mapping

| GADGET PartType | NewDD Type | Description |
|-----------------|------------|-------------|
| PartType0 | GasType | Gas/SPH particles |
| PartType1 | DmType (family=1) | Dark matter (halo) |
| PartType2 | DmType (family=1) | Dark matter (disk, optional) |
| PartType3 | DmType (family=1) | Dark matter (bulge, optional) |
| PartType4 | StarType (family=2) | Star particles |
| PartType5 | SinkType | Black holes/sinks |

### Field Mapping

#### Dark Matter (PartType1 → DmType)

| GADGET Field | NewDD Field | Notes |
|--------------|-------------|-------|
| Coordinates[:, 0] | x | Apply unit conversion |
| Coordinates[:, 1] | y | Apply unit conversion |
| Coordinates[:, 2] | z | Apply unit conversion |
| Velocities[:, 0] | vx | Apply sqrt(a) for peculiar velocity |
| Velocities[:, 1] | vy | |
| Velocities[:, 2] | vz | |
| Masses | mass | From MassTable if uniform |
| ParticleIDs | id | |
| - | levelp | Set to 0 (no AMR in GADGET) |
| - | family | Set to 1 (DM) |

#### Stars (PartType4 → StarType)

| GADGET Field | NewDD Field | Notes |
|--------------|-------------|-------|
| Coordinates | x, y, z | |
| Velocities | vx, vy, vz | |
| Masses | mass | |
| ParticleIDs | id | |
| StellarFormationTime | tp | Scale factor of formation |
| Metallicity | zp | Total metallicity |
| InitialMass | mass0 | If available |
| - | family | Set to 2 (star) |

#### Gas (PartType0 → GasType)

| GADGET Field | NewDD Field | Notes |
|--------------|-------------|-------|
| Coordinates | x, y, z | |
| Velocities | vx, vy, vz | |
| Masses | mass | |
| Density | den | SPH density |
| InternalEnergy | temp | Convert to temperature |
| Metallicity | metallicity | |
| SmoothingLength | cellsize | Use as effective size |
| GFM_Metals | chem[] | If available (Illustris-style) |

#### Black Holes (PartType5 → SinkType)

| GADGET Field | NewDD Field | Notes |
|--------------|-------------|-------|
| Coordinates | x, y, z | |
| Velocities | vx, vy, vz | |
| BH_Mass | mass | True BH mass |
| ParticleIDs | id | |
| BH_Mdot | dMsmbh | Accretion rate |
| - | tbirth | May need to track separately |

---

## Required Modifications

### Files to Modify

#### 1. `Makefile` - Add HDF5 Support

```makefile
# Add these lines
HDF5_INC = -I/path/to/hdf5/include
HDF5_LIB = -L/path/to/hdf5/lib -lhdf5

# Modify OPT to add GADGET flag
OPT = -O3 -DGADGET_HDF5 -DNMEG=20000 ...

# Add HDF5 to includes and libs
INCLUDES = -I./ -I../ $(HDF5_INC)
OSTLIB = -L. -lmyram -lm $(HDF5_LIB)

# Add new object files
TIMEROBJ = ... rd_gadget.o rd_gadget_info.o
```

#### 2. `ramses.h` - Add GADGET-Specific Definitions

Location: After existing structure definitions (~line 227)

```c
#ifdef GADGET_HDF5
/* GADGET-specific structure for header info */
typedef struct GadgetHeaderType {
    int npart[6];           /* Particle count per type */
    double mass[6];         /* Mass table */
    double time;            /* Scale factor or time */
    double redshift;
    double boxsize;
    int nfiles;             /* Files per snapshot */
    double omega0;
    double omegalambda;
    double hubble;
    double unit_length;     /* cm */
    double unit_mass;       /* g */
    double unit_velocity;   /* cm/s */
    int flag_sfr;
    int flag_cooling;
    int flag_stellarage;
    int flag_metals;
} GadgetHeaderType;

/* Function prototypes for GADGET reading */
int rd_gadget_info(RamsesType *ram, char *basename, int ifile);
int rd_gadget_particles(RamsesType *ram, char *basename, int ifile);
int rd_gadget_gas(RamsesType *ram, char *basename, int ifile);
int rd_gadget_bh(RamsesType *ram, char *basename, int ifile);
void gadget_units(RamsesType *ram, GadgetHeaderType *header);
#endif
```

#### 3. `newdd.c` - Add GADGET Main Loop

Location: Replace or add conditional block around line 91-213

The main loop structure needs to change from:
- RAMSES: Loop over CPU domains (icpu = 1 to ncpu)
- GADGET: Loop over snapshot files (ifile = 0 to nfiles-1)

```c
#ifdef GADGET_HDF5
    /* GADGET workflow */
    sprintf(infile, "./snapdir_%03d/snapshot_%03d", istep, istep);
    rd_gadget_info(&ram, infile, 0);

    for(ifile = 0; ifile < ram.nfiles; ifile++) {
        rd_gadget_particles(&ram, infile, ifile);
        rd_gadget_gas(&ram, infile, ifile);
        rd_gadget_bh(&ram, infile, ifile);

        /* Sort and decompose (same as RAMSES) */
        qsort(ram.gas, ram.ngas, sizeof(GasType), gassortx);
        SplitDump(&ram, ram.gas, ram.ngas, GAS, istep, ifile, sinmul, nsplit);
        /* ... etc ... */
    }
#else
    /* Original RAMSES workflow */
    for(icpu = 1; icpu <= ram.ncpu; icpu++) { ... }
#endif
```

#### 4. `rd_info.c` - Add GADGET Metadata Reader Alternative

The `units()` function needs modification to handle GADGET unit system.

---

## New Files to Create

### 1. `rd_gadget.c` - Main GADGET Reader

```c
/*
 * rd_gadget.c - Read GADGET HDF5 snapshot files
 *
 * This file contains functions to read GADGET-format HDF5 simulation data.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "ramses.h"
#include "Memory.h"

#ifdef GADGET_HDF5

/*
 * rd_gadget_info - Read GADGET header/metadata
 *
 * Parameters:
 *   ram      - RamsesType structure to populate
 *   basename - Base path (e.g., "./snapdir_000/snapshot_000")
 *   ifile    - File index (0 for single file, 0-N for multi-file)
 *
 * Returns:
 *   0 on success, -1 on failure
 */
int rd_gadget_info(RamsesType *ram, char *basename, int ifile) {
    char filename[256];
    hid_t file_id, header_id;
    GadgetHeaderType header;

    /* Construct filename */
    if (/* check if multi-file */) {
        sprintf(filename, "%s.%d.hdf5", basename, ifile);
    } else {
        sprintf(filename, "%s.hdf5", basename);
    }

    /* Open HDF5 file */
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        printf("Error: Cannot open %s\n", filename);
        return -1;
    }

    /* Open Header group */
    header_id = H5Gopen2(file_id, "Header", H5P_DEFAULT);

    /* Read header attributes */
    /* NumPart_ThisFile */
    hid_t attr = H5Aopen(header_id, "NumPart_ThisFile", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, header.npart);
    H5Aclose(attr);

    /* MassTable */
    attr = H5Aopen(header_id, "MassTable", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, header.mass);
    H5Aclose(attr);

    /* Time (scale factor) */
    attr = H5Aopen(header_id, "Time", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(attr);

    /* BoxSize */
    attr = H5Aopen(header_id, "BoxSize", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &header.boxsize);
    H5Aclose(attr);

    /* Cosmological parameters */
    attr = H5Aopen(header_id, "Omega0", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &header.omega0);
    H5Aclose(attr);

    attr = H5Aopen(header_id, "OmegaLambda", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &header.omegalambda);
    H5Aclose(attr);

    attr = H5Aopen(header_id, "HubbleParam", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &header.hubble);
    H5Aclose(attr);

    /* Unit system (if available) */
    if (H5Aexists(header_id, "UnitLength_in_cm")) {
        attr = H5Aopen(header_id, "UnitLength_in_cm", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_DOUBLE, &header.unit_length);
        H5Aclose(attr);

        attr = H5Aopen(header_id, "UnitMass_in_g", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_DOUBLE, &header.unit_mass);
        H5Aclose(attr);

        attr = H5Aopen(header_id, "UnitVelocity_in_cm_per_s", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_DOUBLE, &header.unit_velocity);
        H5Aclose(attr);
    } else {
        /* Default GADGET units */
        header.unit_length = 3.085678e21;    /* 1 kpc in cm */
        header.unit_mass = 1.989e43;         /* 10^10 Msun in g */
        header.unit_velocity = 1.0e5;        /* 1 km/s in cm/s */
    }

    H5Gclose(header_id);
    H5Fclose(file_id);

    /* Map to RamsesType structure */
    ram->aexp = header.time;
    ram->omega_m = header.omega0;
    ram->omega_l = header.omegalambda;
    ram->H0 = header.hubble * 100.0;
    ram->boxlen_ini = header.boxsize * header.unit_length / Mpc
                      * header.hubble;  /* Convert to cMpc/h */
    ram->unit_l = header.unit_length;
    ram->unit_d = header.unit_mass / (header.unit_length * header.unit_length
                  * header.unit_length);
    ram->unit_t = header.unit_length / header.unit_velocity;
    ram->nfiles = header.nfiles;
    ram->cosmo = 1;  /* Assume cosmological */

    /* Compute derived units */
    gadget_units(ram, &header);

    return 0;
}

/*
 * rd_gadget_particles - Read DM and star particles
 */
int rd_gadget_particles(RamsesType *ram, char *basename, int ifile) {
    char filename[256];
    hid_t file_id, group_id, dataset_id, dataspace_id;
    hsize_t dims[2];
    int npart_dm, npart_star;

    /* Open file */
    sprintf(filename, "%s.%d.hdf5", basename, ifile);
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    /* Read PartType1 (Dark Matter) */
    if (H5Lexists(file_id, "PartType1", H5P_DEFAULT)) {
        group_id = H5Gopen2(file_id, "PartType1", H5P_DEFAULT);

        /* Get particle count */
        dataset_id = H5Dopen2(group_id, "Coordinates", H5P_DEFAULT);
        dataspace_id = H5Dget_space(dataset_id);
        H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
        npart_dm = dims[0];

        /* Allocate memory */
        double *coords = (double*)Malloc(sizeof(double) * npart_dm * 3, PPTR(coords));
        double *vels = (double*)Malloc(sizeof(double) * npart_dm * 3, PPTR(vels));
        double *masses = (double*)Malloc(sizeof(double) * npart_dm, PPTR(masses));
        long long *ids = (long long*)Malloc(sizeof(long long) * npart_dm, PPTR(ids));

        /* Read datasets */
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, coords);
        H5Dclose(dataset_id);

        dataset_id = H5Dopen2(group_id, "Velocities", H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, vels);
        H5Dclose(dataset_id);

        dataset_id = H5Dopen2(group_id, "ParticleIDs", H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, ids);
        H5Dclose(dataset_id);

        /* Check if masses are in dataset or header */
        if (H5Lexists(group_id, "Masses", H5P_DEFAULT)) {
            dataset_id = H5Dopen2(group_id, "Masses", H5P_DEFAULT);
            H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, masses);
            H5Dclose(dataset_id);
        } else {
            /* Use mass from header MassTable[1] */
            for (int i = 0; i < npart_dm; i++) {
                masses[i] = /* header.mass[1] */;
            }
        }

        H5Gclose(group_id);

        /* Convert to DmType and apply unit conversions */
        ram->particle = (PmType*)Malloc(sizeof(PmType) * npart_dm,
                                        PPTR(ram->particle));
        for (int i = 0; i < npart_dm; i++) {
            ram->particle[i].x = coords[i*3 + 0] * ram->mpcscale_l;
            ram->particle[i].y = coords[i*3 + 1] * ram->mpcscale_l;
            ram->particle[i].z = coords[i*3 + 2] * ram->mpcscale_l;
            ram->particle[i].vx = vels[i*3 + 0] * ram->kmscale_v;
            ram->particle[i].vy = vels[i*3 + 1] * ram->kmscale_v;
            ram->particle[i].vz = vels[i*3 + 2] * ram->kmscale_v;
            ram->particle[i].mass = masses[i] * ram->scale_m;
            ram->particle[i].id = ids[i];
            ram->particle[i].family = 1;  /* DM */
            ram->particle[i].levelp = 0;
        }
        ram->npart = npart_dm;

        Free(coords);
        Free(vels);
        Free(masses);
        Free(ids);
    }

    /* Read PartType4 (Stars) - similar pattern */
    /* ... */

    H5Fclose(file_id);
    return 0;
}

/*
 * rd_gadget_gas - Read gas particles (PartType0)
 */
int rd_gadget_gas(RamsesType *ram, char *basename, int ifile) {
    /* Similar structure to rd_gadget_particles */
    /* Read: Coordinates, Velocities, Masses, Density, InternalEnergy */
    /* Convert InternalEnergy to temperature */
    /* Map to GasType structure */
    return 0;
}

/*
 * rd_gadget_bh - Read black hole particles (PartType5)
 */
int rd_gadget_bh(RamsesType *ram, char *basename, int ifile) {
    /* Read: Coordinates, Velocities, BH_Mass, BH_Mdot, ParticleIDs */
    /* Map to SinkType structure */
    return 0;
}

/*
 * gadget_units - Compute unit conversion factors for GADGET
 */
void gadget_units(RamsesType *ram, GadgetHeaderType *header) {
    double h = header->hubble;
    double a = header->time;  /* scale factor */

    /* Length: internal units to cMpc/h */
    ram->mpcscale_l = header->unit_length / Mpc * h;

    /* Velocity: internal units to km/s */
    /* GADGET stores peculiar velocity * sqrt(a) */
    ram->kmscale_v = header->unit_velocity / 1.0e5 * sqrt(a);

    /* Mass: internal units to Msun/h */
    ram->scale_m = header->unit_mass / Msun * h;

    /* Density */
    ram->scale_d = header->unit_mass /
                   (header->unit_length * header->unit_length * header->unit_length);

    /* Temperature: internal energy to Kelvin */
    /* T = (gamma - 1) * u * mu * mH / kB */
    double gamma = 5.0 / 3.0;
    double mu = 0.588;  /* Mean molecular weight (fully ionized) */
    ram->scale_T2 = (gamma - 1.0) * header->unit_velocity * header->unit_velocity
                    * mu * mH / kB;
}

#endif /* GADGET_HDF5 */
```

### 2. `rd_gadget.h` - Header for GADGET Functions

```c
#ifndef RD_GADGET_H
#define RD_GADGET_H

#ifdef GADGET_HDF5
#include <hdf5.h>
#include "ramses.h"

/* GADGET particle type indices */
#define GADGET_TYPE_GAS   0
#define GADGET_TYPE_HALO  1
#define GADGET_TYPE_DISK  2
#define GADGET_TYPE_BULGE 3
#define GADGET_TYPE_STAR  4
#define GADGET_TYPE_BH    5

/* Function prototypes */
int rd_gadget_info(RamsesType *ram, char *basename, int ifile);
int rd_gadget_particles(RamsesType *ram, char *basename, int ifile);
int rd_gadget_gas(RamsesType *ram, char *basename, int ifile);
int rd_gadget_bh(RamsesType *ram, char *basename, int ifile);
void gadget_units(RamsesType *ram, GadgetHeaderType *header);

/* Helper functions */
int gadget_get_nfiles(char *basename);
int gadget_check_dataset(hid_t group_id, const char *name);

#endif /* GADGET_HDF5 */
#endif /* RD_GADGET_H */
```

---

## Makefile Changes

### Complete Modified Makefile Section

```makefile
# Compiler
CC = mpicc

# HDF5 paths (adjust for your system)
HDF5_DIR = /usr/local/hdf5
# Or use: HDF5_DIR = $(shell h5cc -showconfig | grep "Installation point" | cut -d: -f2 | tr -d ' ')

HDF5_INC = -I$(HDF5_DIR)/include
HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5 -Wl,-rpath,$(HDF5_DIR)/lib

# Choose format: RAMSES or GADGET
# Uncomment ONE of the following:

# For RAMSES format (default):
# FORMAT_OPT =

# For GADGET HDF5 format:
FORMAT_OPT = -DGADGET_HDF5

# Optimization and flags
OPT = -O3 $(FORMAT_OPT) -DNMEG=20000 -DWGROUPSIZE=5 -DUSE_MPI -DNCHEM=9

# Includes
INCLUDES = -I./ -I../ $(HDF5_INC)

# Libraries
OSTLIB = -L. -lmyram -lm $(HDF5_LIB)

# Object files
ifdef GADGET_HDF5
TIMEROBJ = rd_gadget.o utils.o Memory2.o header.o
else
TIMEROBJ = rd_amr.o rd_info.o rd_part.o header.o utils.o Memory2.o \
           find_leaf_gas.o rd_hydro.o rd_sink.o
endif

# ... rest of Makefile remains the same ...
```

### Alternative: Conditional Compilation

```makefile
# Both RAMSES and GADGET in same binary
TIMEROBJ = rd_amr.o rd_info.o rd_part.o header.o utils.o Memory2.o \
           find_leaf_gas.o rd_hydro.o rd_sink.o rd_gadget.o

# Use -DGADGET_HDF5 flag to enable GADGET support
# At runtime, specify input format via command line
```

---

## Implementation Details

### Key Differences to Handle

#### 1. No AMR in GADGET

RAMSES uses AMR with hierarchical cells. GADGET uses SPH particles.

**Solution for GasType**:
- Use `SmoothingLength` as `cellsize`
- Density comes directly from SPH calculation
- No need for `find_leaf_gas()` - all particles are "leaves"

```c
/* In rd_gadget_gas() */
gas[i].cellsize = smoothing_length[i] * ram->mpcscale_l;
gas[i].den = density[i];  /* Already in physical units */
```

#### 2. Velocity Convention

GADGET stores peculiar velocity × √a, while RAMSES stores comoving velocity.

```c
/* Convert GADGET velocity to physical km/s */
double a = ram->aexp;
vel_physical = vel_gadget * sqrt(a) * unit_velocity / 1e5;
```

#### 3. File Structure

RAMSES: One file per CPU per data type
GADGET: One or more files containing all particle types

```c
/* GADGET file loop */
for (ifile = 0; ifile < nfiles; ifile++) {
    /* Each file contains subset of all particle types */
    rd_gadget_particles(&ram, basename, ifile);
    rd_gadget_gas(&ram, basename, ifile);
    rd_gadget_bh(&ram, basename, ifile);

    /* Process and write slabs */
    /* ... */
}
```

#### 4. Mass Table vs Individual Masses

GADGET may store uniform masses in header, not per-particle.

```c
/* Check if mass dataset exists */
if (H5Lexists(group_id, "Masses", H5P_DEFAULT)) {
    /* Read individual masses */
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, ...);
} else {
    /* Use mass from header MassTable */
    for (int i = 0; i < npart; i++) {
        masses[i] = header.mass[particle_type];
    }
}
```

---

## Unit Conversion

### GADGET Default Units

| Quantity | Default Unit | Value in CGS |
|----------|--------------|--------------|
| Length | kpc/h | 3.085678e21 / h cm |
| Mass | 10^10 Msun/h | 1.989e43 / h g |
| Velocity | km/s | 1.0e5 cm/s |
| Time | derived | L/V |
| Density | derived | M/L^3 |

### Conversion to NewDD Units (cMpc/h, Msun/h, km/s)

```c
/* Length: kpc/h -> cMpc/h */
x_cmpc = x_kpc / 1000.0;

/* Mass: 10^10 Msun/h -> Msun/h */
m_msun = m_gadget * 1.0e10;

/* Velocity: already km/s, but correct for sqrt(a) factor */
v_kms = v_gadget * sqrt(a);

/* Temperature from internal energy */
/* u is in (km/s)^2, convert to Kelvin */
double gamma = 5.0/3.0;
double mu = 4.0 / (1.0 + 3.0*XH + 4.0*XH*xe);  /* molecular weight */
T = (gamma - 1.0) * u * 1.0e10 * mu * mH / kB;
```

### Temperature Calculation

GADGET stores thermal energy per unit mass (`InternalEnergy`):

```c
/* Convert internal energy to temperature */
double u = internal_energy[i];  /* (km/s)^2 */
double ne = electron_abundance[i];  /* optional */

/* Mean molecular weight */
double XH = 0.76;  /* Hydrogen mass fraction */
double mu;
if (ne_available) {
    mu = 4.0 / (1.0 + 3.0*XH + 4.0*XH*ne);
} else {
    mu = 0.588;  /* Fully ionized primordial */
}

/* Temperature in Kelvin */
double gamma = 5.0/3.0;
T = (gamma - 1.0) * u * 1.0e10 * mu * mH / kB;
```

---

## Example Implementation

### Minimal Working Example

```c
/* newdd_gadget.c - Simplified main for GADGET */
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include "ramses.h"
#include "Memory.h"

int main(int argc, char **argv) {
    RamsesType ram;
    int istep, nsplit, ifile;
    char basename[256];

    istep = atoi(argv[1]);
    nsplit = atoi(argv[2]);

    Make_Total_Memory();

    /* Construct base filename */
    sprintf(basename, "./snapdir_%03d/snapshot_%03d", istep, istep);

    /* Read header from first file */
    rd_gadget_info(&ram, basename, 0);

    /* Process each file */
    for (ifile = 0; ifile < ram.nfiles; ifile++) {
        printf("Processing file %d of %d\n", ifile+1, ram.nfiles);

        /* Read particles */
        rd_gadget_particles(&ram, basename, ifile);

        /* Read gas */
        rd_gadget_gas(&ram, basename, ifile);

        /* Read black holes */
        rd_gadget_bh(&ram, basename, ifile);

        /* Separate DM and stars, sort, and dump */
        DmType *dm = (DmType*)Malloc(sizeof(DmType)*ram.npart, PPTR(dm));
        int ndm = 0, nstar = 0;

        for (int i = 0; i < ram.npart; i++) {
            if (ram.particle[i].family == 1) {
                dm[ndm++] = ((DmType*)ram.particle)[i];
            }
        }

        qsort(dm, ndm, sizeof(DmType), dmsortx);
        SplitDump(&ram, dm, ndm, DM, istep, ifile, 0, nsplit);

        /* Similar for stars, gas, sinks */
        /* ... */

        Free(dm);
        cleanup_ramses(&ram);
    }

    return 0;
}
```

---

## Testing

### 1. Compile with HDF5 Support

```bash
# Set HDF5 path
export HDF5_DIR=/path/to/hdf5

# Build
make clean
make all
```

### 2. Verify HDF5 Linking

```bash
ldd newdd.exe | grep hdf5
# Should show: libhdf5.so.XXX => /path/to/libhdf5.so.XXX
```

### 3. Test with Small Snapshot

```bash
# Single file test
./newdd.exe 0 4

# Check output
ls -la FoF_Data/NewDD.00000/
```

### 4. Validate Output

```python
import numpy as np

# Read output and verify
dm = np.fromfile('FoF_Data/NewDD.00000/SN.00000.DM.00000.dat',
                 dtype=dm_dtype)
print(f"DM particles: {len(dm)}")
print(f"Position range: {dm['x'].min():.2f} - {dm['x'].max():.2f}")
```

### 5. Compare with Original

If you have the same simulation in both RAMSES and GADGET formats, compare the output statistics (particle counts, mass distributions, etc.).

---

## Summary of Changes

| File | Modification Type | Description |
|------|-------------------|-------------|
| `Makefile` | Modify | Add HDF5 flags and libraries |
| `ramses.h` | Modify | Add GadgetHeaderType, function prototypes |
| `newdd.c` | Modify | Add GADGET main loop (conditional) |
| `rd_gadget.c` | **Create** | GADGET HDF5 reading functions |
| `rd_gadget.h` | **Create** | GADGET function headers |

### Checklist

- [ ] Install HDF5 library with C bindings
- [ ] Modify Makefile with HDF5 paths
- [ ] Create rd_gadget.c with reader functions
- [ ] Add GadgetHeaderType to ramses.h
- [ ] Modify newdd.c with conditional GADGET workflow
- [ ] Test with sample GADGET snapshot
- [ ] Verify unit conversions
- [ ] Validate output format compatibility

---

## References

- [GADGET-2 User Guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf)
- [HDF5 C API Documentation](https://docs.hdfgroup.org/hdf5/v1_14/group___h5.html)
- [Illustris Data Specifications](https://www.illustris-project.org/data/docs/specifications/)
- [AREPO/GADGET-4 Documentation](https://arepo-code.org/documentation)

---

## Contact

For questions about this guide or implementation assistance, please contact the code maintainer.
