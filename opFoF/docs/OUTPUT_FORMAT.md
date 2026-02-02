# opFoF Output Format Specification

This document provides a detailed specification of the binary output files produced by opFoF, including byte-level layouts, data types, and reading examples in multiple languages.

---

## Table of Contents

1. [Overview](#overview)
2. [Halo Catalog Format](#halo-catalog-format)
3. [Member Particle Format](#member-particle-format)
4. [Data Type Definitions](#data-type-definitions)
5. [Byte Order and Alignment](#byte-order-and-alignment)
6. [Reading Examples](#reading-examples)
7. [Format Validation](#format-validation)
8. [Version History](#version-history)

---

## Overview

opFoF produces two output files per snapshot:

| File | Description |
|------|-------------|
| `FoF_halo_cat.<NNNNN>` | Halo properties catalog |
| `FoF_member_particle.<NNNNN>` | Member particle data |

Where `<NNNNN>` is the 5-digit zero-padded snapshot number.

### File Size Estimation

```
Catalog size ≈ 28 bytes (header) + 128 bytes × N_halos

Member file size ≈ 28 bytes (header) + Σ(particle_size × N_particles_per_halo)
```

---

## Halo Catalog Format

### File: `FoF_halo_cat.<NNNNN>`

```
+------------------+
|     HEADER       |  28 bytes
+------------------+
|     HALO 0       |  128 bytes
+------------------+
|     HALO 1       |  128 bytes
+------------------+
|       ...        |
+------------------+
|     HALO N-1     |  128 bytes
+------------------+
```

### Header Section (28 bytes)

| Offset | Size | Type | Name | Description | Units |
|--------|------|------|------|-------------|-------|
| 0 | 4 | float | size | Simulation box size | Mpc/h |
| 4 | 4 | float | hubble | Hubble parameter (H₀/100) | dimensionless |
| 8 | 4 | float | omep | Ωₘ (matter density) | dimensionless |
| 12 | 4 | float | omepb | Ωᵦ (baryon density) | dimensionless |
| 16 | 4 | float | omeplam | Ωₗ (dark energy density) | dimensionless |
| 20 | 4 | float | amax | Maximum scale factor | dimensionless |
| 24 | 4 | float | anow | Current scale factor | dimensionless |

**Total header size: 28 bytes**

### Halo Record (HaloQ Structure, 128 bytes)

| Offset | Size | Type | Name | Description | Units |
|--------|------|------|------|-------------|-------|
| 0 | 8 | size_t | np | Total particle count | count |
| 8 | 8 | size_t | npstar | Star particle count | count |
| 16 | 8 | size_t | npgas | Gas particle count | count |
| 24 | 8 | size_t | npdm | Dark matter particle count | count |
| 32 | 8 | size_t | npsink | Sink particle count | count |
| 40 | 8 | double | x | X position (COM) | Mpc/h |
| 48 | 8 | double | y | Y position (COM) | Mpc/h |
| 56 | 8 | double | z | Z position (COM) | Mpc/h |
| 64 | 8 | double | mass | Total mass | M☉/h |
| 72 | 8 | double | mstar | Stellar mass | M☉/h |
| 80 | 8 | double | mgas | Gas mass | M☉/h |
| 88 | 8 | double | mdm | Dark matter mass | M☉/h |
| 96 | 8 | double | msink | Sink mass | M☉/h |
| 104 | 4 | float | vx | X velocity (COM) | km/s |
| 108 | 4 | float | vy | Y velocity (COM) | km/s |
| 112 | 4 | float | vz | Z velocity (COM) | km/s |
| 116 | 12 | - | padding | Alignment padding | - |

**Total halo record size: 128 bytes**

**Note:** If compiled without `-DXYZDBL`, positions are `float` (4 bytes each), and the structure size changes accordingly.

### Catalog Binary Layout (Hex Dump)

```
Offset    Content
00000000  [Header: 28 bytes]
          XX XX XX XX    size (float)
          XX XX XX XX    hubble (float)
          XX XX XX XX    omep (float)
          XX XX XX XX    omepb (float)
          XX XX XX XX    omeplam (float)
          XX XX XX XX    amax (float)
          XX XX XX XX    anow (float)

0000001C  [Halo 0: 128 bytes]
          XX XX XX XX XX XX XX XX    np (size_t)
          XX XX XX XX XX XX XX XX    npstar (size_t)
          ...

0000009C  [Halo 1: 128 bytes]
          ...
```

---

## Member Particle Format

### File: `FoF_member_particle.<NNNNN>`

```
+------------------+
|     HEADER       |  28 bytes (same as catalog)
+------------------+
| HALO 0 PARTICLES |
|   DM particles   |
|   Gas particles  |
|   Sink particles |
|   Star particles |
+------------------+
| HALO 1 PARTICLES |
|   DM particles   |
|   ...            |
+------------------+
|       ...        |
+------------------+
```

### Particle Record Sizes

| Type | Structure | Size |
|------|-----------|------|
| Dark Matter | DmType | 28 bytes |
| Star | StarType | 32 bytes |
| Gas | GasType | 24 bytes |
| Sink | SinkType | 32 bytes |

### Dark Matter Particle (DmType, 28 bytes)

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 0 | 4 | float | vx | Velocity X (km/s) |
| 4 | 4 | float | vy | Velocity Y (km/s) |
| 8 | 4 | float | vz | Velocity Z (km/s) |
| 12 | 4 | float | mass | Mass (M☉/h) |
| 16 | 8 | long long | id | Particle ID |
| 24 | 4 | int | level | AMR level |

### Star Particle (StarType, 32 bytes)

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 0 | 4 | float | vx | Velocity X (km/s) |
| 4 | 4 | float | vy | Velocity Y (km/s) |
| 8 | 4 | float | vz | Velocity Z (km/s) |
| 12 | 4 | float | mass | Mass (M☉/h) |
| 16 | 8 | long long | id | Particle ID |
| 24 | 4 | float | age | Stellar age (Gyr) |
| 28 | 4 | float | metallicity | Z (metal fraction) |

### Gas Particle (GasType, 24 bytes)

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 0 | 4 | float | vx | Velocity X (km/s) |
| 4 | 4 | float | vy | Velocity Y (km/s) |
| 8 | 4 | float | vz | Velocity Z (km/s) |
| 12 | 4 | float | mass | Mass (M☉/h) |
| 16 | 4 | float | density | Density |
| 20 | 4 | float | temperature | Temperature (K) |

### Sink Particle (SinkType, 32 bytes)

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 0 | 4 | float | vx | Velocity X (km/s) |
| 4 | 4 | float | vy | Velocity Y (km/s) |
| 8 | 4 | float | vz | Velocity Z (km/s) |
| 12 | 4 | float | mass | Mass (M☉/h) |
| 16 | 4 | float | birth_time | Formation scale factor |
| 20 | 4 | float | lx | Angular momentum X |
| 24 | 4 | float | ly | Angular momentum Y |
| 28 | 4 | float | lz | Angular momentum Z |

### Particle Order Within Each Halo

For each halo, particles are written in this order:
1. All DM particles (npdm × 28 bytes)
2. All Gas particles (npgas × 24 bytes)
3. All Sink particles (npsink × 32 bytes)
4. All Star particles (npstar × 32 bytes)

---

## Data Type Definitions

### C Type Mapping

```c
// Standard types
typedef float float32;           // 4 bytes
typedef double float64;          // 8 bytes
typedef int int32;               // 4 bytes
typedef long long int64;         // 8 bytes
typedef unsigned long size_t;    // 8 bytes (64-bit system)

// Position type (depends on compilation)
#ifdef XYZDBL
typedef double POSTYPE;          // 8 bytes
#else
typedef float POSTYPE;           // 4 bytes
#endif
```

### NumPy dtype Specification

```python
# Header dtype
header_dtype = np.dtype([
    ('size', 'f4'),
    ('hubble', 'f4'),
    ('omep', 'f4'),
    ('omepb', 'f4'),
    ('omeplam', 'f4'),
    ('amax', 'f4'),
    ('anow', 'f4')
])

# HaloQ dtype (with XYZDBL)
haloq_dtype = np.dtype([
    ('np', 'u8'),
    ('npstar', 'u8'),
    ('npgas', 'u8'),
    ('npdm', 'u8'),
    ('npsink', 'u8'),
    ('x', 'f8'),
    ('y', 'f8'),
    ('z', 'f8'),
    ('mass', 'f8'),
    ('mstar', 'f8'),
    ('mgas', 'f8'),
    ('mdm', 'f8'),
    ('msink', 'f8'),
    ('vx', 'f4'),
    ('vy', 'f4'),
    ('vz', 'f4'),
    ('_pad', 'V12')  # 12 bytes padding
])

# Particle dtypes
dm_dtype = np.dtype([
    ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
    ('mass', 'f4'), ('id', 'i8'), ('level', 'i4')
])

star_dtype = np.dtype([
    ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
    ('mass', 'f4'), ('id', 'i8'),
    ('age', 'f4'), ('metallicity', 'f4')
])

gas_dtype = np.dtype([
    ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
    ('mass', 'f4'), ('density', 'f4'), ('temperature', 'f4')
])

sink_dtype = np.dtype([
    ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
    ('mass', 'f4'), ('birth_time', 'f4'),
    ('lx', 'f4'), ('ly', 'f4'), ('lz', 'f4')
])
```

---

## Byte Order and Alignment

### Endianness

opFoF writes files in **native byte order** (typically little-endian on x86 systems).

To determine endianness:
```python
import sys
print(sys.byteorder)  # 'little' or 'big'
```

To read on different-endian system:
```python
# Swap byte order if needed
haloq_dtype_swapped = haloq_dtype.newbyteorder('>')  # big-endian
```

### Structure Alignment

All structures are naturally aligned:
- 4-byte types aligned to 4-byte boundaries
- 8-byte types aligned to 8-byte boundaries

Padding is added as needed to maintain alignment.

---

## Reading Examples

### C Example: Complete Reader

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Type definitions
typedef unsigned long size_t;
typedef double POSTYPE;  // Assume XYZDBL

typedef struct {
    float size, hubble, omep, omepb, omeplam, amax, anow;
} Header;

typedef struct {
    size_t np, npstar, npgas, npdm, npsink;
    POSTYPE x, y, z;
    double mass, mstar, mgas, mdm, msink;
    float vx, vy, vz;
    char padding[12];
} HaloQ;

typedef struct {
    float vx, vy, vz, mass;
    long long id;
    int level;
} DmType;

typedef struct {
    float vx, vy, vz, mass;
    long long id;
    float age, metallicity;
} StarType;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <catalog_file>\n", argv[0]);
        return 1;
    }

    FILE *fp = fopen(argv[1], "rb");
    if (!fp) {
        perror("Cannot open file");
        return 1;
    }

    // Read header
    Header header;
    fread(&header, sizeof(Header), 1, fp);

    printf("=== Header ===\n");
    printf("Box size: %.2f Mpc/h\n", header.size);
    printf("Redshift: %.4f\n", 1.0/header.anow - 1.0);
    printf("Omega_m: %.4f\n", header.omep);
    printf("\n");

    // Count halos
    long pos_start = ftell(fp);
    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);
    int nhalo = (file_size - sizeof(Header)) / sizeof(HaloQ);
    fseek(fp, pos_start, SEEK_SET);

    printf("=== Halos: %d ===\n", nhalo);

    // Read halos
    HaloQ *halos = malloc(nhalo * sizeof(HaloQ));
    fread(halos, sizeof(HaloQ), nhalo, fp);
    fclose(fp);

    // Print statistics
    double total_mass = 0;
    size_t total_particles = 0;
    for (int i = 0; i < nhalo; i++) {
        total_mass += halos[i].mass;
        total_particles += halos[i].np;
    }

    printf("Total halos: %d\n", nhalo);
    printf("Total mass: %.4e Msun/h\n", total_mass);
    printf("Total particles: %zu\n", total_particles);

    // Print top 10 most massive
    printf("\n=== Top 10 Massive Halos ===\n");
    printf("%6s %12s %14s %12s %12s %12s\n",
           "Rank", "Np", "Mass", "X", "Y", "Z");

    // Simple bubble sort for top 10
    for (int i = 0; i < nhalo && i < 10; i++) {
        int max_idx = i;
        for (int j = i+1; j < nhalo; j++) {
            if (halos[j].mass > halos[max_idx].mass) {
                max_idx = j;
            }
        }
        if (max_idx != i) {
            HaloQ tmp = halos[i];
            halos[i] = halos[max_idx];
            halos[max_idx] = tmp;
        }
        printf("%6d %12zu %14.4e %12.4f %12.4f %12.4f\n",
               i+1, halos[i].np, halos[i].mass,
               halos[i].x, halos[i].y, halos[i].z);
    }

    free(halos);
    return 0;
}
```

### Python Example: Complete Reader

```python
#!/usr/bin/env python3
"""
fof_reader.py - Complete opFoF output reader
"""
import numpy as np
import struct
from pathlib import Path

class FoFReader:
    """Reader for opFoF output files."""

    # Data types (adjust POSTYPE based on compilation)
    POSTYPE = 'f8'  # 'f8' for double (XYZDBL), 'f4' for float

    def __init__(self, snapshot, data_dir='.'):
        self.snapshot = snapshot
        self.data_dir = Path(data_dir)
        self.catalog_file = self.data_dir / f"FoF_halo_cat.{snapshot:05d}"
        self.member_file = self.data_dir / f"FoF_member_particle.{snapshot:05d}"

        # Define dtypes
        self._define_dtypes()

    def _define_dtypes(self):
        """Define numpy dtypes for binary structures."""
        self.header_dtype = np.dtype([
            ('size', 'f4'),
            ('hubble', 'f4'),
            ('omep', 'f4'),
            ('omepb', 'f4'),
            ('omeplam', 'f4'),
            ('amax', 'f4'),
            ('anow', 'f4')
        ])

        self.haloq_dtype = np.dtype([
            ('np', 'u8'),
            ('npstar', 'u8'),
            ('npgas', 'u8'),
            ('npdm', 'u8'),
            ('npsink', 'u8'),
            ('x', self.POSTYPE),
            ('y', self.POSTYPE),
            ('z', self.POSTYPE),
            ('mass', 'f8'),
            ('mstar', 'f8'),
            ('mgas', 'f8'),
            ('mdm', 'f8'),
            ('msink', 'f8'),
            ('vx', 'f4'),
            ('vy', 'f4'),
            ('vz', 'f4'),
            ('_pad', 'V12')
        ])

        self.dm_dtype = np.dtype([
            ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
            ('mass', 'f4'), ('id', 'i8'), ('level', 'i4')
        ])

        self.star_dtype = np.dtype([
            ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
            ('mass', 'f4'), ('id', 'i8'),
            ('age', 'f4'), ('metallicity', 'f4')
        ])

        self.gas_dtype = np.dtype([
            ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
            ('mass', 'f4'), ('density', 'f4'), ('temperature', 'f4')
        ])

        self.sink_dtype = np.dtype([
            ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
            ('mass', 'f4'), ('birth_time', 'f4'),
            ('lx', 'f4'), ('ly', 'f4'), ('lz', 'f4')
        ])

    def read_catalog(self, min_particles=0):
        """
        Read halo catalog.

        Parameters
        ----------
        min_particles : int
            Minimum particle count to include

        Returns
        -------
        header : dict
            Simulation parameters
        halos : ndarray
            Structured array of halo properties
        """
        with open(self.catalog_file, 'rb') as f:
            # Read header
            header_raw = np.fromfile(f, dtype=self.header_dtype, count=1)[0]
            header = {
                'box_size': float(header_raw['size']),
                'hubble': float(header_raw['hubble']),
                'omega_m': float(header_raw['omep']),
                'omega_b': float(header_raw['omepb']),
                'omega_l': float(header_raw['omeplam']),
                'amax': float(header_raw['amax']),
                'anow': float(header_raw['anow']),
            }
            header['redshift'] = 1.0 / header['anow'] - 1.0

            # Read halos
            halos = np.fromfile(f, dtype=self.haloq_dtype)

        # Filter by particle count
        if min_particles > 0:
            halos = halos[halos['np'] >= min_particles]

        return header, halos

    def read_halo_particles(self, halo_index):
        """
        Read particles for a specific halo.

        Parameters
        ----------
        halo_index : int
            0-based index of halo in catalog

        Returns
        -------
        dict : Particle arrays by type
        """
        # First read catalog to get structure
        header, halos = self.read_catalog()

        if halo_index >= len(halos):
            raise ValueError(f"Halo index {halo_index} out of range")

        # Calculate byte offset to this halo's particles
        offset = self.header_dtype.itemsize
        for i in range(halo_index):
            h = halos[i]
            offset += h['npdm'] * self.dm_dtype.itemsize
            offset += h['npgas'] * self.gas_dtype.itemsize
            offset += h['npsink'] * self.sink_dtype.itemsize
            offset += h['npstar'] * self.star_dtype.itemsize

        target = halos[halo_index]
        particles = {}

        with open(self.member_file, 'rb') as f:
            f.seek(offset)

            if target['npdm'] > 0:
                particles['dm'] = np.fromfile(
                    f, dtype=self.dm_dtype, count=target['npdm'])

            if target['npgas'] > 0:
                particles['gas'] = np.fromfile(
                    f, dtype=self.gas_dtype, count=target['npgas'])

            if target['npsink'] > 0:
                particles['sink'] = np.fromfile(
                    f, dtype=self.sink_dtype, count=target['npsink'])

            if target['npstar'] > 0:
                particles['star'] = np.fromfile(
                    f, dtype=self.star_dtype, count=target['npstar'])

        return particles

    def summary(self):
        """Print summary of catalog."""
        header, halos = self.read_catalog(min_particles=30)

        print("=" * 60)
        print(f"opFoF Catalog: {self.catalog_file.name}")
        print("=" * 60)
        print(f"\nSimulation Parameters:")
        print(f"  Box size:     {header['box_size']:.2f} Mpc/h")
        print(f"  Redshift:     {header['redshift']:.4f}")
        print(f"  Omega_m:      {header['omega_m']:.4f}")
        print(f"  Omega_b:      {header['omega_b']:.4f}")
        print(f"  Omega_L:      {header['omega_l']:.4f}")

        print(f"\nHalo Statistics:")
        print(f"  Total halos:  {len(halos)}")
        print(f"  Total mass:   {halos['mass'].sum():.4e} Msun/h")
        print(f"  Max mass:     {halos['mass'].max():.4e} Msun/h")
        print(f"  Min mass:     {halos['mass'].min():.4e} Msun/h")

        print(f"\nParticle Statistics:")
        print(f"  Total particles: {halos['np'].sum():,}")
        print(f"  DM:              {halos['npdm'].sum():,}")
        print(f"  Stars:           {halos['npstar'].sum():,}")
        print(f"  Gas:             {halos['npgas'].sum():,}")
        print(f"  Sinks:           {halos['npsink'].sum():,}")


# Example usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python fof_reader.py <snapshot_number> [data_dir]")
        sys.exit(1)

    snap = int(sys.argv[1])
    data_dir = sys.argv[2] if len(sys.argv) > 2 else '.'

    reader = FoFReader(snap, data_dir)
    reader.summary()
```

### Fortran Example

```fortran
program read_fof_catalog
    implicit none

    ! Type definitions
    integer, parameter :: sp = 4, dp = 8

    type :: Header
        real(sp) :: size, hubble, omep, omepb, omeplam, amax, anow
    end type Header

    type :: HaloQ
        integer(8) :: np, npstar, npgas, npdm, npsink
        real(dp) :: x, y, z
        real(dp) :: mass, mstar, mgas, mdm, msink
        real(sp) :: vx, vy, vz
        character(12) :: padding
    end type HaloQ

    type(Header) :: hdr
    type(HaloQ), allocatable :: halos(:)
    integer :: unit_num, nhalo, i, ios
    integer(8) :: file_size
    character(256) :: filename

    ! Get filename from command line
    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) then
        print *, "Usage: read_fof_catalog <catalog_file>"
        stop 1
    endif

    ! Open file
    unit_num = 10
    open(unit=unit_num, file=trim(filename), status='old', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file: ", trim(filename)
        stop 1
    endif

    ! Read header
    read(unit_num) hdr

    print *, "=== Header ==="
    print '(A,F10.2)', " Box size: ", hdr%size
    print '(A,F10.4)', " Redshift: ", 1.0/hdr%anow - 1.0
    print '(A,F10.4)', " Omega_m:  ", hdr%omep

    ! Calculate number of halos
    inquire(unit=unit_num, size=file_size)
    nhalo = (file_size - 28) / 128

    print '(A,I10)', " Number of halos: ", nhalo

    ! Allocate and read halos
    allocate(halos(nhalo))
    read(unit_num) halos

    close(unit_num)

    ! Print statistics
    print *, ""
    print *, "=== Top 10 Massive Halos ==="
    do i = 1, min(10, nhalo)
        print '(I6, I12, ES14.4, 3F12.4)', &
            i, halos(i)%np, halos(i)%mass, &
            halos(i)%x, halos(i)%y, halos(i)%z
    enddo

    deallocate(halos)

end program read_fof_catalog
```

---

## Format Validation

### Validation Script

```python
#!/usr/bin/env python3
"""
validate_fof.py - Validate opFoF output files
"""
import numpy as np
import sys

def validate_catalog(filename):
    """Validate catalog file structure and values."""
    errors = []
    warnings = []

    try:
        with open(filename, 'rb') as f:
            # Check header
            header = np.fromfile(f, dtype='f4', count=7)

            if header[0] <= 0:
                errors.append(f"Invalid box size: {header[0]}")
            if header[6] <= 0 or header[6] > 2:
                errors.append(f"Invalid scale factor: {header[6]}")
            if header[2] < 0 or header[2] > 1:
                warnings.append(f"Unusual Omega_m: {header[2]}")

            # Check halos
            halo_dtype = np.dtype([
                ('np', 'u8'), ('npstar', 'u8'), ('npgas', 'u8'),
                ('npdm', 'u8'), ('npsink', 'u8'),
                ('x', 'f8'), ('y', 'f8'), ('z', 'f8'),
                ('mass', 'f8'), ('mstar', 'f8'), ('mgas', 'f8'),
                ('mdm', 'f8'), ('msink', 'f8'),
                ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
                ('_pad', 'V12')
            ])

            halos = np.fromfile(f, dtype=halo_dtype)

            # Validate each halo
            box = header[0]
            for i, h in enumerate(halos):
                # Check particle count consistency
                n_sum = h['npdm'] + h['npstar'] + h['npgas'] + h['npsink']
                if h['np'] != n_sum:
                    errors.append(f"Halo {i}: particle count mismatch "
                                  f"({h['np']} != {n_sum})")

                # Check position within box
                if h['x'] < 0 or h['x'] >= box:
                    errors.append(f"Halo {i}: x position out of bounds")
                if h['y'] < 0 or h['y'] >= box:
                    errors.append(f"Halo {i}: y position out of bounds")
                if h['z'] < 0 or h['z'] >= box:
                    errors.append(f"Halo {i}: z position out of bounds")

                # Check mass consistency
                m_sum = h['mdm'] + h['mstar'] + h['mgas'] + h['msink']
                if abs(h['mass'] - m_sum) / max(h['mass'], 1e-10) > 0.01:
                    warnings.append(f"Halo {i}: mass mismatch "
                                    f"({h['mass']:.3e} != {m_sum:.3e})")

                # Check for negative masses
                if h['mass'] <= 0:
                    errors.append(f"Halo {i}: non-positive mass")

    except Exception as e:
        errors.append(f"Failed to read file: {e}")

    # Print results
    print(f"\n=== Validation Results for {filename} ===\n")

    if errors:
        print(f"ERRORS ({len(errors)}):")
        for e in errors[:10]:  # Limit output
            print(f"  - {e}")
        if len(errors) > 10:
            print(f"  ... and {len(errors)-10} more errors")
    else:
        print("No errors found.")

    if warnings:
        print(f"\nWARNINGS ({len(warnings)}):")
        for w in warnings[:10]:
            print(f"  - {w}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings)-10} more warnings")

    return len(errors) == 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python validate_fof.py <catalog_file>")
        sys.exit(1)

    success = validate_catalog(sys.argv[1])
    sys.exit(0 if success else 1)
```

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial format specification |
| 1.1 | 2024 | Added sink particle support |
| 1.2 | 2025 | Changed particle counts to size_t |
| 1.3 | 2025 | Added padding for alignment |

---

## Notes

1. **Endianness**: Files are written in native byte order. Check and swap if reading on different architecture.

2. **POSTYPE**: If compiled without `-DXYZDBL`, positions use 4-byte floats. Adjust dtypes accordingly.

3. **Padding**: The HaloQ structure includes 12 bytes of padding for alignment. This may change in future versions.

4. **Large Files**: Files may exceed 2GB for large simulations. Ensure your reader supports 64-bit file offsets.
