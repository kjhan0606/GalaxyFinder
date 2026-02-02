# opFoF Usage Examples

This document provides comprehensive examples for using opFoF in various scenarios, from basic usage to advanced workflows.

---

## Table of Contents

1. [Quick Start Examples](#quick-start-examples)
2. [Single Snapshot Analysis](#single-snapshot-analysis)
3. [Multiple Snapshot Processing](#multiple-snapshot-processing)
4. [Reading Output Files](#reading-output-files)
5. [Post-Processing Workflows](#post-processing-workflows)
6. [Batch Processing Scripts](#batch-processing-scripts)
7. [Integration with Other Tools](#integration-with-other-tools)
8. [Visualization Examples](#visualization-examples)

---

## Quick Start Examples

### Minimal Example

```bash
# Navigate to data directory
cd /path/to/simulation/output

# Run FoF on snapshot 100 with 64 data files
mpirun -np 8 /path/to/opfof.exe 100 64

# Output files created:
# - FoF_halo_cat.00100
# - FoF_member_particle.00100
```

### Check Results

```bash
# View halo statistics
/path/to/checkfof FoF_halo_cat.00100

# Calculate mass function
/path/to/fofmassfunc FoF_halo_cat.00100 > mass_function.dat
```

---

## Single Snapshot Analysis

### Basic Usage

```bash
# Syntax: opfof.exe <snapshot_number> <num_files>

# Example: Analyze snapshot 150 with 128 data files
mpirun -np 16 ./opfof.exe 150 128
```

### With Different MPI Configurations

```bash
# Using 4 nodes, 8 tasks per node
mpirun -np 32 -ppn 8 ./opfof.exe 150 128

# Using all cores on single node
mpirun -np $(nproc) ./opfof.exe 150 128

# Explicit host specification
mpirun -np 32 -hosts node01,node02,node03,node04 ./opfof.exe 150 128
```

### Input File Requirements

Ensure these files exist in the current directory:

```
SN.00150.DM.00000.dat
SN.00150.DM.00001.dat
...
SN.00150.DM.00127.dat

SN.00150.STAR.00000.dat
...

SN.00150.GAS.00000.dat
...

SN.00150.SINK.00000.dat
...

SN.00150.00000.info   # Header/info file
```

---

## Multiple Snapshot Processing

### Processing a Range of Snapshots

```bash
# Syntax: opfof.exe <start_snap> <end_snap> <num_files>

# Process snapshots 100 to 200
mpirun -np 32 ./opfof.exe 100 200 128

# This processes: 100, 101, 102, ..., 200
```

### Selective Snapshot Processing

Use a shell loop for non-contiguous snapshots:

```bash
# Process specific snapshots
for snap in 100 150 200 250 300; do
    mpirun -np 16 ./opfof.exe $snap 128
done
```

### Parallel Snapshot Processing

Process multiple snapshots simultaneously (different snapshots on different node groups):

```bash
#!/bin/bash
# parallel_fof.sh

# Array of snapshots to process
snapshots=(100 110 120 130 140 150)

# Process in parallel (background jobs)
for i in ${!snapshots[@]}; do
    snap=${snapshots[$i]}
    # Use different nodes for each snapshot
    mpirun -np 8 ./opfof.exe $snap 128 > log_${snap}.txt 2>&1 &
done

# Wait for all jobs to complete
wait
echo "All snapshots processed"
```

---

## Reading Output Files

### C Example: Reading Halo Catalog

```c
// read_catalog.c
#include <stdio.h>
#include <stdlib.h>
#include "fof.h"

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
    float header[7];
    fread(header, sizeof(float), 7, fp);

    printf("=== Simulation Parameters ===\n");
    printf("Box size: %.2f Mpc/h\n", header[0]);
    printf("Hubble: %.4f\n", header[1]);
    printf("Omega_m: %.4f\n", header[2]);
    printf("Omega_b: %.4f\n", header[3]);
    printf("Omega_L: %.4f\n", header[4]);
    printf("a_max: %.4f\n", header[5]);
    printf("a_now: %.4f\n", header[6]);
    printf("\n");

    // Read halos
    printf("=== Halo List (Np > 100) ===\n");
    printf("%8s %12s %14s %12s %12s %12s\n",
           "Halo", "Np", "Mass", "X", "Y", "Z");
    printf("------------------------------------------------------------\n");

    HaloQ halo;
    int nhalo = 0;
    int nhalo_large = 0;

    while (fread(&halo, sizeof(HaloQ), 1, fp) == 1) {
        nhalo++;
        if (halo.np > 100) {
            nhalo_large++;
            printf("%8d %12zu %14.4e %12.4f %12.4f %12.4f\n",
                   nhalo, halo.np, halo.mass, halo.x, halo.y, halo.z);
        }
    }

    printf("\n");
    printf("Total halos: %d\n", nhalo);
    printf("Halos with Np > 100: %d\n", nhalo_large);

    fclose(fp);
    return 0;
}
```

Compile and run:
```bash
gcc -o read_catalog read_catalog.c -I/path/to/opfof
./read_catalog FoF_halo_cat.00100
```

### Python Example: Reading Halo Catalog

```python
#!/usr/bin/env python3
"""
read_fof_catalog.py - Read opFoF halo catalog
"""
import numpy as np
import struct
import sys

def read_fof_catalog(filename, min_particles=30):
    """
    Read opFoF halo catalog file.

    Parameters
    ----------
    filename : str
        Path to halo catalog file
    min_particles : int
        Minimum particle count to include halo

    Returns
    -------
    dict : Dictionary containing header and halo data
    """
    with open(filename, 'rb') as f:
        # Read header (7 floats)
        header_data = struct.unpack('7f', f.read(28))
        header = {
            'box_size': header_data[0],
            'hubble': header_data[1],
            'omega_m': header_data[2],
            'omega_b': header_data[3],
            'omega_l': header_data[4],
            'amax': header_data[5],
            'anow': header_data[6],
            'redshift': 1.0/header_data[6] - 1.0
        }

        # Define HaloQ structure
        # Adjust POSTYPE based on compilation (double if -DXYZDBL)
        halo_dtype = np.dtype([
            ('np', 'u8'),
            ('npstar', 'u8'),
            ('npgas', 'u8'),
            ('npdm', 'u8'),
            ('npsink', 'u8'),
            ('x', 'f8'),  # double if XYZDBL
            ('y', 'f8'),
            ('z', 'f8'),
            ('mass', 'f8'),
            ('mstar', 'f8'),
            ('mgas', 'f8'),
            ('mdm', 'f8'),
            ('msink', 'f8'),
            ('vx', 'f4'),
            ('vy', 'f4'),
            ('vz', 'f4')
        ])

        # Read all halos
        halos = np.fromfile(f, dtype=halo_dtype)

    # Filter by particle count
    mask = halos['np'] >= min_particles
    halos = halos[mask]

    return {'header': header, 'halos': halos}


def print_statistics(data):
    """Print halo statistics."""
    header = data['header']
    halos = data['halos']

    print("=" * 60)
    print("SIMULATION PARAMETERS")
    print("=" * 60)
    print(f"Box size:    {header['box_size']:.2f} Mpc/h")
    print(f"Redshift:    {header['redshift']:.4f}")
    print(f"Scale factor: {header['anow']:.4f}")
    print(f"Omega_m:     {header['omega_m']:.4f}")
    print(f"Omega_b:     {header['omega_b']:.4f}")
    print(f"Omega_L:     {header['omega_l']:.4f}")
    print()

    print("=" * 60)
    print("HALO STATISTICS")
    print("=" * 60)
    print(f"Total halos: {len(halos)}")
    print(f"Total mass:  {halos['mass'].sum():.4e} Msun/h")
    print()

    print("Mass statistics:")
    print(f"  Min mass:    {halos['mass'].min():.4e} Msun/h")
    print(f"  Max mass:    {halos['mass'].max():.4e} Msun/h")
    print(f"  Mean mass:   {halos['mass'].mean():.4e} Msun/h")
    print(f"  Median mass: {np.median(halos['mass']):.4e} Msun/h")
    print()

    print("Particle count statistics:")
    print(f"  Min Np:    {halos['np'].min()}")
    print(f"  Max Np:    {halos['np'].max()}")
    print(f"  Mean Np:   {halos['np'].mean():.1f}")
    print()

    # Component breakdown
    total_dm = halos['npdm'].sum()
    total_star = halos['npstar'].sum()
    total_gas = halos['npgas'].sum()
    total_sink = halos['npsink'].sum()
    total_all = halos['np'].sum()

    print("Component breakdown:")
    print(f"  Dark matter: {total_dm:>12d} ({100*total_dm/total_all:.1f}%)")
    print(f"  Stars:       {total_star:>12d} ({100*total_star/total_all:.1f}%)")
    print(f"  Gas:         {total_gas:>12d} ({100*total_gas/total_all:.1f}%)")
    print(f"  Sinks:       {total_sink:>12d} ({100*total_sink/total_all:.1f}%)")


def save_ascii_catalog(data, output_file):
    """Save catalog as ASCII file."""
    halos = data['halos']
    header = data['header']

    with open(output_file, 'w') as f:
        f.write(f"# opFoF Halo Catalog\n")
        f.write(f"# Box size: {header['box_size']:.2f} Mpc/h\n")
        f.write(f"# Redshift: {header['redshift']:.4f}\n")
        f.write(f"# Columns: ID Np Mass X Y Z Vx Vy Vz\n")
        f.write(f"# Units: - - Msun/h Mpc/h Mpc/h Mpc/h km/s km/s km/s\n")

        for i, h in enumerate(halos):
            f.write(f"{i:8d} {h['np']:10d} {h['mass']:.6e} "
                   f"{h['x']:.6f} {h['y']:.6f} {h['z']:.6f} "
                   f"{h['vx']:.2f} {h['vy']:.2f} {h['vz']:.2f}\n")

    print(f"Saved ASCII catalog to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <catalog_file> [min_particles]")
        sys.exit(1)

    filename = sys.argv[1]
    min_particles = int(sys.argv[2]) if len(sys.argv) > 2 else 30

    # Read catalog
    data = read_fof_catalog(filename, min_particles)

    # Print statistics
    print_statistics(data)

    # Optional: save ASCII version
    # save_ascii_catalog(data, filename + '.txt')
```

### Python Example: Reading Member Particles

```python
#!/usr/bin/env python3
"""
read_member_particles.py - Read opFoF member particle file
"""
import numpy as np
import struct

def read_member_particles(catalog_file, member_file, halo_index):
    """
    Read particles for a specific halo.

    Parameters
    ----------
    catalog_file : str
        Path to halo catalog
    member_file : str
        Path to member particle file
    halo_index : int
        Index of halo to read (0-based)

    Returns
    -------
    dict : Particle data by type
    """
    # First read catalog to get particle counts
    with open(catalog_file, 'rb') as f:
        f.read(28)  # Skip header

        halo_dtype = np.dtype([
            ('np', 'u8'), ('npstar', 'u8'), ('npgas', 'u8'),
            ('npdm', 'u8'), ('npsink', 'u8'),
            ('x', 'f8'), ('y', 'f8'), ('z', 'f8'),
            ('mass', 'f8'), ('mstar', 'f8'), ('mgas', 'f8'),
            ('mdm', 'f8'), ('msink', 'f8'),
            ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4')
        ])
        halos = np.fromfile(f, dtype=halo_dtype)

    if halo_index >= len(halos):
        raise ValueError(f"Halo index {halo_index} out of range (max: {len(halos)-1})")

    target_halo = halos[halo_index]

    # Calculate offset in member file
    # Sum particles in all previous halos
    offset = 28  # Header size
    for i in range(halo_index):
        h = halos[i]
        offset += h['npdm'] * DM_PARTICLE_SIZE
        offset += h['npgas'] * GAS_PARTICLE_SIZE
        offset += h['npsink'] * SINK_PARTICLE_SIZE
        offset += h['npstar'] * STAR_PARTICLE_SIZE

    # Read particles for this halo
    with open(member_file, 'rb') as f:
        f.seek(offset)

        # Read each component
        particles = {}

        if target_halo['npdm'] > 0:
            dm_dtype = np.dtype([
                ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
                ('mass', 'f4'), ('id', 'i8'), ('level', 'i4')
            ])
            particles['dm'] = np.fromfile(f, dtype=dm_dtype,
                                          count=target_halo['npdm'])

        if target_halo['npgas'] > 0:
            gas_dtype = np.dtype([
                ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
                ('mass', 'f4'), ('density', 'f4'), ('temperature', 'f4')
            ])
            particles['gas'] = np.fromfile(f, dtype=gas_dtype,
                                           count=target_halo['npgas'])

        # ... similar for sink and star

    return particles


# Particle sizes (adjust based on actual structure)
DM_PARTICLE_SIZE = 28    # 3*4 + 4 + 8 + 4
GAS_PARTICLE_SIZE = 24   # 3*4 + 4 + 4 + 4
SINK_PARTICLE_SIZE = 32  # 3*4 + 4 + 4 + 3*4
STAR_PARTICLE_SIZE = 28  # 3*4 + 4 + 8 + 4
```

---

## Post-Processing Workflows

### Calculate Halo Mass Function

```bash
# Using built-in tool
./fofmassfunc FoF_halo_cat.00100 > hmf_z0.dat

# Output format: log10(M), dn/dlogM, N(>M)
```

Python implementation:
```python
#!/usr/bin/env python3
"""
calculate_hmf.py - Calculate halo mass function
"""
import numpy as np

def calculate_hmf(masses, box_size, nbins=50):
    """
    Calculate differential halo mass function.

    Parameters
    ----------
    masses : array
        Halo masses in Msun/h
    box_size : float
        Box size in Mpc/h
    nbins : int
        Number of mass bins

    Returns
    -------
    log_m : array
        Log10 of mass bin centers
    dn_dlogm : array
        dn/dlogM in (Mpc/h)^-3
    """
    volume = box_size**3

    log_masses = np.log10(masses)
    log_m_min = log_masses.min()
    log_m_max = log_masses.max()

    bins = np.linspace(log_m_min, log_m_max, nbins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    dlogm = bins[1] - bins[0]

    counts, _ = np.histogram(log_masses, bins=bins)
    dn_dlogm = counts / volume / dlogm

    return bin_centers, dn_dlogm


# Example usage
from read_fof_catalog import read_fof_catalog

data = read_fof_catalog('FoF_halo_cat.00100', min_particles=50)
masses = data['halos']['mass']
box_size = data['header']['box_size']

log_m, hmf = calculate_hmf(masses, box_size)

# Save to file
np.savetxt('hmf.dat', np.column_stack([log_m, hmf]),
           header='log10(M) dn/dlogM')
```

### Track Halos Across Snapshots

```python
#!/usr/bin/env python3
"""
track_halos.py - Simple halo tracking across snapshots
"""
import numpy as np
from read_fof_catalog import read_fof_catalog

def find_progenitor(halo, halos_prev, search_radius=5.0):
    """
    Find most massive progenitor within search radius.
    """
    dx = halos_prev['x'] - halo['x']
    dy = halos_prev['y'] - halo['y']
    dz = halos_prev['z'] - halo['z']
    dist = np.sqrt(dx**2 + dy**2 + dz**2)

    # Find candidates within search radius
    candidates = np.where(dist < search_radius)[0]

    if len(candidates) == 0:
        return -1

    # Return most massive candidate
    masses = halos_prev['mass'][candidates]
    return candidates[np.argmax(masses)]


def track_most_massive(snapshot_list, data_dir='.'):
    """
    Track most massive halo across snapshots.
    """
    history = []

    for i, snap in enumerate(snapshot_list):
        filename = f"{data_dir}/FoF_halo_cat.{snap:05d}"
        data = read_fof_catalog(filename)
        halos = data['halos']

        if i == 0:
            # Start with most massive halo
            idx = np.argmax(halos['mass'])
        else:
            # Find progenitor of previous halo
            idx = find_progenitor(history[-1], halos)
            if idx < 0:
                print(f"Lost track at snapshot {snap}")
                break

        history.append(halos[idx])
        print(f"Snapshot {snap}: M = {halos[idx]['mass']:.3e} Msun/h")

    return history


# Example usage
snapshots = list(range(100, 201, 10))  # 100, 110, ..., 200
history = track_most_massive(snapshots)
```

---

## Batch Processing Scripts

### SLURM Job Script

```bash
#!/bin/bash
#SBATCH --job-name=opfof
#SBATCH --partition=compute
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --output=opfof_%j.log
#SBATCH --error=opfof_%j.err

# Load modules
module purge
module load intel-mpi/2021.4
module load intel/2021.4

# Set environment
export I_MPI_SHM_LMT=shm

# Change to data directory
cd /path/to/simulation/output

# Run opFoF
OPFOF=/path/to/opfof.exe
SNAP_START=100
SNAP_END=200
NFILES=512

srun $OPFOF $SNAP_START $SNAP_END $NFILES

# Post-processing
echo "Calculating mass functions..."
for snap in $(seq $SNAP_START $SNAP_END); do
    snapstr=$(printf "%05d" $snap)
    /path/to/fofmassfunc FoF_halo_cat.${snapstr} > hmf_${snapstr}.dat
done

echo "Job completed: $(date)"
```

### PBS Job Script

```bash
#!/bin/bash
#PBS -N opfof_job
#PBS -l nodes=8:ppn=8
#PBS -l walltime=24:00:00
#PBS -l mem=128gb
#PBS -j oe
#PBS -o opfof.log

# Load modules
module load openmpi/4.0.5
module load gcc/10.2.0

cd $PBS_O_WORKDIR

# Generate machinefile
cat $PBS_NODEFILE > machines.txt

# Run opFoF
mpirun -np 64 -machinefile machines.txt ./opfof.exe 100 200 512
```

### Array Job for Multiple Snapshots

```bash
#!/bin/bash
#SBATCH --job-name=opfof_array
#SBATCH --array=100-200:10
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=4:00:00
#SBATCH --output=opfof_%A_%a.log

module load intel-mpi

SNAPSHOT=$SLURM_ARRAY_TASK_ID

cd /path/to/data
srun /path/to/opfof.exe $SNAPSHOT 128

echo "Completed snapshot $SNAPSHOT"
```

---

## Integration with Other Tools

### Pipeline: NewDD → opFoF → PGalF

```bash
#!/bin/bash
# full_pipeline.sh - Complete analysis pipeline

SNAPSHOT=100
NFILES=128

# Step 1: Domain decomposition (NewDD)
echo "Running NewDD..."
cd /path/to/newdd
mpirun -np 32 ./NewDD.exe $SNAPSHOT

# Step 2: Halo finding (opFoF)
echo "Running opFoF..."
cd /path/to/output
mpirun -np 32 /path/to/opfof.exe $SNAPSHOT $NFILES

# Step 3: Galaxy finding (PGalF)
echo "Running PGalF..."
mpirun -np 32 /path/to/pgalf.exe $SNAPSHOT

echo "Pipeline complete!"
```

### Export to HDF5

```python
#!/usr/bin/env python3
"""
export_hdf5.py - Export opFoF catalog to HDF5 format
"""
import h5py
import numpy as np
from read_fof_catalog import read_fof_catalog

def export_to_hdf5(input_file, output_file, min_particles=30):
    """Export FoF catalog to HDF5."""

    data = read_fof_catalog(input_file, min_particles)
    header = data['header']
    halos = data['halos']

    with h5py.File(output_file, 'w') as f:
        # Header group
        h = f.create_group('Header')
        h.attrs['BoxSize'] = header['box_size']
        h.attrs['HubbleParam'] = header['hubble']
        h.attrs['Omega0'] = header['omega_m']
        h.attrs['OmegaBaryon'] = header['omega_b']
        h.attrs['OmegaLambda'] = header['omega_l']
        h.attrs['Redshift'] = header['redshift']
        h.attrs['ScaleFactor'] = header['anow']
        h.attrs['NumHalos'] = len(halos)

        # Halo data
        g = f.create_group('Halos')
        g.create_dataset('NumParticles', data=halos['np'])
        g.create_dataset('NumDM', data=halos['npdm'])
        g.create_dataset('NumStar', data=halos['npstar'])
        g.create_dataset('NumGas', data=halos['npgas'])
        g.create_dataset('NumSink', data=halos['npsink'])

        g.create_dataset('Mass', data=halos['mass'])
        g.create_dataset('MassDM', data=halos['mdm'])
        g.create_dataset('MassStar', data=halos['mstar'])
        g.create_dataset('MassGas', data=halos['mgas'])
        g.create_dataset('MassSink', data=halos['msink'])

        pos = np.column_stack([halos['x'], halos['y'], halos['z']])
        g.create_dataset('Position', data=pos)

        vel = np.column_stack([halos['vx'], halos['vy'], halos['vz']])
        g.create_dataset('Velocity', data=vel)

    print(f"Exported {len(halos)} halos to {output_file}")


# Example
export_to_hdf5('FoF_halo_cat.00100', 'halos_z0.hdf5')
```

---

## Visualization Examples

### Plot Mass Function

```python
#!/usr/bin/env python3
"""
plot_hmf.py - Plot halo mass function
"""
import numpy as np
import matplotlib.pyplot as plt
from read_fof_catalog import read_fof_catalog

# Read data
data = read_fof_catalog('FoF_halo_cat.00100', min_particles=50)
masses = data['halos']['mass']
box_size = data['header']['box_size']
z = data['header']['redshift']

# Calculate HMF
volume = box_size**3
log_m = np.log10(masses)
bins = np.linspace(log_m.min(), log_m.max(), 30)
counts, edges = np.histogram(log_m, bins=bins)
centers = 0.5 * (edges[:-1] + edges[1:])
dlogm = edges[1] - edges[0]
hmf = counts / volume / dlogm

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.semilogy(centers, hmf, 'ko-', label=f'opFoF (z={z:.2f})')
ax.set_xlabel(r'$\log_{10}(M / [M_\odot/h])$', fontsize=12)
ax.set_ylabel(r'$dn/d\log M$ [(Mpc/h)$^{-3}$]', fontsize=12)
ax.set_title('Halo Mass Function')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('hmf.png', dpi=150)
plt.show()
```

### Plot Halo Distribution

```python
#!/usr/bin/env python3
"""
plot_halos.py - Plot halo spatial distribution
"""
import numpy as np
import matplotlib.pyplot as plt
from read_fof_catalog import read_fof_catalog

data = read_fof_catalog('FoF_halo_cat.00100', min_particles=100)
halos = data['halos']
box_size = data['header']['box_size']

# Create figure with subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Size proportional to mass
sizes = 10 * (np.log10(halos['mass']) - 10)
sizes = np.clip(sizes, 1, 100)

# X-Y projection
ax = axes[0]
ax.scatter(halos['x'], halos['y'], s=sizes, alpha=0.5, c='blue')
ax.set_xlabel('X [Mpc/h]')
ax.set_ylabel('Y [Mpc/h]')
ax.set_title('X-Y Projection')
ax.set_xlim(0, box_size)
ax.set_ylim(0, box_size)
ax.set_aspect('equal')

# X-Z projection
ax = axes[1]
ax.scatter(halos['x'], halos['z'], s=sizes, alpha=0.5, c='blue')
ax.set_xlabel('X [Mpc/h]')
ax.set_ylabel('Z [Mpc/h]')
ax.set_title('X-Z Projection')
ax.set_xlim(0, box_size)
ax.set_ylim(0, box_size)
ax.set_aspect('equal')

# Y-Z projection
ax = axes[2]
ax.scatter(halos['y'], halos['z'], s=sizes, alpha=0.5, c='blue')
ax.set_xlabel('Y [Mpc/h]')
ax.set_ylabel('Z [Mpc/h]')
ax.set_title('Y-Z Projection')
ax.set_xlim(0, box_size)
ax.set_ylim(0, box_size)
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig('halo_distribution.png', dpi=150)
plt.show()
```

### 3D Visualization

```python
#!/usr/bin/env python3
"""
plot_3d.py - 3D visualization of halos
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from read_fof_catalog import read_fof_catalog

data = read_fof_catalog('FoF_halo_cat.00100', min_particles=100)
halos = data['halos']

# Select massive halos for visualization
mask = halos['mass'] > 1e12
massive = halos[mask]

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Color by mass
colors = np.log10(massive['mass'])
sizes = 20 * (colors - colors.min()) / (colors.max() - colors.min()) + 5

scatter = ax.scatter(massive['x'], massive['y'], massive['z'],
                     c=colors, s=sizes, cmap='viridis', alpha=0.6)

ax.set_xlabel('X [Mpc/h]')
ax.set_ylabel('Y [Mpc/h]')
ax.set_zlabel('Z [Mpc/h]')
ax.set_title(f'Massive Halos (M > 10$^{{12}}$ M$_\\odot$/h)')

cbar = plt.colorbar(scatter, ax=ax, shrink=0.5, pad=0.1)
cbar.set_label(r'$\log_{10}(M / [M_\odot/h])$')

plt.tight_layout()
plt.savefig('halos_3d.png', dpi=150)
plt.show()
```

---

## Tips and Best Practices

1. **Always check input files** before running opFoF
2. **Start with small test runs** to verify setup
3. **Monitor memory usage** during execution
4. **Back up output files** before re-running
5. **Use version control** for analysis scripts
6. **Document parameter choices** in log files
