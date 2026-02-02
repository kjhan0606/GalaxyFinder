# PGalF Usage Examples

This document provides practical examples for running PGalF in various scenarios.

## Table of Contents

1. [Basic Usage](#1-basic-usage)
2. [Job Scripts](#2-job-scripts)
3. [Reading Output Files](#3-reading-output-files)
4. [Analysis Examples](#4-analysis-examples)
5. [Batch Processing](#5-batch-processing)

---

## 1. Basic Usage

### 1.1 Single Snapshot Processing

```bash
# Process snapshot 100 with default settings
cd /path/to/NewGalFinder
./gfind.exe 100

# Expected directory structure:
# ./FoF_Data/FoF.00100/
#   ├── FoF_halo_cat.00100          (input)
#   ├── FoF_member_particle.00100   (input)
#   ├── GALCATALOG.LIST.00100       (output)
#   ├── GALFIND.DATA.00100          (output)
#   └── background_ptl.00100        (output)
```

### 1.2 Parallel Execution

```bash
# Set OpenMP threads
export OMP_NUM_THREADS=8
export OMP_STACKSIZE=256M

# Run with 16 MPI processes
mpirun -np 16 ./gfind.exe 100

# For Intel MPI
mpiexec.hydra -np 16 ./gfind.exe 100

# For SLURM
srun -n 16 ./gfind.exe 100
```

### 1.3 Restart from Checkpoint

If processing was interrupted, restart from specific offsets:

```bash
# Find last processed position from log
grep "offsets=" output.log | tail -1
# Example output: offsets= 123456L 7890123L

# Restart from those offsets
./gfind.exe 100 123456 7890123
```

---

## 2. Job Scripts

### 2.1 SLURM Script (Basic)

```bash
#!/bin/bash
#SBATCH --job-name=pgalf
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --output=pgalf_%j.out
#SBATCH --error=pgalf_%j.err

module load intel-mpi fftw

export OMP_NUM_THREADS=4
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close

cd /path/to/NewGalFinder

srun ./gfind.exe 100
```

### 2.2 SLURM Script (Large Memory)

For processing massive halos:

```bash
#!/bin/bash
#SBATCH --job-name=pgalf_large
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=256G
#SBATCH --partition=largemem
#SBATCH --output=pgalf_%j.out

module load intel-mpi fftw

export OMP_NUM_THREADS=16
export OMP_STACKSIZE=1G
export OMP_PROC_BIND=spread

cd /path/to/NewGalFinder

srun ./gfind.exe 100
```

### 2.3 PBS Script

```bash
#!/bin/bash
#PBS -N pgalf
#PBS -l nodes=2:ppn=32
#PBS -l walltime=24:00:00
#PBS -l mem=256gb
#PBS -o pgalf.out
#PBS -e pgalf.err

module load openmpi fftw

export OMP_NUM_THREADS=8

cd $PBS_O_WORKDIR

mpirun -np 8 ./gfind.exe 100
```

### 2.4 LSF Script

```bash
#!/bin/bash
#BSUB -J pgalf
#BSUB -n 64
#BSUB -R "span[ptile=16]"
#BSUB -W 24:00
#BSUB -o pgalf_%J.out
#BSUB -e pgalf_%J.err

module load intel-mpi fftw

export OMP_NUM_THREADS=4

mpirun ./gfind.exe 100
```

---

## 3. Reading Output Files

### 3.1 Python: Read Galaxy Catalog

```python
import numpy as np
import struct

def read_galcatalog(filename):
    """Read GALFIND.DATA output file"""

    # Define data types matching C structures
    haloinfo_dtype = np.dtype([
        ('nsub', 'i4'),
        ('ndm', 'i4'),
        ('nstar', 'i4'),
        ('nsink', 'i4'),
        ('ngas', 'i4'),
        ('npall', 'i4'),
        ('totm', 'f8'),
        ('mdm', 'f8'),
        ('mgas', 'f8'),
        ('msink', 'f8'),
        ('mstar', 'f8'),
        ('x', 'f8'),
        ('y', 'f8'),
        ('z', 'f8'),
        ('vx', 'f8'),
        ('vy', 'f8'),
        ('vz', 'f8'),
    ])

    subinfo_dtype = np.dtype([
        ('npdm', 'i4'),
        ('npgas', 'i4'),
        ('npsink', 'i4'),
        ('npstar', 'i4'),
        ('npall', 'i4'),
        ('totm', 'f8'),
        ('mdm', 'f8'),
        ('mgas', 'f8'),
        ('msink', 'f8'),
        ('mstar', 'f8'),
        ('x', 'f8'),
        ('y', 'f8'),
        ('z', 'f8'),
        ('vx', 'f8'),
        ('vy', 'f8'),
        ('vz', 'f8'),
    ])

    halos = []
    subhalos = []

    with open(filename, 'rb') as f:
        while True:
            # Read halo header
            data = f.read(haloinfo_dtype.itemsize)
            if len(data) < haloinfo_dtype.itemsize:
                break

            halo = np.frombuffer(data, dtype=haloinfo_dtype)[0]
            halos.append(halo)

            # Read subhalos
            for i in range(halo['nsub']):
                sub_data = f.read(subinfo_dtype.itemsize)
                sub = np.frombuffer(sub_data, dtype=subinfo_dtype)[0]
                subhalos.append(sub)

                # Skip particle data (if present in GALFIND.DATA)
                # Calculate size based on particle counts
                skip_size = (sub['npdm'] * 72 +      # DmType size
                            sub['npgas'] * 120 +     # GasType size
                            sub['npsink'] * 136 +    # SinkType size
                            sub['npstar'] * 72)      # StarType size
                f.seek(skip_size, 1)

    return np.array(halos, dtype=haloinfo_dtype), np.array(subhalos, dtype=subinfo_dtype)

# Usage
halos, subhalos = read_galcatalog('FoF_Data/FoF.00100/GALFIND.DATA.00100')
print(f"Found {len(halos)} halos with {len(subhalos)} total subhalos")
```

### 3.2 Python: Read Catalog Only (LIST file)

```python
import numpy as np

def read_galcatalog_list(filename):
    """Read GALCATALOG.LIST file (no particle data)"""

    haloinfo_dtype = np.dtype([
        ('nsub', 'i4'),
        ('ndm', 'i4'),
        ('nstar', 'i4'),
        ('nsink', 'i4'),
        ('ngas', 'i4'),
        ('npall', 'i4'),
        ('totm', 'f8'),
        ('mdm', 'f8'),
        ('mgas', 'f8'),
        ('msink', 'f8'),
        ('mstar', 'f8'),
        ('x', 'f8'),
        ('y', 'f8'),
        ('z', 'f8'),
        ('vx', 'f8'),
        ('vy', 'f8'),
        ('vz', 'f8'),
    ])

    subinfo_dtype = np.dtype([
        ('npdm', 'i4'),
        ('npgas', 'i4'),
        ('npsink', 'i4'),
        ('npstar', 'i4'),
        ('npall', 'i4'),
        ('totm', 'f8'),
        ('mdm', 'f8'),
        ('mgas', 'f8'),
        ('msink', 'f8'),
        ('mstar', 'f8'),
        ('x', 'f8'),
        ('y', 'f8'),
        ('z', 'f8'),
        ('vx', 'f8'),
        ('vy', 'f8'),
        ('vz', 'f8'),
    ])

    halos = []
    subhalos = []

    with open(filename, 'rb') as f:
        while True:
            data = f.read(haloinfo_dtype.itemsize)
            if len(data) < haloinfo_dtype.itemsize:
                break

            halo = np.frombuffer(data, dtype=haloinfo_dtype)[0]
            halos.append(halo)

            for i in range(halo['nsub']):
                sub_data = f.read(subinfo_dtype.itemsize)
                sub = np.frombuffer(sub_data, dtype=subinfo_dtype)[0]
                subhalos.append(sub)

    return np.array(halos), np.array(subhalos)

# Usage
halos, subs = read_galcatalog_list('FoF_Data/FoF.00100/GALCATALOG.LIST.00100')
```

### 3.3 C: Read Output

```c
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int nsub, ndm, nstar, nsink, ngas, npall;
    double totm, mdm, mgas, msink, mstar;
    double x, y, z, vx, vy, vz;
} HaloInfo;

typedef struct {
    int npdm, npgas, npsink, npstar, npall;
    double totm, mdm, mgas, msink, mstar;
    double x, y, z, vx, vy, vz;
} SubInfo;

int main() {
    FILE *fp = fopen("GALCATALOG.LIST.00100", "rb");
    HaloInfo halo;
    SubInfo sub;
    int nhalo = 0, nsub_total = 0;

    while (fread(&halo, sizeof(HaloInfo), 1, fp) == 1) {
        nhalo++;
        printf("Halo %d: nsub=%d, Mstar=%.2e Msun/h\n",
               nhalo, halo.nsub, halo.mstar);

        for (int i = 0; i < halo.nsub; i++) {
            fread(&sub, sizeof(SubInfo), 1, fp);
            nsub_total++;
            printf("  Sub %d: Mstar=%.2e, pos=(%.2f,%.2f,%.2f)\n",
                   i, sub.mstar, sub.x, sub.y, sub.z);
        }
    }

    fclose(fp);
    printf("\nTotal: %d halos, %d subhalos\n", nhalo, nsub_total);
    return 0;
}
```

---

## 4. Analysis Examples

### 4.1 Stellar Mass Function

```python
import numpy as np
import matplotlib.pyplot as plt

# Read data
_, subhalos = read_galcatalog_list('GALCATALOG.LIST.00100')

# Get stellar masses (already in Msun/h)
mstar = subhalos['mstar']
mstar = mstar[mstar > 0]  # Remove zero-mass entries

# Create mass bins
log_bins = np.linspace(6, 12, 30)
hist, bin_edges = np.histogram(np.log10(mstar), bins=log_bins)

# Convert to mass function (dN/dlogM per volume)
box_volume = 100**3  # (Mpc/h)^3 - adjust for your simulation
dlogM = log_bins[1] - log_bins[0]
phi = hist / box_volume / dlogM

# Plot
plt.figure(figsize=(8, 6))
plt.semilogy(0.5*(bin_edges[:-1]+bin_edges[1:]), phi, 'o-')
plt.xlabel(r'$\log_{10}(M_*/[M_\odot/h])$')
plt.ylabel(r'$\phi$ [dex$^{-1}$ (Mpc/h)$^{-3}$]')
plt.title('Galaxy Stellar Mass Function')
plt.savefig('smf.png', dpi=150)
```

### 4.2 Subhalo Mass-Distance Relation

```python
import numpy as np
import matplotlib.pyplot as plt

halos, subhalos = read_galcatalog_list('GALCATALOG.LIST.00100')

# Calculate distances from host center
sub_idx = 0
distances = []
masses = []

for halo in halos:
    host_x, host_y, host_z = halo['x'], halo['y'], halo['z']

    for i in range(halo['nsub']):
        sub = subhalos[sub_idx]
        dx = sub['x'] - host_x
        dy = sub['y'] - host_y
        dz = sub['z'] - host_z
        dist = np.sqrt(dx**2 + dy**2 + dz**2)

        distances.append(dist)
        masses.append(sub['totm'])
        sub_idx += 1

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(distances, masses, alpha=0.3, s=1)
plt.xlabel('Distance from host center [cMpc/h]')
plt.ylabel('Subhalo mass [Msun/h]')
plt.xscale('log')
plt.yscale('log')
plt.savefig('mass_distance.png', dpi=150)
```

### 4.3 Spatial Distribution

```python
import numpy as np
import matplotlib.pyplot as plt

_, subhalos = read_galcatalog_list('GALCATALOG.LIST.00100')

# 2D projection
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# XY projection
axes[0].scatter(subhalos['x'], subhalos['y'],
                c=np.log10(subhalos['mstar']+1), s=1, alpha=0.5)
axes[0].set_xlabel('X [cMpc/h]')
axes[0].set_ylabel('Y [cMpc/h]')
axes[0].set_title('XY Projection')

# XZ projection
axes[1].scatter(subhalos['x'], subhalos['z'],
                c=np.log10(subhalos['mstar']+1), s=1, alpha=0.5)
axes[1].set_xlabel('X [cMpc/h]')
axes[1].set_ylabel('Z [cMpc/h]')
axes[1].set_title('XZ Projection')

# YZ projection
axes[2].scatter(subhalos['y'], subhalos['z'],
                c=np.log10(subhalos['mstar']+1), s=1, alpha=0.5)
axes[2].set_xlabel('Y [cMpc/h]')
axes[2].set_ylabel('Z [cMpc/h]')
axes[2].set_title('YZ Projection')

plt.tight_layout()
plt.savefig('spatial_distribution.png', dpi=150)
```

### 4.4 Velocity Dispersion Profile

```python
import numpy as np

def velocity_dispersion(subhalos, r_bins):
    """Calculate velocity dispersion in radial bins"""

    # Calculate velocities relative to mean
    vx = subhalos['vx'] - np.mean(subhalos['vx'])
    vy = subhalos['vy'] - np.mean(subhalos['vy'])
    vz = subhalos['vz'] - np.mean(subhalos['vz'])

    # Radial distance from center
    x = subhalos['x'] - np.mean(subhalos['x'])
    y = subhalos['y'] - np.mean(subhalos['y'])
    z = subhalos['z'] - np.mean(subhalos['z'])
    r = np.sqrt(x**2 + y**2 + z**2)

    sigma = []
    for i in range(len(r_bins)-1):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.sum(mask) > 5:
            v2 = vx[mask]**2 + vy[mask]**2 + vz[mask]**2
            sigma.append(np.sqrt(np.mean(v2)))
        else:
            sigma.append(np.nan)

    return np.array(sigma)
```

---

## 5. Batch Processing

### 5.1 Process Multiple Snapshots

```bash
#!/bin/bash
# process_all.sh - Process snapshots 50 to 100

for snap in $(seq 50 100); do
    echo "Processing snapshot $snap..."

    # Format snapshot number
    snapnum=$(printf "%05d" $snap)

    # Check input exists
    if [ ! -f "FoF_Data/FoF.$snapnum/FoF_halo_cat.$snapnum" ]; then
        echo "  Skipping - input not found"
        continue
    fi

    # Run galaxy finder
    mpirun -np 32 ./gfind.exe $snap > logs/snap_$snapnum.log 2>&1

    # Check success
    if [ -f "FoF_Data/FoF.$snapnum/GALFIND.DATA.$snapnum" ]; then
        echo "  Success!"
    else
        echo "  FAILED - check logs/snap_$snapnum.log"
    fi
done
```

### 5.2 SLURM Array Job

```bash
#!/bin/bash
#SBATCH --job-name=pgalf_array
#SBATCH --array=50-100
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --output=logs/pgalf_%a.out

module load intel-mpi fftw

export OMP_NUM_THREADS=4

cd /path/to/NewGalFinder

srun ./gfind.exe $SLURM_ARRAY_TASK_ID
```

### 5.3 Dependency Chain (Process After FoF)

```bash
#!/bin/bash
# submit_pipeline.sh

# Submit FoF job
FOF_JOB=$(sbatch --parsable run_fof.sh)
echo "FoF job: $FOF_JOB"

# Submit galaxy finder dependent on FoF
GAL_JOB=$(sbatch --parsable --dependency=afterok:$FOF_JOB run_gfind.sh)
echo "GalFind job: $GAL_JOB (waits for $FOF_JOB)"

# Submit analysis dependent on galaxy finder
ANAL_JOB=$(sbatch --parsable --dependency=afterok:$GAL_JOB run_analysis.sh)
echo "Analysis job: $ANAL_JOB (waits for $GAL_JOB)"
```

### 5.4 Monitoring Script

```bash
#!/bin/bash
# monitor.sh - Monitor PGalF progress

while true; do
    clear
    echo "=== PGalF Progress Monitor ==="
    echo ""

    # Count completed snapshots
    n_done=$(ls FoF_Data/FoF.*/GALFIND.DATA.* 2>/dev/null | wc -l)
    n_total=$(ls -d FoF_Data/FoF.* 2>/dev/null | wc -l)
    echo "Completed: $n_done / $n_total snapshots"
    echo ""

    # Show running jobs
    echo "Running jobs:"
    squeue -u $USER -n pgalf --format="%.8i %.9P %.20j %.8T %.10M %.6D %R"
    echo ""

    # Show recent log output
    echo "Recent activity:"
    tail -5 logs/*.log 2>/dev/null | tail -20

    sleep 60
done
```

---

## 6. Quick Reference

### Command Summary

| Task | Command |
|------|---------|
| Single snapshot | `./gfind.exe 100` |
| Parallel run | `mpirun -np 32 ./gfind.exe 100` |
| With OpenMP | `OMP_NUM_THREADS=8 mpirun -np 4 ./gfind.exe 100` |
| Restart | `./gfind.exe 100 <offset1> <offset2>` |
| Debug mode | `mpirun -np 2 ./gfind.exe 100 2>&1 \| tee debug.log` |

### Output Files

| File | Contents | Format |
|------|----------|--------|
| `GALCATALOG.LIST.XXXXX` | Halo + subhalo catalog | Binary |
| `GALFIND.DATA.XXXXX` | Catalog + particles | Binary |
| `background_ptl.XXXXX` | Unassigned particles | Binary |

### Environment Variables

| Variable | Purpose | Example |
|----------|---------|---------|
| `OMP_NUM_THREADS` | OpenMP threads | 8 |
| `OMP_STACKSIZE` | Per-thread stack | 256M |
| `OMP_PROC_BIND` | Thread binding | close |
| `LD_LIBRARY_PATH` | Library search | /path/to/fftw/lib |
