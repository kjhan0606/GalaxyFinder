# opFoF Quick Start Guide

Get started with opFoF in 5 minutes. This guide covers the essential steps to compile and run the halo finder.

---

## Prerequisites

- MPI library (Intel MPI, OpenMPI, or MPICH)
- C compiler (Intel ICC or GCC)
- Simulation data in RAMSES format (from NewDD)

---

## Step 1: Compile

```bash
cd opFoF

# Clean any previous build
make clean

# Compile
make this
```

**Expected output:**
```
mpiicc -c opfof.c -o opfof.o ...
mpiicc -c Treewalk.fof.ordered.c -o Treewalk.fof.ordered.o ...
...
mpiicc -o opfof.exe opfof.o ... -L../ -lmyram -lm
```

If successful, `opfof.exe` will be created.

---

## Step 2: Prepare Data

Ensure your simulation data follows this naming convention:

```
SN.00100.DM.00000.dat
SN.00100.DM.00001.dat
...
SN.00100.STAR.00000.dat
SN.00100.GAS.00000.dat
SN.00100.SINK.00000.dat
SN.00100.00000.info
```

---

## Step 3: Run

### Basic Run

```bash
cd /path/to/data

# Run on 8 MPI processes
# Arguments: <snapshot_number> <number_of_files>
mpirun -np 8 /path/to/opfof.exe 100 64
```

### Example with SLURM

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --time=4:00:00

module load intel-mpi
cd /path/to/data
srun /path/to/opfof.exe 100 128
```

---

## Step 4: Check Output

After successful run, two files are created:

```
FoF_halo_cat.00100          # Halo catalog
FoF_member_particle.00100   # Member particles
```

### Quick validation

```bash
/path/to/checkfof FoF_halo_cat.00100
```

### View mass function

```bash
/path/to/fofmassfunc FoF_halo_cat.00100
```

---

## Step 5: Read Results

### Python

```python
import numpy as np

# Read catalog
with open('FoF_halo_cat.00100', 'rb') as f:
    header = np.fromfile(f, dtype='f4', count=7)
    print(f"Box size: {header[0]} Mpc/h")
    print(f"Redshift: {1/header[6] - 1:.2f}")

    # Read halos
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

print(f"Found {len(halos)} halos")
print(f"Most massive: {halos['mass'].max():.2e} Msun/h")
```

---

## Common Issues

### "Cannot allocate memory"

Increase memory pool in `Rules.make`:
```makefile
-DNMEG=32000L  # Increase from 17000L
```
Then recompile: `make clean && make this`

### "File not found"

Check that input files exist and follow naming convention:
```bash
ls SN.00100.*.dat
```

### Slow execution

- Use more MPI ranks
- Ensure number of files is divisible by number of ranks

---

## Next Steps

- Read the full [README.md](../README.md)
- See [EXAMPLES.md](EXAMPLES.md) for detailed usage
- Check [CONFIGURATION.md](CONFIGURATION.md) for tuning options

---

## Quick Reference

| Command | Description |
|---------|-------------|
| `make clean && make this` | Rebuild |
| `mpirun -np N ./opfof.exe SNAP NFILES` | Run single snapshot |
| `mpirun -np N ./opfof.exe START END NFILES` | Run multiple snapshots |
| `./checkfof FILE` | Validate output |
| `./fofmassfunc FILE` | Calculate mass function |
