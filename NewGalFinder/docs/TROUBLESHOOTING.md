# PGalF Troubleshooting Guide

This guide covers common issues, error messages, and their solutions.

## Table of Contents

1. [Compilation Issues](#1-compilation-issues)
2. [Runtime Errors](#2-runtime-errors)
3. [Memory Problems](#3-memory-problems)
4. [MPI Issues](#4-mpi-issues)
5. [Output Problems](#5-output-problems)
6. [Performance Issues](#6-performance-issues)
7. [Scientific Validation](#7-scientific-validation)

---

## 1. Compilation Issues

### 1.1 Missing FFTW3 Library

**Error**:
```
ld: cannot find -lfftw3f
```

**Solution**:
1. Install FFTW3 with single-precision and OpenMP support:
   ```bash
   ./configure --enable-float --enable-openmp --prefix=/path/to/install
   make && make install
   ```

2. Update Makefile:
   ```makefile
   FFTW = /path/to/fftw3
   LIBS = -L$(FFTW)/lib -lfftw3f_omp -lfftw3f -lm
   ```

3. Set library path:
   ```bash
   export LD_LIBRARY_PATH=/path/to/fftw3/lib:$LD_LIBRARY_PATH
   ```

### 1.2 MPI Compiler Not Found

**Error**:
```
mpiicx: command not found
```

**Solution**:
- For Intel MPI:
  ```bash
  source /opt/intel/oneapi/setvars.sh
  ```

- For OpenMPI:
  ```bash
  export PATH=/path/to/openmpi/bin:$PATH
  # Change Makefile: CC = mpicc, FC = mpif90
  ```

### 1.3 OpenMP Not Enabled

**Error**:
```
undefined reference to `omp_get_thread_num'
```

**Solution**:
Ensure OpenMP flag is set:
```makefile
# Intel: -qopenmp
# GCC: -fopenmp
OPT = -O3 -qopenmp ...
```

### 1.4 Incompatible Type Sizes

**Error**:
```
error: 'idtype' has not been declared
```

**Solution**:
Check `ramses.h` for proper type definitions:
```c
#ifdef DOUBLE_PRECISION
typedef long long idtype;
#else
typedef int idtype;
#endif
```

Ensure compile flags match your simulation:
```makefile
OPT = ... -DDOUBLE_PRECISION  # If using 64-bit IDs
```

---

## 2. Runtime Errors

### 2.1 Cannot Open Input File

**Error**:
```
Error opening ./FoF_Data/FoF.00100/FoF_halo_cat.00100
```

**Solution**:
1. Verify directory structure:
   ```bash
   ls -la ./FoF_Data/FoF.00100/
   ```

2. Check file names match snapshot number:
   ```bash
   # Should contain:
   # FoF_halo_cat.00100
   # FoF_member_particle.00100
   ```

3. Verify file permissions:
   ```bash
   chmod 644 ./FoF_Data/FoF.00100/*
   ```

### 2.2 Segmentation Fault at Startup

**Error**:
```
Segmentation fault (core dumped)
```

**Diagnosis**:
```bash
# Run with debugger
mpirun -np 2 xterm -e gdb ./gfind.exe

# Or check stack trace
ulimit -c unlimited
mpirun -np 2 ./gfind.exe 100
gdb ./gfind.exe core.*
```

**Common Causes**:
1. Insufficient stack size:
   ```bash
   ulimit -s unlimited
   export OMP_STACKSIZE=512M
   ```

2. Memory allocation failure (see Section 3)

3. Corrupted input files - regenerate FoF data

### 2.3 MPI Assertion Failure

**Error**:
```
MPI_ABORT was invoked on rank 0
```

**Solution**:
1. Check all ranks have access to input files
2. Verify network connectivity between nodes
3. Ensure sufficient file descriptors:
   ```bash
   ulimit -n 65536
   ```

### 2.4 Infinite Loop / Hang

**Symptoms**:
- Process runs indefinitely without output
- CPU usage at 100% but no progress

**Diagnosis**:
```bash
# Attach debugger to running process
gdb -p <pid>
(gdb) bt  # Show backtrace
```

**Common Causes**:
1. Water-shedding not converging - check `MAXNUMWATERSHEDDING`
2. MPI deadlock - check communication patterns
3. Corrupted particle data - verify FoF output

---

## 3. Memory Problems

### 3.1 Memory Allocation Failure

**Error**:
```
P2 Error initializing mem 90000MB
```

**Solution**:
1. Reduce memory per worker in Makefile:
   ```makefile
   DFLAGS = -DNMEG=30000L  # 30 GB instead of 90 GB
   ```

2. Use more MPI ranks with less memory each

3. Check system limits:
   ```bash
   ulimit -v  # Virtual memory limit
   free -h    # Available memory
   ```

### 3.2 Out of Memory During Processing

**Error**:
```
Cannot allocate memory
```
or
```
std::bad_alloc
```

**Solution**:
1. Reduce grid resolution:
   ```c
   #define TSC_CELL_SIZE 0.008  // Increase from 0.004
   ```

2. Reduce maximum cores:
   ```c
   #define MAXNUMCORE 100000    // Reduce from 1000000
   ```

3. Process fewer halos per node (use more MPI ranks)

### 3.3 Memory Fragmentation

**Symptoms**:
- Memory usage grows over time
- Performance degrades

**Solution**:
The code uses stack-based allocation to prevent this. If issues persist:
1. Restart job periodically
2. Check for memory leaks with Valgrind:
   ```bash
   mpirun -np 1 valgrind --leak-check=full ./gfind.exe 100
   ```

---

## 4. MPI Issues

### 4.1 Message Truncation

**Error**:
```
MPI_ERR_TRUNCATE: message truncated
```

**Solution**:
Large halos may exceed MPI buffer limits. The code handles this with `BIG_MPI_Send/Recv`, but verify:
```c
// In gfind.c, check threshold
if(snp*sizeof(FoFTPtlStruct)>1500000000L)
    BIG_MPI_Send(...);
```

### 4.2 Rank Mismatch

**Error**:
```
MPI error: rank 0 expected message from rank 3, received from rank 5
```

**Solution**:
1. Check for tag conflicts in header.h
2. Ensure all ranks execute same code path
3. Verify no premature exits

### 4.3 Timeout / Slow Communication

**Symptoms**:
- Long pauses between halos
- Network timeouts

**Solution**:
1. Check network bandwidth between nodes
2. Increase MPI timeout:
   ```bash
   export MPICH_NEMESIS_NETMOD_TCP_TIMEOUT=1000000
   ```
3. Use fewer nodes with more cores per node

### 4.4 Master Rank Overwhelmed

**Symptoms**:
- Rank 0 at 100% CPU
- Workers idle waiting

**Solution**:
This is expected for I/O-bound jobs. Optimize by:
1. Using faster storage (SSD, parallel filesystem)
2. Pre-staging data to local disk
3. Reducing output frequency

---

## 5. Output Problems

### 5.1 Empty Output Files

**Problem**:
```bash
ls -la GALFIND.DATA.00100
-rw-r--r-- 1 user group 0 Jan 1 00:00 GALFIND.DATA.00100
```

**Causes**:
1. No halos met minimum particle threshold (30 by default)
2. All halos below stellar mass threshold

**Solution**:
1. Check input data has sufficient particles
2. Lower thresholds:
   ```c
   #define MINCORENMEM 10
   #define MINSTELLARMASS 1.e4
   ```

### 5.2 Missing Subhalos

**Problem**:
Expected more subhalos than found

**Diagnosis**:
1. Check peak threshold isn't too high
2. Verify density calculation:
   ```c
   // Temporarily enable debug output
   #define DEBUG 1
   #define LOG 1
   ```

**Solution**:
```c
#define PEAKTHRESHOLD 1000.L    // Lower threshold
#define MERGINGPEAKLENGTH 2.e-3 // Reduce merging
```

### 5.3 Corrupted Binary Output

**Symptoms**:
- Cannot read output files
- Unexpected values when parsing

**Solution**:
1. Ensure consistent endianness between writing and reading
2. Check structure padding:
   ```c
   printf("sizeof(HaloInfo) = %zu\n", sizeof(HaloInfo));
   printf("sizeof(SubInfo) = %zu\n", sizeof(SubInfo));
   ```
3. Verify dptype matches between writer and reader

### 5.4 Output File Too Large

**Problem**:
Output files are unexpectedly large

**Solution**:
1. Disable particle output (keep only catalog):
   ```c
   // Modify write_data() to skip particle writing
   ```
2. Use compression:
   ```bash
   gzip GALFIND.DATA.00100
   ```

---

## 6. Performance Issues

### 6.1 Slow Density Calculation

**Symptoms**:
- Long time in TSC/smoothing step
- High memory usage during density

**Solution**:
1. Increase cell size:
   ```c
   #define TSC_CELL_SIZE 0.008
   ```
2. Ensure FFTW wisdom is cached:
   ```bash
   export FFTW_WISDOM_ONLY=1
   ```

### 6.2 Poor Parallel Scaling

**Symptoms**:
- Adding more MPI ranks doesn't help
- Workers idle frequently

**Diagnosis**:
```bash
# Check load balance
grep "passing through" output.log | awk '{print $NF}'
```

**Solution**:
1. Large halos dominate - use more OMP threads per rank
2. Pre-sort halos by size in FoF stage
3. Use dynamic scheduling (already implemented)

### 6.3 High Memory Bandwidth Usage

**Symptoms**:
- Performance varies with memory speed
- Multiple threads slow down

**Solution**:
1. Bind threads to cores:
   ```bash
   export OMP_PROC_BIND=close
   export OMP_PLACES=cores
   ```
2. Reduce MAXTHREADS if memory-bound:
   ```c
   #define MAXTHREADS 32
   ```

### 6.4 I/O Bottleneck

**Symptoms**:
- Disk usage at 100%
- Long pauses for file operations

**Solution**:
1. Use parallel filesystem (Lustre, GPFS)
2. Stage data to local SSD:
   ```bash
   cp -r FoF_Data /local/scratch/
   cd /local/scratch/
   ./gfind.exe 100
   ```
3. Increase filesystem striping

---

## 7. Scientific Validation

### 7.1 Unrealistic Galaxy Masses

**Problem**:
Galaxy masses too high or too low

**Checks**:
1. Verify unit conversions:
   ```c
   // In physical_parameters()
   printf("pntmass = %g\n", pntmass);
   printf("potentfact = %g\n", potentfact);
   ```
2. Check input file units match code expectations (Msun/h)

### 7.2 Too Few/Many Subhalos

**Problem**:
Subhalo count doesn't match expectations

**Calibration**:
1. Compare with known simulation (e.g., Illustris subhalo count)
2. Adjust key parameters:
   - `PEAKTHRESHOLD`: Higher = fewer subhalos
   - `MINCORENMEM`: Higher = fewer small subhalos
   - `MERGINGPEAKLENGTH`: Higher = fewer close subhalos

### 7.3 Incorrect Boundedness

**Problem**:
Galaxies losing too many/few particles

**Diagnosis**:
1. Check potential calculation:
   ```c
   // Add debug output in CheckSelfTE()
   printf("KE=%g PE=%g total=%g\n", KE, PE, KE+PE);
   ```

2. Verify velocity units:
   ```c
   printf("r2kineticfact = %g\n", r2kineticfact);
   ```

**Solution**:
```c
#define BOUNDITER 6    // More iterations for convergence
```

### 7.4 Visual Inspection

Always validate results visually:

```python
import numpy as np
import matplotlib.pyplot as plt

# Read subhalo data
data = np.fromfile('GALFIND.DATA.00100', dtype=your_dtype)

# Plot particle distribution
plt.scatter(data['x'], data['y'], s=0.1)
plt.savefig('check.png')
```

---

## 8. Quick Diagnostic Commands

```bash
# Check compilation
make clean && make all 2>&1 | grep -i error

# Test with minimal setup
mpirun -np 2 ./gfind.exe 100

# Monitor memory
watch -n 1 'free -h; ps aux | grep gfind'

# Check MPI communication
mpirun -np 4 ./gfind.exe 100 2>&1 | grep -E "(send|recv|error)"

# Profile performance
mpirun -np 4 valgrind --tool=callgrind ./gfind.exe 100

# Check output integrity
file GALFIND.DATA.00100
hexdump -C GALFIND.DATA.00100 | head
```

---

## 9. Getting Help

If issues persist:

1. **Collect diagnostic information**:
   - Full error message
   - Compilation flags used
   - Input file sizes
   - System specifications (RAM, cores, MPI version)

2. **Minimal reproducible example**:
   - Smallest input that reproduces the issue
   - Exact commands to reproduce

3. **Check log files**:
   ```bash
   ls -la *.log *.err
   cat slurm-*.out
   ```
