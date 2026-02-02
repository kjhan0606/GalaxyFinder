# opFoF Troubleshooting Guide

This document provides solutions for common issues encountered when compiling, running, and analyzing results from opFoF.

---

## Table of Contents

1. [Compilation Issues](#compilation-issues)
2. [Runtime Errors](#runtime-errors)
3. [Memory Issues](#memory-issues)
4. [MPI Problems](#mpi-problems)
5. [Input/Output Issues](#inputoutput-issues)
6. [Performance Problems](#performance-problems)
7. [Result Validation](#result-validation)
8. [Debugging Techniques](#debugging-techniques)
9. [FAQ](#faq)

---

## Compilation Issues

### Error: "mpiicc: command not found"

**Symptom:**
```
make: mpiicc: Command not found
```

**Cause:** Intel MPI compiler not in PATH

**Solutions:**

1. Load Intel MPI module:
   ```bash
   module load intel-mpi
   # or
   source /opt/intel/oneapi/setvars.sh
   ```

2. Use alternative MPI compiler (edit Rules.make):
   ```makefile
   # For OpenMPI/MPICH
   CC = mpicc
   F77 = mpif77

   # Remove Intel-specific flag
   OPT = -g -O2  # Remove -DINTEL
   ```

---

### Error: "undefined reference to `MPI_Init`"

**Symptom:**
```
opfof.o: undefined reference to `MPI_Init'
opfof.o: undefined reference to `MPI_Comm_rank'
...
```

**Cause:** Not using MPI compiler wrapper

**Solution:** Ensure CC uses MPI wrapper:
```makefile
# Wrong
CC = gcc
CC = icc

# Correct
CC = mpicc      # OpenMPI/MPICH
CC = mpiicc     # Intel MPI
```

---

### Error: "cannot find -lmyram"

**Symptom:**
```
/usr/bin/ld: cannot find -lmyram
```

**Cause:** Custom library not found

**Solutions:**

1. Check library exists:
   ```bash
   ls ../libmyram.a
   ```

2. Verify library path in Makefile:
   ```makefile
   LIBS = -L../ -lmyram -lm
   ```

3. Rebuild the library if missing:
   ```bash
   cd ..
   make libmyram.a
   ```

---

### Error: "fof.h: No such file or directory"

**Symptom:**
```
opfof.c:10: fatal error: fof.h: No such file or directory
```

**Cause:** Header file not found

**Solution:** Check include paths:
```makefile
CFLAGS = $(OPT) $(CDFLAGS) $(COMFLAGS) -I. -I..
```

---

### Warning: "implicit declaration of function"

**Symptom:**
```
warning: implicit declaration of function 'FoF_Make_Tree'
```

**Cause:** Missing function prototype

**Solution:** Add proper header includes:
```c
#include "fof.h"
#include "Memory.h"
#include "Time.h"
```

---

### Error: "integer overflow in expression"

**Symptom:**
```
error: integer overflow in expression of type 'int'
```

**Cause:** Large integer without suffix

**Solution:** Use long suffix for large values:
```makefile
# Wrong
-DNMEG=17000

# Correct
-DNMEG=17000L
```

---

## Runtime Errors

### Error: "Segmentation fault"

**Symptom:**
```
[0] Segmentation fault (core dumped)
```

**Possible Causes and Solutions:**

1. **Out of memory**
   ```bash
   # Check memory usage
   ulimit -v unlimited
   # Or increase -DNMEG and recompile
   ```

2. **Invalid input data**
   ```bash
   # Verify input files
   ls SN.00100.*.dat
   file SN.00100.DM.00000.dat
   ```

3. **Stack overflow**
   ```bash
   ulimit -s unlimited
   ```

4. **Debug with gdb**
   ```bash
   mpirun -np 1 xterm -e gdb ./opfof.exe
   ```

---

### Error: "Cannot allocate memory"

**Symptom:**
```
MyMalloc: Cannot allocate 1073741824 bytes
```

**Solutions:**

1. Increase memory pool:
   ```makefile
   -DNMEG=32000L  # Increase from default 17000
   ```

2. Use more MPI ranks to distribute memory:
   ```bash
   mpirun -np 32 ./opfof.exe 100 128  # Instead of -np 8
   ```

3. Check system limits:
   ```bash
   ulimit -v   # Virtual memory limit
   ulimit -m   # Max memory size
   free -h     # Available system memory
   ```

---

### Error: "File not found"

**Symptom:**
```
Cannot open SN.00100.DM.00000.dat
```

**Solutions:**

1. Check file existence and naming:
   ```bash
   ls -la SN.00100.*.dat
   ```

2. Verify naming convention:
   ```
   Expected: SN.NNNNN.TYPE.MMMMM.dat
   Example:  SN.00100.DM.00000.dat
   ```

3. Check permissions:
   ```bash
   chmod +r SN.00100.*.dat
   ```

4. Verify current directory:
   ```bash
   pwd
   ```

---

### Error: "MPI_ABORT was invoked"

**Symptom:**
```
MPI_ABORT was invoked on rank 3
```

**Solutions:**

1. Check all ranks' output:
   ```bash
   mpirun -np 8 ./opfof.exe 100 64 2>&1 | tee output.log
   grep -i error output.log
   ```

2. Run with fewer ranks for debugging:
   ```bash
   mpirun -np 2 ./opfof.exe 100 64
   ```

3. Check for rank-specific issues:
   ```bash
   # Add rank to output
   mpirun -np 8 ./opfof.exe 100 64 2>&1 | grep "Rank"
   ```

---

## Memory Issues

### Issue: Memory Usage Keeps Growing

**Symptom:** Memory usage increases until crash

**Possible Causes:**

1. Memory leak in custom allocator
2. Not freeing temporary buffers
3. Accumulating boundary particles

**Solutions:**

1. Enable memory debugging:
   ```makefile
   CDFLAGS = ... -DMEMORY_DEBUG ...
   ```

2. Check memory at key points:
   ```c
   printf("Memory usage: %zu MB\n", GetMemoryUsage() / (1024*1024));
   ```

---

### Issue: Memory Fragmentation

**Symptom:** Allocation fails even with available memory

**Solution:** The pool allocator should prevent this, but ensure allocations/frees are balanced:

```c
// Always free in reverse order of allocation
ptr3 = MyMalloc(size3);
ptr2 = MyMalloc(size2);
ptr1 = MyMalloc(size1);
// ... use memory ...
MyFree(ptr1);
MyFree(ptr2);
MyFree(ptr3);
```

---

### Issue: OOM Killer Terminates Process

**Symptom:** Process killed without error message

**Solutions:**

1. Check system logs:
   ```bash
   dmesg | grep -i "killed process"
   ```

2. Request more memory in job script:
   ```bash
   #SBATCH --mem=200G
   ```

3. Use more nodes with less memory per node:
   ```bash
   #SBATCH --nodes=8
   #SBATCH --mem-per-cpu=8G
   ```

---

## MPI Problems

### Error: "MPI_Send: Message truncated"

**Symptom:**
```
MPI_Send: Message truncated, error code ...
```

**Cause:** Message larger than MPI buffer

**Solution:** The code uses `BIG_MPI_Send()` for large messages. Verify it's being used:
```c
// In Treewalk.fof.ordered.c
BIG_MPI_Send(buffer, count, MPI_CHAR, dest, tag, comm);
```

---

### Error: "All processes deadlocked"

**Symptom:** Code hangs indefinitely

**Causes and Solutions:**

1. **Mismatched Send/Recv:**
   ```c
   // Ensure each MPI_Send has matching MPI_Recv
   if (myid > 0) MPI_Send(...);
   if (myid < nprocs-1) MPI_Recv(...);
   ```

2. **Use non-blocking calls:**
   ```c
   MPI_Isend(..., &request);
   // ... do other work ...
   MPI_Wait(&request, &status);
   ```

3. **Check with timeout:**
   ```bash
   mpirun --timeout 300 -np 8 ./opfof.exe 100 64
   ```

---

### Error: "No route to host"

**Symptom:**
```
[node1] No route to host connecting to node2
```

**Solutions:**

1. Check network connectivity:
   ```bash
   ping node2
   ssh node2 hostname
   ```

2. Use correct MPI fabric:
   ```bash
   export I_MPI_FABRICS=shm:tcp  # For Ethernet
   export I_MPI_FABRICS=shm:ofi  # For InfiniBand
   ```

---

### Issue: Poor MPI Scaling

**Symptom:** Adding more ranks doesn't improve performance

**Solutions:**

1. Check load balance:
   ```bash
   # Add timing to each rank
   if (myid == 0) printf("Rank 0 time: %.2f s\n", elapsed);
   ```

2. Profile communication:
   ```bash
   # With Intel MPI
   export I_MPI_STATS=10
   mpirun -np 8 ./opfof.exe 100 64
   ```

3. Reduce communication:
   - Increase chunk size for boundary exchange
   - Batch multiple small messages

---

## Input/Output Issues

### Error: "Corrupted input file"

**Symptom:** Unexpected values when reading

**Diagnosis:**
```bash
# Check file size
ls -la SN.00100.DM.00000.dat

# Check header
od -A d -t f4 SN.00100.DM.00000.dat | head

# Compare with expected size
# Expected = header_size + n_particles × particle_size
```

**Solutions:**

1. Re-download or regenerate data
2. Check for truncated files
3. Verify byte order (endianness)

---

### Error: "Output file too large"

**Symptom:** Cannot write complete output

**Solutions:**

1. Ensure large file support:
   ```makefile
   -D_LARGE_FILES -D_FILE_OFFSET_BITS=64
   ```

2. Check filesystem limits:
   ```bash
   df -h .
   ```

3. Use different filesystem:
   ```bash
   cd /scratch/large_storage
   mpirun ... /path/to/opfof.exe ...
   ```

---

### Issue: Slow I/O Performance

**Symptom:** Most time spent reading/writing

**Solutions:**

1. Use parallel filesystem:
   ```bash
   # On Lustre
   lfs setstripe -c 8 .
   ```

2. Aggregate small files:
   ```bash
   # Combine multiple small reads
   ```

3. Use buffered I/O:
   ```c
   setvbuf(fp, buffer, _IOFBF, 1024*1024);
   ```

---

## Performance Problems

### Issue: Code Runs Slowly

**Diagnosis:**
```bash
# Profile with perf
perf stat mpirun -np 1 ./opfof.exe 100 64

# Time different sections
# Add timing calls in code
```

**Common Causes and Solutions:**

1. **Unoptimized build:**
   ```makefile
   OPT = -O3 -xHost -ipo  # Enable optimizations
   ```

2. **Poor data locality:**
   - Sort particles by position before processing
   - Use cache-friendly access patterns

3. **Excessive tree depth:**
   - Adjust `NODE_HAVE_PARTICLE`
   - Balance tree after construction

---

### Issue: Load Imbalance

**Symptom:** Some ranks finish much earlier

**Diagnosis:**
```c
double my_time = elapsed_time();
double max_time, min_time;
MPI_Reduce(&my_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
MPI_Reduce(&my_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
if (myid == 0) printf("Imbalance: %.2f\n", max_time / min_time);
```

**Solutions:**

1. Dynamic load balancing
2. Redistribute particles
3. Use adaptive domain decomposition

---

## Result Validation

### Checking Output Integrity

```bash
# Basic check
./checkfof FoF_halo_cat.00100

# Expected output:
# Header: OK
# Halos: 12345
# Total particles: 67890123
```

### Comparing with Expected Values

```python
# compare_results.py
import numpy as np
from read_fof_catalog import read_fof_catalog

data = read_fof_catalog('FoF_halo_cat.00100')
halos = data['halos']

# Check basic statistics
print(f"Number of halos: {len(halos)}")
print(f"Total mass: {halos['mass'].sum():.3e} Msun/h")
print(f"Max halo mass: {halos['mass'].max():.3e} Msun/h")

# Sanity checks
assert len(halos) > 0, "No halos found!"
assert halos['mass'].min() > 0, "Negative mass found!"
assert halos['np'].min() >= 30, "Halos below minimum size!"
```

### Validating Halo Properties

```python
# Check center of mass is inside box
box = data['header']['box_size']
assert np.all(halos['x'] >= 0) and np.all(halos['x'] < box)
assert np.all(halos['y'] >= 0) and np.all(halos['y'] < box)
assert np.all(halos['z'] >= 0) and np.all(halos['z'] < box)

# Check component masses sum to total
for h in halos:
    m_sum = h['mdm'] + h['mstar'] + h['mgas'] + h['msink']
    assert abs(h['mass'] - m_sum) / h['mass'] < 0.01, "Mass mismatch!"

# Check particle counts
for h in halos:
    n_sum = h['npdm'] + h['npstar'] + h['npgas'] + h['npsink']
    assert h['np'] == n_sum, "Particle count mismatch!"
```

---

## Debugging Techniques

### Enable Debug Output

```makefile
# In Rules.make
OPT = -DINTEL -g -O0 -DDEBUG -DVERBOSE
```

### Use GDB with MPI

```bash
# Single rank debugging
mpirun -np 1 gdb ./opfof.exe

# Multi-rank with xterm
mpirun -np 4 xterm -e gdb ./opfof.exe

# Attach to running process
gdb -p <pid>
```

### Add Checkpoints

```c
// Add to critical sections
void checkpoint(const char *msg) {
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    printf("[Rank %d] Checkpoint: %s\n", myid, msg);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
}

// Usage
checkpoint("Before tree building");
FoF_Make_Tree();
checkpoint("After tree building");
```

### Memory Debugging with Valgrind

```bash
# Single process
valgrind --leak-check=full ./opfof.exe 100 64

# With MPI (limited utility)
mpirun -np 2 valgrind ./opfof.exe 100 64
```

### Core Dump Analysis

```bash
# Enable core dumps
ulimit -c unlimited

# Run until crash
mpirun -np 8 ./opfof.exe 100 64

# Analyze core
gdb ./opfof.exe core.12345
(gdb) bt    # Backtrace
(gdb) info locals
```

---

## FAQ

### Q: How many MPI ranks should I use?

**A:** Guidelines:
- One rank per 1-10 million particles
- Match ranks to available files (nfiles divisible by nranks)
- Consider memory per rank (need sufficient for particles + tree + buffers)

```bash
# Example: 1 billion particles, 512 files
# 64-512 ranks is reasonable
mpirun -np 128 ./opfof.exe 100 512
```

---

### Q: Why are some halos missing?

**A:** Possible reasons:
1. Below minimum particle threshold (default: 30)
2. Split across domain boundaries (check boundary handling)
3. Input data incomplete

Check:
```python
# Count halos above different thresholds
for thresh in [10, 30, 50, 100]:
    n = np.sum(halos['np'] >= thresh)
    print(f"Halos with Np >= {thresh}: {n}")
```

---

### Q: Why does the mass function look wrong?

**A:** Common issues:
1. Wrong linking length
2. Minimum halo size too high
3. Box size incorrect in calculation
4. Units mismatch

Verify:
```python
# Check cosmology
print(f"Box: {header['box_size']} Mpc/h")
print(f"Omega_m: {header['omega_m']}")

# Check linking length in code
# b = 0.2 × (m / (2.7755e11 × Ω_m))^(1/3)
```

---

### Q: How do I reduce memory usage?

**A:** Options:
1. Use single precision (`-DXYZDBL` removed)
2. Process fewer files at once
3. Use more MPI ranks
4. Reduce tree depth (`NODE_HAVE_PARTICLE = 16`)

---

### Q: Can I restart a failed run?

**A:** Currently no built-in checkpointing. Workarounds:
1. Process snapshots individually
2. Skip completed snapshots:
   ```bash
   for snap in {100..200}; do
       if [ ! -f "FoF_halo_cat.${snap}" ]; then
           mpirun -np 8 ./opfof.exe $snap 128
       fi
   done
   ```

---

### Q: How do I know if results are correct?

**A:** Validation steps:
1. Run `checkfof` on output
2. Compare mass function with theory
3. Check total mass conservation
4. Verify halo positions are within box
5. Compare with other FoF codes on same data

---

## Getting Help

If issues persist:

1. **Check logs carefully** for first error
2. **Reproduce with minimal case** (fewer ranks, smaller data)
3. **Collect system information**:
   ```bash
   uname -a
   mpirun --version
   cat /proc/meminfo
   ```
4. **Create detailed bug report** with:
   - Exact error message
   - Steps to reproduce
   - Configuration (Rules.make)
   - Input data description
