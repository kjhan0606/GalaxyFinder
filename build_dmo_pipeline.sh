#!/bin/sh
set -eu

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
mpi_cc=${MPI_CC:-mpicc}

common_flags="-g -O2 -DNBODY -DLONGINT -DOUTPUT_PARTICLE_POTENTIAL -DNENER=0 -DNPRE=8"
newdd_flags="$common_flags -DNMEG=20000 -DWGROUPSIZE=5 -DUSE_MPI -DDEBUG=1 -DLOG=1"
opfof_flags="$common_flags -DWGROUPSIZE=8 -DNMEG=17000L -DINCLUDE_TREE_FORCE -D_LARGE_FILES -DSAVESLICE -DPMSEEDFORCE"

make -C "$repo_dir/NewDD" clean
make -C "$repo_dir/NewDD" all CC="$mpi_cc" OPT="$newdd_flags"

make -C "$repo_dir/opFoF" clean
make -C "$repo_dir/opFoF" this \
	CC="$mpi_cc" \
	OPT="-g -O2" \
	CDFLAGS="$opfof_flags" \
	RAMLIBS="-L../NewDD -lmyram"

printf '%s\n' "Built NewDD/newdd.exe and opFoF/opfof.exe with $mpi_cc"
