# opFoF variant: in-memory component labels

This subdirectory builds a shared library that labels one supplied set of
particles. It is the linker used when a stellar catalogue is already in
memory and the same particles are grouped at many linking lengths. It does
not read a RAMSES snapshot and it does not replace `opfof.exe`.

`opfof_sb_labels` takes contiguous `x`, `y`, `z`, and `link02` arrays and
returns one component label per particle. Two particles are joined when their
separation is strictly smaller than the average of their linking lengths.
The walk lives in this directory's `Treewalk.fof.ordered.c`. Relative to the
parent tree walk, the pair test compares squared distances, and the link loop
stops on the last accepted particle instead of reading one past it.

Build with the Intel MPI compiler:

```sh
cd opFoF/variant
make
```

`libopfof_sb.so` is a local build product and is not part of the source tree.
The entry point is `opfof_sb_labels` in `sb_adapter.c`.
