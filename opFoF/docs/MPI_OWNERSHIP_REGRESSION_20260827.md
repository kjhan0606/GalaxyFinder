# MPI halo-ownership regression (2026-08-27)

## Contract

An opFoF catalogue is MPI-rank converged only when all tested decompositions
produce the same sets of particle IDs per halo. Matching total halo counts is
not sufficient. The additional hard gates are:

- no duplicate member ID;
- member-file length equals the sum of `HaloQ.np`;
- one final owner for every resolved global component; and
- identical membership digest across MPI decompositions.

## Distributed algorithm

1. Each rank runs FoF on all NewDD slabs in its contiguous domain.
2. Adjacent ranks exchange only particles within the global maximum linking
   length of their shared face.
3. Cross-face component edges are identified with the same variable-link
   pair criterion used by the main FoF walk.
4. Component labels propagate between adjacent ranks until a global
   no-change reduction is reached.
5. Components below 20 particles and absent from the boundary graph are
   discarded before member transport.
6. Remaining member records move one rank per hop around an adjacent-rank
   ring. A record leaves the ring at its label owner. Receives are posted
   before sends and large payloads are chunked.
7. Owner ranks write in rank order. No boundary-particle gather on rank 0 and
   no `MPI_Alltoall` or `MPI_Alltoallv` is used.

## Regression inputs and results

Both tests used the same 5 Mpc/h, 256-cubed, z=12 matched snapshots and 64
NewDD slabs. Runs were executed on `grammar-debug` with Intel MPI 2021.17.

| representation | ranks | haloes | members | duplicate IDs | maximum halo | membership SHA-256 |
|---|---:|---:|---:|---:|---:|---|
| one fluid | 8, 16, 32, 64 | 6,204 | 465,841 | 0 | 3,309 | `f93a9a2e794b9de1dcb4d99de6e40e307caaec8439de67576a05c2c6a9d7d306` |
| two fluid, variable link | 8, 16, 32, 64 | 8,634 | 1,040,984 | 0 | 7,311 | `abf28d313b154aa8939d7fb6876479b67a24541a33e35e453742dfc9804d1123` |

The digest is constructed by sorting member IDs within each halo, hashing
each halo, sorting the halo hashes, and hashing their concatenation. It is
therefore insensitive to catalogue and member ordering but sensitive to any
fragmentation, merger, omission, or duplicate ownership.

The former one-hop implementation split a 679-particle one-fluid halo into
571- and 108-particle catalogued objects at 32 ranks. The adjacent-rank label
path restores the 8/16-rank membership exactly at both 32 and 64 ranks.
