#ifndef OPFOF_DISTRIBUTED_H
#define OPFOF_DISTRIBUTED_H

#include "fof.h"

/* Resolve FoF components across the one-dimensional MPI slab decomposition
 * and write each globally owned halo exactly once. */
void DistributedMergeAndWrite(
		FoFTPtlStruct *ptl, size_t np, size_t nhalo,
		FoFTStruct **tree, long long *tree_capacity, Box box,
		POSTYPE boundary_link, POSTYPE domain_low, POSTYPE domain_high,
		POSTYPE input_fof_link, int periodic,
		POSTYPE box_x, POSTYPE box_y, POSTYPE box_z,
		const char *halofile, const char *memberfile);

#endif
