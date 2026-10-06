#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "../ramses.h"
#include "fof.h"

/*
 * The upstream translation unit also contains legacy catalog writers, whose
 * unused code references this property routine. The adapter calls only the
 * tree construction/linking API and deliberately does not expose those
 * writers; this inert compatibility symbol prevents loading their dead-code
 * dependency into the Python process.
 */
__attribute__((visibility("hidden")))
HaloQ haloproperty(particle *member, size_t nmem)
{
    HaloQ unused = {0};
    (void)member;
    (void)nmem;
    return unused;
}

/*
 * In-memory entry point for the existing opFoF tree walker.
 * GALFIND star records are decoded by the Python caller; only the four
 * FoF-relevant arrays cross this ABI, so the 120-byte GALFIND record is never
 * mistaken for this build's StarType layout.
 */
int opfof_sb_labels(size_t n, const double *x, const double *y,
                    const double *z, const double *link,
                    uint64_t *labels)
{
    size_t i, ncomp = 0, tree_capacity;
    double xmin, xmax, ymin, ymax, zmin, zmax, span, width;
    FoFTPtlStruct *ptl = NULL;
    FoFTStruct *tree = NULL;
    particle *linked = NULL;
    Box box;

    if (n == 0 || x == NULL || y == NULL || z == NULL ||
        link == NULL || labels == NULL) return 1;
    if (n > SIZE_MAX / sizeof(*ptl) || n > SIZE_MAX / sizeof(*linked)) return 2;

    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];
    for (i = 0; i < n; ++i) {
        if (!isfinite(x[i]) || !isfinite(y[i]) || !isfinite(z[i]) ||
            !isfinite(link[i]) || link[i] < 0.0) return 3;
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
        if (y[i] < ymin) ymin = y[i];
        if (y[i] > ymax) ymax = y[i];
        if (z[i] < zmin) zmin = z[i];
        if (z[i] > zmax) zmax = z[i];
    }

    span = fmax(xmax - xmin, fmax(ymax - ymin, zmax - zmin));
    if (span < 4.0 * MINCELLSIZE) span = 4.0 * MINCELLSIZE;
    width = exp2(ceil(log2(span * (1.0 + 1.0e-10))));
    box.width = width;
    box.x = 0.5 * (xmin + xmax - width);
    box.y = 0.5 * (ymin + ymax - width);
    box.z = 0.5 * (zmin + zmax - width);

    ptl = calloc(n, sizeof(*ptl));
    linked = malloc(n * sizeof(*linked));
    tree_capacity = (n / NODE_HAVE_PARTICLE) * 10 + 1024;
    if (tree_capacity > SIZE_MAX / sizeof(*tree)) goto allocation_error;
    tree = calloc(tree_capacity, sizeof(*tree));
    if (ptl == NULL || linked == NULL || tree == NULL) goto allocation_error;

    for (i = 0; i < n; ++i) {
        ptl[i].type = TYPE_STAR;
        ptl[i].included = NO;
        ptl[i].x = x[i];
        ptl[i].y = y[i];
        ptl[i].z = z[i];
        ptl[i].link02 = link[i];
        ptl[i].haloindx = SIZE_MAX;
        ptl[i].sibling = (i + 1 < n) ? (void *)&ptl[i + 1] : NULL;
    }

    FoF_Make_Tree(tree, ptl, n, box);
    for (i = 0; i < n; ++i) {
        if (ptl[i].included == NO) {
            new_fof_link(&ptl[i], 0.0, tree, ptl, linked, ncomp);
            ++ncomp;
        }
    }
    for (i = 0; i < n; ++i) labels[i] = (uint64_t)ptl[i].haloindx;

    free(tree);
    free(linked);
    free(ptl);
    return 0;

allocation_error:
    free(tree);
    free(linked);
    free(ptl);
    return 4;
}
