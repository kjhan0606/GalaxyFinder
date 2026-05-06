# PGalF Algorithm Details

This document provides in-depth technical details of the galaxy finding algorithms implemented in PGalF.

## Table of Contents

1. [Density Estimation](#1-density-estimation)
2. [Peak Finding](#2-peak-finding)
3. [Water-Shedding Algorithm](#3-water-shedding-algorithm)
4. [Boundedness Determination](#4-boundedness-determination)
5. [Tidal Radius Calculation](#5-tidal-radius-calculation)
6. [Final Membership Assignment](#6-final-membership-assignment)
7. [Periodic Boundary Handling](#7-periodic-boundary-handling)
8. [Dark-Galaxy Detection (Optional)](#8-dark-galaxy-detection-optional)

---

## 1. Density Estimation

### 1.1 Triangular-Shaped Cloud (TSC) Interpolation

The TSC method distributes each particle's mass to nearby grid cells using a piecewise quadratic kernel:

```
W(x) = {
    0.75 - x²           for |x| ≤ 0.5
    0.5 × (1.5 - |x|)²  for 0.5 < |x| ≤ 1.5
    0                   otherwise
}
```

**Implementation** (`tsc_omp2.c`):

```c
// For particle at position (px, py, pz):
for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
        for (k = -1; k <= 1; k++) {
            ix = (int)(px/cellsize) + i;
            iy = (int)(py/cellsize) + j;
            iz = (int)(pz/cellsize) + k;

            wx = TSC_weight((px - ix*cellsize) / cellsize);
            wy = TSC_weight((py - iy*cellsize) / cellsize);
            wz = TSC_weight((pz - iz*cellsize) / cellsize);

            grid[ix][iy][iz] += mass * wx * wy * wz;
        }
    }
}
```

### 1.2 Gaussian Smoothing

After TSC interpolation, the density field is smoothed using FFT-based Gaussian convolution:

```
ρ_smooth(r) = ∫ ρ(r') × G(|r - r'|) d³r'

G(r) = (2πσ²)^(-3/2) × exp(-r²/(2σ²))
```

Where σ = `Gaussian_Smoothing_Length` (default: 0.008 cMpc/h)

**Implementation** (`gsmooth.c`):

```c
// In Fourier space, Gaussian becomes:
// G(k) = exp(-k²σ²/2)

fftwf_execute(forward_plan);   // ρ(r) → ρ(k)

for (all k-modes) {
    k2 = kx² + ky² + kz²;
    rho_k[i] *= exp(-k2 * sigma² / 2);
}

fftwf_execute(backward_plan);  // ρ(k) → ρ_smooth(r)
```

### 1.3 Particle Density Assignment

After smoothing, particle densities are read back from the grid:

```c
for (each particle i) {
    // Trilinear interpolation from grid to particle
    wp[i].den = interpolate_density(grid, p[i].x, p[i].y, p[i].z);
}
```

---

## 2. Peak Finding

### 2.1 k-d Tree Construction

For efficient nearest neighbor queries, a k-d tree is built from particle positions:

**Implementation** (`nnost.c`):

```c
typedef struct NNTreeNode {
    int split_dim;        // 0=x, 1=y, 2=z
    float split_value;
    struct NNTreeNode *left, *right;
    int particle_indices[MAX_LEAF_SIZE];
    int n_particles;
} NNTreeNode;

NNTreeNode* build_kd_tree(particles, start, end, depth) {
    if (end - start <= MAX_LEAF_SIZE) {
        return create_leaf(particles, start, end);
    }

    dim = depth % 3;
    median = find_median(particles, start, end, dim);

    node->split_dim = dim;
    node->split_value = median;
    node->left = build_kd_tree(particles, start, mid, depth+1);
    node->right = build_kd_tree(particles, mid, end, depth+1);

    return node;
}
```

### 2.2 Nearest Neighbor Search

Find the `NUMNEIGHBOR` (32) nearest particles:

```c
void find_k_nearest(particle p, int k, int *neighbors) {
    // Priority queue to track k nearest
    PriorityQueue pq;

    search_tree(root, p, k, &pq);

    // Return indices of k nearest particles
    copy_from_queue(pq, neighbors);
}
```

### 2.3 Peak Identification

A particle is identified as a density peak if:

1. Its density exceeds `PEAKTHRESHOLD` (2000 h² Msun/ckpc³)
2. It has the highest density among its k nearest neighbors

**Implementation** (`subhaloden.mod6.c`):

```c
for (i = 0; i < np; i++) {
    if (wp[i].den < PEAKTHRESHOLD) continue;

    find_k_nearest(p[i], NUMNEIGHBOR, neighbors);

    is_peak = 1;
    for (j = 0; j < NUMNEIGHBOR; j++) {
        if (wp[neighbors[j]].den > wp[i].den) {
            is_peak = 0;
            break;
        }
    }

    if (is_peak) {
        SET_PEAK(i);
        core[ncore].peak = i;
        core[ncore].density = wp[i].den;
        ncore++;
    }
}
```

---

## 3. Water-Shedding Algorithm

### 3.1 Concept

The water-shedding algorithm assigns each particle to a peak by following the steepest density gradient uphill, analogous to water flowing downhill and collecting in basins.

### 3.2 Algorithm

```c
int find_watershed_peak(int particle_idx) {
    int current = particle_idx;

    while (!IS_PEAK(current)) {
        // Find neighbor with highest density
        find_k_nearest(p[current], NUMNEIGHBOR, neighbors);

        int max_neighbor = current;
        float max_den = wp[current].den;

        for (j = 0; j < NUMNEIGHBOR; j++) {
            if (wp[neighbors[j]].den > max_den) {
                max_den = wp[neighbors[j]].den;
                max_neighbor = neighbors[j];
            }
        }

        if (max_neighbor == current) {
            // Local maximum found (new peak or saddle)
            break;
        }

        current = max_neighbor;
        SET_VISITED(current);
    }

    return find_peak_for(current);  // Return peak index
}
```

### 3.3 Recursive Path Compression

To speed up subsequent lookups, paths are compressed:

```c
int watershed_with_compression(int idx) {
    if (IS_PEAK(idx)) return idx;
    if (wp[idx].haloid >= 0) return wp[idx].haloid;  // Already assigned

    int peak = watershed_with_compression(get_max_density_neighbor(idx));
    wp[idx].haloid = peak;  // Cache the result

    return peak;
}
```

### 3.4 OpenMP Parallelization

The water-shedding is parallelized across particles:

```c
#pragma omp parallel for schedule(dynamic, 1024)
for (i = 0; i < np; i++) {
    int thread_id = omp_get_thread_num();

    // Use thread-local visited flags
    int peak = watershed_parallel(i, thread_id);
    wp[i].haloid = peak;
}
```

---

## 4. Boundedness Determination

### 4.1 Energy Criterion

A particle is considered bound if its total energy in the center-of-mass frame is negative:

```
E_total = KE + PE < 0  ⟹  Bound
```

Where:
- KE = ½m|v - v_CoM|² (kinetic energy)
- PE = -Gm × Σ(m_j / r_ij) (potential energy)

### 4.2 Iterative Unbinding

The boundedness test is applied iteratively to remove unbound particles:

```c
for (iter = 0; iter < BOUNDITER; iter++) {
    // Calculate center-of-mass position and velocity
    calc_center_of_mass(core, &x_cm, &y_cm, &z_cm, &vx_cm, &vy_cm, &vz_cm);

    int n_removed = 0;

    for (each particle i in core) {
        // Kinetic energy in CoM frame
        double dvx = vx[i] - vx_cm;
        double dvy = vy[i] - vy_cm;
        double dvz = vz[i] - vz_cm;
        double KE = 0.5 * mass[i] * (dvx² + dvy² + dvz²) * r2kineticfact²;

        // Potential energy (tree-based calculation)
        double PE = calc_potential(i, core_particles) * potentfact;

        if (KE + PE > 0) {
            // Particle is unbound
            UNSET_CORE(i);
            SET_SHELL(i);
            n_removed++;
        }
    }

    if (n_removed == 0) break;  // Converged
}
```

### 4.3 Potential Calculation

The gravitational potential is calculated using a spline-softened kernel:

```c
double calc_potential(int i, int *members, int nmem) {
    double potential = 0;

    for (j = 0; j < nmem; j++) {
        if (i == members[j]) continue;

        double dx = x[i] - x[members[j]];
        double dy = y[i] - y[members[j]];
        double dz = z[i] - z[members[j]];
        double r = sqrt(dx² + dy² + dz²);

        // Spline softening to avoid singularities
        if (r < epsilon) {
            potential += mass[j] * spline_potential(r/epsilon) / epsilon;
        } else {
            potential += mass[j] / r;
        }
    }

    return -G * potential;
}
```

The spline kernel (`force_spline.mod3.c`) provides smooth forces at small separations.

---

## 5. Tidal Radius Calculation

### 5.1 NFW Profile Model

The tidal radius is calculated assuming an NFW density profile for the host halo:

```
ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
```

Where:
- r_s = R_vir / c (scale radius)
- c = concentration parameter
- R_vir = virial radius

### 5.2 Tidal Stripping Criterion

A subhalo is tidally stripped when the host's tidal force exceeds the subhalo's self-gravity:

```
F_tidal > F_self
GM(r)/r² × (2r_t/r) > Gm/r_t²
```

Solving for the tidal radius r_t:

```
r_t = r × (m / (2M(r)))^(1/3)
```

### 5.3 Lookup Table Implementation

To avoid expensive calculations, a 3D lookup table is precomputed:

```c
// mkRtidal.c
float nfw_rtidal(float mM, float dRv, float c) {
    // Inputs:
    //   mM  = m/M (subhalo mass / host mass ratio)
    //   dRv = d/R_vir (distance / virial radius)
    //   c   = concentration parameter

    // Table dimensions: 128 × 128 × 32
    int i = (log10(mM) + 4.0) / 4.0 * Nm;   // Mass ratio index
    int j = (log10(dRv) + 3.0) / 4.0 * Nd;  // Distance index
    int k = (c - 1.0) / 32.0 * Nc;          // Concentration index

    // Trilinear interpolation
    return interpolate_3d(tidal_table, i, j, k);
}
```

### 5.4 Core Application

```c
void GetTidalRCenterCore(Coretype *core, int ncore,
                          float host_mass, float host_Rvir, float host_c) {
    for (i = 0; i < ncore; i++) {
        // Distance from core to host center
        float dx = core[i].cx - host_cx;
        float dy = core[i].cy - host_cy;
        float dz = core[i].cz - host_cz;
        float dist = sqrt(dx² + dy² + dz²);

        // Mass ratio
        float mM = core[i].mass / host_mass;

        // Distance ratio
        float dRv = dist / host_Rvir;

        // Get tidal radius
        core[i].Rtidal = dist * pow(10, nfw_rtidal(mM, dRv, host_c));

        // Clamp to maximum
        if (core[i].Rtidal > MAX_TIDAL_R) {
            core[i].Rtidal = MAX_TIDAL_R;
        }
    }
}
```

---

## 6. Final Membership Assignment

### 6.1 Shell Division

Non-core particles are divided into `NSHELLDIVIDE` (10) iso-density shells:

```c
// Sort particles by density
qsort(shell_particles, n_shell, sizeof(int), compare_density);

// Divide into shells
int shell_size = n_shell / NSHELLDIVIDE;
for (i = 0; i < NSHELLDIVIDE; i++) {
    shell_start[i] = i * shell_size;
    shell_end[i] = (i + 1) * shell_size;
}
```

### 6.2 FoF Linking for Membership

Particles in each shell are linked to cores using FoF:

```c
void assign_shell_membership(int shell_idx) {
    float link_length = FOFLINK4MEMBERSHIP;  // 0.005 cMpc/h

    for (each particle i in shell) {
        if (wp[i].haloid >= 0) continue;  // Already assigned

        // Find nearby core particles
        int nearest_core = -1;
        float min_dist = link_length;

        for (each core c) {
            for (each particle j in core c) {
                float dist = distance(i, j);
                if (dist < min_dist) {
                    min_dist = dist;
                    nearest_core = c;
                }
            }
        }

        if (nearest_core >= 0) {
            wp[i].haloid = nearest_core;
        }
    }
}
```

### 6.3 Tree-Based Optimization

The FoF linking uses a Barnes-Hut tree for efficiency:

```c
// ost.c
int new_fof_link(particle *p, float link_length,
                  FoFTStruct *tree, FoFTPtlStruct *particles) {

    if (tree->type == TYPE_TREE) {
        // Check if tree node can contain linked particles
        if (tree->minlink02 > link_length²) return 0;

        // Recursively search daughter nodes
        return new_fof_link(p, link_length, tree->daughter, particles);
    }
    else {
        // Leaf node - check particle
        float dx = p->x - particles->x;
        float dy = p->y - particles->y;
        float dz = p->z - particles->z;
        float r2 = dx² + dy² + dz²;

        if (r2 < link_length²) {
            // Link particles
            return 1;
        }
        return 0;
    }
}
```

---

## 7. Periodic Boundary Handling

### 7.1 Detection

Halos crossing periodic boundaries are detected by checking coordinate ranges:

```c
void mklocalize(particles, np, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax) {
    // Check if halo spans the boundary
    if (xmin < boundary_buffer && xmax > size - boundary_buffer) {
        // Halo crosses X boundary
        handle_periodic_x(particles, np);
    }
    // Similar for Y and Z
}
```

### 7.2 Coordinate Shifting

Particles are shifted to create a continuous distribution:

```c
void handle_periodic_x(particles, np) {
    // Find the largest gap in X coordinates
    sort_by_x(particles, np);

    float max_gap = 0;
    int gap_idx = 0;

    for (i = 0; i < np - 1; i++) {
        float gap = particles[i+1].x - particles[i].x;
        if (gap > max_gap) {
            max_gap = gap;
            gap_idx = i;
        }
    }

    // Shift particles on the low side across the boundary
    for (i = 0; i <= gap_idx; i++) {
        particles[i].x += size;
    }
}
```

### 7.3 Final Coordinate Normalization

After processing, coordinates are normalized back to [0, size]:

```c
for (i = 0; i < np; i++) {
    while (p[i].x >= size) p[i].x -= size;
    while (p[i].x < 0)     p[i].x += size;
    // Similar for y, z
}
```

---

## 8. Dark-Galaxy Detection (Optional)

Enabled by compiling with `-DDARK_GAL`. The goal is to identify *dark
galaxies* — gravitationally bound DM-only subhalos that host no resolved
stellar component — using exactly the same watershed pipeline that finds
luminous galaxies, with no separate dual-pass code path.

### 8.1 Motivation

A pure DM-only second pass (build a DM density grid, find DM peaks, dedup
against stellar peaks, append) was implemented and rejected because the
shared watershed shell loop bins by *stellar* density: a dark core sitting at
a DM peak has near-zero stellar density there, so the shell test
`wp[i].den > dthreshold` fails and the dark core grows zero members. Forcing
dark cores into the last shell created an asymmetric pipeline that was
fragile to tune and that broke the equality assumption between core types.

The unified weighted-density formulation below removes this asymmetry: a
single density field, a single peak set, a single watershed, and a single
shell loop. Dark vs. luminous classification happens *post-hoc* on the
already-formed cores.

### 8.2 Unified Density Field

```
ρ_total(x) = ρ_star(x) + w_DM × ρ_DM(x)        (DM_DENSITY_WEIGHT = w_DM)
```

Each component is built independently:

```c
// nnost.c :: lagFindCoreParam, ptype == TYPE_STAR_DM
assign_density_TSC(bp, np, denGrid,    nx,ny,nz, ..., TYPE_STAR);
gaussian_Smoothing(denGrid,    nx,ny,nz, cellsize, Gaussian_Smoothing_Length);

if (DM_DENSITY_WEIGHT != 0) {
    assign_density_TSC(bp, np, denGrid_dm, nx,ny,nz, ..., TYPE_DM);
    gaussian_Smoothing(denGrid_dm, nx,ny,nz, cellsize, DM_GAUSSIAN_SMOOTHING_LENGTH);
    for (ic = 0; ic < ncells; ic++)
        denGrid[ic] += DM_DENSITY_WEIGHT * denGrid_dm[ic];
}
```

Note:

- Stars and DM are smoothed **at different scales**. DM is dynamically
  hotter; smoothing it on the stellar scale produces noisy peaks at the
  granularity of individual heavy DM particles.
- `w_DM = 0.1` keeps stellar peaks dominant at galaxy centers while still
  letting pure-DM concentrations rise above `PEAKTHRESHOLD`. Numerically:
  - star-only peak (ρ_star = 100, ρ_DM = 0) → ρ_total = 100
  - galactic core (ρ_star = 100, ρ_DM = 1000) → ρ_total = 200
  - DM-only subhalo (ρ_star = 0, ρ_DM = 1000) → ρ_total = 100
- When `w_DM = 0` the DM TSC, DM smoothing, and the combine loop are all
  skipped — the path collapses to the stellar-only case.

### 8.3 FoF-Halo Gate

`subhalo_den()` first decides whether to use the grid path or fall back to
SPH-density peak finding. Under DARK_GAL the gate is OR-combined:

```c
int do_grid =
       ((nstar > NUMNEIGHBOR) && (mstar >= MINSTELLARMASS))   // luminous halo
    || (mdm   >= MINDMMASS);                                  // DM-rich halo
```

`MINDMMASS` defaults to `10 × MINSTELLARMASS`; halos below both thresholds
are unlikely to host substructure worth resolving and use the cheaper SPH
path.

### 8.4 Linked-List for Peak Resolution

The peak finder picks the densest particle in each peak cell as the core
center. Under `TYPE_STAR_DM` the linked-list grid includes both stars and
DM, so a peak in a DM-dominated cell resolves to a DM particle and a peak in
a stellar cell resolves to a star — which is exactly what the post-hoc
classification then reads.

### 8.5 Post-hoc Dark Classification

After `FindCoreDensity` populates `core[i].numstar` and `core[i].starmass`
(both filtered by `bp[j].type == TYPE_STAR`, so they are accurate even though
the watershed ran on combined density), each core is flagged:

```c
// subhaloden.mod6.c :: subhalo_den(), under #ifdef DARK_GAL
for (i = 0; i < numcore; i++) {
    if (core[i].numstar  < MINCORENMEM
     || core[i].starmass < MINSTELLARMASS) {
        core[i].is_dark = 1;
    }
}
```

`is_dark` cores then participate in the shell loop, boundedness test, and
final membership assignment on equal footing with luminous cores. Downstream
catalogs can be split by the `is_dark` flag stored in `Coretype`.

### 8.6 What Was Removed

The earlier dual-pass implementation in `subhalo_den()` — a separate DM TSC
pass calling `lagFindDarkCore`, dedup against stellar peaks via
`STAR_DM_DEDUP_LENGTH`, and append-with-`is_dark=1` — has been removed. The
`lagFindDarkCore` wrapper itself is retained in `nnost.c` for callers that
want a stand-alone DM-only peak finder, but `subhalo_den()` no longer
invokes it.

### 8.7 Code Touchpoints

| File | Change |
|------|--------|
| `params.h` | `DM_DENSITY_WEIGHT`, `DM_GAUSSIAN_SMOOTHING_LENGTH`, `DM_TSC_CELL_SIZE`, `MINDMMASS` (and legacy `DM_PEAKTHRESHOLD`, `DM_MERGINGPEAKLENGTH`, `DM_MINCORENMEM`, `STAR_DM_DEDUP_LENGTH`). |
| `tree.h` | `enum { ..., TYPE_STAR_DM = 6 }`; `Coretype.is_dark`. |
| `nnost.c` | `lagFindCoreParam` parameterizes thresholds and dispatches the unified path on `TYPE_STAR_DM`; new wrapper `lagFindTotalCore`; `lagFindDarkCore` retained for stand-alone use. |
| `subhaloden.mod6.c` | OR-combined FoF gate; `lagFindStellarCore` → `lagFindTotalCore` under DARK_GAL; post-hoc `is_dark` classification; `MergingPeak` accepts a per-call FoF-link length. |

### 8.8 Tuning Notes

- `DM_DENSITY_WEIGHT`: lower if pure-DM peaks crowd out luminous galaxies in
  cluster cores; raise if known DM-only subhalos are missed. The 0–1 range
  is the design envelope; values > ~0.3 will start to dominate over stars.
- `DM_GAUSSIAN_SMOOTHING_LENGTH`: tied to the DM mass resolution. Rule of
  thumb: roughly the inter-particle separation for a 50–100 particle clump
  at the lightest DM species in the simulation.
- `MINDMMASS`: serves only the FoF-halo gate, not per-core selection. Cores
  below this DM mass are still found; they are simply not allowed to
  *trigger* the grid path on a stellar-empty halo.

---

## Algorithm Complexity Summary

| Step | Time Complexity | Space Complexity |
|------|-----------------|------------------|
| TSC Interpolation | O(N) | O(N_grid³) |
| Gaussian Smoothing | O(N_grid³ log N_grid) | O(N_grid³) |
| k-d Tree Build | O(N log N) | O(N) |
| Peak Finding | O(N × k) | O(N) |
| Water-shedding | O(N × k × h) | O(N) |
| Boundedness | O(N_core × iter) | O(N_core) |
| Final FoF | O(N log N) | O(N) |

Where:
- N = total particles
- N_grid = grid cells per dimension
- k = number of neighbors (32)
- h = average path length to peak
- iter = BOUNDITER (4)
- N_core = particles per core

---

## Tuning Guidelines

### For High-Resolution Simulations

```c
#define TSC_CELL_SIZE 0.002           // Finer grid
#define Gaussian_Smoothing_Length 0.004  // Smaller smoothing
#define MINCORENMEM 50                // More particles required
#define PEAKTHRESHOLD 5000            // Higher density threshold
```

### For Low-Resolution Simulations

```c
#define TSC_CELL_SIZE 0.008           // Coarser grid
#define Gaussian_Smoothing_Length 0.016  // Larger smoothing
#define MINCORENMEM 20                // Fewer particles required
#define PEAKTHRESHOLD 1000            // Lower density threshold
```

### For Cluster-Scale Halos

```c
#define MERGINGPEAKLENGTH 0.01        // Larger merge distance
#define NSHELLDIVIDE 20               // More shells
#define BOUNDITER 6                   // More iterations
```
