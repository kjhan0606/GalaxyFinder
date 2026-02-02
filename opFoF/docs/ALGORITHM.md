# opFoF Algorithm Documentation

This document provides a comprehensive description of the algorithms implemented in opFoF, including the Friends-of-Friends method, tree construction, parallel domain decomposition, and halo property calculations.

---

## Table of Contents

1. [Friends-of-Friends Algorithm](#friends-of-friends-algorithm)
2. [Oct-Tree Construction](#oct-tree-construction)
3. [Tree Walking and Neighbor Search](#tree-walking-and-neighbor-search)
4. [Particle Linking Process](#particle-linking-process)
5. [Periodic Boundary Conditions](#periodic-boundary-conditions)
6. [MPI Domain Decomposition](#mpi-domain-decomposition)
7. [Boundary Halo Management](#boundary-halo-management)
8. [Halo Property Calculation](#halo-property-calculation)
9. [Complexity Analysis](#complexity-analysis)
10. [Algorithm Variants](#algorithm-variants)

---

## Friends-of-Friends Algorithm

### Theoretical Foundation

The Friends-of-Friends (FoF) algorithm, introduced by Davis et al. (1985), is a percolation algorithm for identifying gravitationally bound structures in N-body simulations. The fundamental principle is simple:

> Two particles belong to the same group (halo) if they are separated by less than the **linking length** `b`.

This creates a transitive relationship: if particle A is linked to B, and B is linked to C, then A, B, and C all belong to the same halo, even if A and C are far apart.

### Linking Length Definition

The linking length is typically expressed as a fraction of the mean inter-particle separation:

```
b = b₀ × d̄
```

Where:
- `b₀` = dimensionless linking parameter (typically 0.2)
- `d̄` = mean inter-particle separation

The mean inter-particle separation is calculated from:

```
d̄ = (mp / ρ̄)^(1/3)
```

Where:
- `mp` = particle mass
- `ρ̄` = mean matter density = ρ_crit × Ω_m

### Physical Interpretation

The choice of `b₀ = 0.2` has physical significance:

1. **Overdensity Connection**: Groups identified with b₀ = 0.2 correspond approximately to regions with mean interior overdensity δ ≈ 180 (close to the virial overdensity)

2. **Percolation Threshold**: For b₀ < 0.2, halos are under-linked (fragmented). For b₀ > 0.2, halos become over-linked (merged)

3. **Universal Applicability**: This value works across different cosmologies and redshifts

### Mathematical Formulation

For N particles with positions **r**ᵢ, the FoF algorithm defines groups G_k such that:

```
∀ i,j ∈ G_k : ∃ sequence (i = n₀, n₁, ..., n_m = j) where |r_{n_l} - r_{n_{l+1}}| < b
```

This means any two particles in the same group are connected by a chain of links, each shorter than b.

---

## Oct-Tree Construction

### Purpose

Direct particle-particle distance calculations would require O(N²) operations. The oct-tree reduces this to O(N log N) by organizing particles hierarchically in space.

### Tree Structure

```c
typedef struct FoFTStruct {
    POSTYPE dist;                    // Maximum extent of node
    POSTYPE monox, monoy, monoz;     // Geometric center
    void *daughter;                  // First child node
    void *sibling;                   // Next sibling node
    int Nparticle;                   // Particle count in subtree
} FoFTStruct;
```

### Construction Algorithm

```
FUNCTION FoF_Make_Tree(particles[], N):
    root = CreateNode(particles, N)
    BuildTree(root)
    RETURN root

FUNCTION BuildTree(node):
    IF node.Nparticle <= NODE_HAVE_PARTICLE (=8):
        // Leaf node: store particle pointers
        node.daughter = particle_list
        RETURN

    // Calculate node bounds
    ComputeBoundingBox(node)

    // Partition particles into 8 octants
    FOR each octant (i,j,k) in {0,1}³:
        subset = particles in octant
        IF subset not empty:
            child = CreateNode(subset)
            AddChild(node, child)
            BuildTree(child)  // Recursive
```

### Octant Division

The space is divided into 8 octants based on the center point:

```
Octant 0: x < cx, y < cy, z < cz  (---)
Octant 1: x >= cx, y < cy, z < cz (+--)
Octant 2: x < cx, y >= cy, z < cz (-+-)
Octant 3: x >= cx, y >= cy, z < cz (++-)
Octant 4: x < cx, y < cy, z >= cz (--+)
Octant 5: x >= cx, y < cy, z >= cz (+-+)
Octant 6: x < cx, y >= cy, z >= cz (-++)
Octant 7: x >= cx, y >= cy, z >= cz (+++)
```

### Node Properties

Each node maintains:

| Property | Description | Use |
|----------|-------------|-----|
| `dist` | Maximum distance from center to any contained particle | Pruning criterion |
| `monox,y,z` | Geometric center of contained particles | Distance calculations |
| `Nparticle` | Total particles in subtree | Statistics |

---

## Tree Walking and Neighbor Search

### Opening Criterion

When searching for neighbors of particle p within distance b, we can prune entire subtrees if they cannot contain any neighbors:

```
FUNCTION CanPrune(node, particle_p, link_length_b):
    d = distance(node.center, particle_p.position)
    RETURN d - node.dist > b
```

If the closest possible particle in the node (distance d - dist) is farther than b, skip the entire subtree.

### Tree Walk Algorithm

```
FUNCTION TreeWalk(node, particle_p, link_length_b, neighbors[]):
    IF CanPrune(node, particle_p, b):
        RETURN  // Skip this subtree

    IF IsLeaf(node):
        FOR each particle q in node:
            IF distance(p, q) < b AND q not linked:
                neighbors.append(q)
        RETURN

    // Recursively walk children
    FOR each child in node.children:
        TreeWalk(child, particle_p, b, neighbors)
```

### Distance Calculation

For non-periodic boundaries:
```c
dx = p.x - q.x;
dy = p.y - q.y;
dz = p.z - q.z;
dist2 = dx*dx + dy*dy + dz*dz;
```

For periodic boundaries (box size L):
```c
dx = p.x - q.x;
if (dx > L/2) dx -= L;
if (dx < -L/2) dx += L;
// Similar for dy, dz
dist2 = dx*dx + dy*dy + dz*dz;
```

---

## Particle Linking Process

### Main Linking Algorithm

The core FoF algorithm implemented in `new_fof_link()`:

```
FUNCTION new_fof_link(particles[], N, link_length):
    halo_index = 0

    FOR each particle p in particles:
        IF p.haloindx != UNLINKED:
            CONTINUE  // Already assigned

        // Start new halo
        halo_index++
        linked[] = [p]
        p.haloindx = halo_index

        // Grow halo by finding all connected particles
        i = 0
        WHILE i < linked.size():
            current = linked[i]
            neighbors = TreeWalk(tree, current, link_length)

            FOR each neighbor n in neighbors:
                IF n.haloindx == UNLINKED:
                    n.haloindx = halo_index
                    linked.append(n)

            i++

        // linked[] now contains all particles in this halo
        ProcessHalo(linked[], halo_index)

    RETURN halo_index  // Total number of halos
```

### Linked Array Management

The `linked[]` array serves as a queue for breadth-first search:

```
Initial:    [seed_particle]
            ^ process this

After 1st:  [seed, neighbor1, neighbor2, neighbor3]
                  ^ process this

After 2nd:  [seed, n1, n2, n3, n4, n5]
                      ^ process this
...
Final:      [all particles in halo]
```

### Halo Index Assignment

Each particle stores its halo membership:

```c
typedef struct FoFTPtlStruct {
    ...
    size_t haloindx;  // 0 = unlinked, >0 = halo ID
} FoFTPtlStruct;
```

---

## Periodic Boundary Conditions

### The Wrapping Problem

In periodic simulations, a halo near the box edge may have members that wrap around:

```
Box [0, L]:
         |.....*****|          (particles near right edge)
|***....            |          (particles near left edge, same halo)
```

### Periodic Distance Calculation

```c
FUNCTION periodic_distance(p, q, BoxSize):
    dx = p.x - q.x
    dy = p.y - q.y
    dz = p.z - q.z

    // Minimum image convention
    IF dx > BoxSize/2:  dx -= BoxSize
    IF dx < -BoxSize/2: dx += BoxSize
    IF dy > BoxSize/2:  dy -= BoxSize
    IF dy < -BoxSize/2: dy += BoxSize
    IF dz > BoxSize/2:  dz -= BoxSize
    IF dz < -BoxSize/2: dz += BoxSize

    RETURN sqrt(dx² + dy² + dz²)
```

### Periodic Tree Walk

The `pnew_fof_link()` function handles periodic boundaries:

```
FUNCTION pnew_fof_link(particles[], N, link_length, BoxSize):
    // Same as new_fof_link but with periodic distance
    // Must also check ghost images near boundaries

    FOR each particle p:
        // Check real neighbors
        neighbors = TreeWalk_Periodic(tree, p, link_length, BoxSize)

        // For particles near boundaries, also check wrapped positions
        IF p.x < link_length OR p.x > BoxSize - link_length:
            // Check across x-boundary
            ...
```

### Center of Mass for Wrapped Halos

Special handling required for halos spanning the periodic boundary:

```c
FUNCTION compute_com_periodic(particles[], BoxSize):
    // Use angular coordinates to avoid discontinuity
    theta_x = 2π × x / BoxSize

    // Compute mean angle
    sum_cos = Σ cos(theta_x)
    sum_sin = Σ sin(theta_x)
    mean_theta = atan2(sum_sin, sum_cos)

    // Convert back to position
    x_com = mean_theta × BoxSize / 2π
    IF x_com < 0: x_com += BoxSize
```

---

## MPI Domain Decomposition

### Slab-Based Decomposition

The simulation volume is divided into Z-slabs, one per MPI rank:

```
Rank 0: z ∈ [0, L/nprocs)
Rank 1: z ∈ [L/nprocs, 2L/nprocs)
...
Rank n-1: z ∈ [(n-1)L/nprocs, L)
```

### File Distribution

Input files are distributed among ranks:

```c
initfile = (nfile / nid) * myid;
finalfile = (nfile / nid) * (myid + 1);

// Rank myid reads files [initfile, finalfile)
```

### Domain Boundaries

Each rank maintains:

```c
float zmin;  // Lower Z boundary
float zmax;  // Upper Z boundary
```

### Particle Assignment

```
FUNCTION assign_particles_to_ranks(particle):
    rank = floor(particle.z / BoxSize * nprocs)
    IF rank == nprocs: rank = nprocs - 1  // Handle boundary
    RETURN rank
```

---

## Boundary Halo Management

### The Boundary Problem

Halos near domain boundaries may have particles on multiple ranks:

```
Rank 0      |  Rank 1
    ***     |  **
      **    |  *
        *   |
--------zmax|zmin--------
```

### Boundary Detection

```c
typedef struct HaloBound {
    size_t nmem;      // Number of particles
    POSTYPE zmin;     // Minimum z of halo
    POSTYPE zmax;     // Maximum z of halo
    int boundflag;    // Boundary status
} HaloBound;

// boundflag values:
// 0 = Internal halo (no boundary contact)
// 1 = Touches lower boundary
// 2 = Touches upper boundary
// 3 = Touches both boundaries
```

### Boundary Particle Exchange

The `WriteBottomFaceContact()` and `ReadBottomFaceContact()` functions handle this:

```
ALGORITHM: Boundary Halo Assembly

1. Each rank identifies boundary halos
   FOR each halo:
       IF halo.zmin < zmin + link_length:
           Mark as lower boundary halo
       IF halo.zmax > zmax - link_length:
           Mark as upper boundary halo

2. Exchange boundary particles
   Rank i sends lower boundary particles to Rank i-1
   Rank i receives upper boundary particles from Rank i+1

3. Merge boundary halos
   Match particles across boundaries
   Combine matching halos

4. Final output
   Lower rank writes merged boundary halos
```

### Large Message Handling

For halos with millions of particles, standard MPI may fail. opFoF uses chunked transfers:

```c
FUNCTION BIG_MPI_Send(buffer, count, datatype, dest, tag, comm):
    chunk_size = 1GB / sizeof(datatype)
    offset = 0

    WHILE offset < count:
        this_chunk = min(chunk_size, count - offset)
        MPI_Send(&buffer[offset], this_chunk, datatype, dest, tag, comm)
        offset += this_chunk
```

---

## Halo Property Calculation

### Center of Mass Position

```c
FUNCTION compute_com_position(halo):
    total_mass = 0
    com_x = com_y = com_z = 0

    FOR each particle p in halo:
        m = p.mass
        total_mass += m
        com_x += m * p.x
        com_y += m * p.y
        com_z += m * p.z

    com_x /= total_mass
    com_y /= total_mass
    com_z /= total_mass

    RETURN (com_x, com_y, com_z)
```

### Center of Mass Velocity

```c
FUNCTION compute_com_velocity(halo):
    total_mass = 0
    com_vx = com_vy = com_vz = 0

    FOR each particle p in halo:
        m = p.mass
        total_mass += m
        com_vx += m * p.vx
        com_vy += m * p.vy
        com_vz += m * p.vz

    com_vx /= total_mass
    com_vy /= total_mass
    com_vz /= total_mass

    RETURN (com_vx, com_vy, com_vz)
```

### Component Mass Calculation

```c
FUNCTION compute_component_masses(halo):
    m_dm = m_star = m_gas = m_sink = 0
    n_dm = n_star = n_gas = n_sink = 0

    FOR each particle p in halo:
        SWITCH p.type:
            CASE DM (2):
                m_dm += p.mass
                n_dm++
            CASE STAR (1):
                m_star += p.mass
                n_star++
            CASE GAS (3):
                m_gas += p.mass
                n_gas++
            CASE SINK (4):
                m_sink += p.mass
                n_sink++

    RETURN (m_dm, m_star, m_gas, m_sink, n_dm, n_star, n_gas, n_sink)
```

### Velocity Dispersion (Optional)

```c
FUNCTION compute_velocity_dispersion(halo, com_velocity):
    sigma2 = 0
    N = 0

    FOR each particle p in halo:
        dvx = p.vx - com_velocity.vx
        dvy = p.vy - com_velocity.vy
        dvz = p.vz - com_velocity.vz
        sigma2 += dvx² + dvy² + dvz²
        N++

    sigma = sqrt(sigma2 / N / 3)  // 1D velocity dispersion
    RETURN sigma
```

---

## Complexity Analysis

### Time Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Tree Construction | O(N log N) | Recursive partitioning |
| Tree Walk (per particle) | O(log N) average | Depends on clustering |
| Full FoF Linking | O(N log N) | Amortized over all particles |
| Halo Property Calculation | O(N) | Linear scan of particles |

### Space Complexity

| Component | Memory | Notes |
|-----------|--------|-------|
| Particle Array | O(N) | ~100 bytes per particle |
| Tree Nodes | O(N) | ~48 bytes per node |
| Linked Array | O(N_max_halo) | Largest halo size |
| Boundary Buffer | O(N_boundary) | Boundary particles |

### Parallel Scaling

**Strong Scaling** (fixed problem size):
- Ideal: T(P) = T(1) / P
- Actual: T(P) = T(1) / P + T_comm + T_boundary
- Efficiency degrades when communication dominates

**Weak Scaling** (fixed particles per processor):
- Ideal: T(P) = T(1)
- Actual: T(P) ≈ T(1) + O(log P) for global operations

---

## Algorithm Variants

### Standard FoF (Implemented)

- Links all particles within linking length
- No distinction between particle types for linking
- Type-aware for property calculation

### 6D FoF (Not Implemented)

Links particles based on both position and velocity:
```
d6D² = dx² + dy² + dz² + (dvx² + dvy² + dvz²)/σ_v²
```

### Hierarchical FoF (Not Implemented)

Runs FoF at multiple linking lengths to identify substructure:
```
b₁ = 0.2 (main halos)
b₂ = 0.1 (subhalos)
b₃ = 0.05 (sub-subhalos)
```

### Adaptive FoF (Not Implemented)

Uses local mean separation instead of global:
```
b(r) = b₀ × d̄(r)
```

Where d̄(r) varies with local density.

---

## References

1. Davis, M., Efstathiou, G., Frenk, C. S., & White, S. D. M. (1985). "The evolution of large-scale structure in a universe dominated by cold dark matter." ApJ, 292, 371-394.

2. Springel, V., White, S. D. M., Tormen, G., & Kauffmann, G. (2001). "Populating a cluster of galaxies - I. Results at z=0." MNRAS, 328, 726-750.

3. More, S., Kravtsov, A. V., Dalal, N., & Gottlöber, S. (2011). "The overdensity and masses of the friends-of-friends halos and universality of halo mass function." ApJS, 195, 4.

4. Knebe, A., et al. (2011). "Haloes gone MAD: The Halo-Finder Comparison Project." MNRAS, 415, 2293-2318.
