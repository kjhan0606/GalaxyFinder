//*******************
// The number of interation to determine 
// the boundedness of member particles
#define BOUNDITER 4
//*******************
//


//*******************
// fractional error to measure 
// the core density in ratio
#define COREDENRESOLUTION (1.e-3) 
//*******************


//*******************
// lowest stellar-density peak
// in unit of h^2 Msun/ckpc^3
#define PEAKTHRESHOLD 50.L
//*******************

//--------------------------------
// please tune these six parameters 
// for higher resolution simulations 
//
//*******************
// the minimum number of star/dm
// particles to identify a core
#define MINCORENMEM 100
//*******************

//*******************
// the number of neighbors to build 
// the neighbor network 
// the bigger the better
// it should be no larger than MAX_NUM_NEAR in tree.h
#define NUMNEIGHBOR 15
//*******************

//*******************
// the number of iso-den division of non-core
// particles
#define NSHELLDIVIDE 10
//*******************

//*******************
// the separation limit of peaks 
// to merge in cMpc/h 
#define MERGINGPEAKLENGTH 3.e-3
//*******************

//*******************
// the minimum stellar mass for a core 
#define MINCORESTARMASS -1 
//*******************

//*******************
// the minimun stellar mass of the FoF halo 
// for galaxy finding with stellar density 
#define MINSTELLARMASS 1.e5
//*******************
//--------------------------------

//*******************
// Maximum number of cores 
#define MAXNUMCORE 1000000
//*******************


//*******************
// (obsolete) The number of nearby stars 
// to measure stellar density
#define NUMNEARDEN 10
//*******************

//*******************
// (obsolete) The number of stellar neighbors
// for the core detection to find density peaks 
#define NUMSTELLARNEIGHBORS 30 
//*******************


//*******************
// (obsolete) minimum smoothing length 
// for stellar density in cMpc/h
#define MIN_CONST_R_SMOOTHING 0.010
//*******************

//*******************
// The cellsize for TSC of the stellar density 
// in unit of cMpc/h
#define TSC_CELL_SIZE 0.003
//*******************

//*******************
// Gaussian Smoothing Length 
// for stellar density in unit of cMpc/h. 
// This should be not smaller 
// than 2*TSC_CELL_SIZE.
// Note that the smoothing may 
// lower the stellar density 
// and PEAKTHRESHOLD should be 
// lowered accordingly.
#define Gaussian_Smoothing_Length 0.006
//*******************
//
//*******************
// Nucleus radius (cMpc/h) for the self-bound rest frame of a satellite
// buried in a heavier core. Only stars that already belong to that
// satellite and lie within this radius of the watershed core center
// are used. The central galaxy keeps the all-candidate stellar COM.
// Fewer than NUCLEUS_MIN_STARS falls back to the watershed core
// velocity, not to the mixed candidate COM.
#define NUCLEUS_RADIUS 0.012
#define NUCLEUS_MIN_STARS 5
//*******************
//
//*******************
// mininum density
#define DENFLOOR 1.
//*******************
//
//
//
// the (even) number of cell boundary buffer. 
#define NCELLBUFF 10
//*******************
//
//
//
//
//
// the size of deep linked lists (omp parallelized)
// if the number of particles linked to a cell is smaller
// than DEEPSIZE, the omp parallelization is located outside the loop.
// In the other cases, the omp parallelization goes down the loop.
#define DEEPSIZE 1024
//*******************
//
//
// The lower Number of particles for the OMP-parallized FoF
#define NOMPFoF 500000
//*******************
//
//
//
// the linking length to finalize the membership (obsolete)
#define FOFLINK4MEMBERSHIP 0.005
//*******************
//
// the maximum number of linking in water shedding to find core density
// It should be sufficiently larger than 
// the maximum number of core particles.
#define MAXNUMWATERSHEDDING 100000000L
//*******************
//
//
// the maximum number of threads
#define MAXTHREADS 64

//*******************
// H-maxima / persistence-based pruning of watershed peaks.
//
// For each core c produced by FindCoreDensity:
//   persistence(c) = peak_density(c) - saddle_density(c)
// where saddle_density(c) is the density at which c's watershed basin
// first touches a peer peak (i.e., core[c].coredensity from the bisection).
// If persistence(c)/peak_density(c) < PERSISTENCE_TAU, c is folded into
// the higher-density peer it first touched.
//
// Suppresses over-segmentation in crowded fields (BCG / cluster centres)
// where small density fluctuations spawn spurious sub-cores. The density
// field is already gaussian_Smoothing'd at Gaussian_Smoothing_Length, so
// this is a second, topology-aware merge that doesn't blur real galaxies.
//
// Set 0 (or negative) to disable. Recommended: 0.20 - 0.40 for cluster
// runs; lower values merge less aggressively. Use 0 for byte-equivalence
// with the pre-pruning code path.
#define PERSISTENCE_TAU 0.30f
//*******************

//*******************
// Dormant-core skip optimization.
// After a core gets nmem==0 in DORMANT_EMPTY_STREAK consecutive shell
// iterations of the same FoF halo, it is marked dormant and skipped in
// the remaining shells of that halo. Dormant flag is reset per halo.
// Set to 0 to disable the skip (every core is examined in every shell).
// Recommended: 3-5. Lower values are more aggressive.
#define DORMANT_EMPTY_STREAK 0
//*******************

//*******************
// Sparse-shell fast-path optimization.
// In the shell-loop, when a shell has very few particles but many enclosed
// cores (e.g. outer iso-density bins of the BCG halo: np=3 with cores=6917),
// the per-core setup cost (halonmem/halomass O(np) scan, qsort, tidal-radius
// scan) dominates. If np_shell < FAST_SHELL_RATIO * nc_shell, take the
// fast path: assign all shell particles to the single most-massive core in
// the shell (skipping AdGetTidalRCenterCore + per-core boundedness). Bound
// particles get corrected in later boundedness iterations; mass loss to the
// max-mass core is negligible (< 1e-5 of halo mass). The last shell is
// never fast-pathed because it holds all unassigned leftovers via
// SaveRemainingParticles2LastShell. Set 0 to disable.
// Recommended: 0.30 (catches ~7% of shells in BCG runs).
#define FAST_SHELL_RATIO 0.10f
//*******************

//*******************
// Dark-galaxy (DM-only subhalo) finding parameters.
// Always built in; set DM_DENSITY_WEIGHT = 0 to disable (the DM TSC
// + Gaussian-smoothing path is then short-circuited and the run is
// byte-equivalent to the stellar-only finder).
//*******************

// the threshold for DM-density peaks
// in unit of h^2 Msun/ckpc^3.
// DM density per cell is typically much higher than stellar
// in absolute units, so set this larger than PEAKTHRESHOLD.
#define DM_PEAKTHRESHOLD 1.e3

// Gaussian smoothing length for the DM density field, cMpc/h.
// DM is dynamically hotter than stars, so smooth more aggressively.
#define DM_GAUSSIAN_SMOOTHING_LENGTH 0.012

// TSC cell size for the DM density field, cMpc/h.
// Default to the same grid as stars; set larger if memory is tight.
#define DM_TSC_CELL_SIZE TSC_CELL_SIZE

// merging length for DM-only peaks among themselves, cMpc/h.
#define DM_MERGINGPEAKLENGTH 5.e-3

// minimum number of DM particles required to keep a dark core
#define DM_MINCORENMEM 50

// dedup distance: a DM peak within this radius (cMpc/h) of any
// stellar peak is considered the DM core of an already-found
// luminous galaxy and is dropped.
#define STAR_DM_DEDUP_LENGTH 5.e-3

// minimum DM mass within a FoF halo to attempt dark-galaxy finding.
// Roughly 10x MINSTELLARMASS (since halos with so little DM are
// unlikely to host DM-only subhalos worth resolving).
#define MINDMMASS (10.*MINSTELLARMASS)

// weight of DM density when building the unified star+DM density field.
// Stars get weight 1; DM is downweighted because DM is much more numerous
// and dynamically hotter, so without downweighting every halo center would
// look like a peak. 0.1 keeps stellar peaks dominant while still letting
// pure-DM subhalos rise above PEAKTHRESHOLD.
#define DM_DENSITY_WEIGHT 0.0f
//*******************
