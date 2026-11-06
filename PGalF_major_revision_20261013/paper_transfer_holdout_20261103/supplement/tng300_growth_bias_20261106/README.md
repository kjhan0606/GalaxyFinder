# TNG300-1 stellar-growth analysis

This directory contains the 36-branch TNG300-1 analysis of stellar-mass,
half-light-radius, and half-mass-radius growth. The analysis compares
isophotal particle selections, stellar FoF groups, projected and spherical
apertures, and the stars bound to the central SubLink subhalo. Snapshots 50,
67, and 84 correspond to physical redshifts near 1.0, 0.5, and 0.2.

## Included inputs and results

The complete local review package includes the cutouts and branch files below.
The Overleaf review copy carries the analysis code and derived result tables,
but omits the binary HDF5 cutouts and compiled linking library. Reproducing
the image measurements from Overleaf therefore requires the local cutouts (or
a fresh TNG API download) and a compatible library build.

- `lineage_cohort.csv` and `lineage_epochs.csv` define the locked 36-branch
  sample and the three target epochs.
- `tree_epoch_membership/` contains the star-particle cutouts and exact
  central/satellite membership flags used by the renderer. The original
  SubLink branch files are `mpb_rank_*.hdf5`.
- `growth_scores_r1audit_20261106/` contains the versioned per-branch
  measurements with isophote detection, fallback, region-area, and image-edge
  fields. `growth_r1r4_audit_20261106/` contains the per-view and per-epoch
  detection tables, full-cohort and detection-qualified growth comparisons,
  mass-recovery ratios, and size-growth summaries.
- `growth_scores/fsps_cache/` contains the stellar-population flux grids used
  by the renderer. The grids are local inputs; the renderer does not need a
  user-specific home directory.
- `assembly_population_20261106/` contains the exact-ID central-mass ledger,
  population-weighting controls, radial age/metallicity/M-to-L profiles, and
  branch-paired summary script and tables. The analysis uses the same locked
  36 branches and the existing parent-FoF cutouts.
- `source_dependencies/` contains the Python measurement modules, the stellar
  linking library, and the C source and build files for that library.

## Reproducing the R1--R4 tables and figure

Run the commands from this directory. Use new output names because the scripts
preserve existing results and refuse to overwrite files.

```sh
python render_growth_r1audit.py --workers 2 \
  --output growth_scores_r1audit_rebuild
python analyze_growth_r1r4_20261106.py \
  --scores growth_scores_r1audit_rebuild \
  --output growth_r1r4_audit_rebuild
python plot_growth_r1r4_20261106.py \
  --summary growth_r1r4_audit_rebuild/r1r3_growth_bias_summary.csv \
  --output ../../figures/tng300_tree_growth_bias_r1r4d_rebuild.pdf
```

The renderer writes 378 rows for each of the 36 branches. The analysis script
uses 20,000 branch bootstrap resamples with seed 20261106. The plot has two
vertically joined panels with separate axes and compact logarithmic-ratio
axis titles. The orange diamond in the isophotal category shows the $g$-band
result for the 29 branches with central-pixel detections in every view at
both endpoint snapshots; circles show the full 36-branch cohort. The
manuscript uses `../../figures/tng300_tree_growth_bias_r1r4d.pdf`.

The scripts were run with Python 3.13.11, NumPy 2.3.3, SciPy 1.18.0, h5py
3.15.1, Astropy 7.2.0, Matplotlib 3.10.8, and Numba 0.63.1. The cached FSPS
grids remove the need for FSPS during the rerun above. If a cache is absent,
install python-fsps 0.4.7 with its FSPS stellar libraries and set `SPS_HOME` to
the local FSPS data directory before running the renderer.

To reproduce the supplied assembly and stellar-population summaries from
the saved per-view tables, run from this directory:

```sh
python assembly_population_20261106/analyze_assembly_population.py
python assembly_population_20261106/summarize_assembly_population.py
```

The first command reads the included particle cutouts and writes the
particle-ID ledger, selection capture fractions, controlled radii, and
radial profiles. The second command pairs measurements by branch and
projection, then writes the compact tables. Each 95 per cent interval
resamples the 36 descendant branches 20,000 times. Summary-statistic seeds
are deterministic and keyed by metric, with base seed 20261106. The results
README documents the interpretation limits and the new output files.

The bundled `libopfof_sb.so` was built with the Intel MPI toolchain used for
the published calculations. On a system with a compatible MPI installation,
rebuild the library from the included source with:

```sh
make -C source_dependencies/opfof_source/variant
cp source_dependencies/opfof_source/variant/libopfof_sb.so \
  source_dependencies/libopfof_sb.so
```

The scripts `fetch_sublink_branches.py` and `fetch_tree_epoch_membership.py`
document the TNG API retrieval. They are not needed to reproduce the supplied
measurements because the branch files and particle cutouts are included. API
retrieval requires a user-held TNG API key. Re-fetching snapshot 84 cutouts
also requires the original parent FoF cutouts, supplied to the fetch script
with `--snapshot84-parent-dir`.

The exact-ID mass ledger and population-weighting controls use the locked
branch sample and included particle cutouts. To recreate these products in a
fresh directory from this directory, run:

```sh
python assembly_population_20261106/analyze_assembly_population.py \
  --output assembly_population_20261106/results_rebuild
python assembly_population_20261106/summarize_assembly_population.py \
  --input assembly_population_20261106/results_rebuild \
  --output assembly_population_20261106/results_rebuild
```

The first command recomputes the particle-level ledger, population controls,
and radial profiles; the second creates paired branch summaries. The compact
summary and its scope limits are described in
`assembly_population_20261106/results/README.md`.

To verify the complete supplement after extracting the paper package, run
`sha256sum -c SHA256SUMS` from the `supplement/` directory.
