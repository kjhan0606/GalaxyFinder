# TNG300 assembly and population-weighting analysis

This analysis follows the locked 36-branch TNG300-1 sample at snapshots 50, 67,
and 84. It uses the parent-FoF stellar cutouts included in
`../../tree_epoch_membership/`. No new simulation run or network fetch was
required.

## Analysis 1: particle-ID mass ledger

For each interval, the ledger matches `ParticleIDs` and decomposes the central
Subfind stellar-mass change into the current-mass change of stars assigned to
the central subhalo at both endpoints, pre-existing stars assigned to the
central subhalo only at the later endpoint, stars formed between endpoints
and assigned to the later central subhalo, and the earlier mass of stars no
longer assigned to the central subhalo. The ledger closes to a maximum absolute
fractional residual of `9.668e-16` of the starting central mass.

The earlier `RemoveSatelliteGalaxy` flag identifies stars in resolved
non-central subhaloes with at least 100 stellar particles. The flag is boolean,
so a source subhalo ID is unavailable. Pre-existing stars missing from the
earlier parent-FoF cutout have unknown prior locations. The analysis does not
label stars formed between outputs as in situ because the birth host is not
provided. The current cutouts contain no particle velocities or potential.

The six selection masks from the paper renderer were recreated at the physical
observer redshift and compared against the saved per-particle score tables.
Maximum absolute differences are `0.000e+00` dex in selected mass and
`0.000e+00` dex in projected half-light radius across `3888` branch, snapshot,
band, projection, and method records.

## Analysis 2: stellar-population size controls

The analysis measures projected half-light radii for central-subhalo stars
using the existing dust-free SSP flux grids. The full calculation retains each
particle's age and total metal mass fraction. Three counterfactual controls
replace the particle ages, metallicities, or both with the central population's
mass-weighted mean log value at each output. The controls preserve particle
positions and formed masses. The Shapley split of the log-radius difference
shares the age-metallicity interaction equally between the two terms.

The radii use full particle fluxes at particle positions, as in the released
central-reference measurement. The calculation excludes image PSF and
isophotal selection and therefore isolates population weighting. The physical
observer-redshift condition changes the observer-frame filter between outputs.
The fixed-redshift condition keeps the observer redshift at 0.197284.

Radial profiles use projected radius in units of the projected current-mass
half-mass radius. Inner and outer profiles combine bins at $R<R_{50,\rm mass}$
and $R>R_{50,\rm mass}$, respectively. The compact summary script pairs
endpoints and observer conditions by branch and projection before reducing
projections to branch medians. It uses 20,000 branch bootstrap resamples; the
seed is deterministic for each metric and keyed to base seed 20261106.

## Summary results

From snapshot 50 ($z=0.997$) to 84 ($z=0.197$), median central stellar mass
increases by `0.869` times its starting value. The largest median ledger term
is `0.789` times the starting mass in pre-existing stars newly assigned to the
central subhalo. A median `88.1` per cent of that mass had the earlier
resolved-satellite flag. The flag is boolean and supplies no donor subhalo ID;
this result is a catalogue-membership accounting, not a measured merger
history.

At physical observer redshifts, central-subhalo half-light-radius growth is
`0.020` dex in $V$ and `-0.003` dex in $g$, compared with `0.386` dex in both
bands when age and metallicity are held spatially uniform. The paired
full-population minus uniform-population effect is `-0.255` dex in $V$ and
`-0.324` dex in $g$. The outer population has lower mass-to-light ratio than
the inner population, especially at snapshot 50. These are results for the
adopted dust-free SSP model and one TNG300-1 volume.

## Outputs

- `assembly_mass_budget_per_branch.csv` and `assembly_mass_budget_summary.csv`
- `selection_metrics_per_view.csv`
- `selection_population_capture_per_view.csv` and
  `selection_population_capture_summary.csv`
- `central_size_metrics_per_view.csv`, `central_size_growth_per_branch.csv`,
  and `central_size_growth_summary.csv`
- `radial_population_profiles_per_view.csv` and
  `radial_population_profiles_summary.csv`
- `key_results.csv` and `central_size_growth_paired_per_view.csv`
- `radial_gradient_summary.csv` and `radial_gradients_per_view.csv`
- `assembly_mass_budget.png` and `size_growth_population_controls.png`
- `analysis_manifest.json`

The TNG data documentation defines `GFM_StellarFormationTime` as the stellar
formation scale factor, `GFM_InitialMass` as birth mass, and
`GFM_Metallicity` as total metal mass fraction, not [Fe/H]. See the
[TNG data specifications](https://www.illustris-project.org/data/docs/specifications/).

## Reproducing the products

From `supplement/tng300_growth_bias_20261106/`, write a fresh particle-level
output directory and derive the paired summaries with:

```sh
python assembly_population_20261106/analyze_assembly_population.py \
  --output assembly_population_20261106/results_rebuild
python assembly_population_20261106/summarize_assembly_population.py \
  --input assembly_population_20261106/results_rebuild \
  --output assembly_population_20261106/results_rebuild
```

The particle-level script reproduces the membership ledger, rendering
controls, and radial profiles. The summary script pairs endpoints and
observer conditions by projection within each branch before summarizing. It
uses 20,000 bootstrap resamples per statistic and deterministic metric-keyed
seeds based on 20261106. The supplied tables can be inspected without rerunning
the particle-level analysis.
