Referee controls, 2026-10-13

Scope: 12 distinct central galaxies at z=0.190, cached b_star=0.10 and
0.20 memberships, and 15-node image curves. No new FoF runs or simulations.
Raw catalogue calculations ran on Lageunha with two single-thread workers.
All twelve entries completed. The main manuscript is a calibration paper.

Portable commands from this directory (NumPy, pandas, SciPy, matplotlib):
  python scripts/summarize_controls.py
  python scripts/plot_controls.py
The summary uses the released baseline and recentered curves. The plot uses
outer_cut.tsv and population_controls.tsv. Source reconstruction scripts
retain the original project paths and require unreleased raw catalogues,
full membership caches, and the earlier reader/library setup.

results/outer_cut.tsv: actual and outer-truncation half-mass radii.
results/population_controls.tsv: paired age and metallicity, radial support,
coverage, fixed-support bin merging, young-star cuts, joint radius/age cells,
and ellipsoidal bins. coverage always divides by the original excluded mass;
eligible_coverage divides by the mass surviving a control's cuts. Fixed-support
rows therefore have eligible_coverage=1. Reported radial limits were corrected
to use the actual fixed support after collection; weighted means are unchanged.
results/components.tsv: retained, excluded, matched-excluded and unmatched-
excluded properties, including fractions assigned to central/small/background
catalogue categories. These fractions do not identify physical birth sites.
results/absolute_properties.tsv: both b values, current masses in Msun,
proper radii in kpc, observed flux densities in microJy, AB magnitudes, ages
in Gyr, metallicity in solar units (Zsun=0.0134), and velocities in km/s.
Group projected half-light radii use density-peak centres and stellar flux
weights before PSF convolution. Isophotal equivalent radii use primary images.
results/normalized_centre_offsets.tsv: 3D offsets / 3D half-mass radii and
view-median projected offsets / primary isophotal radii.
results/recentred_targets.tsv: shrinking-centred isophotal sizes and half-light
radii measured in convolved images. Do not treat different centre/PSF choices
as a pure membership comparison.
results/control_summary.json: halo medians, scatter, sign counts, and conditional
96.1% order-statistic median intervals (third/tenth ordered values for n=12).
centres/: saved centre coordinates and recorded iteration traces reused here.

Interpretation: the measured half-mass radius decrease is reproduced by outer
truncation and is not independent evidence of structural contraction. Radial
population differences survive young-star cuts but vary with age matching
and geometry. No stellar-origin or ICL identification is established.
Full 76-entry recentering, realistic sky/dust/2D masks, independent simulations,
and stellar histories remain outside the completed controls.
