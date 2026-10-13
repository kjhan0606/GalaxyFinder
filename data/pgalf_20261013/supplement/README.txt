Canonical final-analysis tables

photometric_curves.tsv: 15 linking nodes, all 76 entries, six band-view
combinations, primary and threshold/wing/brightest-pixel variants. Primary
rows contain positional, PSF aperture, and ±1/±2 arcsec boundary-band flux
sums and ratios. Non-primary auxiliary fields are undefined.
photometric_halo_metrics.tsv and photometric_view_selection.tsv: fitted
parameters and endpoint flags, with all-entry fixed-parameter scores.
new_summary.json: canonical halo medians and 16th/84th percentiles.
photometric_strata.tsv: lowz≤0.23, highz≥0.5, and equivalent radius split
at 3 arcsec. The radius is the median across views for each entry.
fine_validation_*.tsv: all-entry held-out evaluation with positional-neighbour
exclusion. best_common_b_iou means max_b(median_view score), distinct from
fitted_view_iou=median_view(max_b score). Shortfalls use the common-b reference.

stellar_population_metrics.tsv: current masses, luminosity-weighted formation
ages, signed mass ratios and age differences. Reference age uses full flux
for positional selection and q*f for the PSF aperture. Fractional age error
is group mean age/reference mean age−1. The summary takes absolute values
by view, then medians by halo and sample. Formation ages precede SSP age
clipping; the flux weights retain the same SSP assumptions as the images.
stellar_mass_linking_quantiles.json: actual per-entry mass and linking
quantiles in comoving, proper, and finest-cell units.
fine_mass_correlations.tsv: correlations with parent group mass. Nominal
probabilities do not account for unidentified repeated progenitors.
fine_matched_pair_comparison.tsv: ten near-output positional pairs. Epoch
and catalogue settings change together, so the differences are not an
isolated catalogue-setting experiment.

conserved_sampling_curves.tsv: nested fractions 1/2,1/4,1/8 with retained masses
and fluxes multiplied by 2, 4, and 8, respectively. Totals are preserved in expectation.
sampling_total_ratios.tsv records realized mass and flux ratios. Fixed masks
and rebuilt images are separate. Endpoint-censored fitted distributions do
not replace all-entry performance. The sampling experiment does not establish hydrodynamical convergence.
wide_boundary_summary.tsv: every set is restricted to the ±1/±2 arcsec band.
A halo statistic remains undefined if any view is undefined. n_defined records
coverage for each band, width, parameter, and ratio.
failure_view_catalogue.tsv and failure_halo_catalogue.tsv: every low-score
view and halo at b=0.1 or the held-out parameter, with separate failure flags
for both criteria. Group mass, area, offset, seed inclusion, and the number
of bright components are recorded. failure_flags_summary.tsv counts common-
parameter descriptive flags alongside held-out view failure counts. No merger cause is inferred.
conditional_uncertainty.json: 2000 conditional draws of 66 groups, keeping the
ten known positional pairs together. The intervals are not independent-volume
confidence limits. Original catalogue identifiers do not establish a full tree.

failure_drawing.npz/json: bounded corrected image cutouts and sampled group
positions. cluster_illustration.npz: corrected contour with reusable star
positions and mass drawing arrays. rest_frame.tsv supplies the unchanged
absolute-flux, unredshifted control. See scripts/README.txt for reproduction.

Older table and drawing products remain as provenance only. Final manuscript
statistics and figures use the canonical products above. Full particle
catalogues and full per-halo images/group caches are not included.
