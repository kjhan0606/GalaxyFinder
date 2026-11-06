# CV0--19 score-unseen CAMELS holdout

## Frozen analysis plan (2026-11-03)

This analysis extends the PGalF fixed-parameter transfer test to CAMELS Swift-EAGLE
boxes CV0--CV19. These boxes do not overlap the earlier CV20--CV26 analysis. The
holdout scores have not been inspected as of the plan lock. The centre rule was
selected during earlier work on CV20--CV26 and TNG300-1, so this is a score-unseen
holdout under a locked analysis, not a fully preregistered blind experiment.

### Sample

- Snapshot 090 at simulation redshift zero; it is rendered at the same mock
  observed redshift, 0.19043796228481114, used in the earlier transfer analysis.
- Select FoF hosts with catalogue `GroupMass` in `[1000, 10000)` in units of
  `1e10 Msun/h` (that is, `[1e13, 1e14) Msun/h`).
- The central galaxy is the subhalo with the largest stellar-particle count and
  must contain at least 100 stellar particles.
- The sample is frozen by the generated manifest. Expected size: 130 hosts in
  20 CV boxes. No host is removed based on image or metric outcomes.

### Locked mock image and target contour

Use the existing z=0.1904379623 FSPS grid, V and SDSS g limits, image scale,
PSF, stellar selection, catalogue-ID membership, and linking-parameter grid
from the earlier CAMELS transfer calculation. Use the catalogue central
subhalo's potential minimum (`SubhaloPos`) as the image centre and select the
connected isophotal region containing that centre. This centre choice was
already made from earlier analyses; it will not be retuned using CV0--CV19.
Each host is measured along the same three orthogonal projections. The target
definition does not identify galaxy membership dynamically; it defines the
image region used in this structural comparison.

For Swift-EAGLE, the particle membership input is the catalogue `/IDs/ID` block.
This is the prescribed CAMELS catalogue route for suites whose snapshot
particles are not laid out in SUBFIND order. All 130 catalogue star-ID lists
resolve to the expected snapshot stars. In five hosts, the snapshot
`FOFGroupIDs` labels have zero intersection with the selected catalogue group;
those labels are retained as a diagnostic and do not replace catalogue-ID
membership. The CV8 and CV11 group catalogues were SHA-256 checked against the
current public CAMELS files and match exactly. The official CAMELS documentation
also records a historical Swift-EAGLE group-catalogue bug and says the public
catalogues have since been corrected: [Subfind catalogue notes](https://camels.readthedocs.io/en/latest/subfind.html).

### Comparisons

For every host and band, compare:

1. `b_star=0.20`, the original linking-parameter baseline;
2. a common fixed `b_star=0.10` in both bands;
3. the low-redshift Darwin calibration used for the transfer analysis
   (`b_star=0.12` in V and `0.10` in g);
4. a 30 proper-kpc projected aperture.

No best-fit value will be selected from CV0--CV19. The fine grid is evaluated
only to generate the stated fixed comparisons and diagnose whether the fixed
values lie on plausible portions of the measured curves; any per-host optimum
is not used as transfer evidence.

The paper's all-output median calibration is `b_star=0.08` in V and `0.10` in
g. The `0.12/0.10` pair used here is the already reported low-redshift Darwin
calibration, which matches this mock image redshift; the common `0.10` run is
retained as a separate control.

### Metrics and failure rules

First reduce the three projections to one median per host and band. Report
light-weighted IoU, completeness, precision, contamination (`1 - precision`),
signed `log10(L_selected/L_contour)` and `log10(R50_selected/R50_contour)` biases,
and the 90th percentile of the absolute host-level luminosity and R50 errors.
Also report the fraction of hosts with median IoU below 0.90.

A projected target contour is invalid if the centre selection falls back to a
different isophotal component, contains no pixels, has zero contour luminosity,
or touches the image edge. Invalid views are assigned IoU and completeness of
zero; precision and photometric/size errors are undefined. A host is counted
once as a large failure for a method if its host-median luminosity or R50 differs
by more than a factor of two, median contamination exceeds 0.20, or any of its
three projected target contours is invalid. The number of hosts with any
orientation-specific large failure is reported separately.
The 90th-percentile absolute luminosity and R50 errors are conditional on
valid target contours; invalid contours remain in the failure fractions.

Hosts are the statistical units. CV boxes are resampled first and hosts within
each sampled box second (2,000 bootstrap draws); the three projections are
never treated as independent hosts. Report box-level results as well as the
pooled host summaries. Paired method contrasts use the same hosts and the same
nested box-host resampling.

The full per-host membership cross-check is in `membership_audit.csv`; its
summary is `membership_audit_summary.json`. The per-host FoF and image metrics
are stored in `results/` and the final host-level summaries are stored under
`summary_cv019_20261103_v2_potential_centered/`.
The verified numerical findings and their interpretation are in `RESULTS.md`.

### Provenance and limitations

The calculation evaluates image-structure transfer in one Swift-EAGLE model at
one simulation output. It does not test redshift evolution, automatic centre
identification from noisy data, dust, sky noise, or dynamic membership, and does
not measure a physical intracluster-light fraction. The centre rule is inherited
from earlier work and was not blinded; the new CV boxes keep the metric outcomes
unseen until the locked calculation is complete.
