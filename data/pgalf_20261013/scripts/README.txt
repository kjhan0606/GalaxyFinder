Reproducing the PGalF revision

Run from the package root with the versions in supplement/calculation_config.json.

python scripts/primary_diagnostics.py
python scripts/plot_scores.py
python scripts/plot_failures.py
python scripts/plot_cluster.py

primary_diagnostics.py calculates the consistent-distance calibration, all-76
omitted-output validation, output-specific gains/failures and fixed-b luminosity
errors. Training is not filtered by a fitted interior maximum. It also reports
endpoint counts and all-original-entry scores for archived fixed-mask thinning.
plot_scores.py draws primary consistent-distance scores, selected parameters,
completeness/precision and total-versus-boundary curves. The unredshifted control
figure overlays corresponding primary image curves.
Cluster and failure illustrations retain archived image normalization and are
labelled as construction/selection controls in the manuscript.

The archived control summaries can be regenerated separately with
python scripts/calibration_diagnostics.py
python scripts/luminosity_diagnostics.py
python scripts/revision_diagnostics.py
Those commands preserve historical definitions and do not replace primary files.

With python-fsps 0.4.7 and the recorded libraries under SPS_HOME,
python scripts/audit_photometry.py
reproduces the seven archived grids and distance factors. The primary image uses
the distance factor before convolution and isophote selection, never after scoring.

measurement_methods.py contains the numerical image and score definitions.
source_archive includes production readers, linking/image routines, C sources,
headers/build recipe, original simulation namelist and PSF reference source.
Absolute paths in the producing environment require configuration elsewhere.
The binary hash is recorded, but the compiled library and particle catalogues
are not included. No claim of full membership reproduction from tables is made.
Source archives preserve original whitespace and historical labels.

audit_targets.py derives target_membership.tsv from existing membership caches
and the catalogue ordering. audit_particles.py in the source archive audits
SSP-grid clipping from original particles on Lageunha. The scripts require the
original project paths and caches and are not runnable using only released tables.
