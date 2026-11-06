# PGalF manuscript review release

Release ID: PGalF-MNRAS-editorial-20261107-r5

This is the single entry point for the current local manuscript and its
supplementary package. The manuscript and derived tables are synchronized to Overleaf; source
code is maintained in the GalaxyFinder GitHub repository. The complete
local data package has not been deposited publicly.

## Verify and rebuild the manuscript

From this directory, run:

    ./reproduce_current_release.sh

The script verifies RELEASE_SHA256SUMS and the supplementary
SHA256SUMS, then builds main.pdf and
supplement/extended_notes.pdf with latexmk. The checksums identify the
frozen files in this local review release. Rebuilding with a different TeX
installation may change PDF metadata and therefore the PDF hashes.

## Current analyses and data

The main manuscript uses the corrected Darwin-cosmology luminosity distances.
Current analysis inputs, results, and scripts are grouped by analysis:

- The six-output Darwin calibration and corrected overlap curves are in
  supplement/; the overlap and selection tables list the retained sample.
- The 130-host CV0--CV19 validation sample, analysis protocol, measurements,
  and scripts are
  in supplement/cv0_19_holdout_20261103/.
- The TNG300-1 central-subhalo growth analysis, 36-branch sample list, rendered
  inputs, overlap measurements, and scripts are in
  supplement/tng300_growth_bias_20261106/. The final growth figure is
  figures/tng300_tree_growth_bias_r1r4d.pdf; its orange diamond marks the
  detection-qualified $g$-band result. The exact-particle-ID stellar-mass budget and
  stellar-population size controls are in
  supplement/tng300_growth_bias_20261106/assembly_population_20261106/.
  Their branch-paired summaries and radial gradients are in that folder's
  results/ directory. The top-level analysis_manifest.json points to the
  detailed manifest in the results/ directory.
- The TNG300-1 40-host image-transfer validation is in
  supplement/tng300_subhalopos_r50_transfer_20261103/.
- Within-output validation and the equal-output weighting sensitivity are in
  supplement/within_output_validation_z019_20261106/ and
  supplement/output_weighting_sensitivity_20261106/.
- The corrected-distance PSF-aperture overlap, five boundary-band widths,
  and central-pixel detection flags are in
  supplement/astra_review_response_20261106/results/. The calculations and
  figure sources are in the matching review-response directory.
- The branch-level TNG300-1 growth-error and mass-recovery figure is
  supplement/astra_review_response_20261106/branch_level_growth_recovery.pdf;
  its source is plot_branch_growth_recovery.py in the same directory.

Each analysis README gives the exact command to reproduce its supplied
tables. The TNG growth README gives the renderer, analysis, and figure
commands. The supplied tables support inspection and regeneration of the
reported summary statistics without rerunning the image pipeline.

## Archived products

main.tex and main.pdf define the current manuscript. The corrected
distance results and figures are primary. The archived-normalization cluster
image, failure maps, boundary-distance table, and no-dimming control now
appear only in the supplementary notes and are identified as historical
illustrations. Legacy plotting scripts such as plot_scores.py,
plot_failures.py, and plot_cluster.py reproduce archived-normalization
figures; they are not the entry point for current primary figures.

## Data access limits

The original Darwin particle catalogues are not included in this local
package. Reconstructing the exact stellar membership and rerunning the
primary Darwin image measurements requires those catalogues and the
archived stellar-population grids. Requests for access to the original
Darwin inputs should be directed to the corresponding author. The local TNG
analysis package includes the particle cutouts and SubLink branch data used
for the reported sample. The Overleaf review copy carries the analysis code
and derived tables but omits the binary TNG cutouts and compiled linking
library; those files are present in the local release. Rerunning the image
measurements requires the local cutouts (or a fresh TNG API download with a
user-held API key) and a compatible library build.

This author-review release has not been deposited publicly; no public URL or
DOI is available yet. The paper's Data availability statement records this
status explicitly. A persistent public link will be needed for journal
submission.

Regenerate the current diagnostic figures with `python scripts/render_editorial_figures.py`.
The editorial renderer uses saved measurements, 1.2-pt axes, inward ticks,
and mathematical axis labels. The former high-redshift bar panel is now
a table of detection counts and mean overlaps in the main manuscript.
