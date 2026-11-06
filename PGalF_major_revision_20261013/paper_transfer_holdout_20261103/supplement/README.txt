For the current manuscript release PGalF-MNRAS-review-20261106-r3, use the
single entry point in ../REPRODUCIBILITY.md. This file indexes supplementary
data and records calculation provenance.

The 2026-11-06 Astra response calculations are in
astra_review_response_20261106/. They report corrected-distance PSF-aperture
overlap, a 15-node linking-parameter grid, boundary-band scores at five
angular half-widths, and corrected central-pixel detection flags. The paired
figure sources also set axis ticks inward and show ticks on both sides of
joined panels. See the response document in that directory for the comment-
by-comment record. The release is local to the co-authors and has no public
repository URL or DOI yet.

Supplementary data for the stellar linking-length calibration

Scope and identification
The tables describe 76 output entries, with three views and two bands per entry.
The 75 entries with interior maxima in every view form the baseline matched
sample. Join tables by (step, nstar), with los and band for individual curves.
The pair (step, catalogue_row) identifies the original catalogue entry. Row
numbers start at zero and are not merger-tree identifiers. Output 00967 uses
the stellar-link-off top100/0 catalogue. Other outputs reassemble the parent
stars from the corresponding stellar-link-on catalogue. No particle catalogue
or membership was modified for this supplement.

Files
halo_metadata.tsv supplies the catalogue row and byte offset, group and stellar
masses in h^-1 solar masses, parent and satellite counts, and published interior
flags. halo_positions.tsv supplies the halo-header centres in h^-1 comoving Mpc
and original catalogue velocities without a new velocity-unit interpretation.
The centres allow the positional matches to be checked without loading particles.
lineage_dm_id_overlaps.tsv records exact dark-matter particle-ID intersections.
The accepted links require the later halo to be the best-overlap descendant
and to contain more than half of the earlier halo's dark-matter particles.
lineage_position_candidates.tsv records reciprocal-nearest centre matches;
these are provisional and are used only for the full-sample sensitivity check.

sightline_baseline.tsv, sightline_variants.tsv, and sightline_refined.tsv are
unchanged copies of the deposited tables at data commit 5e3cc51. The baseline
classification is authoritative. In particular, 00967:17494810 remains censored
along x and y in both bands. Recomputing argmax from rounded scores can change
ties and must not silently replace the published classification. The PSF counts
13 and 22 compare with published positional lengths. Comparing rounded refined
POS scores with rounded PSF scores can instead give 12 and 21. The different
counts concern a tie convention, not newly measured trajectories.

rest_frame.tsv preserves the original tagged records for the four rest-frame
controls. These are archived-normalization qualitative checks, not the corrected-distance primary calibration. The first column is step; the record schemas appear in the header.
Boundary areas remain angular areas, although the rest-frame brightness uses
flux per proper square parsec and unredshifted absolute magnitudes.

boundary_distance.tsv records fluxes near the adopted boundary at the published
best length of each view. L_miss is group light outside the region; L_extra is
region light outside the group. le1, le2, and le5 denote distances no greater
than 1, 2, and 5 arcsec on the appropriate side. Distances are Euclidean distances
between pixel centres, evaluated at stellar positions. The ring instead uses
8-connected morphological dilation and erosion. L_miss_off is exterior group
light outside the stamp. empty and border_pixels record region diagnostics.
The original published and recomputed overlaps are both retained. Halo summaries
use the median over three views, then the median over fully interior halos.

boundary_ring.tsv preserves the earlier exploratory comparison. The corrected-distance
comparison restricted to the boundary band is in boundary_ring_corrected.tsv,
for each corrected best length, 0.10, and 0.20. The numerator is L(G & I & R)
and the denominator L((G | I) & R), with stellar group G, isophotal selection I,
and boundary band R. Completeness and precision use the same restricted sums.
The region area, published best length, full-image score, and endpoint flag are
preserved from the original table. Undefined ratios are nan rather than zero.
The earlier table compared G with R and included central group light outside
the band. The corrected-distance selected-parameter median overlaps in
boundary_ring_corrected.tsv are 0.4614 in V and 0.4760 in g; the associated
completeness and precision summaries are described in extended_notes.tex.

image_diagnostics.tsv preserves three tagged record types. Step precedes each
record. Fields after the tag are:
MASS: z nstar n_parent particle_mass_median particle_mass_p16 particle_mass_p84
      ell_b0.10_median ell_b0.10_p16 ell_b0.10_p84 group_mass stellar_mass
      Masses are h^-1 solar masses; ell is h^-1 comoving Mpc.
CMP: z nstar los band component_count adopted_area seed_component_area
     same_component fallback centre_offstamp seed_offstamp centre_seed_arcsec
PHO: z nstar los band positional_inside positional_total convolved_inside
     convolved_total satellite_positional_inside satellite_positional_total
     satellite_convolved_inside satellite_convolved_total
PHO baseline totals already exclude satellites. Do not subtract the satellite
columns from the baseline fluxes. No parent-light fraction is inferred here.

mass_correlations.tsv lists Spearman coefficients and nominal two-sided
asymptotic p-values recomputed with scipy.stats.spearmanr. Values use published
halo medians and group masses. Ties, small samples, repeated systems, and changing
mass cuts restrict interpretation; these are not independent tests of evolution.

sensitivity_summary.tsv reports paired halo-median changes relative to the
baseline. N counts halos interior on all views in both comparisons. Percentiles
use linear interpolation. n_good_baseline_sightlines counts individual interior
views with baseline overlap >= 0.90 across all 76 entries, including qualifying
views of the excluded halo. n_edge and n_below_0.90_interior are disjoint counts.
The fixed-image tests use the retained subset for both linking and scoring.
Retained particle masses and luminosities are unchanged. The seed is selected
again from the retained particles. The original image boundary stays fixed.
The rebuilt-image tests change both the image and the linked subset.
The random subsets are nested, with seed sequence [20261007, nstar].

Photometry
The seven fsps_*.npz files are byte-for-byte copies of the cached input grids.
Each contains age (Gyr), logz = log10(Z/0.0142), flux, and z. Flux has dimensions
(age, metallicity, band), with bands V and SDSS g, and equals 10^(-0.4*m) for one
solar mass formed. At z > 0, m is the FSPS apparent AB magnitude. At z=0, m is
the absolute AB magnitude. Stars use bilinear interpolation of log flux in
log age and log metallicity, clipping outside the grid. Multiply by formed mass
in solar masses, converting the catalogue h^-1 mass first. The image uses flux
per angular pixel area; no extra luminosity-distance or Tolman factor is applied
on top of the redshifted apparent magnitudes. For the rest-frame maps, divide
absolute relative flux by proper area in pc^2, then add
2.5*log10((206265/10)^2) to -2.5*log10(flux_per_pc2).

filter_v.tsv and filter_sdss_g.tsv archive the installed FSPS transmissions.
The producing script sets Chabrier IMF, zcontinuous=1, sfh=0, dust2=0,
and add_neb_emission=False. AB is the default magnitude system.
The available python-fsps installation is version 0.4.7. The flux caches did
not record the generating version or compiled stellar-library configuration.
On 2026-10-13, all 32*12*2 values in every one of the seven archived grids were
reproduced exactly with the available installation and the recorded parameters.
photometry_audit.json records zero differences, input hashes, and distance factors.
fsps_build.json records the compiled library names (mist, miles, and DL07),
the extension hash, and the Python version on Lageunha. The attenuation parameters
dust1 and dust2 are zero.
The default AGB circumstellar dust model remains enabled. The build record
also contains all initial population parameters. Grid ages, metallicities,
and redshifts are supplied separately by audit_photometry.py.
The packaged Fortran luminosity-distance source uses Omega_m=0.27,
Omega_Lambda=0.73, and H0=72 km/s/Mpc. Angular sizes in the analysis use the
simulation cosmology Omega_m=0.315 and H0=67.4 km/s/Mpc. No correction between
the two distance conventions is made in the producing script. The generic
Python docstring quotes different cosmological values; the source constants
are reported here. Historical binary settings cannot be recovered from cache
metadata. The archived grids preserve the actual numerical normalization.
The isophotes and scores were recalculated for all 76 entries on Lageunha.
The flux multipliers range from 0.8873 to 0.9589, corresponding to 0.1298 to
0.0456 magnitudes of dimming. The same 75 entries remain in the calibration.
The median parameter becomes 0.08 in V and remains 0.10 in g. No finer-grid,
PSF-fraction, or sampling experiment was repeated with the corrected distances.
Those experiments continue to refer to the explicitly documented archived flux
normalization. The main paper and abstract report the distance comparison.

PSF score
q_i is the fraction of the normalized, truncated kernel falling within the
baseline region, sampled at the stellar pixel. q_i=0 for off-stamp stars.
L_group=sum_G f_i, L_inside=sum_all q_i*f_i, L_both=sum_G q_i*f_i.
IoU=L_both/(L_group+L_inside-L_both), completeness=L_both/L_inside,
precision=L_both/L_group. The group denominator retains whole-star flux.
The archived TSV columns recall and purity have the same numerical definitions
as completeness and precision in the revised manuscript. Column names remain
fixed so existing readers of the deposited tables continue to work.
The weight does not enter the linking operation.

Provenance and limits
provenance.json distinguishes recovered installation metadata from historical
cache metadata. SHA256SUMS hashes the supplement, manuscript, figures, and scripts.
Run sha256sum -c SHA256SUMS from the supplement directory. Paths beginning with
../ refer to the manuscript root. The merger-versus-catalogue
origin of the large image-centre offset remains undetermined. Supplementary
files supplement the public deposit; preparing these files does not update the
GitHub repository. The 2026-10-12 supplement reused existing calculations.

Editorial revision on 2026-10-12
-------------------------------
extended_notes.tex and extended_notes.pdf contain the watershed and tidal
recipe, mass and positional-match diagnostics, complete bootstrap and PSF
comparisons, and image-construction and particle-sampling tables moved from
the main text. Compile main.tex first. From the project root, compile
supplement/extended_notes.tex with xr able to read main.aux. References
imported from the main paper carry the prefix main-. The supplied PDF can
be read without compilation.

The scripts directory contains drawing code using the archived baseline
and rest-frame tables. Published b_best nodes and endpoint flags are
retained. A common redshift palette is used in every score plot.
The failure_panels directory stores six cropped diagnostic raster panels
from the completed calculations. The crops retain the original colours.
Their extents and core coordinates are recorded in scripts/plot_failures.py.
The original save_figure function in sb_revise.py used gray_r from 20 to
28 mag arcsec^-2. The shared brightness scale follows that normalization.
New marker overlays and display windows change presentation only.

No particle catalogues, links, flux grids, or sensitivity curves were
recomputed for the editorial revision. The five main-method equations and
the two retained numerical table bodies match the preceding version.

Scientific revision on 2026-10-13
-------------------------------
All 76 particle calculations completed on Lageunha. Every original score triple
and group count reproduced before completion was recorded. The seven photometric
grids reproduced exactly. The corrected boundary values replace the earlier ring
table. Original full-image calibration tables remain unchanged and retain their
archived flux normalization. The new distance comparison is reported separately.

recomputed_curves.tsv records full-precision scores, group/inside/intersection
light sums, and their boundary restrictions at every original grid point for
both baseline and consistent_distance images. distance_sensitivity.tsv pairs
the published original choice with the maximum of each corrected curve.
selected_parameters.tsv preserves raw recomputed maxima as an additional
diagnostic. Original published choices and flags remain authoritative for the
baseline comparisons, including any ties. revision_summary.json gives sample
counts, boundary statistics, distance changes, and catalogue checks, regenerated
by scripts/revision_diagnostics.py.

catalogue_audit.tsv records each selected galaxy, brightest galaxies in V and g,
most massive galaxy, star counts, and reconstructed current stellar masses.
catalogue_galaxies.tsv contains all 58995 catalogue galaxy records. galaxy_index
identifies a galaxy within the halo block, not across outputs. The x, y, and z
columns are catalogue galaxy centres in h^-1 comoving Mpc. current_mass is in
h^-1 solar masses, and flux_V and flux_g are relative apparent AB fluxes. The
selected galaxy is most massive in all 76 halo entries and brightest in 52
entries in each band. The sets of 52 entries need not coincide.

validation_folds.tsv records each omitted output, the number of calibration
entries, additional positional exclusions, and the selected common b.
validation_predictions.tsv retains every halo entry, including the originally
excluded halo, and the baseline interior flag. validation_summary.json gives
16th, 50th, and 84th percentiles for the 75 retained entries and a separate
summary for all 76 entries. The calibration
maximizes the median overlap after taking medians over views. Calibration
excludes the held-out output and remaining halos within 0.5 h^-1 cMpc of any
held-out halo, with periodic separations. More distant progenitors cannot be
identified without merger trees. The statistical summaries count output entries.

radial_luminosity.tsv contains recomputed full-precision sums for b=0.10 and 0.20.
The four annuli have edges at 0, 0.5, 1.0, 1.5, and infinity in units of the
equivalent isophotal radius. L_group, L_inside, and L_both are full stellar
relative fluxes selected by projected position, without PSF aperture weights.
luminosity_summary.json gives signed and absolute luminosity differences and
exterior-light fractions. Absolute values are taken in each view before halo
and sample medians. The original seven-significant-digit RAD sums reproduce
within their printed precision, with maximum fractional difference 4.930e-7.
radial_reproduction.json records the comparison. Results use the original
75-halo selection and archived flux normalization.

calculation_config.json records explicit image parameters, source hashes,
central-galaxy selection, and available package versions. See scripts/README.txt
for commands and the limits of reproduction without original particle data.
scripts/source_archive contains the production particle readers, image and
measurement routines, and stellar FoF C sources with the required headers and
build recipe. Original absolute input and output paths are preserved and require
configuration outside the producing environment. No compiled binary or particle
catalogue is bundled. Source hashes and the production binary hash are recorded
separately in calculation_config.json.
cluster_illustration.npz stores only the arrays needed for the cluster drawing,
including sampled stellar positions, the isophotal region, mass image, and
galaxy ellipses. Removing unused full-resolution stellar arrays reduces the
archive to about 1 MB. plot_cluster.py reproduces the original rendered pixels.
The public deposit has not been updated for the 2026-10-13 revision.


Revision 2026-10-30 corrected-distance and external-validation files
The corrected-distance version of recomputed_curves.tsv and selected_parameters.tsv is the primary Darwin calibration. sample_manifest.tsv records the 76 Darwin entries and the distinct calibration and artificial-redshift inclusion flags. validation_corrected_folds.tsv, validation_corrected_predictions.tsv, and validation_corrected_summary.json report the leave-one-output transfer check.

The centered CAMELS Swift-EAGLE transfer products are camels_host_metrics.csv, camels_box_summary.csv, and camels_paired_contrasts.csv. Artificial-redshift detection fractions and overlap means are in artificial_redshift_summary.csv. These tables distinguish entry-level descriptive values from independent-system inference. The raw artificial-redshift measurements, corrected boundary-ring table, and CAMELS measurement-source archive are supplied with this revision.


CAMELS catalogue identification and resolution
The independent transfer sample uses CAMELS Swift-EAGLE first-generation L25n256 CV boxes CV20--CV26, snapshot 090 at z=0, and the paired public files snapshot_090.hdf5 and groups_090.hdf5. The group catalogue header records Omega_m=0.3, h=0.6711, and a 25 h^-1 cMpc box. The snapshot Header InitialMassTable gives an initial gas-element mass of 1.88725e7 Msun; current star particle masses are variable. The catalogue Parameters group records SofteningMaxPhysType4=5.84 in h^-1 pkpc. The analysis identifies the central subhalo from the catalogue and uses its SubhaloPos (minimum-potential position) as the common image/aperture/radius origin.

The combined raw measurements are in `artificial_redshift_metrics_v6.tsv.gz` (27,000 view/parameter records). To regenerate the summary, run `python scripts/summarize_artificial_redshift.py supplement/artificial_redshift_metrics_v6.tsv.gz supplement/artificial_redshift_summary.csv`. The script also accepts a results directory containing the per-entry metrics.tsv products from the proper-size, Darwin-distance-normalized run. The script gives every view equal weight and reports the detection fraction, mean overlap among detections, and mean overlap across all views; views without a detected isophote have zero overlap in the all-view mean. It also checks that the all-view mean equals the detection fraction times the mean among detected views. The main figures are regenerated with `python plot_revised_validation.py`.


Direct fixed-b versus aperture comparison
camels_b010_vs_aperture.csv compares b=0.10 with the projected 30 pkpc aperture on identical hosts, with positive differences denoting b=0.10 minus aperture. Both host and box-then-host 95% intervals use 5,000 bootstrap draws. Run python scripts/compare_camels_b010_aperture.py from the manuscript root to regenerate the table from camels_host_metrics.csv. The g-band overlap interval includes zero, although b=0.10 reduces luminosity and size errors and increases completeness in both bands. The aperture includes less extra-isophotal light.

Corrected-distance fixed-parameter image checks
camels_image_checks/ contains nine one-factor image conditions for the same 49 CAMELS hosts. Brightness limits, pixel scale, PSF wing, a flat positive sky residual, and three Gaussian pixel-noise amplitudes change the isophotal selection while b=0.10, b=0.20, and a projected 30 pkpc aperture remain fixed. The nominal calibration images remain noiseless. Viewing directions and noise realisations are summarized within host before inference. The supplied compressed projection table has 13,230 records after removal of the repeated baseline from the noise extension. Run python scripts/summarize_camels_image_checks.py --output /tmp/camels-image-reproduction from the paper directory to reproduce the five result tables without original snapshots. Full conditions, inference, exact generating sources, and validation are in camels_image_checks/README.md and export_validation.json. The image changes perturb boundary selection, while luminosity and radius retain noiseless particle fluxes. The checks do not measure noisy aperture-photometry errors or test automatic centre finding.

Lineage-grouped cross-validation
Run python scripts/validation_lineage_groups.py from the paper directory to
reproduce lineage_cv_components.tsv, lineage_cv_folds.tsv,
lineage_cv_predictions.tsv, and lineage_cv_summary.json. The exact-ID result
uses 53 entries at outputs 00178, 00475, and 00786; the full-sample result
adds positional candidates for uncovered output intervals and is explicitly
provisional. Both hold out each whole connected group, remove any training
group with a member within 0.5 h^-1 cMpc of the held-out group, and give each
training group equal weight when selecting b. Bootstrap intervals resample
groups, not output entries. The positional result is conditional on the
candidate matches and is not a merger-tree validation.

Exact-ID paired artificial-redshift image controls
population_redshift_exact_lineage_components_major_20261030.tsv contains
component-level metrics for the 53 entries at outputs 00178, 00475, and 00786.
The 20 major-descendant exact dark-matter particle-ID links define 33
components. The table reduces the three projections at each output and then
reduces repeated outputs within each component. Run
`python scripts/summarize_exact_lineage_redshift_modes.py` from the paper
directory to regenerate
population_redshift_exact_lineage_paired_modes_major_20261030.tsv and its
metadata file. The paired table compares each control with the full image at
fixed b=0.10 and bootstraps the 33 components 5,000 times. A positive
difference means a larger value in the comparison image. For the
non-detection fraction, a negative difference means fewer views without a
detected isophote in the comparison image. The 33 components are an exact-ID
grouping, not a complete census of independent halo lineages.


Transfer-validation revision on 2026-11-03
----------------------------------------
The manuscript now reports a score-unseen holdout of 130 CAMELS Swift-EAGLE
hosts in boxes CV0--CV19. The frozen protocol, sample manifest, catalogue-ID
membership audit, per-host metrics, box summaries, paired contrasts, and
calculation scripts appear in cv0_19_holdout_20261103/. The holdout tests
transfer across boxes within the Swift-EAGLE model. Scores remained unopened
until the sample and analysis protocol were fixed, although the centre rule
was selected in earlier work.

The TNG300-1 transfer results now use the central subhalo SubhaloPos as the
image centre and projected half-light radii as the size measure. Host metrics,
per-projection records, position manifest, API provenance, run status, and
source scripts appear in tng300_subhalopos_r50_transfer_20261103/. The earlier
root-level tng300_*.csv files preserve the SubhaloCM-centred measurements as a
centre-sensitivity comparison. The new TNG analysis is post hoc because earlier
TNG scores had already been inspected. Both transfer analyses require original
simulation catalogues or cutouts for full reruns.

Within-output validation and output-weighting sensitivity on 2026-11-06
-----------------------------------------------------------------------
The within-output test holds out one of the 12 retained z=0.1904 haloes and
selects b from the other 11 using the median of the three view-specific
maximizing values for each training halo. Every fold selects b=0.12 in V and
b=0.10 in g. Fold parameters and per-view and per-host scores are in
within_output_validation_z019_20261106/. The folds share training haloes and
one output and are not independent-volume validation.

The output-weighting sensitivity uses the corrected-distance curves. Each
output contributes its median-over-halo IoU curve with equal weight. The
result selects b=0.10 in both bands, compared with entry-level medians of
b=0.08 in V and b=0.10 in g. The script, output-level table, summary, and
input/output hashes are in output_weighting_sensitivity_20261106/.

TNG300-1 merger-tree growth analysis, 2026-11-06
-------------------------------------------------
The R1--R4 revision adds per-view isophote-detection flags and full versus
detection-qualified growth comparisons. Per-epoch mass recovery, paired
growth errors, observer-condition-dependent half-light growth, and projected
half-mass growth are in tng300_growth_bias_20261106/growth_r1r4_audit_20261106/.
The 36-branch input cutouts, source dependencies, and rerun instructions are
in tng300_growth_bias_20261106/README.md. This separate TNG package includes
derived star-particle cutouts and an in-memory stellar-linking library. The
earlier statement that no particle catalogue or compiled binary is bundled
describes the Darwin source archive and does not apply to the TNG growth data.

An exact ParticleID ledger decomposes central-subhalo mass changes, and
controlled dust-free SSP renderings measure the contribution of spatial age
and total-metallicity gradients to half-light-radius growth. The ledger,
population controls, radial profiles, paired summary script, CSV tables, and
interpretation limits are in
tng300_growth_bias_20261106/assembly_population_20261106/ and its results/
subdirectory. From snapshot 50 to 84, pre-existing stars newly assigned to
the central are the largest median ledger term; most had the resolved-
satellite mask at snapshot 50. This mask does not identify a donor or a
physical merger. The modelled population gradients lower paired median
half-light-radius growth by 0.255 dex in V and 0.324 dex in g at physical
observer redshifts. These results use one TNG300-1 volume and dust-free SSP
weights.
