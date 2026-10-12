PGalF stellar linking-length calibration — revision 2026-10-13

Primary measurement
The primary photometry uses the simulation cosmology for angular size and
luminosity distance. Multiply archived apparent fluxes by the per-output factor
in photometry_audit.json BEFORE constructing images and selecting regions.
The fixed nine-node grid is 0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20.
Do not combine the archived finer grid with corrected-distance selections.

Join tables by (step, nstar), adding los and band for individual views. Entries
identify halos at a particular output, not merger-tree identities. Three views
and two bands give 456 curves for 76 entries. Exactly 75 entries have interior
maxima in every view; all 76 remain in fixed-parameter validation.
The catalogue at output 00967 used stellar linking disabled. Other outputs used
stellar linking enabled. Reassembling parents does not remove differences in
satellite assignment between catalogues.

Primary files
recomputed_curves.tsv contains full-precision positional scores and light sums
at every original node. variant=consistent_distance is PRIMARY; variant=baseline
is the archived-distance control. The three light sums are L_group, L_inside,
and L_both. The boundary sums restrict every stellar set to the same pixel band.
primary_view_selection.tsv supplies corrected maxima, endpoint flags and areas.
primary_output_summary.tsv and primary_summary.json give the primary statistics.
validation_corrected_folds.tsv and validation_corrected_predictions.tsv record
validation on the primary grid. Training includes all other-output entries after
excluding periodic neighbours within 0.5 h^-1 cMpc of a held-out entry. Training
has no fitted-interior filter. Evaluation always includes all 76 entries.
validation_output_summary.tsv reports gains and failures at every output.
A positional exclusion cannot guarantee progenitor independence.

Targets and catalogues
halo_metadata.tsv retains original row identifiers, byte offsets, group/stellar
masses, parent/satellite counts, and archived fitted-interior flags. Row indices
start at zero. Masses are h^-1 Msun. halo_positions.tsv supplies catalogue centres
in h^-1 cMpc. catalogue_audit.tsv reconstructs the original parent counts and
masses. catalogue_galaxies.tsv records all 58995 galaxies with galaxy_index,
n_stars, current_mass, selected, centres, and archived relative AB fluxes.
The selected central has the largest star-particle count and largest stellar
mass in every entry, and the largest flux in 52 entries in each band.
target_membership.tsv uses the cached seed index and catalogue particle ordering
to establish the seed's galaxy. All 76 seeds belong to the selected central.
The seed is a number-density maximum under a fixed-neighbour radius definition.
A mass-density peak or mass-weighted image centre need not coincide with the seed.

Archived controls
sightline_baseline.tsv, sightline_variants.tsv and sightline_refined.tsv preserve
the originally deposited records at commit 5e3cc51. Legacy column names recall
and purity mean completeness and precision respectively. Preserving the raw
records is intentional. New tables and figures use completeness and precision.
POS and PSF in the refined table use fifteen nodes; thin2/thin4/thin8 use nine.
PSF scores assign each star's convolved flux fraction inside the region.
All refined/PSF/sampling/sky/threshold/satellite controls use ARCHIVED distances.
The retained sample is used for BOTH linking and score sums in fixed-image
thinning; individual masses and luminosities are unchanged. The seed is chosen
again in the retained sample. No full-luminosity recovery is inferred.
thinning_all_entries.tsv reports upper/lower endpoint counts and fixed-b scores
without conditioning on a modified interior maximum. At factor 8, 45 of the
75 original entries per band have upper-endpoint maxima and only 30 remain.
sensitivity_summary.tsv preserves conditional paired changes and the original
view count conventions. Published nodes/flags are retained for historical controls.
validation_folds.tsv, validation_predictions.tsv and validation_summary.json
are ARCHIVED-distance finer-grid validation controls with interior-filtered
training. Use validation_corrected_* for the main manuscript.

Boundary and radial measurements
boundary_ring.tsv is the corrected FORMULA applied to ARCHIVED images. It
replaces the erroneous legacy comparison of full group G with ring R.
The valid formula is L(G & I & R)/L((G | I) & R), with inside mask I.
The primary boundary values are in consistent_distance rows of
recomputed_curves.tsv and primary_summary.json, not boundary_ring.tsv.
All boundary widths are one 0.2-arcsec pixel layer on each side, with corner
contacts included. Zero denominators produce NaN, not zero. Halo medians remain
undefined if any view ratio is undefined. Quantiles omit undefined halo values.
boundary_distance.tsv stores ARCHIVED image light near isophotes at published
best b. L_miss denotes exterior group light; L_extra denotes interior non-member
light. le1/le2/le5 give Euclidean pixel-centre separations in arcseconds.
radial_luminosity.tsv and luminosity_summary.json use ARCHIVED normalization.
The annuli have edges 0, 0.5, 1, 1.5 and infinity in units of equivalent radius.
The primary absolute luminosity errors are regenerated from corrected light
sums by primary_diagnostics.py. Positional luminosity is not PSF aperture flux.
Absolute errors are taken per view before halo and sample medians.

Photometry and image diagnostics
fsps_*.npz are immutable archived SSP grids, with age, logz, flux, and z.
Flux is per solar mass formed; catalogue formed mass is divided by h first.
Bilinear log-flux interpolation uses log age and log metallicity, clipping
outside the grid. All 5376 archived grid values reproduce exactly with the
recorded FSPS installation. fsps_build.json and calculation_config.json give
libraries, hashes, cosmology, image parameters, and population settings.
Interstellar attenuation and IGM absorption are absent. AGB circumstellar dust
remains enabled. filter_v.tsv and filter_sdss_g.tsv provide transmissions.
ssp_clipping.tsv and ssp_clipping_summary.json count particles outside the age
and metallicity grids and measure their clipped flux fractions. The audit uses
all 76 non-satellite input lists, not unique stars across epochs. No error versus
a finer/younger SSP grid is inferred from the clipped-light fraction.
rest_frame.tsv records controls without dimming and with unredshifted spectra.
The control changes sampled wavelength as well as dimming. Angular image areas
are retained even though the control brightness uses flux per proper pc^2.
image_diagnostics.tsv contains archived tagged records MASS, CMP and PHO.
MASS records particle-mass and linking-distance quantiles.
CMP records components, areas, selection fallback and centre-seed offsets.
PHO records positional and convolved fluxes for baseline and satellite stars.
Baseline PHO totals already exclude satellites; do not subtract satellites again.
mass_correlations.tsv preserves ARCHIVED selected-b mass correlations, not the
primary corrected calibration. Nominal p-values do not establish evolution.

Reproduction and access
scripts/README.txt gives commands. primary_diagnostics.py and plot_scores.py
need no particle catalogues. The source archive contains production readers,
image/measurement routines, stellar FoF C sources and headers, a build recipe,
the simulation namelist, and the external PSF source/version.
Raw particle catalogues, cached membership and compiled FoF library are absent.
Reconstructing images or membership requires those original inputs and adjusted
absolute paths. Contact the corresponding author for catalogue-access requests.
cluster_illustration.npz and failure_panels store ARCHIVED image drawing inputs;
cluster/failure figures are identified as controls. Primary score and boundary
figures use consistent-distance tables. SHA256SUMS verifies the supplied files.
The Overleaf manuscript package also contains extended_notes.tex and PDF.
