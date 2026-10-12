Calculation source archive, 2026-10-13

The directory preserves the calculation sources and stellar FoF implementation
with their relative source layout. File hashes are in
supplement/calculation_config.json. The archive is supplied for inspection and
reproduction with the original particle inputs. No raw catalogue, flux cache,
compiled library, or historical measurement file is bundled here.

Darwin1/sb_vs_b.py reads GALFIND binary star records, reconstructs parent lists,
removes satellites, calculates stellar light, constructs images, and calls the
stellar FoF library. The function named brightest selects the catalogue galaxy
with the largest star-particle count. The function name is historical.

Darwin1/analysis_sb/revision/sb_revise.py contains image variations and the
corrected boundary score. The boundary score restricts both original selections
to the band. Archived rev_*.txt files were generated with the earlier incorrect
band comparison and are not outputs of the corrected function.
sb_fine.py contains finer-grid, PSF-fraction, and fixed-image sampling operations.
sb_dist.py measures separations and light fractions on either side of the contour.

major_revision/scripts preserves the exact particle driver and frozen helpers
used for the 2026-10-13 recalculation, plus the aggregation script. The driver
reproduces every original curve and group count, then measures the corrected
boundary score and images with a consistent luminosity distance. The recorded
task layout uses inputs/supplement for metadata and flux grids, and results for
photometry_audit.json and per-halo output. Absolute raw-catalogue and linker
paths must be configured before reproduction outside the original environment.
The portable scripts/audit_photometry.py supplies the distance factors.
Do not run the archived driver against existing output files
without first selecting a new output directory.

opFoF/variant contains the C adapter, tree walker, and build recipe. The required
project headers are included in opFoF. The tree walker includes MPI interfaces.
Build in opFoF/variant with an MPI C compiler using make CC=mpicc, or the compiler
appropriate to the local MPI installation. The historical shared library hash
is recorded separately in calculation_config.json. The binary is not distributed.
The in-memory entry point labels particles by chains of links using the mean
of the two particle linking distances. Python selects the component containing
the density seed after receiving labels.

The Python sources retain absolute paths from the producing environment.
Reproduction elsewhere requires setting ROOT, REV, SPS_HOME, the catalogue paths,
and the output directory, and providing the binary catalogue files, their offset
tables, and the original measurement tables used by the reproduction checks.
Raw stellar coordinates must be unwrapped about the selected central galaxy
before calling the in-memory linker. No periodic wrapping occurs inside the C
entry point. Large particle calculations require a compute server.

For the manuscript figures and statistical summaries, use the portable drawing
and diagnostic scripts one directory above. The archived tables suffice for
the drawing and diagnostic commands, and no catalogue read or C compilation is required.

Additional audit sources in the Astra revision
The source opFoF/ramses2read.c documents the b=0.20 species normalization.
Darwin1/namelist_00967.txt preserves the actual simulation settings.
references/legacy_halos.py records the external DECam halo model at commit
ac5c2f982c358bdfd6a50b46b02884e3cdb7f6e5. The derivative profile uses the
g/r amplitudes and 7-to-8 arcsec transition from that source.
