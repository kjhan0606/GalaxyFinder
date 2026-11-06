CAMELS Swift-EAGLE transfer calculations used by this manuscript

These files were copied byte-for-byte from
/gpfs/kjhan/Others/Jaehyun/2026.Apr/PGalF_blind_transfer_20261029/.
SHA-256 hashes are recorded in the parent supplement/SHA256SUMS manifest.
The sample_manifest.csv contains the 49 selected CV-box/group rows.
run_camels_blind_batch_centered.py records the calculation tag
blind20261030_v4_potential_centered and catalogue-ID membership source.

Measurement definitions checked against measure_camels_halo.py:
- Each view projects stellar positions onto a two-dimensional plane. The
  30 pkpc comparison is a circular projected aperture (projected radius <=
  30 pkpc), not a three-dimensional sphere.
- R50 is the projected radius at which the cumulative full stellar flux of
  the selected set first reaches half its total. FoF, aperture, and isophotal
  selections each use their own selected-set flux total.
- The luminosity error is abs(log10(L_selected/L_contour)); the radius error
  is abs(log10(R50_selected/R50_contour)). The reported host value is the
  median over the three projections; the paper reports the 90th percentile
  across hosts.
- The PSF-convolved image is used to find the isophotal region. Flux sums and
  projected R50 use full particle fluxes at projected particle positions,
  without PSF-aperture fractions.

To reproduce the table summaries from the original per-halo metrics, run
summarize_camels_blind_batch.py with --manifest sample_manifest.csv,
--tag blind20261030_v4_potential_centered, --target-variant potential,
--results pointing to the original metrics directory, and a new empty
--outdir. The exact historical paths are embedded in the original measuring
script and must be configured for another environment. Rebuilding the per-halo
metrics requires the public CAMELS snapshot and group-catalogue HDF5 files,
the Darwin SSP grid, and the PGalF image/measurement source archive. Those
large inputs are not included here.
