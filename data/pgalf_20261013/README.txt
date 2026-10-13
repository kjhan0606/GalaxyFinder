PGalF final-analysis data and code release, 2026-10-13

The canonical analysis uses simulation distances for both angular size and
stellar flux, a 15-node linking grid, positional and PSF aperture comparisons,
and inverse-probability reweighted sampling. Start with supplement/README.txt
and scripts/README.txt. Original nine-node and other older products remain
for provenance and are not inputs to the final headline statistics.

The released products contain aggregate curves for all 76 output entries,
validation predictions, boundary bands, stellar-mass and mean-age comparisons,
image controls, failure flags, conditional resampling results, and compact
corrected drawing inputs. Full particle catalogues, full membership caches,
full per-halo images, manuscript TeX and manuscript PDFs are absent.

Check hashes from the release root with
sha256sum -c SHA256SUMS

Population independence, hydrodynamical resolution convergence, kinematic
accuracy, realistic survey noise, dust attenuation, and IGM absorption are
not established by the supplied idealized-image calibration. Known positional
pairs are grouped only for conditional resampling. Sampling reweighting
preserves mass and light in expectation, rather than exactly per draw.

Raw reconstruction requires the original catalogues, cache layout and compiled
linking library. Source files and library hash document the producing setup.
Contact the corresponding author for catalogue access.
