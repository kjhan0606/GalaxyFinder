# Per-sightline curves

These tables are the observer-frame scores behind the stellar linking-length
calibration. There are 76 halos and three orthogonal sightlines, in the V
band and in SDSS g.

`sightline_baseline.tsv` is the adopted measurement. Each row is one
sightline at one linking length. The nine values of `b` are 0.02, 0.03,
0.05, 0.08, 0.10, 0.12, 0.15, 0.18 and 0.20. `iou`, `recall` and `purity`
use the positional whole-star flux. `area_arcsec2` is the adopted isophote.
`b_best` is the first maximum on this grid. `edge=1` means that maximum
sits on a grid end. An endpoint is not a matched linking length.

`sightline_variants.tsv` gives the overlap at the same nine values of `b`
after each change of the image or the isophote. The tag `base` repeats the
adopted overlap at higher precision. The other tags change one threshold,
add a flat sky, remove the outer wing, replace the adopted isophote by the
component that contains the density-peak seed, drop the nearest-component
fallback, return satellites to the image, or rebuild the image from a
nested subset of the stars.

`sightline_refined.tsv` gives overlap, recall and purity on the refined
grid and for the thinning that holds the adopted isophote fixed. `POS` and
`PSF` use fifteen nodes, the nine above plus 0.06, 0.07, 0.09, 0.11, 0.13
and 0.14. The PSF weight is the fraction of the kernel inside the isophote
and does not enter the link. `thin2`, `thin4` and `thin8` use the original
nine nodes.

The last comment line of each file names the columns.
