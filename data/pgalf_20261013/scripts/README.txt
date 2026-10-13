Reproducing the final PGalF analysis

Run from the package root. The final table-based commands require Python,
NumPy, SciPy, pandas, and matplotlib. Producing versions and numerical
settings are recorded in supplement/opus_calculation_config.json.

python scripts/aggregate_new.py --data supplement --out regenerated
python scripts/conditional_uncertainty.py --data supplement --out regenerated
python scripts/plot_revision.py --data supplement --out figures
python scripts/plot_cluster.py

The first command reproduces fine-grid selections, all-entry common-b scores,
output-held-out validation, radius/redshift strata, sampling summaries,
wide-band ratios, failure flags, mass correlations, and stellar-population
summaries from the released flat tables. The bootstrap uses the ten known
positional pairs and 66 groups of entries. Unidentified repeated progenitors
are not removed. Aggregate products use positional or PSF weights explicitly.

plot_revision.py regenerates the corrected parameter, score, completeness-
precision, failure-map, photometry/sampling, and unredshifted-control figures.
Failure maps use bounded cutouts. plot_cluster.py uses compact saved drawing
arrays with the corrected contour. The unchanged mass map and stellar point
geometry do not depend on luminosity-distance normalization. A uniform flux
factor also leaves the relative flux-weighted sampling probabilities unchanged.

The source archive retains original readers, compiled-library build sources,
and producing scripts. Full particle reconstruction requires the original
catalogues/cache layout, configured project paths, and an appropriate compute
host. The producing source and library hashes accompany the configuration.
Image and membership reconstruction is not possible from aggregate tables.
The release excludes full per-halo images, full group caches, and raw catalogues.

Older diagnostic and plotting scripts reproduce prior products and are not
the canonical entry points for the final manuscript. The fixed Git history
preserves the previous analysis. Use the commands above for final figures.
