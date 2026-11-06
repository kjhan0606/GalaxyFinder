# Fixed-parameter CAMELS image checks

The sample contains the same 49 hosts from seven Swift-EAGLE CV boxes as
the nominal transfer test. The stellar linking parameters remain fixed at
0.10 and 0.20. The projected 30 proper kpc aperture uses the same simulated
potential-minimum centre. No parameter is fitted to the external sample.

Nine conditions comprise the nominal image, limits brighter or fainter by
0.5 mag, 0.4 arcsec pixels, no outer PSF wing, a positive flat sky residual
of 29 AB mag/arcsec^2, and Gaussian pixel noise with rms equal to the
threshold pixel flux divided by 5, 3, or 1. Only one image choice changes
per condition. Noiseless conditions have one realisation. Noisy conditions
have three realisations, with paired standard-normal noise patterns across
amplitudes. The observing redshift is 0.1904379623 for the fixed z=0
structures. Stellar ages, masks, and particle photometry remain fixed.

The perturbations change the isophotal target. Luminosity and half-light
radius use noiseless particle fluxes at projected positions, rather than
noisy aperture sums. The calculation tests boundary sensitivity and does
not test the photometric uncertainty of an observational pipeline. The
calculation omits dust, correlated sky errors, Poisson counting statistics,
and image-based centre finding. High-redshift intergalactic absorption is
not tested by the low-redshift image checks.

## Files

- `projection_metrics.csv.gz` contains 13,230 projection-level records.
  The repeated baseline of the stronger-noise extension is excluded after
  exact numerical comparison. Each unique host/condition/realisation/band/
  method has three views. Failed central targets have zero overlap and
  completeness. Precision and luminosity/radius errors are undefined for
  failed targets. The present sample has no failed target in any condition.
- `host_realisation_metrics.csv` gives three-view medians for each
  realisation. Failure fractions use the mean of the three views.
- `host_metrics.csv` averages the realisation summaries within each host.
  Viewing directions and noise realisations are repeated measurements.
- `sensitivity_summary.csv` gives equal-host means, conditional error
  quantiles, and separate failure summaries for 54 condition/band/method
  combinations. Error quantiles describe measurable targets only.
- `paired_condition_changes.csv` gives changes relative to the nominal
  image for each selection. `paired_method_contrasts.csv` compares fixed
  0.10 with 0.20 or the aperture. Both use the same 5,000 host and
  box-then-host bootstrap draws, with seed 20261030.
- `fixed_b_image_sensitivity.pdf` compares mean overlap and 90th-percentile
  luminosity errors. Its shared legend is inside panel (a), in the blank strip
  above the curves.
- `configuration_*.json`, `metadata_*.json`, `export_validation.json`, and
  `measurement_sources/` retain the frozen design and generating sources.
  The strong-noise extension inherited a descriptive string mentioning
  noise divisors 5 or 3. The extension's active parameter `[1.0]`, raw
  metrics, and generating source define the actual divisor 1.

The column `contamination` is the selected-light fraction outside the
isophote, equal to one minus precision. Columns `logL` and `logR50` are
base-10 logarithms of selected-to-isophotal luminosity and radius ratios.
Positive values indicate a larger selected luminosity or radius. The
`abs_` columns are absolute logarithmic errors in dex. A noisy host
summary averages three view-median log errors before the host-level
quantile is evaluated. `b_star` is an archived data-column name for the
stellar linking parameter, displayed as b_link in the figures.

## Reproduce the tables

From the paper directory, with numpy and pandas installed, run

```bash
OPENBLAS_NUM_THREADS=1 python scripts/summarize_camels_image_checks.py --output /tmp/camels-image-reproduction
```

Use a new output directory. The script reads the compressed projection
table and needs no original snapshots, yt, or stellar-population library.
Five exported tables, including every paired confidence interval, were
checked against the generating analyses and a compressed-CSV reanalysis
at relative tolerance 1e-10 and absolute tolerance 1e-12. The archived
measurement sources retain the original analysis environment and paths.
Rerendering images additionally requires the original snapshots,
photometric grids, and helper modules.

The image conditions were designed after the nominal transfer results.
The original fixed-parameter transfer remains the calibration-independent
test. The additional conditions map sensitivity without retuning the
linking parameter. A larger overlap or a smaller error in a perturbed
image describes a changed observational target, rather than a physical
improvement to the stellar group.
