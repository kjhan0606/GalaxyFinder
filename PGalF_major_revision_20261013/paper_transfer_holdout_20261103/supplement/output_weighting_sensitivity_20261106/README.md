# Output-weighting sensitivity

This analysis checks whether the Darwin calibration changes when each of the
six simulation outputs contributes the same weight. It uses the corrected
luminosity-distance curves and the 75 halo-output entries retained by the
paper's all-view, both-band interior-maximum criterion.

For each entry, the script takes the median IoU over the three orthogonal
projections at each tested linking parameter. The paper's reported parameter
is the median of the entry-level median-over-view best parameters. For the
equal-output sensitivity, the script first takes the median IoU curve within
each output, then averages the six output curves with equal weights. The
selected parameter maximizes that curve; ties select the smaller grid value.
The table also records each output's median entry-level best parameter and
its output-specific curve maxima.

The analysis is descriptive. The outputs come from one simulation volume and
include repeated systems, so output weighting does not provide an independent
sample, a redshift trend, or a confidence interval.

## Inputs and reproduction

Inputs are `../halo_metadata.tsv` and `../recomputed_curves.tsv`. The script
uses only rows labelled `consistent_distance` in the curve table. Run:

```bash
python3 analyse_output_weighting.py
```

The script writes `output_weighting_by_output.csv`,
`output_weighting_summary.json`, and `SHA256SUMS`. From this directory, verify
the recorded inputs and outputs with:

```bash
sha256sum -c SHA256SUMS
```
