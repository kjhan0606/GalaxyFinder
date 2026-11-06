# Within-output validation at (z=0.1904)

This test isolates prediction within one Darwin output. The sample contains
12 halo-output entries with interior overlap maxima in all three projections
and both photometric bands. The test excludes the (z=0.227) catalogue,
which uses a different satellite-selection setting.

For each held-out host, the script selects a linking parameter from the other
11 hosts. The script finds each training host's best parameter separately in
the three projections, takes the median across projections, and then takes
the median across training hosts. The script scores the held-out host at the
selected parameter in every projection. Fixed parameters (b=0.10) and
(b=0.20) provide comparison rows. The image and overlap scores use the
corrected Darwin-cosmology luminosity distances.

The 12 leave-one-out folds share most training hosts and come from one output
in one simulation volume. The results are descriptive and do not provide an
independent-volume estimate or a confidence interval.

## Inputs and reproduction

Inputs are `../halo_metadata.tsv` and `../recomputed_curves.tsv`. Run:

```bash
python3 analyse_within_output_cv.py
```

The script writes per-fold parameters, per-host summaries, per-view scores,
and `summary.json`. From this directory, verify the recorded files with:

```bash
sha256sum -c SHA256SUMS
```
