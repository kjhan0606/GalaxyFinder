# Satellite metallicity–kinematics study: manuscript plan

This is a planning document, not a results manuscript. The PGalF surface-
brightness calibration paper has a separate scope. Do not insert exploratory
statistics from the legacy ASCII catalogue as final results.

## Scientific question

At fixed satellite stellar mass and host, is stellar birth metallicity
associated with (1) the satellite's internal stellar velocity dispersion and
(2) the host-scale orbital velocity dispersion of its satellite population?
For the internal relation, also test the NewGalFinder-assigned subhalo
dark-matter-mass control and whether a residual is associated with stellar
size. For the host-scale comparison, contrast the orbital dispersions of
stellar-mass-matched, metal-rich and metal-poor satellites within each host.
Keep these two velocity measures separate. A statistical relation alone will
not establish a physical cause. Host-mass cuts use the available FoF sum of
dark-matter and star particles; they are not full host masses because the
NewDD input deliberately omits gas and sink particles.

The metallicity--kinematics connection itself is not new. Using SDSS data,
Spindler & Wake (2017) compared central and satellite galaxies across
stellar-mass, stellar-metallicity, and core-velocity-dispersion relations;
they found almost no central--satellite stellar-metallicity difference after
separating star-forming and quiescent galaxies (\cite{SpindlerWake2017}). More
directly, Fattahi et al. (2018) examined $[\mathrm{Fe/H}]$ versus $V_{1/2}$ for
Local Group dwarfs. They found that tidal stripping can flatten the relation
at low velocity and that inferred progenitor corrections restore a tighter,
more monotonic trend (\cite{FattahiEtAl2018}). The present study must be
positioned as a conditional test in a different simulation regime, not as the
first metallicity--dispersion study. Separately, EAGLE work found stellar-
metallicity enhancement in satellites relative to equal-mass field galaxies
and examined environmental gas-removal and enrichment processes
(\cite{BaheEtAl2017}); this motivates controlling for host environment but
does not by itself establish a residual relation with internal stellar
dispersion. NIHAO XV compared the stellar mass--metallicity and stellar
mass--velocity-dispersion relations of simulated satellites and field dwarfs
separately, finding broad overlap between the populations and lower satellite
dispersions at fixed stellar mass (\cite{BuckEtAl2019}). It did not test the
within-host metallicity association at fixed stellar and present-day subhalo
mass targeted here. The present study's potential contribution is to quantify
whether a metallicity association with **internal stellar** dispersion
remains at fixed stellar mass, host, and finder-assigned subhalo mass, and
whether metal-rich and metal-poor satellites in the same host have different
orbital velocity dispersions after stellar-mass matching. This host-scale
comparison is distinct from the internal-dispersion analysis. Any claimed
novelty must be limited to what those controls and the audited measurements
actually establish.

The recent Auriga Streams III analysis shows that the observed small scatter
in the Local Group satellite mass--metallicity relation can coexist with
substantial tidal stellar-mass loss (\cite{RileyEtAl2026}). This further
motivates treating the mass--metallicity relation as context rather than as a
test that rules out tidal processing; it does not directly measure the
conditional metallicity--dispersion association studied here.

Interpret the present-day subhalo-mass control as a distinct conditional
question, not as a causal isolation of metallicity's effect. Tidal stripping
can lower a satellite's internal velocity dispersion while leaving its mean
stellar metallicity nearly unchanged; it can therefore create or reshape a
metallicity--dispersion trend. Because stripping also changes the present-day
subhalo mass, conditioning on that mass can remove part of this physical
pathway. Report the host-and-stellar-mass-controlled relation alongside the
additional subhalo-mass-controlled relation. Without a validated bound-mass
history or pre-infall measurements, describe the latter as the residual
association at fixed present-day finder-assigned mass, and do not treat its
attenuation as evidence that tides or other physical processes are absent.

Keep the kinematic definitions distinct when comparing with Fattahi et al.:
their $V_{1/2}=\sqrt{3}\,\sigma_{\rm los}$ is based on luminosity-weighted
line-of-sight stellar velocities, whereas this analysis uses a mass-weighted
one-dimensional dispersion of stellar particles, globally and within the
three-dimensional stellar half-mass radius.

Keep internal stellar velocity dispersion distinct from the orbital velocity
dispersion of satellite galaxies within their host. The internal measure is
the mass-weighted one-dimensional dispersion of member stars; report the
global value and the value within the three-dimensional stellar half-mass
radius separately. The host-scale check compares the one-dimensional orbital
dispersion, calculated from the three velocity-component variances of galaxy
centres, for the upper and lower metallicity halves of each host after
one-to-one stellar-mass matching. Treat it as a separate outcome, not as a
measurement of the host potential profile.

## Proposed manuscript structure

1. **Introduction** — Establish why satellite stellar kinematics and metallicity
   can respond differently to stellar mass, host environment, and tidal mass
   loss. Position the conditional relation relative to the existing literature.
   State both kinematic questions and explain that internal stellar motions
   and host-scale satellite orbits are measured separately.
2. **Simulation and outputs** — Describe the cosmological hydrodynamical run,
   the four selected outputs, resolution, and relevant simulation units. The
   selected outputs are snapshots 972, 977, 980, and 984, at redshifts
   `0.222581`, `0.217836`, `0.213243`, and `0.208698`, respectively, computed
   from the scale factors in their RAMSES `info` files. Their total span is
   `Δz = 0.013883`, so they support catalogue replication, not a
   redshift-evolution claim.
3. **Satellite catalogue** — Document NewDD input, FoF host identification,
   NewGalFinder assignments, central/satellite definition, halo and particle
   counts, mass cuts, and repeated star-ID audit. State that FoF host masses
   and finder-assigned subhalo particle masses omit gas and sinks in this
   extraction, and that subhalo dark-matter mass is a finder-assigned particle
   sum unless an independent boundness validation is established.
4. **Measurements and statistical methods** — Define stellar mass, the
   current-mass-weighted mean birth-metallicity proxy, each internal-dispersion
   aperture, host-orbital comparison, host and subhalo mass proxies, sample
   cuts, host fixed effects, matched comparisons, and host-level bootstrap
   intervals. State that the internal regressions use ordinary least squares
   with equal weight per satellite. Do not label the stored metallicity as
   `[Fe/H]` unless its absolute solar normalization is verified.
5. **Results** — First report catalogue validation and sample sizes. Then show
   the unadjusted relation, the relation after controlling for stellar mass
   and host, and the primary result after adding subhalo dark-matter mass.
   Present effect sizes in dex and percent, uncertainty intervals, and
   between-galaxy scatter. Follow with aperture, mass-cut, host-bin, age, and
   size sensitivities. Report the separate high-versus-low-metallicity orbital
   dispersion comparison with its host count, stellar-mass matching balance,
   effect size, and host-bootstrap interval. Convert the NewGalFinder birth-time
   field to proper age using the output cosmology and the RAMSES part2birth
   method; report age-controlled results only after checking valid-time
   coverage and the age distribution in all four catalogues. Recent SFR and
   activity classes remain out of the primary analysis.
6. **Discussion** — Compare the conditional result with stellar mass–metallicity,
   size–dispersion, and tidal-stripping expectations. Discuss a physical
   mechanism only if the audited measurements distinguish its predictions;
   otherwise describe the correlation as conditional and leave its cause open.
7. **Conclusions** — Answer whether the relation is reproducible, independent
   of the tested mass controls, and large enough to matter relative to the
   intrinsic scatter. State the limits imposed by the metallicity proxy,
   snapshot spacing, and catalogue construction.

## Planned figures and tables

- **Figure 1:** satellite sample in stellar mass and host mass, with the
  applied cuts and per-output counts.
- **Figure 2:** stellar metallicity versus internal stellar dispersion, showing
  raw and mass-conditioned relations; separate global and half-mass-radius
  dispersions.
- **Figure 3:** percentage change in internal stellar dispersion for a
  0.3-dex increase in the birth-metallicity proxy, with whole-host bootstrap
  intervals. Show host-plus-stellar-mass control beside the additional
  finder-assigned dark-matter-mass control for each mass threshold and
  dispersion aperture. Treat the four outputs as separate checks, not a
  redshift-evolution sequence.
- **Table 3 or a compact supplemental panel:** the within-host orbital
  dispersion ratio for stellar-mass-matched high- and low-metallicity
  satellites. Include a main-text panel only if the validated result is robust
  and materially changes the physical interpretation.
- **Figure 4:** stellar-size diagnostic, clearly marked as a secondary test
  rather than causal evidence. Omit the planned activity split unless recent
  SFRs are derived and validated for the primary catalogue.
- **Table 1:** simulation and output properties.
- **Table 2:** catalogue validation and sample counts at every cut.
- **Table 3:** primary coefficients, bootstrap intervals, effect sizes, and
  residual scatter for the prespecified kinematic measures.

## Final MNRAS submission checks

- Prepare the submission manuscript with the MNRAS LaTeX class. Keep the
  abstract as one paragraph of at most 250 words and select one to six key
  words from the journal's approved list. Provisional choices are `galaxies:
  abundances`, `galaxies: clusters: general`, `galaxies: kinematics and
  dynamics`, `galaxies: evolution`, and `methods: numerical`; confirm the final
  set against the completed results.
- Number sections and place the unnumbered Acknowledgements, Data Availability,
  and References after the Conclusions. Cite every figure and table in order.
- Use vector artwork for quantitative plots; label every axis and unit, check
  line weights and lettering at final print size, and use colours that remain
  distinguishable for colour-blind readers.
- Include concise alt text for every figure directly below its caption in the
  manuscript file. Make sure the captions and surrounding text explain the
  plotted data without requiring the reader to infer undocumented details.
- Supply machine-readable catalogue tables and column descriptions as online
  supporting material if they are part of the published results.

These checks follow the [MNRAS Instructions to Authors](https://academic.oup.com/mnras/pages/General_Instructions).

## Evidence required before drafting claims

- FoF and NewGalFinder complete successfully for the selected outputs; verify
  catalogue structure, host assignments, finite values, particle counts, and
  repeated star IDs.
- Recompute every reported relation from the audited NewGalFinder extraction,
  not from the legacy `Mdm` column or its undocumented metallicity field.
- Confirm the internal-dispersion result under fixed stellar mass, host
  identity, and finder-assigned subhalo mass; quantify sensitivity to aperture,
  resolution cuts, host-mass boundaries shifted by `+/-0.2 dex`, and repeated
  systems. Independently validate the within-host orbital comparison, including
  the usable-host count and stellar-mass balance of matched pairs.
  Because the FoF input omits gas, report host cuts as dark-matter-plus-star
  particle-mass cuts rather than full halo-mass thresholds.
- Report effect size and residual scatter alongside statistical significance.
  Do not call a small non-zero slope physically important on the basis of a
  confidence interval alone.
- Restrict multi-output claims to replication or consistency because the four
  redshifts are closely spaced.
