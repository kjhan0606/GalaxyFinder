# Stellar metallicity and internal kinematics of satellite galaxies

**Working manuscript.** The results, tables, abstract, and conclusions below are
not submission-ready until the validated NewGalFinder catalogues and analyses
for all four outputs are complete. Bracketed items are deliberate placeholders;
no exploratory legacy-catalogue values are carried into this draft.

## Abstract

[To be written after catalogue validation and the four-output relation analysis.
State the simulation sample, the current-mass-weighted mean stellar
birth-metallicity mass fraction, the global and half-mass-radius internal
one-dimensional dispersions, the host and stellar-mass controls, the
finder-assigned subhalo-mass and stellar-age controls, the measured effect size
and uncertainty, and the conclusion about physical importance.]

## 1. Introduction

Stellar metallicity records the chemical enrichment of the gas from which a
galaxy's stars formed. A galaxy's internal stellar velocity dispersion
measures the spread of stellar motions. Both depend on stellar mass, and both can also
change in a dense environment. At fixed stellar mass, satellite galaxies can
be more metal-rich than field galaxies, possibly because environmental gas
loss reduces the supply of fresh, metal-poor gas \citep{BaheEtAl2017}. Tidal
stripping can reduce a satellite's dark-matter mass and internal stellar
velocity dispersion \citep{FattahiEtAl2018}. These effects mean that a
metallicity--dispersion relation must be tested after accounting for stellar
mass and environment.

Previous studies have examined metallicity and kinematics with different
samples and questions. In SDSS, Gallazzi et al. (2021) found that satellites
below $10^{10.5}\,M_\odot$ become more metal-rich as host halo mass
increases at fixed stellar mass. Much of the satellite--central difference
reflects the higher passive fraction among satellites
\citep{GallazziEtAl2021}.
Using SDSS galaxies, Spindler \& Wake (2017) studied
how stellar metallicity varies with core velocity dispersion in central and
satellite galaxies. After separating star-forming and quiescent systems, they
found little metallicity difference between centrals and satellites
\citep{SpindlerWake2017}. In
the NIHAO XV simulations, Buck et al. (2019) found that satellite and
field dwarfs follow nearly the same stellar mass--metallicity and stellar
mass--velocity-dispersion relations, while satellite and nearby dwarfs have
lower velocity dispersions overall \citep{BuckEtAl2019}.
Fattahi et al. (2018) directly examined the metallicity--velocity relation of
Local Group dwarfs using $V_{1/2}$, estimated from luminosity-weighted
line-of-sight stellar dispersion. They found that stripped satellites flatten
the relation at low velocities. Their correction for tidal mass loss assumes
that stripping leaves each satellite's mean stellar metallicity unchanged;
with that assumption, the inferred progenitors follow a tighter relation
\citep{FattahiEtAl2018}. These studies connect stellar metallicity with galaxy
kinematics and document environmental differences. We ask a complementary
question: among satellites in the same cosmological host, is metallicity
associated with internal stellar dispersion after controlling for stellar
mass and present-day, finder-assigned subhalo dark-matter mass?

Riley et al. (2026) also show that substantial tidal stellar-mass loss can
coexist with the small scatter in the Local Group satellite mass--metallicity
relation; agreement with that relation alone therefore does not rule out
tidal processing \citep{RileyEtAl2026}.

Here we ask whether a satellite's mean stellar birth metallicity, stored as
the total metal mass fraction $z_p$, is associated with its internal stellar
velocity dispersion at fixed stellar mass and host. We then add the
dark-matter mass assigned by the galaxy finder
as a separate control. This mass is measured at the same epoch as the
dispersion, so the additional model asks whether metallicity remains
associated with dispersion at fixed present-day finder-assigned mass. Because
tidal stripping can change both quantities, a weaker relation in this model
would not show that stripping has no role. We measure dispersion across all
member stars and within the three-dimensional stellar half-mass radius, and
report the effect size, uncertainty, and residual scatter.

We analyse four closely spaced outputs from the Darwin-1 cosmological hydrodynamical
simulation. Their small redshift span lets us check whether the result is
reproduced across outputs, but not whether it evolves with redshift. The
stellar-particle field $z_p$ stores total metallicity as a mass fraction;
it is not an iron abundance, so we do not convert it to
$[\mathrm{Fe/H}]$. We use the
galaxy finder's dark-matter mass as a finder-defined present-day control and
do not independently recompute its particle membership. The analysis measures
conditional associations; it cannot by itself establish that metallicity
causes a change in velocity dispersion or identify the physical process
behind any association.

## 2. Simulation and methods

### Simulation outputs and galaxy catalogue

We analyse four closely spaced outputs from the Darwin-1 simulation. Table 1 lists the
scale factor, redshift, and RAMSES coarse-step number from each output's
metadata. Their redshift span is only $0.01388$, allowing us to check whether
the measured relation is reproduced across outputs but not to measure redshift
evolution. The outputs
come from a $44.171\,h^{-1}\,\mathrm{Mpc}$ comoving box with
$H_0=67.4\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$, $\Omega_m=0.315$,
$\Omega_\Lambda=0.685$, and $\Omega_b=0.045$.
The RAMSES metadata reports configured minimum and maximum AMR levels of 10
and 17, respectively, in each output. We identify friends-of-friends (FoF) host
systems from the dark-matter and stellar particles supplied to NewDD. Gas and
sink particles are omitted from this extraction. Within each FoF system, we
classify the member with the largest stellar mass as the central galaxy and
the remaining members as satellites. The FoF host masses therefore represent
dark-matter-plus-star particle masses, rather than total masses including gas
and sinks.

NewGalFinder assigns particles to density basins, then applies an iterative
tidal-radius cut and a self-binding energy test to satellite candidates. It
also applies a stellar linking step to the remaining member stars. We use the
resulting catalogue assignments without independently rerunning these
membership procedures. We therefore call the subhalo mass the
*finder-assigned dark-matter mass*; it is a finder-defined present-day mass,
not an independently measured bound mass. We also retain the
finder-assigned dark-matter-plus-star mass as a separate control. Because
NewDD uses the absolute value of a signed particle ID as its tracking key, we
audit particle identities using $|\mathrm{ID}|$. Within each FoF host, we
exclude galaxies with repeated keys or keys shared with another galaxy. We
reject an output if any key occurs in more than one FoF host.

**Table 1. Metadata for the four simulation outputs.**

| Snapshot | Scale factor | Redshift | RAMSES coarse step |
|---:|---:|---:|---:|
| 00972 | 0.817942 | 0.222581 | 82328 |
| 00977 | 0.821128 | 0.217836 | 82514 |
| 00980 | 0.824237 | 0.213243 | 82693 |
| 00984 | 0.827336 | 0.208698 | 82870 |

### Stellar mass, metallicity, and internal velocity dispersion

We calculate each satellite's stellar mass by converting the NewGalFinder
stellar mass from units of $M_\odot/h$ to $M_\odot$ using the snapshot's
Hubble parameter. NewGalFinder stores stellar-particle velocities in
\si{km\,s^{-1}}, positions in comoving $\mathrm{Mpc}/h$, and particle masses
in $M_\odot/h$; we use the stored physical velocities directly. We calculate
the three-dimensional stellar half-mass radius from the periodic positions
and convert it to physical kpc using the output scale factor and Hubble
parameter. For each galaxy, we calculate the current-stellar-mass-weighted
mean of the stellar-particle total metallicity mass fraction $z_p$, recorded
at star formation. We do not convert $z_p$ to $[\mathrm{Fe/H}]$, because
$z_p$ is total metallicity rather than iron abundance. A number-weighted mean
is retained as a sensitivity measure.

NewGalFinder also retains each stellar particle's RAMSES birth-time field. We
convert this time coordinate to proper age by interpolating the output-specific
Friedmann expansion history, following the method used by the RAMSES
\texttt{part2birth} utility \citep{Teyssier2002}. For each galaxy, we
calculate the current-stellar-mass-weighted mean age and record the fraction
of its stellar particles with a valid birth time. We use mean stellar age as
an additional control in the internal-dispersion analysis.

Our primary kinematic measurements are the one-dimensional stellar velocity
dispersions. We calculate the global value from all member stars with finite
velocities and positive current masses. We calculate a second value using
only stars inside the three-dimensional, stellar-mass-weighted half-mass
radius. We determine that radius from the periodic positions of stars with
finite positions and positive current masses, centred on their stellar-mass-
weighted centroid, and take the radius enclosing half of their total current
stellar mass. For a selected set of stars $S$, the mass-weighted mean velocity
of component $j$ is $\bar v_j=\sum_{i\in S}m_i v_{ij}/\sum_{i\in S}m_i$;
the reported one-dimensional dispersion is
\[
\sigma_{\star,\mathrm{1D}}=\left[\frac{1}{3}\sum_{j=x,y,z}
\frac{\sum_{i\in S}m_i(v_{ij}-\bar v_j)^2}{\sum_{i\in S}m_i}
\right]^{1/2}.
\]
For the global measurement, $S$ contains all member stars with finite
velocities and positive current masses; for the aperture measurement, $S$
contains only those stars also lying within the three-dimensional half-mass
radius. We report these internal stellar dispersions separately from the
host-scale orbital motion of satellite galaxies.

As a separate host-scale test, we split the satellites in each selected FoF
host at the median birth metallicity (at least eight satellites with valid
centre velocities), then pair the lower and upper halves by stellar mass
(difference at most $0.2$ dex; at least four pairs). For each half, we
calculate the one-dimensional orbital dispersion from the matched
galaxy-centre velocities,
$\sigma_{\rm orb,1D}=[(s_x^2+s_y^2+s_z^2)/3]^{1/2}$, with equal weight per
satellite. We report the median across hosts of
$\log_{10}(\sigma_{\rm orb,highZ}/\sigma_{\rm orb,lowZ})$ and its 95 per cent
interval from 2,000 whole-host bootstrap resamples. This comparison is
distinct from each galaxy's internal dispersion; by itself it neither measures
the host potential profile nor identifies a cause.

### Sample cuts and statistical comparison

For the main analysis, we require at least 100 member stars with finite
velocities and positive current masses, matching the sample used for the
mass-weighted global dispersion; we also require positive, finite galaxy
stellar mass and positive, finite mean stellar birth metallicity. The
half-mass-radius measurement
additionally requires at least 30 stars with finite velocities and positive
current masses inside that radius. We analyse stellar-mass
thresholds of $10^8$ and $10^9\,M_\odot$. We select hosts with
$10^{13}\le M_{\rm host}<10^{14}\,M_\odot$, where $M_{\rm host}$ is the FoF
dark-matter-plus-star particle-mass sum described above. We test sensitivity
by shifting each host-mass boundary up and down by $0.2$ dex while holding the
other boundary fixed. These cuts are not full total-halo-mass thresholds.

The current catalogue does not provide a validated recent star-formation rate,
so we do not divide the sample by star-formation activity.

For each output and each kinematic measure, we fit the following relation by
ordinary least squares, giving each satellite equal weight:
\[
\log_{10}\sigma_i = \alpha_{h(i)} + \beta_Z\log_{10}z_{p,i}
  + \beta_\star\log_{10}M_{\star,i} + \epsilon_i,
\]
where $\alpha_{h(i)}$ is a host-specific intercept. The coefficient $\beta_Z$
therefore measures the metallicity association at fixed stellar mass and
host. We then add the logarithm of the finder-assigned subhalo dark-matter
mass as a separate control; a further model uses the finder-assigned
dark-matter-plus-star mass instead. We report the fitted coefficient,
uncertainty interval, sample and host counts, and residual scatter. To make
the magnitude interpretable, we convert each slope to the predicted
percentage change in dispersion across a 0.3-dex difference in the
mean stellar birth metallicity. The model with present-day subhalo mass answers a distinct
conditional question: since tidal stripping can change that mass and the
stellar velocity dispersion \citep{FattahiEtAl2018}, attenuation after this
control is not interpreted as evidence that stripping has no role.

As a population-age sensitivity check, we fit another model that adds the
logarithm of mean stellar age to the host- and stellar-mass-controlled model.
We report its whole-host bootstrap interval and residual scatter. This tests
whether the metallicity association remains among satellites with similar mean
stellar ages; it does not identify a causal role for age.

We report host-clustered standard errors from the fixed-effects fits and
percentile 95 per cent confidence intervals from 2,000 bootstrap resamples
of whole host systems. We analyse the four outputs separately; repeated
galaxies across nearby outputs are treated as repeated
measurements and are not pooled as independent systems. We compare the global
and half-mass-radius results, current-mass- and number-weighted metallicity,
stellar-mass thresholds, and the alternative subhalo-mass controls as
prespecified robustness checks. As a structural sensitivity analysis, we also
add the logarithm of the stellar half-mass radius to the host- and stellar-
mass-controlled model. This asks whether the metallicity association remains
among galaxies of similar size; because size may itself respond to the same
processes as dispersion, we do not interpret this adjustment as a causal
mediation test. We report its whole-host bootstrap interval and residual
scatter alongside the other models. We assess physical importance from effect
size and residual scatter, not from a non-zero coefficient alone.

## 3. Results

### 3.1 Catalogue validation and sample

[Insert per-output host and galaxy counts; stellar-particle ID audit; catalogue
stellar-mass versus particle-mass consistency; and counts after each mass,
particle-number, and host-mass cut. Reject any output that fails the global ID
check.]

### 3.2 Metallicity and internal stellar dispersion

[Show the raw relation, then the equal-weight host-fixed-effects fit controlling
for stellar mass. Report the slope, host-bootstrap interval, the predicted
percentage change in dispersion for a 0.3-dex increase in mean stellar birth
metallicity, and residual scatter. Repeat the fit with the number-weighted
mean birth metallicity and report its effect size and host-bootstrap interval.]

### 3.3 Age, structural, and subhalo-mass controls and robustness

[Report the age-adjusted fit and the size-adjusted fit, then the fits adding
finder-assigned subhalo dark-matter mass and, separately, finder-assigned
dark-matter-plus-star mass. For each, give the effect size, whole-host
bootstrap interval, and residual scatter. Compare global and half-mass-radius
dispersions, both stellar-mass thresholds, host-mass boundaries shifted one at
a time by +/-0.2 dex, and the four outputs. Treat the outputs as replication
checks, not as a redshift-evolution sequence.]

### 3.4 Satellite orbital dispersion within hosts

[Report the hosts with mass-matched pairs, exclusions for nonpositive
dispersions, the total and typical matched-pair counts, the median
high-to-low-metallicity orbital-dispersion ratio, its whole-host bootstrap
interval, and the residual stellar-mass imbalance. Keep this host-scale
comparison separate from the internal stellar-dispersion fits.]

## 4. Discussion

[Interpret only the associations supported by the audited measurements. Compare
the effect size with the residual galaxy-to-galaxy scatter and existing work.
Make literature comparisons qualitative: Local Group studies commonly use
stellar $[\mathrm{Fe/H}]$ and line-of-sight dispersions, whereas we use total
birth-metallicity mass fraction and three-dimensional particle dispersions in
more massive hosts.
Discuss tidal stripping and other possible mechanisms as interpretations, not
causal conclusions, unless additional evidence distinguishes them. State the
limits of using total birth-metallicity mass fraction rather than an iron
abundance, the finder-assigned present-day mass,
the DM-plus-star-only host definition, and the closely spaced outputs.]

## 5. Conclusions

[Answer whether the metallicity association is reproduced across outputs,
remains after the specified mass controls, and is large enough to matter
relative to the residual scatter. Do not infer physical importance from a
non-zero coefficient alone.]

## Data availability

[Complete after confirming the simulation-data access conditions.]

## Acknowledgements

K.J.H. provided the MPI+OpenMP codes used for halo finding and galaxy finding.
OpenAI Codex (GPT-6) assisted with developing the satellite-catalogue
extraction and validation, statistical-analysis, and figure-generation
scripts, and with drafting and editing the manuscript. Claude Opus 5.5 and
GPT-6 Astra were used to obtain referee-style critiques of manuscript drafts.
K.J.H. directed these AI-assisted tasks. The authors are responsible for
checking the AI-assisted outputs and for the final analysis and manuscript.

## References

Bibliographic entries are maintained in `references.bib`; format them for the
selected journal after the target journal is confirmed.
