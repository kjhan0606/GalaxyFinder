# Satellite sample and measurements: methods draft

*Working text for the new satellite-galaxy study. This document records the
implemented measurement and analysis definitions; catalogue counts and all
results must be filled from the audited NewGalFinder outputs.*

## Simulation outputs and galaxy catalogue

We analyse four closely spaced outputs from the Darwin-1 cosmological hydrodynamical simulation: snapshot 00972 at $z=0.22258$, 00977 at
$z=0.21784$, 00980 at $z=0.21324$, and 00984 at $z=0.20870$. Their redshift span
is only $0.01388$, allowing us to check whether the measured relation is
reproduced across outputs but not to measure redshift evolution. The outputs
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

## Stellar mass, metallicity, and internal velocity dispersion

We calculate each satellite's stellar mass by converting the NewGalFinder
stellar mass from units of $M_\odot/h$ to $M_\odot$ using the snapshot's
Hubble parameter. NewGalFinder stores stellar-particle velocities in
\si{km\,s^{-1}}, positions in comoving $\mathrm{Mpc}/h$, and particle masses
in $M_\odot/h$; we use the stored physical velocities directly. We calculate
the three-dimensional stellar half-mass radius from the periodic positions
and convert it to physical kpc using the output scale factor and Hubble
parameter. For each galaxy, we calculate the current-stellar-mass-weighted
mean of the RAMSES stellar birth-metallicity field $z_p$. We treat
$z_p$ as a simulation metallicity proxy and do not convert it to
$[\mathrm{Fe/H}]$ without a verified solar normalization. A number-weighted
mean is retained as a sensitivity measure.

NewGalFinder also retains each stellar particle's RAMSES birth-time field.
We convert this time coordinate to proper age by interpolating the
output-specific Friedmann expansion history, following the method used by the
RAMSES \texttt{part2birth} utility \citep{Teyssier2002}. For each galaxy, we
calculate the current-stellar-mass-weighted mean age and record the fraction
of its stellar particles with a valid birth time.

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

## Sample cuts and statistical comparison

For the main analysis, we require at least 100 member stars with finite
velocities and positive current masses, matching the sample used for the
mass-weighted global dispersion; we also require positive, finite galaxy
stellar mass and metallicity proxy. The half-mass-radius measurement
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
metallicity proxy. The model with present-day subhalo mass answers a distinct
conditional question: since tidal stripping can change that mass and the
stellar velocity dispersion \citep{FattahiEtAl2018}, attenuation after this
control is not interpreted as evidence that stripping has no role.

As a population-age sensitivity check, we fit another model that adds the
logarithm of mean stellar age to the host- and stellar-mass-controlled model.
We report its whole-host bootstrap interval and residual scatter. This tests
whether the metallicity association remains among satellites with similar mean
stellar ages; it does not identify a causal role for age.

As a structural sensitivity analysis, we add the logarithm of the stellar
half-mass radius to the host- and stellar-mass-controlled model. This tests
whether the metallicity association remains among galaxies of similar size.
Because size may respond to the same processes as dispersion, we treat this
as a conditional comparison rather than a causal mediation test, and report
its whole-host bootstrap interval and residual scatter.

We report host-clustered standard errors from the fixed-effects fits and
percentile 95 per cent confidence intervals from 2,000 bootstrap resamples
of whole host systems. We analyse the four outputs separately; repeated
galaxies across nearby outputs are treated as repeated
measurements and are not pooled as independent systems. We compare the global
and half-mass-radius results, current-mass- and number-weighted metallicity,
stellar-mass thresholds, and the alternative subhalo-mass controls as
prespecified robustness checks. We assess physical importance from effect
size and residual scatter, not from a non-zero coefficient alone.

## Results to insert after catalogue validation

Replace this note only after the FoF and NewGalFinder chains have completed
and the extracted catalogues have passed the particle-ID, particle-count,
finite-value, mass-conversion, and host-assignment audits. Insert per-output
catalogue counts, fitted coefficients and intervals, effect sizes, residual
scatter, and robustness results. Do not use the legacy ASCII catalogue for
the primary claims.
