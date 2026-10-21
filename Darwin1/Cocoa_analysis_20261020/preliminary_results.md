# Preliminary satellite metallicity–kinematics result

Status: exploratory result from one existing ASCII catalogue, not a paper result.

## Kinematic scope and FoF update (2026-10-21 18:43 KST)

The primary internal-stellar-dispersion analysis and the host-scale orbital
comparison are now specified separately in the methods and manuscript plan.
For the orbital comparison, the implemented analysis splits satellites at
each host's metallicity median, matches the high- and low-metallicity halves
one-to-one within 0.2 dex in stellar mass, requires at least four pairs per
host, and bootstraps whole hosts. Its summary now records the hosts with
mass-matched pairs, hosts excluded for nonpositive dispersion, and total,
median, and minimum pair counts. It will be reported separately from the
per-galaxy internal dispersion; no result is available until the validated
catalogues are extracted.

At this checkpoint, 00972, 00977, and 00980 remain `RUNNING`; 00984 is
`PENDING (Resources)`, and downstream jobs remain dependency-pending. The
latest MPI CPU counters were 2-13:58:11 (00972), 2-23:59:33 (00977), and
2-14:20:15 (00980), with maximum resident memory near 111 GB for 00972/00980
and 95 GB for 00977. The 00972 and 00980 files/logs last advanced at 18:37
KST; the 00977 output/log timestamps remain 17:19/17:34 KST. Slurm still
reports all three jobs `RUNNING`. The 00977 task has elapsed 21 h 42 min of
its 48 h limit, with an active MPI step; its FoF files are not complete. The
00984 task remains `PENDING (Resources)`, and all downstream catalogue and
analysis jobs remain dependency-pending. No task has exited successfully.

## Darwin-1 FoF liveness update (2026-10-21 18:29 KST)

The 00972, 00977, and 00980 FoF array tasks remain `RUNNING`; 00984 remains
`PENDING (Resources)`, and all NewGalFinder, extraction, relation, host-mass,
and figure tasks remain dependency-pending. Snapshot 00972's halo catalogue
grew to 37,398,836 bytes at 18:25 KST, and 00980's grew to 39,335,412 bytes
at 18:30 KST; both logs advanced with their outputs. The 00977
FoF files and log have not advanced since 17:19 and 17:34 KST, but its active
MPI step's accumulated CPU time increased from 2-23:02:56 to 2-23:15:28 over
the latest check; this supports ongoing computation, not a confirmed hang.
No FoF task has exited successfully, so no downstream catalogue or physical
relation is ready to interpret.

The scheduler detail for 00984 confirms that it requests 10 CPU ranks on two
nodes and no GPU; it remains resource-pending, with a projected start of
2026-10-22 21:00 KST. The older 00977 task was submitted with 60 CPUs and
768 GB total reserved memory and has used about 21 h 25 min of its 48 h limit.
Its MPI step is still active, so it was left running.

## Darwin-1 FoF and dependent-analysis status (2026-10-21 18:23 KST)

Slurm still reports FoF tasks for 00972 (430610_0), 00977 (430418_1), and
00980 (430611_2) as `RUNNING`; 00984 (430622_3) is `PENDING (Resources)`.
The 00972 and 00980 MPI steps report about 2.5 days of aggregate CPU time and
roughly 111 GB maximum resident memory each. Their logs advanced at 17:51 and
18:15 KST, respectively. The 00977 log last advanced at 17:34 KST; its job is
still live, so the unchanged log alone does not establish a failure. The FoF
files remain intermediate until successful job exits. NewGalFinder, satellite
extraction, relation analysis, host-mass sensitivity, and figure tasks are
still dependency-pending. No primary satellite catalogue or validated
metallicity–dispersion result is available yet; do not interpret the partial
FoF files as completed catalogues.

## Darwin-1 FoF liveness check (2026-10-21 14:14 KST)

Slurm still reports 00972, 00977, and 00980 as `RUNNING`; their live MPI CPU
counters increased to 1-15:52:03, 2-09:10:03, and 1-16:06:31. Their halo and
member-particle files remain unchanged since 13:54:29, 13:49:34, and 13:47:59,
respectively. The latest logs also remain unchanged since 14:02:08, 13:49:34,
and 14:09:11. FoF 00984 remains priority-pending, and its downstream finder,
catalogue, analysis, and figure jobs remain dependency-pending. No primary
satellite list is ready.

## Relation-figure pipeline checked (2026-10-21 14:11 KST)

Reviewed `make_satellite_relation_figure.py` and its Slurm wrapper against the
analysis output fields and manuscript plan. The figure requires all four
snapshot relation CSVs, checks that each requested model's host-bootstrap
interval is finite and ordered around its estimate, and plots the host-plus-
stellar-mass model alongside the dark-matter and DM-plus-star controls. Its
four panels cover both stellar-mass cuts and both global and half-mass-radius
dispersions; the four nearby outputs are categorical checks rather than a
redshift trend. Expected field names match `analyse_satellite_relations.py`.
The figure job remains dependency-pending until those CSVs exist; no code
change was needed.

## Darwin-1 FoF progress (2026-10-21 14:10 KST)

The three FoF tasks remain `RUNNING` and their MPI-step CPU counters continue
to rise. No halo/member output file has changed since the previous checkpoint.
The 00980 log advanced at 14:09:11: in its current boundary-check phase, P2
and P5 have completed their checks; their transfer lines have not yet appeared,
and other ranks have not logged completion for this phase. The 00972 log
remains at 14:02:08 with P2's boundary check complete, while 00977 remains at
13:49:34 with P5 still missing from the last root-consolidation transfer set.
00984 remains priority-pending; all galaxy-finder and analysis jobs remain
dependency-pending. No primary satellite catalogue is ready.

## Darwin-1 FoF liveness check (2026-10-21 14:06 KST)

FoF 00972, 00977, and 00980 remain `RUNNING` in Slurm; their live MPI-step CPU
counters have increased to 1-15:15:04, 2-08:45:16, and 1-15:29:21,
respectively. No FoF file or stdout log has advanced since the previous
recorded checkpoints (00972: 13:54:29 files / 14:02:08 log; 00977: 13:49:34;
00980: 13:47:59 files / 13:55:03 log). Snapshot 00984 remains priority-pending
with no FoF outputs. The NewGalFinder, extraction, relation, sensitivity, and
figure tasks remain dependency-pending; no primary satellite catalogue is
available.

## Darwin-1 FoF work-stage detail (2026-10-21 14:04 KST)

The 00972 log updated at 14:02:08. All ten ranks had completed tree building
for the current shard, and P2 then completed its boundary check on 19.6 million
particles. The same shard includes substantially larger slabs on P8 (142.5
million particles), P7 (93.9 million), and P9 (74.6 million); their boundary
checks and transfers have not yet been recorded in the current phase. This
explains why files may remain unchanged while the live job continues scanning
large slabs. The 00977 and 00980 logs remain at 13:49:34 and 13:55:03. All
three FoF jobs are still live; 00984 and their downstream jobs remain queued.
No primary galaxy catalogue has been produced.

## Darwin-1 FoF status (2026-10-21 14:02 KST)

All three active FoF jobs remain `RUNNING`; MPI-step CPU counters increased to
2-08:29:46 (00977), 1-14:51:56 (00972), and 1-15:06:07 (00980). FoF files
remain at the previous sizes, with latest updates at 13:54:29, 13:49:34, and
13:47:59, respectively. Snapshot 00984 still has no FoF output and remains
priority-pending. NewGalFinder, extraction, relation, sensitivity, and figure
jobs remain dependency-pending. There is no completed primary catalogue to
analyse yet; the live FoF jobs were left untouched.

## Darwin-1 FoF liveness check (2026-10-21 14:00 KST)

Slurm continues to report FoF 00972, 00977, and 00980 as `RUNNING`. Their
MPI-step CPU counters increased again (1-14:46:16, 2-08:25:58, and
1-15:00:27), while the FoF output files remained unchanged from their last
updates at 13:54:29, 13:49:34, and 13:47:59. The latest logs remain at
13:57:23 for 00972, 13:49:34 for 00977, and 13:55:03 for 00980. The 00984
FoF task is still priority-pending and all downstream catalogue and analysis
jobs remain dependency-pending. No primary satellite catalogue is ready.

## Darwin-1 extraction validation reviewed (2026-10-21 13:58 KST)

Code review of `extract_satellite_properties.py` confirms that the stellar-ID
audit uses NewDD's `abs(ID)` key, bit-packs the full signed-32-bit key range,
counts repeats within and across galaxies in each FoF host, and aborts the
snapshot if a key reappears in a different FoF host. It writes to a temporary
file and only publishes the CSV after the entire catalogue has been read, so a
failed audit cannot leave a partial file presented as complete. The existing
extract/analysis chain therefore has an explicit catalogue-integrity gate; no
code change was indicated by this review.

## Darwin-1 FoF progress (2026-10-21 13:58 KST)

The 00972, 00977, and 00980 FoF tasks remain `RUNNING`, with increasing MPI
CPU counters; 00984 is still priority-pending. The 00972 log advanced to
13:57:23, when all ten ranks had passed tree building for the next shard.
The 00980 log remains at 13:55:03, and 00977 at 13:49:34; neither has logged
the next stage since then. FoF output files have not changed since 13:54:29
(00972), 13:49:34 (00977), and 13:47:59 (00980). No successful FoF exit or
primary satellite catalogue is available, so all dependent analysis jobs
remain pending.

## Darwin-1 FoF progress (2026-10-21 13:57 KST)

The three FoF jobs remain `RUNNING`; their CPU counters are still increasing,
and none of the FoF files has changed since the previous check. Logs show
ongoing work: for 00972, ranks P2, P3, P4, P6, P7, and P9 have passed tree
building for the current shard by 13:56:36; for 00980, P2 has completed its
boundary check while other ranks were still in tree building at 13:55:03.
The 00977 log remains at 13:49:34 and still lacks P5's transfer marker in the
current root-consolidation phase. Slurm continues to report 00984 as
priority-pending and all downstream catalogue/analysis/figure jobs as
dependency-pending. No primary satellite catalogue is ready.

## Darwin-1 FoF progress (2026-10-21 13:54 KST)

Slurm still reports 00972, 00977, and 00980 as `RUNNING`; their MPI-step CPU
counters continue increasing. Snapshot 00972's halo and member-particle files
grew at 13:54:29, and its log updated at 13:54:39. The log shows every rank
sending its isolated-halo data in the preceding batch, followed by reads from
the next particle batch, so FoF processing has advanced rather than ended.
Snapshots 00977 and 00980 have no newer log or file updates than 13:49:34 and
13:50:26, respectively, although both Slurm tasks remain live. 00984 remains
priority-pending. No FoF task has yet completed successfully, so dependent
NewGalFinder, extraction, and analysis jobs are still waiting and no primary
satellite catalogue is ready.

## Darwin-1 FoF status (2026-10-21 13:53 KST)

Slurm still reports 00972, 00977, and 00980 as `RUNNING`, with MPI-step CPU
counters at 1-14:10:04, 2-08:01:42, and 1-14:24:04, respectively. The FoF
files have not changed since 13:34 (00972), 13:49 (00977), and 13:47 (00980).
The latest available 00977 consolidation log still lacks P5's transfer marker;
00980's latest log shows multiple ranks building trees. 00984 has no FoF
catalogue files and remains priority-pending. All downstream tasks remain
dependency-pending, so no primary satellite catalogue or relation output is
available for physical interpretation.

## Darwin-1 FoF progress (2026-10-21 13:52 KST)

All three active FoF tasks remain `RUNNING`, and their MPI-step CPU counters
increased since the previous check. No FoF file size changed during that
interval. The 00977 log remains at 13:49:34; P5 is still the only rank without
a transfer marker in the current consolidation phase. The 00980 log advanced
to 13:50:26 and now reports ranks P3, P4, P6, P7, P8, and P9 passing tree
building for its current shard. The 00972 log remains at 13:34:12. These are
active processing stages, not catalogue completion. FoF 00984 remains
priority-pending, and the downstream galaxy-finding and relation-analysis
arrays remain dependency-pending.

## Darwin-1 FoF liveness check (2026-10-21 13:51 KST)

The three active FoF tasks remain `RUNNING`: 00977 has elapsed 16 h 51 min,
and 00972 and 00980 have each elapsed 7 h 41 min. Their live MPI-step CPU
counters are increasing (2-07:53:29, 1-13:57:50, and 1-14:11:47,
respectively), so unchanged output timestamps do not indicate a stopped job.
The latest file updates remain 13:49:34 for 00977, 13:47:59 for 00980, and
13:34:12 for 00972. In the 00977 consolidation phase, P5 is still the only
rank without a transfer marker. The 00984 FoF task remains priority-pending;
all downstream galaxy-finder, extraction, analysis, and figure tasks remain
dependency-pending. No primary catalogue is ready for analysis.

## Darwin-1 FoF progress (2026-10-21 13:49 KST)

Slurm continues to report FoF 00972 (`430610_0`), 00977 (`430418_1`), and
00980 (`430611_2`) as `RUNNING`; FoF 00984 (`430622_3`) remains pending for
priority. The 00977 halo and member-particle files grew again at 13:49:34 to
57,212,740 and 151,835,899,824 bytes. Its log now records P9 completing the
current root-consolidation transfer; P5 is the remaining rank not yet recorded
as sent, and no root-completion marker is present. The 00980 files last grew
at 13:47:59, and its log advanced at 13:49:45; its MPI CPU counter continues
to rise. The 00972 files and log remain unchanged since 13:34, while its live
MPI CPU counter also continues to rise. All three products remain intermediate
until their FoF jobs exit successfully. NewGalFinder, extraction, relation,
host-mass-sensitivity, and figure tasks remain dependency-pending; no primary
satellite catalogue or relation result exists yet.

## Darwin-1 FoF progress (2026-10-21 13:48 KST)

Slurm still reports the 00972 (`430610_0`), 00977 (`430418_1`), and 00980
(`430611_2`) FoF tasks as `RUNNING`; 00984 (`430622_3`) is `PENDING` for
priority. The 00980 FoF halo and member-particle files grew at 13:47:59, and
its log advanced at 13:48:24 into the next tree-building stage. It is not a
completed catalogue. The 00972 files/log are unchanged since 13:34, and the
00977 files/log since 13:40; both tasks remain live. NewGalFinder, extraction,
relation analysis, host-mass sensitivity, and figure jobs remain dependency-
pending. No primary satellite-properties or relation CSV exists yet.

## Darwin-1 FoF run rechecked (2026-10-21 13:46 KST)

The active Slurm query confirms that snapshots 00972 (`430610_0`), 00977
(`430418_1`), and 00980 (`430611_2`) are still running; 00984 (`430622_3`)
is pending for priority. In the latest 00977 consolidation phase, ranks P1,
P2, P3, P4, P6, P7, and P8 have sent their isolated-halo data to P0. P5 and
P9 have not yet logged that transfer, and the root completion marker is absent.
The 00977 log and FoF files have not changed since 13:40 KST. The 00972 and
00980 logs and files likewise remain unchanged since 13:34 and 13:15 KST,
respectively. Their output files under `FoF_Data/FoF.00972`, `FoF.00977`, and
`FoF_Data/FoF.00980` are still intermediate because Slurm has not reported
successful completion. No primary satellite catalogue or relation table is
available, and the dependent NewGalFinder and analysis tasks remain queued.

## Primary-analysis definitions audited (2026-10-21)

The queued pipeline and methods draft agree on the planned comparison: identify
the most stellar-massive NewGalFinder object in each FoF group as its central,
retain the other objects as satellites, reject satellites with repeated stellar
particle IDs, and apply the predeclared stellar- and host-mass cuts. The
primary metallicity is the current-stellar-mass-weighted mean birth-metallicity
field; the kinematic outcomes are the global and three-dimensional half-mass-
radius mass-weighted one-dimensional stellar dispersions. The regression
compares log dispersion with log metallicity while controlling for stellar
mass and host identity, then adds finder-assigned subhalo dark-matter mass as a
separate conditional model. The outputs are analysed separately because their
redshifts are closely spaced. The paper plan correctly avoids treating the
conditional association as causal and notes that host mass omits gas and sink
particles. No method change is needed before primary catalogues arrive.

## Primary Darwin-1 catalogue run (2026-10-21 13:44 KST)

The primary sample uses snapshots 00972, 00977, 00980, and 00984; the exploratory
ASCII-catalogue numbers below are not substituted for those results. FoF jobs
for 00972, 00977, and 00980 are still running. Snapshot 00984's FoF job remains
pending on the `a10` partition with scheduler reason `Priority`; Slurm now
projects a start at 2026-10-22 21:00:08 KST, subject to rescheduling.
The scheduler identifies the 00977 array task as JobId 430419 and reports it
RUNNING on three nodes. Its stdout log has not been flushed since 12:47 KST,
but the live MPI-step CPU counter advanced to 2-07:31:44 by 13:44 KST,
confirming that the calculation is still active. Jobs 00972 and 00980 also
remain RUNNING, with live task CPU counters of 1-13:25:19 and 1-13:39:13,
respectively. Non-empty FoF halo
and member-particle files are present for 00972,
00977, and 00980, but their Slurm tasks are still RUNNING; these files are not
treated as complete products until the jobs exit successfully. The matched
NewGalFinder, extraction, and relation-analysis jobs remain dependency-pending.
No primary GALFIND catalogue, satellite-property table, or relation result is
available yet.

At this check, FoF 00977 had run for 16 h 44 min of its 48-h limit. FoF 00972
and 00980 had each run for 7 h 34 min of their 24-h limits.

I attempted to extend the running 00972 and 00980 jobs to the current script's
48-h limit, but Slurm denied both updates (`Access/permission denied`). Their
24-h limits remain unchanged, with the scheduler end time at 2026-10-22
06:10:13 KST.

FoF 00972 last showed fresh output at 13:34 KST: its halo catalogue grew to
22,817,612 bytes and its member-particle file to 32,953,460,728 bytes. The
latest log shows rank 0 writing an intermediate batch while other ranks save
their isolated halos. These files and the log had not changed again by 13:38
KST, while the MPI CPU counter continued to rise. Slurm still reports the job
RUNNING, so these appended files remain intermediate rather than completed
FoF catalogues.

FoF 00977 produced a newer intermediate batch at 13:40 KST. Its halo catalogue
grew from 54,853,444 to 55,820,100 bytes, and its member-particle file grew
from 141,569,775,536 to 145,877,555,120 bytes. Rank 0 was writing the batch;
seven other ranks had logged their saves while the remaining ranks were still
in that phase. No newer files were present by 13:43 KST, although MPI CPU
counters continued to increase. The job remains RUNNING, so these files are
also intermediate.

The pending 00984 FoF task had been submitted to `a100` despite requesting no
GPU. I moved that pending task (430622_3) to `a10`, matching the partition in
the current CPU-only FoF script, and extended its time limit from 24 to 48 h.
Its allocation remains 2 nodes, 10 CPU tasks, and 448 GiB total (224 GiB per
node); the current script requests 3 nodes and 256 GiB per node for 48 h. The
submitted per-node memory fits the 512-GiB A10 nodes. Slurm projects a start
after the current A10 FoF reservation ends, but the estimate can change.
This avoids scheduling the CPU-only task on A100. NewGalFinder still requires
A100 memory and waits on successful FoF completion; the active 00972 and 00980
FoF jobs remain on A100 until they finish. The 00984 task may also wait for
the currently occupied A10 nodes to free up.

The opFoF source confirms that the 100 input shards are divided across the 10
MPI ranks; each rank processes its assigned shards in sequence, appending
intermediate halo/member data and exchanging boundary particles. Final-halo
writing occurs only on each rank's last shard. Thus, existing non-empty FoF
files are expected intermediate output, and the advancing MPI CPU counters
with unchanged file timestamps do not indicate that the jobs have finished.
The latest flushed per-rank shard messages are 00007 for all ranks in 00977
and 00003 for all ranks in 00972 and 00980 (local shard indices within each
rank's ten-shard block). These log positions are stale snapshots rather than
current progress markers: the corresponding stdout files were last modified
at 12:47, 11:33, and 13:15 KST, respectively, and no completion ETA can be
inferred from them.

## Input and extraction

- Simulation: Darwin-1. Historical input path recorded for this exploratory catalogue: `cocoa_01_run_20261019/catalogue/sn0972/galaxy_catalogue_rcut05_972.txt`. The directory label is a storage name, not the simulation name.
- The first line gives `a = 0.817941556`, or `z = 0.22258124`.
- The file has 112,848 galaxy rows. The final `host flag` column was read as
  1 for the central and 0 for satellites, yielding 16,500 satellite rows.
- Validation in hosts with `M200c >= 10^13 Msun` found exactly one flag-1
  central in each of 53 hosts, 7,692 flag-0 satellites, no other flag values,
  and no repeated galaxy IDs in the input catalogue.
- Extracted galaxy and host IDs, host `M200c`, stellar mass, catalogue stellar
  metallicity `Z(star)`, three stellar velocity-dispersion measures, rotation
  speed, mass-weighted age, and stellar half-mass radius.
- Extraction and analysis scripts are `extract_existing_satellites.py` and
  `analyse_satellite_relations.py`. Machine-readable outputs are the
  `satellite_properties_ascii_v2_00972.csv` table and
  `preliminary_satellite_physics_v2_00972.csv` summaries.

## Main provisional comparison

Restricting to hosts with `10^13 <= M200c/Msun < 10^14`, positive catalogue
metallicity, and satellites with `Mstar >= 10^8 Msun` gives 3,628 satellites
in 51 hosts. At `Mstar >= 10^9 Msun`, there are 1,343 satellites in the same
51 hosts.

- The unadjusted Spearman coefficient between `log Z(star)` and the global
  internal `log sigma` is +0.637 for the `10^8 Msun` cut.
- A host-fixed-effect regression controlling linearly for `log Mstar` gives
  slope `-0.201 +/- 0.009` dex in `log sigma` per dex in `log Z(star)`.
- Resampling the 51 host haloes gives a 95 per cent interval
  `[-0.219,-0.184]` for that slope.
- Controlling for host and `0.3 dex` stellar-mass bins gives
  `-0.162 +/- 0.009` using 3,434 satellites in 51 hosts. The aperture
  dispersion within the stellar half-mass radius gives a consistent slope
  `-0.174 +/- 0.009` in the continuous-mass model.
- The host-bootstrap interval for the binned-mass slope is
  `[-0.179,-0.143]`.
- At `Mstar >= 10^9 Msun`, the corresponding global-dispersion slopes are
  `-0.136 +/- 0.017` (continuous mass control) and
  `-0.112 +/- 0.023` (0.3-dex mass bins).
- The fixed-host, fixed-stellar-mass relation between metallicity and
  three-dimensional stellar half-mass radius is also negative
  (`-0.700 +/- 0.020` for the `10^8 Msun` cut). Mass-weighted age is lower
  at higher metallicity (`-0.190 +/- 0.015`). The global rotation-to-dispersion
  ratio has no significant trend at this mass cut.

Thus the raw positive metallicity–dispersion trend reverses after controlling
for host and stellar mass. This is a provisional residual relation, not yet a
physical interpretation.

## Matched host-mass definition for comparison with NewGalFinder

The initial legacy-catalogue analysis used host `M200c`, whereas the new
NewGalFinder extraction supplies FoF total mass. To avoid comparing those
different host-mass cuts directly, the legacy catalogue was reanalysed using
its `host_mtot_msun` field, which is the closest available match to the FoF
total-mass proxy. For `10^13 <= Mhost/Msun < 10^14`, this selects 60 hosts.
At `Mstar >= 10^8 Msun`, the fixed-host, fixed-stellar-mass slope is
`-0.209 +/- 0.010`, with host-bootstrap 95 per cent interval
`[-0.230,-0.191]`; the 0.3-dex mass-bin result is `-0.168 +/- 0.011`, with
interval `[-0.190,-0.148]`. The mass-bin result at `Mstar >= 10^9 Msun` is
`-0.114 +/- 0.024`, with interval `[-0.162,-0.071]`. The matched orbital
velocity comparison remains consistent with zero (`+0.0017 dex`, interval
`[-0.0300,+0.0269]`). These results support using catalogue total host mass
for the cross-catalogue comparison, but the legacy and FoF mass definitions
still need explicit verification before a paper claim.

## Age and size diagnostics in the legacy sample

For `Mstar >= 10^8 Msun`, the continuous fixed-host model for global internal
dispersion gives a `log sigma`--`log Z` slope of `-0.209` when controlling
stellar mass, `-0.195` when also controlling mass-weighted stellar age, and
`-0.352` when controlling three-dimensional half-mass radius. Including both
age and radius gives `-0.352`. At `Mstar >= 10^9 Msun`, these slopes are
`-0.145`, `-0.109`, `-0.394`, and `-0.372`, respectively. These are
diagnostic multivariable fits from one snapshot and do not establish a
mechanism; the marked shift after adding radius shows that internal structure
is an important covariate and should be included in the final analysis where
the catalogue definitions permit it. The calculations are in
`preliminary_satellite_relations_age_size_v2_00972.csv`.

## Star-formation activity split in the legacy sample

The original ASCII catalogue provides the 100-Myr-averaged SFR. Using
`sSFR = SFR_100Myr/Mstar`, galaxies were classed as star-forming for
`sSFR >= 10^-10.5 yr^-1`, transitional for
`10^-11.5 < sSFR < 10^-10.5 yr^-1`, and quiescent for
`sSFR <= 10^-11.5 yr^-1`. The group slopes below use the global internal
dispersion, control continuously for stellar mass and with host fixed effects,
and obtain 95 per cent intervals by resampling whole host haloes.

- At `Mstar >= 10^8 Msun`, the slopes are `-0.143`
  (`[-0.196,-0.095]`, 615 galaxies) for star-forming satellites, `-0.221`
  (`[-0.257,-0.181]`, 466) for transitional satellites, and `-0.215`
  (`[-0.236,-0.196]`, 2,558) for quiescent satellites.
- At `Mstar >= 10^9 Msun`, the star-forming slope is `-0.032`
  (`[-0.085,+0.016]`, 396), consistent with zero; the transitional slope is
  `-0.154` (`[-0.221,-0.095]`, 258), and the quiescent slope is `-0.191`
  (`[-0.249,-0.143]`, 690).

This one-snapshot diagnostic suggests the residual relation depends on recent
star-formation activity, especially above `10^9 Msun`; it does not establish
that activity causes the trend. The thresholds follow a commonly used
three-way sSFR split, but need checking against this simulation's SFR
calibration and redshift before interpretation. Results:
`preliminary_satellite_relations_sfr_v2_00972.csv`.

### Direct test of the activity-dependent slope

I fitted a pooled host-fixed-effects model with separate stellar-mass slopes
for star-forming and quiescent satellites, then bootstrapped whole hosts 2,000
times. This directly tests the difference between the two metallicity slopes.
For `Mstar >= 10^9 Msun`, `10^13 <= Mhost < 10^14 Msun`, and the baseline
sSFR cuts (`-10.5` and `-11.5` in `log10(yr^-1)`), the quiescent-minus-
star-forming slope is `-0.189` with host-bootstrap 95 per cent interval
`[-0.258,-0.121]`. The difference remains negative when both boundaries are
shifted by 0.5 dex: `-0.219 [-0.290,-0.141]` for cuts `(-10,-11)` and
`-0.173 [-0.233,-0.111]` for `(-10.7,-11.7)`. Thus, in this single output,
the slope difference is not an artefact of the exact activity thresholds
within this range. This is still a conditional statistical result, not a
causal or multi-epoch physical conclusion; validate the SFR calibration and
repeat on the four NewGalFinder outputs before using it in the manuscript.
With the baseline activity cuts, the slope difference is similar for all
three internal-dispersion definitions. At `Mstar >= 10^8 Msun`, it is
`-0.086 [-0.135,-0.036]` for global dispersion, `-0.091 [-0.144,-0.036]`
for the mass-range estimate, and `-0.082 [-0.131,-0.032]` within the
half-mass radius. At `Mstar >= 10^9 Msun`, the corresponding differences are
`-0.189 [-0.258,-0.121]`, `-0.194 [-0.263,-0.124]`, and
`-0.186 [-0.246,-0.123]`. In the higher-mass sample, all three definitions
also retain a negative difference under both shifted activity cuts. In the
lower-mass sample, the half-mass-radius result for cuts `(-10,-11)` has an
interval that includes zero; the other two definitions do not. The complete
set of 18 direct interaction fits is in
`preliminary_activity_interaction_sensitivity_v3_00972.csv`. This strengthens
the one-snapshot statistical result but still needs replication in the new
catalogues and validation of the SFR and dispersion definitions. Earlier
interaction files are superseded; the initial one also excluded zero-SFR
objects and must not be used.

## Current subhalo dark-matter-mass diagnostic

The legacy ASCII table labels column 42 `Mdm (Msun)`. I extracted it and
repeated the fixed-host, fixed-stellar-mass fits while also controlling for
`log Mdm`, with 2,000 whole-host bootstrap resamples. For `Mstar >= 10^8 Msun`,
the baseline metallicity slopes of `-0.209`, `-0.222`, and `-0.182` (global,
mass-range, and half-mass-radius dispersions) fall to `-0.037 [-0.057,-0.017]`,
`-0.031 [-0.052,-0.012]`, and `-0.027 [-0.049,-0.007]`. At
`Mstar >= 10^9 Msun`, the controlled slopes are `-0.035 [-0.073,+0.001]`,
`-0.031 [-0.069,+0.007]`, and `-0.023 [-0.061,+0.012]`; all three intervals
include zero. This strong attenuation makes present-day subhalo mass a key
covariate to examine, but it does not establish a causal pathway. The source
IDL routine that defines the catalogue's `Mdm` column is unavailable, so its
aperture and exact physical meaning remain unverified. Treat this as a
diagnostic only. Data and bootstrap results are in
`satellite_properties_ascii_sfr_mdm_00972.csv`,
`preliminary_satellite_relations_sfr_mdm_00972.csv`, and
`preliminary_mdm_control_host_bootstrap_00972.csv`.

In effect-size units, the controlled global-dispersion slope at the
`10^8 Msun` cut corresponds to about a 2.5 per cent lower dispersion for a
0.3-dex increase in the metallicity proxy (host-bootstrap interval: about
1.2--3.9 per cent lower). At the `10^9 Msun` cut, the corresponding interval
is about 4.9 per cent lower to 0.05 per cent higher, consistent with no
residual change. These are transformations of the legacy-ASCII diagnostic,
not validated results; they illustrate why the new catalogue analysis must
judge physical importance by effect size as well as by the slope interval.

I also fitted an activity interaction model that allows both the stellar-mass
and `Mdm` slopes to vary between star-forming and quiescent satellites. At
`Mstar >= 10^9 Msun`, the quiescent-minus-star-forming metallicity-slope
difference remains negative after this control: `-0.146 [-0.228,-0.068]`
for global dispersion, `-0.142 [-0.232,-0.061]` for the mass-range estimate,
and `-0.141 [-0.217,-0.069]` within the half-mass radius. At the lower mass
cut, two of the three intervals include zero. This indicates that the
activity-dependent result is not removed by the catalogue `Mdm` value in the
higher-mass sample, although the uncertain definition of `Mdm`, the SFR
calibration, and the single snapshot prevent a physical interpretation.
Results are in `preliminary_activity_mdm_interaction_00972.csv`.

## Related published work checked

- Spindler & Wake (2017) compare SDSS central and satellite size, stellar mass,
  stellar metallicity, and core velocity dispersion. Their question is
  environment-dependent scaling relations; metallicity differences are
  small after separating star-forming and quiescent galaxies. This is close
  context but is not the same satellite-only residual relation tested here:
  https://doi.org/10.1093/mnras/stx427
- Buck et al. (2019, NIHAO XV) measure the stellar 1D dispersion within the
  3D half-light radius and stellar metallicity in simulated satellites and
  field dwarfs. They present the metallicity--stellar-mass and
  dispersion--stellar-mass relations separately, finding that the
  populations largely overlap and that satellites tend to have lower
  dispersions at fixed stellar mass. They do not test the residual
  metallicity--dispersion association at fixed host and stellar mass studied
  here: https://doi.org/10.1093/mnras/sty2913
- Wang, Wang & Chen (2023) study host-mass dependence of satellite
  gas-phase metallicity in EAGLE and IllustrisTNG, including suppression of
  gas accretion and stripping. Their metallicity is gas-phase, and their
  analysis is not an internal stellar-kinematics test:
  https://doi.org/10.3847/1538-4357/acd633
- Lacerna et al. (2017) examine stellar and gas metallicities of group and
  cluster satellites in EAGLE and the origin of satellite metallicity
  enhancement, mainly at `Mstar > 10^10 Msun`:
  https://doi.org/10.1093/mnras/stw2329
- Fattahi et al. (2018, APOSTLE) directly examine the dwarf-galaxy
  `[Fe/H]--V_{1/2}` relation. They find that tidally stripped satellites can
  lower and flatten the low-velocity end, because stripped satellites have
  lower velocities than field dwarfs at comparable metallicity. Their
  `V_{1/2}` is defined as `sqrt(3)` times the luminosity-weighted line-of-sight
  stellar dispersion, so this is directly relevant to stellar kinematics,
  although it is not the same estimator as our mass-weighted particle
  dispersion. Their
  progenitor reconstruction assumes that tides leave mean stellar metallicity
  unchanged:
  https://doi.org/10.1093/mnras/sty408
- Frings et al. (2017) follow hydrodynamical satellite analogues after infall.
  They report that satellites can lose dark matter and lower their internal
  stellar velocity dispersion while retaining nearly the same mean stellar
  metallicity; their paper discusses the mass--metallicity and
  size--dispersion relations separately rather than testing a conditional
  metallicity--dispersion relation:
  https://doi.org/10.1093/mnras/stx2171
- Errani et al. (2022) use N-body models to show that tidally limited stellar
  remnants do not follow one universal structural track; their sizes and
  dispersions depend on the progenitor stellar distribution, while strongly
  stripped remnants trace the surviving dark halo:
  https://doi.org/10.1093/mnras/stac476

These papers motivate a careful comparison but are not yet a complete
novelty review. The present analysis must not claim novelty until the
cross-snapshot catalogues, definitions, and targeted literature search are
complete.

## Host orbital velocities

Within each host, satellites were divided at the host median metallicity,
then paired one-to-one across the two groups by nearest stellar mass with a
0.1-dex caliper. Hosts needed at least four pairs. At the `10^8 Msun` cut,
44 hosts remain; the median `log10(sigma_high-Z/sigma_low-Z)` is `+0.0035 dex`
with a host-bootstrap 95 per cent interval `[-0.0261,+0.0370]`. The median
stellar-mass offset of the matched groups is `0.048 dex`. This does not show
a detectable orbital-dispersion difference. The `M200c >= 10^14 Msun` sample
has only two hosts and cannot support a cluster-only conclusion.

## Required checks before interpretation or manuscript edits

1. Reproduce the trends using the newly generated NewGalFinder catalogues for
   snapshots 972, 977, 980, and 984, and compare the overlapping 972 galaxies
   against this legacy catalogue.
2. The ASCII header identifies the measured columns as `Z(star)`,
   `SFR_100Myr`, global `Vsig`, `Vsig(M1/3-2/3)`, and
   `Vsig(r<R1/2_3D)`. Confirm the exact weighting, apertures, and formulae,
   along with host flags and catalogue selection, from the source code. The
   source IDL routines are symlink targets outside this workspace and are not
   readable from the login node. The new binary parser checks repeated stellar
   IDs within each galaxy and across galaxies in each FoF host; affected
   galaxies are excluded from the new-catalogue relation analysis. The legacy
   ASCII catalogue contains no stellar particle IDs, so that duplicate check
   is unavailable for its preliminary result.
3. Test alternate stellar-mass thresholds, host-mass bins, metallicity
   definitions, and robust uncertainty estimates based on host resampling.
4. Only then decide whether to include a concise satellite result in the
   existing surface-brightness calibration paper or keep it as a separate
   study. Do not combine satellite orbital dispersion with internal stellar
   dispersion under the single term “velocity dispersion”.

## Current processing queue

The four working directories `output_00972`, `output_00977`, `output_00980`,
and `output_00984` each report `ncpu=6000` and contain 6,000 AMR, gravity,
particle, and hydro files. The numbered first and last files (`out00001` and
`out06000`) were checked for all four file types in all four snapshots. This
allowed NewDD to be submitted without waiting for the separate full-tree copy:

- `430263` was cancelled after its 972 output showed a zero-byte STAR file:
  `NEWDD_DM_ONLY` had assigned every particle to the dark matter family. That
  output is preserved as `NewDD.00972.partial_dm_only_430263`.
- `430282` and `430294` both failed while parsing optional particle records.
  A direct record-length inspection of four snapshot-972 CPU files showed that
  the input contains three chemistry arrays after birth time and metallicity,
  followed by the initial-mass array. The NewDD reader had incorrectly tried
  to read nine chemistry arrays. It now reads the three input arrays while
  retaining the nine-element in-memory struct layout required by the existing
  NewDD/NewGalFinder binaries. The reader also skips AMR and hydro mesh data
  and returns immediately for empty CPU domains. The executable has been
  explicitly relinked with these parser changes. No usable output was produced
  by either failed attempt.
- `430305`: current NewDD array for snapshots 972, 977, 980, and 984, using 25
  CPU ranks across two A10 nodes, 300 GB per node, and no GPU request. Its
  batch script requires a non-empty STAR output.
- NewDD task `430305_0` completed snapshot 972 successfully in 12:40 with
  exit code 0. All ranks reported nonzero star counts; the first DM and STAR
  slab files are 2.3 GB and 4.1 GB. Task `430305_1` is processing snapshot
  977. Parser errors were absent from the completed 972 log.
  A direct read of five 168-byte STAR records confirmed family=2, finite
  positions and masses, positive initial masses, and nonzero stellar
  metallicities in the expected field. This checks the particle-record
  mapping, not the catalogue-level abundance calibration.
  All 100 DM and 100 STAR slabs are present and non-empty, with every file
  size divisible by the 168-byte particle-record size. Summing the 6,000 CPU
  records in the NewDD log gives 1,314,008,696 stars and 1,073,741,824 DM
  particles; these totals exactly match the slab counts. The combined 972
  NewDD particle files occupy about 401 GB.
- NewDD task `430305_1` completed snapshot 977 successfully in 12:26. All
  6,000 CPU records sum to 1,318,136,826 stars and 1,073,741,824 DM particles;
  all 100 DM and 100 STAR slabs are non-empty and have sizes divisible by 168.
  Task `430305_2` is now processing snapshot 980.
- The first FoF attempt (`430347_0`) was cancelled after its old executable
  reported file-size errors and impossible particle coordinates outside the
  44.17 cMpc/h box. The partial directory is preserved as
  `FoF.00972.partial_stale_opfof_430347`. The executable predated the current
  `ramses2read.c`; opFoF has now been rebuilt against that reader, and the
  previous binary is preserved as
  `opFoF/opfof.exe.pre_ramses2read_sync_20261020`. The first retry
  (`430395`) was cancelled because Slurm had already purged NewDD task 0, so
  its `aftercorr` dependency could not resolve. The failed FoF output remains
  preserved and is not used.
- The current chain runs snapshot 972 as array index 0 and snapshots 977--984
  as indices 1--3, preserving per-snapshot progression: `430417` (FoF index 0,
  no dependency because verified NewDD task 0 has aged out of Slurm), `430418`
  (FoF indices 1--3 aftercorr `430305`); `430422` (NewGalFinder index 0
  aftercorr `430417`), `430423` (NewGalFinder indices 1--3 aftercorr `430418`);
  `430426` and `430427` (catalogue extraction after the corresponding
  NewGalFinder jobs); and `430430` and `430431` (relation analysis after the
  corresponding extraction jobs). Each stage refuses to overwrite an
  existing product. The rebuilt FoF reader has not yet been validated on a
  complete snapshot.
  The box-size override `NEWDD_BOXLEN_CMPCH=44.171267363` remains in place and
  agrees across the four outputs and the valid NewDD header for snapshot 967.
  Earlier partial outputs `NewDD.00972.partial_from_430252_0`,
  `NewDD.00977.partial_from_430237_1`, and
  `NewDD.00972.partial_dm_only_430263` remain preserved and unusable.
- `429983`: separate full-tree transfer, still running independently.

The binary parser's particle layout and velocity units were cross-checked
against `NewGalFinder/ramses.h` (`x,y,z` in cMpc/h, velocities in km/s, stellar
masses in Msun/h). RAMSES documentation identifies the particle field `zp` as
the stellar particle's birth metallicity. The extractor therefore measures
the current-mass-weighted mean of stellar birth-metallicity values; it should
not be described as an element-specific abundance such as `[Fe/H]` without
additional conversion. The local header documents the stored `zp` only as
being in RAMSES units, so the exact physical conversion and the legacy
`Z(star)` convention remain to be confirmed before comparing absolute values.
A constant unit conversion would not change a log--log slope, but a definition
or weighting difference could. RAMSES particle-field documentation:
https://ramses-organisation.readthedocs.io/en/latest/dev_docs/lecture/particles.html

The satellite selection follows NewGalFinder's central definition. Its
`src/binding.c::central_stellar_core` selects the core with the largest
stellar mass; the extractor applies the same maximum-stellar-mass rule within
each FoF host, then excludes that core from the satellite sample. This was
verified against the source rather than inferred from subhalo ordering.

The queued extraction now also exports each subhalo's assigned dark-matter,
gas, sink, stellar, and total masses. The queued analysis will fit the
metallicity--internal-dispersion relation with subhalo total mass and dark
matter mass as additional controls, alongside host and stellar mass. This
will test whether the residual trend tracks the present-day subhalo mass
budget. It also stores both mass-weighted and number-weighted stellar
metallicity, and will compare their slopes. The particle-count threshold is
applied to stars with valid velocities, which are the particles entering the
dispersion estimate. The binary extraction does not yet provide stellar ages or radii;
those controls are currently available only in the legacy ASCII catalogue.

## Pipeline update (2026-10-20 20:53 KST)

- NewDD tasks `430305_0`, `_1`, and `_2` have completed snapshots 972, 977,
  and 980. For snapshot 980, the 6,000 CPU log records sum to 1,322,162,847
  stars and 1,073,741,824 DM particles. All 100 DM and 100 STAR slabs are
  non-empty, and every file size is divisible by the 168-byte particle-record
  size; their record totals match the log exactly. Task `_3` (snapshot 984)
  remains pending with Slurm reason `JobArrayTaskLimit`.
- Rebuilt opFoF task `430417_0` has started snapshot 972. Its initial log shows
  particle coordinates within the expected 44.171 cMpc/h box and plausible
  per-slab ranges; the previous impossible-coordinate failure has not
  reappeared. It reports missing GAS and SINK files as empty, consistent with
  the NewDD particle-only output. The `RamsesType` warning concerns failure to
  read the binary parameter block at matching struct size; the reader still
  parses the ASCII header. Treat FoF as unvalidated until it completes and
  its output and membership counts are checked.
- No NewGalFinder satellite catalogue exists yet. The relation analysis must
  still be run on the new catalogues before the legacy single-snapshot result
  can support any physical claim.

## Pipeline update (later on 2026-10-20)

- Snapshot 980 is complete: all 6,000 NewDD log records sum to 1,322,162,847
  stars and 1,073,741,824 DM particles. All 100 DM and 100 STAR slabs are
  non-empty; each file size is divisible by 168 bytes, and the summed slab
  record counts match the log.
- FoF task `430417_0` is processing snapshot 972. Its input-coordinate ranges
  so far are inside the 44.171 cMpc/h box, unlike the failed stale-executable
  run. Its non-empty GAS/SINK files are absent and are handled as empty; treat
  the FoF catalogue as unverified until the job finishes and its catalogue and
  member counts are audited.
- NewDD task `430305_3` was held with `JobArrayTaskLimit` despite tasks 0--2
  being complete, so that pending task and its three dependent snapshot-984
  jobs were cancelled. Snapshot 984 has been re-submitted as a single-item
  chain: NewDD `430473_3` -> FoF `430474_3` -> NewGalFinder `430476_3` ->
  extraction `430477_3` -> relation analysis `430479_3`. No snapshot-984
  output existed when this replacement chain was submitted.
- A second safety audit found that `FoF_Garbage/Garb.00972` was a symlink to
  the original run's auxiliary directory. opFoF opens its first
  `BottomFaceContactHalo` output in write/truncate mode. FoF task `430417_0`
  was cancelled before reaching that write: all ten original contact files
  still have their pre-run sizes and Oct 17--18 modification times. Its
  partial `FoF_Data/FoF.00972` directory was preserved under a
  `.partial_cancelled_before_auxpath_fix_430417` name, and the symlink itself
  was preserved as `Garb.00972.preserved_original_symlink_20261020`; a new
  empty local `Garb.00972` directory is in place. The FoF batch script now
  refuses symlinked or non-empty auxiliary-output directories. Snapshot 972
  was re-submitted as FoF `430490_0` -> NewGalFinder `430491_0` -> extraction
  `430492_0` -> relation analysis `430493_0`. Snapshot 977 task `430418_1`
  is running with its own real, initially empty auxiliary-output directory.

## Pipeline update (2026-10-20 21:13 KST)

- FoF task `430418_1` for snapshot 977 is still running on three nodes. Its
  log records that ranks P1, P6, P2, P3, and P4 passed halo-boundary checks;
  the final catalogue and member files have not yet been written. Slurm
  reports the job active, with about 35 GB maximum resident memory. The FoF
  request is `cpu=60,mem=768G,node=3` and requests no GPU.
- Snapshot 972 FoF `430490_0` is pending node resources; snapshot 984 NewDD
  `430473_3` is pending priority. The remaining 977/980 array and downstream
  catalogue jobs are waiting on their configured array limits or dependencies.
  No satellite catalogue or new metallicity--dispersion result is available
  yet.
- A source review confirmed that the extractor preserves the stored stellar
  `zp` values and computes mass-weighted and number-weighted means, while the
  analysis uses host fixed effects and stellar-mass controls. The RAMSES
  header still labels particle `zp` only as “RAMSES units”; its physical
  abundance convention must be established before assigning absolute
  metallicities or interpreting cross-epoch differences.

At 21:14 KST, task `430418_1` was still active and its log had advanced: rank
P8 also passed the halo-boundary check. The two final FoF output files remained
at their initial sizes (28-byte halo header and empty member file), so the
catalogue is not ready for downstream processing yet.

At 21:17 KST, rank P5 had also passed the check. Slurm still reported active
CPU use; no FoF catalogue records had been written yet.

At 21:18 KST, rank P9 passed the same check, bringing the count to eight of
ten MPI ranks. The halo and member output files were still at their initial
sizes; the job remained active.

## Metallicity-field definition and extraction schema (2026-10-20 21:22 KST)

- RAMSES particle documentation describes `zp` as a star particle's birth
  metallicity ([official RAMSES particle documentation](https://ramses-organisation.readthedocs.io/en/latest/dev_docs/lecture/particles.html)).
  The extracted mass-weighted statistic uses each particle's
  current mass, so it is a current-mass-weighted mean of birth metallicities,
  not an element-specific abundance such as `[Fe/H]`. The local code still
  labels the stored value only as being in RAMSES units; its absolute
  conversion remains unresolved.
- Before any new satellite-property catalogue existed, the extractor and
  analysis scripts were updated to use the explicit field name
  `mean_birth_metallicity_current_mass_weighted`. A source-reference check confirmed that
  both scripts now use the same name. This avoids leaving the mass definition
  implicit in the generated data products.
- A small read of the first 10,000 stellar records in the snapshot-972 STAR
  slab found finite `zp` values spanning `1.34e-5` to `1.46e-2`. The sum of
  the stored `chem[1:]` components tracks `zp` very tightly in this sample
  (log-space correlation `0.9999`, median ratio `0.746`), but the component
  mapping is not established. This is a consistency check only; it does not
  determine the absolute `zp` conversion or justify treating it as `[Fe/H]`.
- At this update, FoF task `430418_1` was still running and its dependent
  NewGalFinder task remained pending. No satellite catalogue or new
  metallicity--dispersion result was available.

At 21:32 KST, the relation-analysis output schema was also updated to record
the metallicity definition used in each result row. Legacy `Z(star)` results
will be labelled as having unverified weighting/convention; new-catalogue
results will be labelled as current-mass-weighted RAMSES stellar birth
metallicity in RAMSES units. A source-reference check confirmed that both
input branches populate this field before any analysis output is generated.

At 21:24 KST, the live log showed why the task is taking longer: ranks P0 and
P7 hold the two largest slab payloads (41.2 and 54.5 million particles).
Both had completed tree construction but had not yet printed their FoF halo
counts, while the other eight ranks had completed that pass. This points to
the particle-linking pass on the two largest slabs as the current workload;
Slurm continued to report the job running with CPU time increasing.

## Cluster resource check (2026-10-20 21:29 KST)

`sinfo -s` reports zero idle nodes in every listed compute partition and no
CPU-only partition. FoF `430418_1` remains active; the earlier Slurm TRES audit
showed only `cpu=60,mem=768G,node=3` and no GPU request. Thus the queued
snapshot-972 FoF job is waiting for node availability, not a GPU quota, while
the snapshot-984 NewDD job is waiting on priority.

## Satellite-size control and subhalo-mass provenance (2026-10-20 21:46 KST)

- The 00972 preliminary catalogue's `Mdm` definition is still unavailable, so
  its strong attenuation after controlling for `Mdm` cannot yet be interpreted
  as a verified physical control. In the new NewGalFinder source,
  `catalogue.c` sums the masses of all dark-matter particles assigned to each
  subhalo. The new extraction therefore provides a clearly defined subhalo
  dark-matter mass for an independent check.
- The new extractor previously emitted no stellar size, despite the 00972
  legacy analysis showing that size is an important covariate. It now records
  the current-stellar-mass-weighted 3D half-mass radius in physical kpc,
  measured around the stellar mass centroid with periodic box wrapping. The
  calculation unwraps coordinates around a stellar-particle anchor, avoiding
  dependence on the catalogue's total-mass centre near a periodic boundary. It
  also records mass-weighted 1D stellar velocity dispersion inside that same
  3D half-mass radius, alongside the all-member dispersion. The relation
  analysis now maps the size into its radius-control fit and tests the
  aperture-based dispersion separately. The size convention is explicit and
  differs from the legacy catalogue's radius until that legacy
  aperture/centering definition is recovered.
- The extractor records how many stellar velocities enter the half-radius
  estimate. For the new catalogue, the relation analysis requires at least 30
  such particles and applies that same sample cut to all dispersion
  definitions, so aperture and all-member slopes use the same satellites.
  The threshold can be varied on reruns; the legacy catalogue has no aperture
  particle count and is unchanged by this cut.
- No test suite was run. Static source inspection confirmed the extractor
  writer and analysis reader use matching size, aperture-dispersion, and
  aperture-count field names. The extractor is
  dependency-pending, so no generated catalogue was overwritten or produced
  from this change yet.

At 21:39 KST, FoF task `430418_1` remained active at 38:58 elapsed time. Its
MPI step reported 2:09:35 aggregate CPU time and about 35.5 GB maximum RSS;
ranks P0 and P7 had still not completed their linking pass, and both 00977 FoF
output files remained at their initial sizes. Snapshot 972 FoF was still
pending node resources; no GPU is requested by either FoF job.

At 21:44 KST, task `430418_1` was still running, and output writing had begun.
The log shows rank P0 completed its halo-boundary check and began preparing
5.17 million isolated halos; ranks P1, P2, P3, P4, P5, P6, P8, and P9 had sent
their isolated-halo data to P0. Rank P7 had not yet reported its halo count, so the
snapshot catalogue is not complete. At 21:42, the output files had grown to
180,252 bytes and 11,616,256 bytes, respectively. Do not use either file until
Slurm reports successful task completion and the batch script's output checks
pass.

At 21:46 KST, Slurm still reported `430418_1` running. The MPI step's
aggregate CPU time had risen to 2:32:47 while maximum RSS remained about
35.5 GB. The 00977 output files had not grown since 21:42, and rank P7 still
had no halo-count line in the log. This is ongoing CPU work rather than a
completed catalogue; downstream jobs remain dependency-pending.

At 21:49 KST, the same task was still active and aggregate CPU time continued
to rise (2:43:03 at 21:49), while output sizes and the last P7 log milestone
were unchanged. Stderr contains only the `RamsesType` size warning and notices
that empty GAS/SINK component files are absent. The NewDD text header records
`BoxSize=44.171267`, `H0=67.400002`, and `a=0.821128`; the FoF startup log
reports the same box size. Source inspection shows `NewDD::write_head` writes
the ASCII parameter header and `opFoF::read_head` parses it. The warning path
skips the optional binary-struct read; its diagnostic about precision loss is
not itself evidence that the parsed ASCII box size or cosmology was wrong.
The FoF output still requires normal completion and content validation.

## Targeted literature check (2026-10-20)

The broad claim that stellar metallicity is related to satellite velocity is
not new. Fattahi et al. directly discuss the Local Group dwarf `[Fe/H]`–
`V_{1/2}` relation and argue that tidal stripping can lower a satellite's
velocity at roughly unchanged mean metallicity, flattening the low-velocity
end ([MNRAS 476, 3816](https://academic.oup.com/mnras/article/476/3/3816/4875945)).
Frings et al. follow satellite analogues after infall and find that stellar
mass–metallicity can remain nearly unchanged while stellar velocity
dispersion and size evolve under tides ([MNRAS 472, 3378](https://academic.oup.com/mnras/article/472/3/3378/4093087)).
Spindler & Wake compare stellar mass, metallicity, size, and central velocity
dispersion for observed centrals and satellites ([MNRAS 468, 333](https://academic.oup.com/mnras/article/468/1/333/3038254)).

Therefore, any contribution here must be narrower and tested directly: a
conditional relation in resolved satellites of massive hosts, after matching
or controlling stellar mass, host, present subhalo mass, and stellar size,
with its evolution across the available outputs. The current one-output
legacy result is not enough to establish that contribution. The remaining
novelty review must also determine whether prior work has already tested this
same conditional relation in cluster satellites at multiple redshifts.

## Snapshot 00977 FoF progress check (2026-10-20 21:56 KST)

Slurm still reports array task `430418_1` running. `sstat`/`sacct` did not
return current per-process CPU usage, so a single short, read-only `srun` on
the allocated nodes inspected the running ranks: all ten `opfof.exe 977 100`
MPI processes were present and each was using about 99--100% of one CPU. The
unchanged output files therefore do not indicate a stalled job; the ranks are
still doing CPU work. Keep the task running and wait for rank P7 and the batch
script's final output checks before using snapshot 00977 or releasing its
downstream steps.

## FoF resource-use audit (2026-10-20 21:58 KST)

The active snapshot-00977 allocation reserves 60 CPUs (10 MPI ranks x 6 CPUs
per rank) across three nodes. A read-only process inspection showed ten
`opfof.exe` ranks each using approximately one CPU, and source inspection found
no OpenMP pragmas in the `opFoF` C sources. This executable is therefore
operating as a single-threaded MPI rank per process; the current request appears
to reserve substantially more CPU capacity than the ranks use. Do not alter or
restart the active task. For any future re-submission, consider requesting one
CPU per MPI rank and reassessing whether 10 ranks can fit on two nodes; retain
adequate per-node memory because the earlier step reported a maximum RSS near
35.5 GB. Existing pending tasks were submitted with the old resource request,
so changing the batch script alone will not change their allocations.

At 21:59 KST, a fresh Slurm poll confirmed task `430418_1` was still running
(59:50 elapsed). A lightweight process inspection on its allocated nodes showed
all ten MPI ranks in running state at about 99--100% CPU; their CPU times had
advanced with elapsed time. The log and partial output sizes were unchanged
since 21:42, consistent with rank P0 continuing its large local-halo pass while
waiting for rank P7. This is a verified live computation, not a safe reason to
restart. Snapshot 00972 remains pending resources and all dependent catalogue
steps remain queued.

## FoF queue request adjustment (2026-10-20 22:05 KST)

The pending 00972 FoF job requested 3 nodes, 60 CPUs, and 256 GiB per node.
With `syn02` excluded, `syn01` had only 20 idle CPUs while the job required at
least 24 CPUs on every node; the other a10 nodes had 29--38 idle CPUs. The
application uses ten single-threaded MPI ranks. I therefore submitted
replacement FoF tasks with 2 nodes, 10 total CPUs (one per rank), 5 ranks per
node, and 224 GiB per node; the 00972 inputs were present and the output
location was unused. The 00972 galaxy-finder task now depends on the new FoF
job. The old pending FoF task was cancelled only after that dependency update.

I applied the same change to snapshot 00980 after confirming all 100 DM shard
files were present and its output location unused. Its galaxy-finder task was
rewired before cancelling the old array task. Extraction and analysis remain
attached to their existing galaxy-finder array tasks. Both replacement FoF
jobs request no GPU. Their time limits were reduced from 48 to 24 hours to
improve backfill eligibility, but Slurm still projects their earliest start as
2026-10-22 21:11 KST, with priority as the current pending reason. The active
00977 task was not modified; its MPI ranks remain CPU-active. Thus the new
requests reduce resource reservation, but the cluster queue still controls
when they start.

Snapshot 00984's NewDD conversion remains a separate pending resource request
(`430473_3`: 2 nodes, 25 MPI ranks, 300 GiB per node; no GPU requested). I left
that request unchanged because no measured peak-RSS evidence for this conversion
was available to justify lowering its memory limit.

## Snapshot 00977 FoF milestone (2026-10-20 22:05 KST)

A later log poll showed new progress after the previous check: all ten MPI ranks
reported `passed Tree building`, including P7. The two partial FoF files grew
from 180,252/11,616,256 bytes to 204,508/12,818,568 bytes at 22:03:45, and the
log continued updating through 22:04:50. Slurm still reports the task running,
so these outputs remain incomplete and downstream jobs must continue waiting
for the batch script's successful completion checks.

## Snapshot 00984 NewDD memory adjustment and 00977 live check (2026-10-20 22:10 KST)

Completed NewDD array tasks for snapshots 00972, 00977, and 00980 provide
measured peak RSS values of 5,989,388 K, 5,184,864 K, and 7,935,492 K,
respectively, with successful 12--17 minute runtimes. On that evidence, I
replaced the pending 00984 NewDD request (300 GiB/node) with job `430621_3`,
requesting 128 GiB/node for the same 25 ranks on two nodes and a 6-hour limit.
The 128-GiB request allows headroom above 13 ranks at the prior maximum RSS.
The 00984 FoF task was also replaced by `430622_3` (10 ranks, one CPU/rank,
two nodes, 224 GiB/node, 24 hours). Its dependency points to the new NewDD
job; the galaxy-finder dependency was rewired before cancelling the old FoF
job. The extraction and analysis dependencies remain attached to that
same galaxy-finder task.

A live inspection at 22:10 found all ten 00977 FoF MPI ranks still running at
about 99% CPU, with resident memory around 5--16 GiB per rank. The job remains
active; no further output or log update has appeared since 22:04:50, so the
partial catalogue remains unusable. Scheduler refreshes currently show
`430610_0` pending for resources, `430611_2` pending priority, and `430621_3`
pending priority. Their reduced requests have not yet produced a start; the
cluster queue remains the timing constraint.

At 22:10 KST, snapshot 00977 FoF remained active (1:09:30 elapsed). The log
advanced beyond tree construction: rank P1 reported 5,567,027 halos from
19,909,588 particles and completed its halo-boundary check. The main catalogue
files had not yet changed beyond the partial sizes recorded earlier, so the
batch-level completion checks are still required before downstream jobs run.
The replacement 00972, 00980, and 00984 tasks remain queued; 00984's NewDD
request is now pending priority rather than the previous 300-GiB resource
request.

A refreshed queue snapshot at 22:10 KST still showed 00977 running at 1:10:00
elapsed. Its latest log entry remained rank P1's completed boundary check; the
FoF files are still partial. The new 00972 FoF task is pending for resources,
00980 FoF and 00984 NewDD are pending priority, and 00984 NewDD now has a
projected start of 2026-10-21 17:59 KST. All downstream jobs remain correctly
held on their respective success dependencies.

## Snapshot 00984 NewDD active run (2026-10-20 22:17 KST)

The 00984 NewDD conversion was moved from the pending two-node a10 request to
job `430626_3` on one h200 node (25 MPI ranks, 128 GiB, 6 hours, no GPU). This
uses the measured RSS of earlier NewDD snapshots and fits the h200 node's
available CPU and memory. The old `430621_3` task was cancelled after the
FoF dependency was rewired to `430626_3`; downstream 00984 FoF remains job
`430622_3`. At 2:18 elapsed, Slurm showed `430626_3` running on `syn104`, with
all particle shards being processed and 00984 DM/star outputs already being
written. A process sample showed 25 CPU-active ranks at roughly 0.2--0.7 GiB
resident memory each (their virtual address spaces include 20 GiB buffers).

The 00972 and 00980 replacement FoF jobs now target a100 rather than a10;
neither requests GPUs. Both remain pending priority. Snapshot 00977 remains
active; its latest new log milestones are boundary checks completed by P1, P2,
P3, and P6, while its catalogue files are still partial.

## 00977 / 00984 processing progress (2026-10-20 22:18 KST)

At 22:17 KST, `430626_3` (NewDD for 00984) had run for three minutes on
`syn104`; its log was progressing through particle sorting. The first output
shards had grown to 441 MB (DM) and 489 MB (STAR), so the one-node CPU run is
actively producing data. The dependent FoF job `430622_3` was moved to the
a100 partition, matching the available memory/CPU path used for the other
reduced-resource FoF tasks; no GPU was requested.

Snapshot 00977 FoF remains active and incomplete. Since the previous check,
rank P4 joined P1, P2, P3, and P6 in reporting its halo count and completing
the boundary check. The catalogue files have not advanced past their partial
sizes, so they still await later ranks and final batch validation.

## Latest active-job verification (2026-10-20 22:19 KST)

Snapshot 00984 NewDD job `430626_3` remained active at 4:09 elapsed, with
particle-sorting stages still appearing in its log. The first 00984 DM and
STAR shard files grew to 621 MB and 602 MB, respectively. The 00977 FoF job
`430418_1` remained active at 1:18 elapsed; a fresh node-side process sample
showed all ten MPI ranks running at about 99% CPU, with roughly 4.7--17.3 GB
RSS each. Its last new log output was rank P4's boundary check at 22:16:33,
and the two catalogue files had not grown since 22:03:45. These FoF outputs
remain partial. Replacement FoF jobs for 00972 and 00980 are pending priority
on a100; the 00984 FoF task waits for NewDD success.

At 22:19 KST, NewDD job `430626_3` was still running (5:01 elapsed) and its
sorting log continued to advance. The first 00984 DM and STAR shards had grown
to 824 MB and 962 MB, respectively. FoF task `430418_1` remained live at
1:19:07, with the latest recorded halo-boundary completion from P4 at 22:16:33;
its catalogue files remained incomplete. No downstream extraction was released.

## opFoF NewDD-header correction and current queue (2026-10-20 22:23 KST)

Source inspection confirmed NewDD writes `.info` as an ASCII header only. The
opFoF reader had nevertheless compared the uninitialized `ramses_sizeof` field
against `sizeof(simpar)` and conditionally attempted a binary read. I removed
that invalid branch in both supported NewDD invocation modes and zero-initialized
the parsed `RamsesType` before reading the text header. The previous executable
was preserved as `opFoF/opfof.exe.pre_ascii_header_fix_20261020`; `make -j4`
rebuilt `opfof.exe` successfully at 22:22 KST. The compiler reported 17 format
warnings in existing `%d`/`size_t` calls; these are unrelated to the header fix.
No tests were run. The active 00977 MPI processes continue the already-loaded
old executable; queued FoF tasks will use the rebuilt executable when launched.

At the same check, NewDD 00984 job `430626_3` remained active at 8:41 elapsed;
its log had entered stage 9 and the first DM/STAR shards had reached 1.73/3.44
GB. The 00977 FoF task also remained active, and its output files were still
partial. The reduced-resource 00972/00980 FoF jobs remained pending priority.

## Active pipeline check (2026-10-20 22:23 KST)

The current analysis directory is `Darwin1/Cocoa_analysis_20261020`. Slurm
reported 00977 FoF task `430418_1` and 00984 NewDD task `430626_3` still
running; 00972/00980 FoF replacements `430610_0` and `430611_2` were pending
priority, while 00984 FoF `430622_3` remained dependency-pending. The 00977
FoF catalogue files were unchanged at 200 KB and 13 MB and remain partial.
The 00984 NewDD log had advanced through stage 8 and reported stellar sorting
on several MPI ranks. These are progress indicators only: no output is yet
accepted as complete, and extraction/analysis have not begun. The batch scripts
require successful program exit plus non-empty upstream files before advancing;
extraction and analysis also refuse to overwrite existing outputs.

A follow-up check at 22:25 KST confirmed that the 00977 FoF output files are
still partial, but the job is not idle: four opFoF ranks observed on `syn01`
were each using about 99 per cent CPU after 1:25 elapsed. Its log has not
advanced beyond 22:16:33, so the unchanged catalogue sizes are consistent with
ongoing rank-local computation rather than a confirmed exit or failure. The
00984 NewDD log and shard files were actively growing at 22:25 KST (about
1.94 GB DM and 4.13 GB STAR for the first shard); the job remains unvalidated
until all shards and the batch completion checks pass.

## Live-job check (2026-10-20 22:26 KST)

The same scheduler handles remain live: 00977 FoF `430418_1` and 00984
NewDD `430626_3` are running; 00972/00980 FoF replacements remain pending
priority; 00984 FoF remains dependency-pending. All ten 00977 opFoF MPI ranks
were present across `syn01`, `syn03`, and `syn04`, each using approximately
99 per cent CPU at 1:25 elapsed. Thus the long interval without log or output
file growth is ongoing compute, not a verified hang. The 00984 NewDD log and
its first DM shard continued growing through 22:25:52; the job has not reached
its success marker yet.

## Progress update (2026-10-20 22:26 KST)

Both active jobs are still running. The 00977 FoF log advanced at 22:26:16:
rank P5 completed its halo-boundary check, joining P1, P2, P3, P4, and P6;
the halo and member catalogues remain partial until all ranks finish. The
00984 NewDD log advanced through additional particle sorting and opening
`part_00984` records. At 22:26:29 its first DM and STAR shards were about
2.01 GB and 4.14 GB and still growing. No extraction or relation analysis has
started because neither upstream task has passed its batch completion checks.

## Static measurement-definition audit (2026-10-20)

While the large catalogue jobs run, I checked the extraction definitions against
NewGalFinder's source. Its `PmType` stores positions and velocities as
`cMpc/h` and `km/s`, particle mass as `Msun/h`, and stellar `zp` in RAMSES
units. The extractor's 168-byte stellar record layout matches the declared
fields through `zp`, nine chemical abundances, and initial mass. NewGalFinder
accumulates subhalo stellar mass from each assigned star's current mass. The
new analysis therefore treats `zp` as a relative birth-metallicity proxy; its
absolute solar normalization and equivalence to `[Fe/H]` are not established,
so no `[Fe/H]` claim is justified. It measures both all-member and
three-dimensional-half-mass-radius, mass-weighted 1D stellar dispersions and
will require particle-count and duplicate-ID cuts. These definitions must be
stated explicitly when interpreting final slopes.

## Progress update (2026-10-20 22:28 KST)

The 00977 FoF job `430418_1` and 00984 NewDD job `430626_3` remain live;
00984 FoF `430622_3` is still dependency-pending. The 00977 log advanced
through rank P5's boundary check and rank P0's member sorting, but the halo and
member catalogues have not been finalized. NewDD 00984 has created the expected
100 numbered `.info`, DM, and STAR shard names (300 files total, excluding
`.` and `..`) and continues reading and sorting input particle records. Shard
existence is not completion evidence: sizes and successful exit still need
verification for every shard.

## Manuscript-scope note (2026-10-20)

The PGalF handover at the project root defines the existing MNRAS manuscript as
a calibration of stellar FoF linking length against surface-brightness
contours and explicitly excludes unrelated physical-result papers. The
satellite metallicity--dispersion study is not part of that handover's stated
scope. I will establish and document the satellite result first, then keep its
manuscript treatment separate from PGalF unless a later explicit scope decision
says otherwise; the existing Overleaf manuscript is not being altered on the
basis of preliminary satellite statistics.

## Progress update (2026-10-20 22:30 KST)

Slurm still reports 00977 FoF `430418_1` and 00984 NewDD `430626_3` running;
00984 FoF `430622_3` remains after-OK dependent, and 00972/00980 FoF tasks are
pending priority. NewDD continues stage-2 input reads and particle sorting; its
first DM and STAR shards grew to about 2.20 GB and 4.20 GB by 22:30:11. The
00977 FoF log and two catalogue files have not changed since 22:26, so these
remain incomplete; the last log marker shows rank P5 finished its boundary
check. The task has not been declared stalled or restarted because its handle
is live and the MPI ranks were actively computing on all allocated nodes at
the preceding process check.

## Progress and reader-path check (2026-10-20 22:30 KST)

The same Slurm jobs remain live. NewDD 00984 advanced through further
particle sorting; its first DM and STAR shards reached about 2.25 GB and
4.23 GB by 22:30:42. The latest log shows AMR path names, so I checked the
active `NEWDD_SKIP_AMR` branch in source: NewDD prints those paths but does not
call `rd_amr`, and it jumps past hydro reading while retaining particle
records. The FoF 00977 catalogues remain at their partial sizes and its last
log boundary-check marker is P5; all rank-local progress will need final
completion checks before these files can be consumed.

## Sink-particle omission audit (2026-10-20)

The active `NEWDD_SKIP_AMR` path omits the separate RAMSES sink file, and
opFoF reports missing SINK shards for 00977. Before deciding whether this
invalidates the satellite analysis, I read the sink count and mass records
from the four existing `sink_NNNNN.out` files and applied NewDD's cosmological
mass conversion (`units()` in `NewDD/rd_info.c`; box length is one code unit).
Counts are 98,338, 98,211, 98,147, and 98,091 at outputs 00972, 00977, 00980,
and 00984. Their total sink masses are respectively `3.198`, `3.201`, `3.203`,
and `3.206 x 10^6 Msun/h`; the largest individual sink is about `9.3 x 10^4
Msun/h`. This is far below the `10^13 Msun` host-mass threshold and below
0.1 per cent of a `10^8 Msun` stellar-mass-cut satellite. It supports treating
the omission as negligible for the mass scales in this analysis, while
retaining it as a catalogue limitation and checking the final sample's minimum
stellar masses. No sink particles were added to or removed from simulation
outputs.

The 00977 FoF stderr also contains the old `RamsesType`-size warning and
messages that absent SINK/GAS component shards are treated as empty. The active
job predates the ASCII-header fix. The NewDD `.info` source is ASCII and its
header values are parsed as doubles; its writer prints six digits after the
decimal. The warning says binary precision was lost, but does not itself show
that the parsed values or particle records were corrupted. The completed
catalogue still requires structural checks (counts, finite values, duplicate
IDs, and agreement with NewGalFinder output) before scientific use; redo only
if those checks reveal a discrepancy.

## NewDD 00984 complete and shard validation (2026-10-20 22:34 KST)

NewDD job `430626_3` finished with Slurm state `COMPLETED`, exit code `0`,
and printed `NewDD complete snapshot=00984`. I checked all 100 numbered
outputs: all `.info` headers are non-empty; all DM and STAR shard files exist
and are non-empty (zero missing or empty files). The dependent 00984 FoF job
`430622_3` has been released by its dependency and is pending scheduler
priority. No galaxy extraction can start until FoF completes successfully.
The separate 00977 FoF job `430418_1` remains active; its log has advanced to
rank P8's completed halo-boundary check, while its catalogues are still
partial.

## Pipeline status (2026-10-20 22:37 KST)

The 00984 NewDD output has passed the all-shard check and its job exit is `0`.
Its FoF job `430622_3` is dependency-ready and pending only scheduler
priority. The 00977 FoF job `430418_1` remains live at 1:37 elapsed, with no
new output files since the partial catalogues were opened. A current node check
found all ten MPI ranks across `syn01`, `syn03`, and `syn04` using about 99 per
cent CPU, with several ranks' RSS slightly higher than the previous sample.
This supports ongoing computation; no restart or duplicate submission was
made. Rank P8 is the latest logged boundary-check completion. The older
header warning remains under review, so the completed catalogue will receive
post-run validation before acceptance.

## Pipeline status (2026-10-20 22:38 KST)

NewDD 00984 remains successful and fully shard-validated. Its FoF task
`430622_3` is ready but awaiting scheduler priority. FoF task `430418_1` for
00977 remains active: rank P9 has now reported its halo count and passed the
boundary check, joining P1--P8. The two catalogue files are still partial and
no successful completion marker has appeared. All ten MPI ranks were still
using about 99 per cent CPU at the latest node check, so the task remains in
progress rather than being restarted. Galaxy finding and satellite extraction
remain downstream of each successful FoF catalogue.

## Input readiness check (2026-10-20 22:39 KST)

I checked all four snapshots' NewDD directories, not just 00984: for each of
00972, 00977, 00980, and 00984, all 100 numbered `.info`, DM, and STAR shard
files exist and are non-empty. This confirms that the remaining FoF jobs have
complete input shards available. Snapshot 00977 FoF `430418_1` remains active;
all ten ranks were still at about 99 per cent CPU at 1:39 elapsed, with rank P0
not yet reporting a boundary-check completion. The 00984 FoF task and 00972/
00980 replacements are pending Slurm priority.

## Snapshot redshift baseline (2026-10-20 22:43 KST)

I read the `aexp` values directly from the four simulation `info_NNNNN.txt`
files and computed `z = 1/aexp - 1`. Outputs 00972, 00977, 00980, and 00984
have redshifts 0.222581, 0.217836, 0.213243, and 0.208698, respectively.
These outputs span only `0.0139` in redshift, so this set cannot support a
meaningful redshift-evolution claim by itself. Treat them as a closely spaced
late-time sequence for catalogue construction and checks; any evolution
analysis requires additional snapshots over a substantially wider redshift
range.

At 22:42 KST, FoF job `430418_1` was still `RUNNING` at 1:42:42 elapsed. Its
log had not advanced beyond the completed boundary checks for ranks P1--P9
since 22:38, and the P0 boundary check was still absent. No completion marker
had appeared, so the partial 00977 FoF catalogues remain unusable downstream.
The other three FoF jobs remained pending Slurm priority. No job was restarted
or duplicated.

## Subhalo-mass-conditioned relation check (2026-10-20 22:45 KST)

I compared the existing 00972 host-fixed models for satellites with
`Mstar >= 1e8 Msun` and host mass `1e13--1e14 Msun`. With stellar mass as the
only satellite-mass control, the slope of `log10(sigma_star)` on
`log10(Z_star)` is `-0.209` for global 1D dispersion (60 hosts; host-bootstrap
95% interval `[-0.230, -0.191]`). Adding the legacy catalogue's `Mdm` column
(whose aperture and exact definition are unknown) changes the slope to
`-0.037` (3636 satellites, 60 hosts;
95% interval `[-0.057, -0.017]`). The analogous slopes change from `-0.222`
to `-0.031` for the mass-range dispersion and from `-0.182` to `-0.027` for
the dispersion within the stellar half-mass radius; all three conditional
intervals remain slightly below zero. Thus the mass-controlled negative trend
is attenuated by roughly 82--85 per cent when the legacy `Mdm` proxy is
included. Since the source IDL definition is unavailable, this attenuation
cannot be attributed to a well-defined present-day subhalo mass. It is not
evidence for a large independent metallicity--dispersion
relation; at most, the current table suggests a small residual conditional
trend.

As a diagnostic of the fit, after subtracting each host's mean, correlations
among `logZ`, `logMstar`, `logMdm`, and `log(sigma)` are respectively
`corr(logZ, logMstar)=0.807`, `corr(logMstar, logMdm)=0.777`, and
`corr(logMstar, log(sigma))=0.923`. The within-host VIF for `logZ` against
`logMstar` and `logMdm` is 4.57. This covariance warrants care but does not
alone invalidate the conditional fit.

These numbers are still exploratory. The input
`satellite_properties_ascii_sfr_mdm_00972.csv` has the legacy ASCII schema:
its metallicity weighting/convention is unverified and it contains no
particle-ID duplicate audit. The newer extractor reports the birth-metallicity
weighting explicitly and audits repeated star IDs, so this relation must be
recomputed from that validated catalogue before making a paper claim. The
existing analysis therefore motivates a specific next check: test whether the
attenuation persists when the explicitly defined NewGalFinder dark-matter mass
replaces the legacy `Mdm` proxy, and whether any small residual survives in the
audited catalogue and across completed outputs.

## FoF wait check (2026-10-20 22:47 KST)

A second Slurm check one minute after the 22:45 status found 00977 FoF job
`430418_1` still `RUNNING` at 1:46:59; 00984, 00972, and 00980 FoF jobs
remain pending priority. The 00977 stdout still ends at rank P9's boundary
check and shows no P0 completion or final marker. The allocated job has not
reported a terminal state, so it is still being treated as live; no restart or
second submission was made. The `sstat` query did not return step statistics
for this array task, so no current CPU-utilization claim is made.

## Conditional effect-size check (2026-10-20 22:49 KST)

For the legacy 00972 global-dispersion model (`Mstar >= 1e8 Msun`, host mass
`1e13--1e14 Msun`), I residualized metallicity and dispersion against host,
stellar mass, and the legacy `Mdm` proxy (aperture and definition unknown).
The residual metallicity--dispersion
correlation is `r=-0.0755` (partial `R^2=0.0057`). A one-standard-deviation
change in metallicity left after the mass controls corresponds to only
`0.0041 dex`, or about a `0.9 per cent` change in dispersion; the model
residual scatter in dispersion is `0.0537 dex`. Under the same linear model,
a `0.3 dex` metallicity difference corresponds to about a `2.5 per cent`
dispersion difference. The bootstrap interval excludes zero, but this effect
size is small compared with the remaining galaxy-to-galaxy scatter. Thus the
current exploratory evidence supports statistical detectability more than a
large physical effect; assess scientific importance only after recomputing
from the audited NewGalFinder catalogue, whose catalogue-reported subhalo
dark-matter mass will provide a defined comparison quantity.

A targeted check on allocated node `syn01` at 22:49 found four active opFoF
MPI processes at about 99 per cent CPU after 1:48:55 elapsed. This confirms the
00977 job is still computing despite the unchanged log. The expected 00972 and
00977 `GALFIND.DATA` files and `satellite_properties_00972.csv` are absent, so
no updated catalogue can yet be extracted or analysed.

## FoF aggregation and downstream queue (2026-10-20 22:52 KST)

FoF job `430418_1` for 00977 remains `RUNNING` (1:52:21 elapsed). The log has
advanced: P0 passed its halo-boundary check and started preparing to write
6,850,503 isolated halos after the distributed ranks sent their halo data.
The halo catalogue and member-particle file are present at about 7.6 MB and
5.17 GB, with 22:48 modification times; because there is no final completion
marker yet, they are still partial and must not be read downstream.

The user queue already contains dependency-pending NewGalFinder jobs
(`430423_1`, `430423_2`, `430491_0`, `430476_3`), extraction jobs
(`430427_1`, `430427_2`, `430477_3`, `430492_0`), and analysis jobs
(`430431_1`, `430431_2`, `430479_3`, `430493_0`). I did not submit duplicates.
The other three FoF tasks remain pending priority; no 00972/977 audited
satellite catalogue exists yet.

## Dependency audit and active FoF check (2026-10-20 22:55 KST)

I inspected the actual Slurm dependency expressions. The 00977 NewGalFinder
array task `430423_1` is correctly waiting on the matching successful FoF task
`430418_1`; extraction `430427_1` waits on NewGalFinder, and analysis
`430431_1` waits on extraction. The other snapshots also have chains attached
to their replacement FoF jobs: 00972 (`430610_0` -> `430491_0`), 00980
(`430611_2` -> `430423_2`), and 00984 (`430622_3` -> `430476_3`), followed by
their extraction and analysis tasks. No duplicate jobs were submitted.

FoF 00977 is still running and its output files have not changed since 22:48,
but a targeted check on `syn01` at 22:55 found four opFoF processes each using
about 99 per cent CPU at 1:55:16 elapsed. P0's halo aggregation is therefore
still actively computing; the job remains live and its partial files remain
unusable until the success marker and terminal Slurm state appear.

## Meaning of the NewGalFinder dark-matter field (2026-10-20 22:58 KST)

I checked the source used to interpret the new catalogue's dark-matter mass.
In `NewGalFinder/gfind.c`, `GetSubHaloInfo` sums the masses of DM particles
whose `rptl2halonum` equals the galaxy ID; the extractor converts the
catalogue's `Msun/h` value to `Msun` using the snapshot Hubble parameter.
This defines the reported quantity operationally as the total DM mass assigned
to that galaxy by NewGalFinder, without an aperture cut in the accumulation
loop. It is distinct from the legacy IDL column, whose definition is unknown.
The finder membership/binding rules still need to be described accurately
before calling it a gravitationally bound mass.

## FoF process recheck (2026-10-20 23:00 KST)

A follow-up process listing on `syn01` confirms that the 00977 opFoF MPI job
has not stalled: four `opfof.exe` processes remain active at 99.3--99.6 per
cent CPU after 1:59:42 elapsed. One earlier narrow process-name query returned
no matching lines; the full user-process listing resolved that transient
observation. Slurm still reports the task as `RUNNING`, and the 00977 outputs
remain partial. No restart or duplicate job was submitted.

## FoF live-process verification (2026-10-20 23:02 KST)

The 00977 FoF task remains `RUNNING` at about 2:02 elapsed. A second targeted
process check on `syn01` found the same four opFoF ranks using 99.4--99.7 per
cent CPU. The halo catalogue and member-file sizes are unchanged from 22:48,
and the log remains at P0's isolated-halo write stage. This is active CPU
work, not an observed terminal failure; keep the pending NewGalFinder chain
intact and do not consume the partial files.

## FoF output-writer diagnostics (2026-10-20 23:03 KST)

To distinguish a long write stage from a hung process, I inspected P0 process
`1762664` on `syn01`. It is in Linux `R` (running) state, uses about 11.9 GB
RSS, and has both 00977 FoF output files open for writing. Its cumulative
`wchar` is about 6.34 GB. GPFS reports zero `write_bytes` in `/proc/io`, so that
counter is not useful here; the active CPU state and open output descriptors
show that P0 remains in its output/aggregation stage. I will compare cumulative
I/O and file sizes on a later check before judging whether it is progressing.

## FoF CPU-progress sample (2026-10-20 23:07 KST)

I sampled P0's CPU and output counters 61 seconds apart. Its CPU time rose from
2:05:41 to 2:07:07 (about 86 seconds of CPU over the 61-second interval,
consistent with one fully used core), while `wchar`, `syscw`, and both output
file sizes remained unchanged. Slurm still reports `430418_1` as `RUNNING`.
This confirms sustained CPU work in the large isolated-halo loop, but gives no
completion-time estimate; no restart has been made.

## FoF writer cost from source (2026-10-20 23:08 KST)

The active binary is built with `Treewalk.fof.ordered.o`. Its P0 writer loops
through the full local halo-bound array; for each accepted group, `FWRITE`
allocates a temporary particle buffer, scans the members by particle type,
and writes the rearranged records. For 6.85 million local halo entries this
can create a long single-core phase before the next log line. This source path
is consistent with the observed high P0 CPU and no new stdout marker; it does
not prove when the loop will finish. Keep waiting on the live process and
preserve the partial catalogue.

## FoF status check (2026-10-20 23:10 KST)

At 23:10, Slurm still showed 00977 FoF `430418_1` as `RUNNING` at 2:09:59.
P0 remained in `R` state at 99.6 per cent CPU with about 11.9 GB RSS. Its
cumulative `wchar` and `syscw`, plus both FoF file sizes, matched the preceding
sample; the output files are still partial. The other three FoF tasks remain
pending priority, while the matching downstream jobs remain dependency-pending.
No files were passed to NewGalFinder and no job was restarted.

## FoF second CPU sample (2026-10-20 23:11 KST)

At 23:11, P0's CPU time had advanced to 2:10:50 and the process was still
running at 99.6 per cent CPU. The 00977 files remained 7,552,732 bytes and
5,169,371,272 bytes; cumulative `wchar` and `syscw` were unchanged. Slurm showed
`430418_1` still running at 2:11:11. This is continued single-core execution
inside the isolated-halo traversal, with no additional output progress yet.

## FoF output growth and Tree-building progress (2026-10-20 23:12 KST)

The 00977 FoF job `430418_1` is still `RUNNING`, but this check shows new
progress. Both output files grew since the prior sample: the halo catalogue
from 7,552,732 to 8,361,148 bytes and the member file from 5,169,371,272 to
5,830,589,688 bytes. P0's cumulative `wchar` rose from 6,343,109,984 to
8,944,413,762 bytes and `syscw` from 164,750 to 172,797. Its CPU remains near
99.6 per cent and RSS is about 12.3 GB. The latest log reports Tree-building
completion for P0--P6 and P9. The long isolated-halo stage has therefore
advanced into subsequent FoF work. Outputs are still incomplete: no final
success marker, NewGalFinder catalogue, or satellite-property catalogue is
present. The existing dependency chain remains pending; no restart or duplicate
submission was made.

## FoF tree-build completion across ranks (2026-10-20 23:13 KST)

A later 00977 log check shows all ten MPI ranks P0--P9 have now reported
`passed Tree building`. The FoF job itself remains `RUNNING` at 2:13:24, and
no final completion marker or downstream GALFIND file exists yet. The two FoF
output files reached 8,361,148 bytes and 5,830,589,688 bytes at 23:11; they are
still partial until the job exits successfully. The dependent NewGalFinder,
extraction, and analysis stages remain queued; no duplicate work was added.

## FoF post-tree link-assignment stage (2026-10-20 23:14 KST)

All ten ranks have passed `FoF_Make_Tree`. The active source then loops over
local particles not yet included and assigns FoF links (`new_fof_link` or
`pnew_fof_link`) before halo-boundary checking and the next write stage. At
23:14, P0 remained at 99.6 per cent CPU and its RSS had risen to about 13.3 GB;
`wchar` and `syscw` were unchanged since the prior sample, and the log had not
yet printed per-rank halo totals. This is consistent with computation in the
post-tree particle-link loop. The FoF job remains live; its files are still
partial.

## FoF post-tree loop remains active (2026-10-20 23:15 KST)

One minute after all ranks passed Tree building, 00977 FoF remained `RUNNING`
at 2:15:26. P0's CPU time advanced to 2:15:02, with 99.6 per cent CPU and
RSS about 13.4 GB; its cumulative output counters and the FoF files were
unchanged, and the log had not yet reported halo totals. This is consistent
with the source's post-tree per-particle linking loop. The three other FoF
jobs still await priority, and the 00977 downstream chain remains dependency-
pending.

## First post-tree rank completes halo checks (2026-10-20 23:16 KST)

The 00977 log has advanced beyond Tree building: rank P2 completed FoF link
assignment for 5,177,392 halos from 15,534,395 particles and passed its halo
boundary check. The task remains `RUNNING`; the latest check has not yet shown
P2's isolated-halo save/transfer marker or a final job marker. The outputs and
satellite catalogues remain incomplete, and the dependent NewGalFinder chain
is still waiting.

## FoF rank activity after P2 boundary check (2026-10-20 23:18 KST)

The 00977 job remains live at 2:17:48 and the FoF outputs have not changed
since 23:11. On `syn01`, all four allocated opFoF processes (P0--P3) were in
`R` state at 99.3--99.6 per cent CPU, with RSS values from 4.2 to 13.5 GB.
P2 has completed its local halo-boundary check, while P0 has not yet logged its
own halo total; the ranks are still performing substantial local work. No
GALFIND or satellite-property output exists, and the dependency chain remains
pending.

## Second rank completes post-tree FoF (2026-10-20 23:19 KST)

The 00977 FoF log has advanced again: rank P3 completed FoF linking for
5,961,065 halos from 24,434,175 particles and passed its halo-boundary check,
after P2's earlier completion. The job is still `RUNNING` at 2:19:07; neither
output file has changed since 23:11 and no global completion marker is present.
The dependency-pending NewGalFinder/extraction/analysis jobs remain in place.

## Third rank completes post-tree FoF (2026-10-20 23:20 KST)

Rank P1 has now also completed FoF linking and passed its boundary check:
5,511,013 halos from 23,557,896 particles. P2 and P3 had already completed
this stage. `430418_1` remains `RUNNING` at 2:19:48; the FoF files still have
not grown beyond their 23:11 sizes, and the other ranks have not yet reported
halo totals. No final marker or GALFIND output exists; downstream jobs remain
blocked on the successful FoF dependency.

## Fourth rank completes post-tree FoF (2026-10-20 23:20 KST)

The latest 00977 log reports rank P6 completed link assignment for 6,164,772
halos from 23,357,879 particles and passed its boundary check. P1, P2, P3,
and P6 have now reported halo totals and passed this check; P0 and the remaining
ranks are still working. On `syn01`, the four local opFoF ranks remain near
99 per cent CPU, with P0 RSS about 13.8 GB. Output sizes are unchanged from
23:11 and no final or NewGalFinder output exists yet.

## Fifth rank completes post-tree FoF (2026-10-20 23:24 KST)

The 00977 log now reports rank P4 completed link assignment and passed its
boundary check: 6,670,232 halos from 31,459,202 particles. The job remains
`RUNNING` at 2:23:51. P0, P1, P2, P3, P4, and P6 have reported totals; the
other ranks are still working. Both FoF files remain at their 23:11 partial
sizes, and the already-submitted NewGalFinder chain remains dependency-blocked.

## 00977 FoF still computing (2026-10-20 23:26 KST)

A scheduler and process recheck confirms `430418_1` is still `RUNNING` (2:26:08).
The four opFoF ranks on `syn01` remain CPU-active at about 99 per cent; the log
has not advanced beyond P4 and the two partial output files retain their 23:11
sizes. The queued downstream jobs are still correctly waiting on FoF completion.

## Galaxy-finder resource request audit (2026-10-20 23:28 KST)

Slurm lists no CPU-only partition on `syntax`; the available compute partitions
are all GPU-equipped. The pending NewGalFinder jobs request CPUs and memory on
`a100` but have no `--gres=gpu` request (`ReqTRES=cpu,mem,node,billing`, with no
`gres/gpu`). Their current scheduler reason is the FoF dependency, not a GPU
quota reason. The configured request is two A100 nodes, 100 CPUs, and 600 GB
per node, so this CPU-only MPI job still needs Slurm and currently targets
GPU-equipped, high-memory nodes; running it on the login node is not suitable.

## Sixth rank completes post-tree FoF (2026-10-20 23:29 KST)

Rank P5 has completed FoF linking and passed its boundary check, with 5,771,243
halos from 34,144,265 particles. P0, P1, P2, P3, P4, P5, and P6 have now
reported totals; P7, P8, and P9 have not. Job `430418_1` remains `RUNNING` at
2:29:14. Its catalogue and member files are still the partial files last
written at 23:11, and the dependency-queued downstream jobs have not started.

## Correction and active-rank check (2026-10-20 23:31 KST)

Correction to the 23:29 entry: only six ranks (P1--P6) have reported halo
counts; P0 has not, along with P7--P9. A process sample across `syn01`,
`syn03`, and `syn04` finds all ten `opfof.exe` MPI ranks still running at about
99.5 per cent CPU. The log and partial output files have not advanced since
23:29/23:11, respectively. This confirms active computation rather than an
idle or exited process, but completion time remains unknown.

## 00977 FoF late-rank workload visible (2026-10-20 23:32 KST)

The log has advanced with the spatial ranges for two late slabs: P9 reports
60,142,363 particles and P8 reports 114,890,669. All ten ranks have passed
Tree building, but only P1--P6 have completed local FoF linking and boundary
checks. `430418_1` remains `RUNNING` at 2:32:27; outputs remain partial. This
shows the late ranks are processing substantially larger particle loads, which
is consistent with their longer completion time.

## FoF shard-loop progress from source and live log (2026-10-20 23:38 KST)

I checked the active `opfof.c` loop: for 100 input shards and 10 MPI ranks,
each rank processes a consecutive set of 10 files, one at a time. The live log
is on the third file per rank (`.00002`, `.00012`, ..., `.00092`): P1--P6
have completed the third file's halo-boundary checks; P0 and P7--P9 are still
processing it. Thus 26 of 100 files have completed by file count, with four
more in flight; this is not a wall-time fraction because particle loads vary
strongly by shard. `WriteFinalHalo` is only called on each rank's last file,
so no FoF output can be considered final yet. The job remains live and its
already-submitted downstream jobs correctly wait for success.

## Bootstrap uncertainty added for mass-controlled relations (2026-10-20 23:42 KST)

Static review of the new-catalogue analysis found that it reported cluster-
robust standard errors for slopes controlling for finder-assigned subhalo
mass, but omitted the host-level bootstrap intervals used for the primary
stellar-mass-only slope. I added a block-bootstrap calculation that resamples
whole hosts and computes intervals for the metallicity coefficient while
controlling jointly for stellar mass and subhalo total mass, and while
controlling for stellar mass and subhalo dark-matter mass. The results will be
written to `satellite_relations_*.csv` as separate low/high 95 per cent
columns. The already-queued analysis jobs call this script at runtime and will
use the updated code after their catalogue dependencies complete. This does
not change any legacy-catalogue result or constitute a new measured result.

## Third FoF shard-group advances (2026-10-20 23:46 KST)

Rank P9 has now completed its third input file (`SN.00977.00092.dat`) and
passed the halo-boundary check, reporting 8,982,772 halos from 60,142,363
particles. P1--P6 and P9 are complete for this third shard group; P0, P7, and
P8 have not yet reported their halo totals. The 00977 job remains `RUNNING` at
2:46:11, and downstream NewGalFinder jobs remain dependency-pending. The FoF
catalogue files are still not final.

## Effect-size and residual-scatter reporting added (2026-10-20 23:48 KST)

The relation table will now also report the predicted percent change in
internal stellar dispersion for a 0.3-dex metallicity difference, with the
host-bootstrap interval transformed to that scale. It will report the
within-host residual scatter in `log10(sigma)` after controlling for host and
stellar mass, and after additionally controlling for subhalo total mass or
subhalo dark-matter mass. This makes the statistical slope interpretable as a
physical effect size relative to the remaining galaxy-to-galaxy scatter. These
are output definitions only; no new catalogue measurements have been made.

## Bibliography metadata audit and corrections (2026-10-20)

Created `references.bib` with six papers whose metadata were checked against
Oxford Academic journal records. Two author attributions in the older
exploratory-literature notes above need correction: DOI `10.1093/mnras/stw2329`
is Bah\'e et al. (2017), *The origin of the enhanced metallicity of satellite
galaxies*, not Lacerna et al.; DOI `10.1093/mnras/sty2913` is Buck et al.
(2019), *NIHAO XV: the environmental impact of the host galaxy on galactic
satellite and field dwarf galaxies*, not Obreja et al. The 2023 Wang, Wang &
Chen citation remains excluded from the bibliography because its record has
not yet been verified against a primary journal source. These bibliographic
corrections do not change the legacy-catalogue exploratory measurements.

## Direct prior work on metallicity and internal kinematics (2026-10-20)

A targeted literature search confirms two direct precedents. Spindler & Wake
(2017) compare SDSS stellar metallicities and central stellar velocity
dispersions for central and satellite galaxies, including star-forming and
quiescent subsamples (Oxford Academic,
https://academic.oup.com/mnras/article/468/1/333/3038254). Fattahi et al.
(2018) explicitly study the dwarf-galaxy `[Fe/H]`--`V_{1/2}` relation and show
how tidally stripped satellites can flatten its low-velocity end (Oxford
Academic, https://academic.oup.com/mnras/article/476/3/3816/4875945). The
present paper therefore cannot claim that the metallicity--dispersion
connection is previously unstudied. Its potential contribution is a
conditional analysis of a different simulated satellite population with
stellar-mass, host, and finder-assigned subhalo-mass controls, plus a
cross-output replication test. That distinction remains provisional until
the catalogue is complete and the literature review is broadened.

## 00977 FoF live recheck (2026-10-20 23:55 KST)

Slurm confirms `430418_1` remains `RUNNING` at 2:55:18 on `syn[01,03-04]`.
The live log still ends with P9 completing the third shard group and P0/P7/P8
not yet reporting that group's FoF totals. The final halo and member files
remain at 8,361,148 and 5,830,589,688 bytes, last written at 23:11; they are
partial and must not be read by NewGalFinder or extraction. The queued
NewGalFinder/extraction/analysis stages remain dependency-pending. This check
confirms the producer job is still live; it does not provide a completion
estimate or a new measurement.

## Stellar-mass consistency diagnostic added to the catalogue pipeline (2026-10-20)

Static review found that extraction used NewGalFinder's stellar mass but did
not independently compare it with the sum of current masses in the serialized
stellar-particle records. The extractor now writes that particle-mass sum and
its ratio to the finder stellar mass. The relation-analysis output summarizes
the number of objects with this check, the median ratio and absolute fractional
difference, and the fractions within 1 and 5 per cent. This is a diagnostic
only; it does not remove objects or define a pass threshold. Queued extraction
and analysis jobs load these scripts when they run, so the changes are ready
for the existing dependency chain. No calculations or tests were run in this
edit; the check will become measured evidence only after successful job output
and catalogue audit.

The NewGalFinder source `src/catalogue.c` confirms that each `SubInfo.mstar`
is accumulated from the same `raw[i].p.star.mass` values copied into its
stellar-particle payload. Therefore, a ratio near unity is expected up to
floating-point summation and invalid-record handling; any larger mismatch will
need investigation before using the catalogue mass as a validated quantity.

## 00977 FoF remains active; no final files yet (2026-10-20 23:58 KST)

A fresh Slurm query shows `430418_1` still `RUNNING` at 2:58:02 on
`syn[01,03-04]`. The live output has not advanced beyond the third shard-group
FoF totals, and both output files retain their 23:11 partial sizes. The
NewGalFinder, extraction, and analysis jobs remain dependency-pending. No
catalogue extraction was started from these partial FoF files.

## Atomic publication of extraction and analysis outputs (2026-10-21)

The extraction and summary scripts previously wrote directly to final CSV paths.
An interrupted run could leave a truncated file that the no-overwrite guards
would then refuse to replace. Both scripts now write to unique same-directory
partial files and publish the final path with an exclusive hard link only after
writing succeeds. Existing final outputs remain protected; failed-run partials
remain separately identifiable. This edit was not executed or tested, and the
queued jobs will use it when they reach these scripts.

## FoF and downstream pipeline state (2026-10-21 00:01 KST)

Slurm still reports `430418_1` as `RUNNING` after 3:01:35 on three `a10`
partition nodes. Its allocation is 60 CPUs and 768 GB RAM with no GPU in
`ReqTRES` or `AllocTRES`. The log and partial FoF file sizes have not changed
since the 23:58 check. NewGalFinder and extraction arrays remain pending on the
FoF dependency, so no complete 00977 satellite catalogue exists yet.

## 00977 FoF resource activity and prior-attempt audit (2026-10-21 00:04 KST)

The active Slurm step `430419.0` is accumulating CPU time (`AveCPU=10:09:09`)
and reports about 401 GB average RSS and 64 GB maximum RSS per task in `sstat`;
this supports that the job is actively processing despite no recent log or
output-file growth. Its FoF files remain partial, so no downstream reader was
started. `sacct` also confirms earlier snapshot-972 attempts `430347_0` and
`430417_0` did not produce successful FoF outputs; the 430347 log records
incorrect file-size checks on STAR shards 00060 and 00070, and 430417 was
cancelled. Those stale partial directories remain separate from the active
00977 output and were not used.

## NewDD inputs audited for all four selected outputs (2026-10-21 00:06 KST)

Each selected snapshot (00972, 00977, 00980, 00984) has a successful
`NewDD complete` marker in its Slurm log. For every snapshot, all 100 DM and
100 STAR shard files are present, and every file size is divisible by 168
bytes, matching the configured double-precision, nine-element `StarType` and
`DmType` record size. This resolves the old 00972 shard-size concern for the
current on-disk NewDD inputs; it does not rehabilitate the cancelled FoF
attempts or prove the FoF outputs complete. The active 00977 job and queued
FoF jobs still need to finish before catalogue extraction.

## 00977 FoF has begun assembling rank outputs (2026-10-21 00:07 KST)

The live log has advanced: rank P0 completed its local boundary check and is
preparing to write 6,101,967 halos; P1--P6 and P9 report that their isolated
halo data have been sent to and saved by P0. P7 and P8 have not yet reported
completion. The FoF output files grew to 13,948,092 and 16,968,670,456 bytes
at 00:05:56, so they are being written but are still not final. Slurm still
lists the job as running, and downstream jobs correctly remain blocked on its
successful completion.

## Correction: FoF and subhalo masses omit gas and sinks (2026-10-21 00:10 KST)

Source inspection of `newdd_array.sbatch` and `NewDD/newdd.c` shows that
`NEWDD_SKIP_AMR=1` jumps directly to the particle dump, skipping AMR, hydro/gas,
and sink reads. The NewDD files for the selected outputs therefore contain DM
and star particles only. Consequently, `opFoF` host mass and NewGalFinder
assigned subhalo total mass are DM+star particle sums, not full physical total
masses. I corrected the extractor column names, analysis metadata and control
labels, and the manuscript plan. The earlier exploratory note calling the
FoF host mass a matched total-mass proxy is not valid as written; comparisons
to the legacy catalogue's host total mass are not one-to-one and must be
revisited. The analysis will report host selection as DM+star mass cuts and
include boundary sensitivity. This does not change the current 00977 FoF run,
but queued extraction/analysis read the updated Python scripts at runtime.

## 00977 FoF aggregation still running (2026-10-21 00:10 KST)

The FoF job is still `RUNNING` at 3:10:21; its MPI step has accumulated
10:31:08 of CPU time and about 400 GB average RSS, so the P0 assembly stage is
active rather than abandoned. The log still shows P0 preparing the 6.1-million
halo local file while P7/P8 have not reported their sends. The partial FoF
files remain unchanged since 00:05:56, and no GALFIND output exists. No
consumer was started.

## Analysis-method review and FoF wait (2026-10-21 00:15 KST)

Static review confirmed that the simulation catalogue estimates internal stellar
velocity dispersion from mass-weighted particle velocities and separately
measures the dispersion inside the 3D stellar half-mass radius. The primary
within-host model controls for stellar mass, and the analysis excludes
satellites with duplicated stellar-particle IDs. I found and corrected a stale
comment in the legacy-catalogue adapter that incorrectly implied its host mass
was equivalent to the NewDD FoF DM+star mass; the legacy mass definition remains
unverified and must not be presented as a matched control. No tests or
calculations were run during this review.

The 00977 FoF Slurm task `430418_1` remained `RUNNING` at 00:14 KST, but its
MPI log and output file sizes had not advanced since 00:05:56. The latest log
shows P0 has received isolated-halo data from P1–P6 and P9; no final catalogue
or GALFIND output exists yet. Existing dependent jobs remain pending, and no
replacement jobs were submitted.

## Rank-deficient fits are now treated as unmeasurable (2026-10-21 00:16 KST)

Static review found that the primary host-fixed-effects fit could return a zero
cluster error when only one host contributed, and bootstrap resamples with
singular predictor matrices could either crash or return a non-identifiable
coefficient. The analysis now returns missing fit values when the original fit
is rank-deficient or has fewer than two hosts, and skips singular bootstrap
draws before forming percentile intervals. No tests or scientific calculations
were run; the edit is awaiting use by the already queued analysis jobs.

## 00977 FoF aggregation wait traced to missing worker completion messages (2026-10-21 00:18 KST)

The Slurm task `430418_1` is still `RUNNING`, while its FoF files and output
log remain unchanged since 00:05:56. The latest aggregation log records P1--P6
and P9 completing their isolated-halo sends; P7 and P8 have no corresponding
completion lines. Static inspection of `opFoF/Treewalk.fof.ordered.c` shows
that rank 0 waits to receive the end-of-stream marker from every worker before
closing the FoF files. This is consistent with the incomplete catalogue being
held open while waiting on P7/P8, but does not prove whether those ranks are
computing, blocked in MPI, or stalled on I/O. No live job was cancelled or
restarted, and no dependent jobs were duplicated.

## Live rank inspection: FoF workers are actively computing (2026-10-21 00:20 KST)

A lightweight `srun --overlap` inspection of the existing allocation (underlying
job 430419, array element 430418_1) found all ten `opfof.exe` processes across
syn01, syn03, and syn04 in runnable state, each using about 99 per cent of one
CPU at the sampling instant. This shows that the task is computationally active;
the missing P7/P8 completion lines do not indicate an idle MPI deadlock. The
most likely interpretation is that workers are still processing the final
large slab/halo set, though the exact rank-to-process mapping was not checked.
No new compute job was submitted and no running process was interrupted.

## 00977 FoF rank 7 completed; aggregation continues (2026-10-21 00:21 KST)

The FoF log advanced at 00:20:16: rank P7 checked its halo boundary, sent its
isolated-halo data to rank 0, and reported completion. The partial output files
grew to 14,726,332 bytes (halo catalogue) and 18,003,090,680 bytes (member
particles), confirming that the master is continuing the aggregation. Slurm
still reports task `430418_1` as `RUNNING`; NewGalFinder task `430423_1` remains
pending on its dependency. The outputs are still partial and were not read by
downstream jobs.

## 00977 FoF remains live after P7; no further output in the next minute (2026-10-21 00:22 KST)

A one-minute recheck found task `430418_1` still `RUNNING`, with NewGalFinder
still dependency-pending. The stdout and both FoF output files had not changed
since P7's 00:20:16 completion message. The catalogue remains partial; this
short observation interval does not establish a second failure or justify a
restart.

## 00977 FoF processes continue consuming CPU; outputs not yet updated (2026-10-21 00:23 KST)

A second lightweight inspection of the existing allocation, about one minute
later, again found all ten `opfof.exe` ranks runnable at roughly 99 per cent of
one CPU. Their cumulative CPU times had advanced to about 3 h 22 min, while the
stdout and FoF output mtimes remained at 00:20:16. This confirms sustained
CPU-bound work despite the lack of new log lines; the task remains live and the
partial catalogue is still not safe for downstream reading.

## 00977 FoF CPU activity persists without file growth (2026-10-21 00:25 KST)

A later process sample again found all ten MPI ranks runnable near 99 per cent
CPU, with cumulative CPU time increasing by about two minutes over the observation
interval. The stdout and FoF file sizes still matched the 00:20:16 values. The
job is CPU-active rather than idle, but the current log/output cadence is too
slow to estimate completion; downstream NewGalFinder remains gated on success.

## Rank-resolved counters show CPU-bound FoF work (2026-10-21 00:27 KST)

Mapped `PMI_RANK` for all ten `opfof.exe` processes on the existing allocation.
Ranks 0--9 all remained runnable at about 99 per cent CPU. Comparing `/proc`
counters 60 seconds apart, their CPU times advanced by roughly 80 seconds,
while each rank's `rchar` and `wchar` stayed unchanged; the FoF log and output
sizes also remained unchanged. This points to a long CPU-bound calculation in
the current stage rather than ongoing file I/O. The processes are live; the
partial outputs must not be consumed.

## FoF execution is CPU-only and the binary has no OpenMP build flags (2026-10-21 00:30 KST)

The active Slurm allocation remains CPU/memory-only (`cpu=60`, `mem=768G`,
three nodes; no GPU TRES in the prior `scontrol show job` record). Static
inspection of the compiled program's source/build rules found no OpenMP
pragmas or OpenMP compiler/link flags, consistent with the ten MPI ranks each
using about one CPU. The existing allocation was not altered; any resource
change would require safely replacing queued jobs and their dependent job IDs.
The current outputs remain partial, so downstream work is still pending.

## Symbolized rank stacks locate the long-running FoF stage (2026-10-21 00:33 KST)

A brief GDB attach/detach on PMI rank 8, with the executable's debug symbols
loaded, located it in `new_fof_link()` at `Treewalk.fof.ordered.c:919`, checking
particle coordinates while constructing a FoF group. This corrects the earlier
unsymbolized-address guess that it was in `FoF_divide_node()`. A symbolized
sample of rank 0 located its active instruction in Intel MPI's libfabric CQ
error/progress path (`MPIDI_OFI_handle_cq_error`), consistent with the master
polling while workers finish. These samples indicate a slow FoF-linking
straggler rather than a demonstrated deadlock. GDB detached cleanly; the job
was not interrupted. The output files remain partial.

## Correct-binary stack trace identifies the 00977 bottleneck (2026-10-21 00:37 KST)

The executable actually loaded by the running ranks is the preserved
`opfof.exe.pre_ascii_header_fix_20261020` (280,704 bytes, 20:29 build); the path
`opfof.exe` was replaced at 22:22 while the job retained the old inode. The
previous stack interpretation based on the replacement executable is therefore
superseded. With the exact loaded binary, rank 8 is in `pnew_fof_link()` at
`Treewalk.fof.ordered.c:812`; its current group has `ncount=23,365,419` and
`now=21,710,300`. Rank 0 is in `PMPI_Probe(source=MPI_ANY_SOURCE, tag=987)` inside
`WriteIsolatedHalo`, waiting for the worker's completion signal. This directly
explains the partial files: rank 8 is still traversing a very large FoF group,
while the master waits for it. Both GDB sessions detached cleanly. The run was
not stopped, and no downstream job was duplicated.

## 00977 large FoF component is advancing (2026-10-21 00:38 KST)

Using the exact loaded binary, a second GDB sample of rank 8 about one minute
after the first showed `pnew_fof_link()` still processing halo index 5, with
`ncount` increasing from 23,365,419 to 23,678,779 and `now` from 21,710,300 to
21,732,777. This is measurable progress through a very large connected group,
not a frozen counter. Rank 0 remains in `PMPI_Probe` awaiting a worker message;
Slurm still shows the FoF task running and NewGalFinder dependency-pending. The
FoF output files remain partial and unchanged, so no catalogue extraction has
started.

## Rank 8 continues the large-component traversal (2026-10-21 00:39 KST)

A new stack sample using the exact loaded binary again found rank 8 in
`pnew_fof_link()`, still at local halo index 5. Its counters rose to
`ncount=23,777,680` and `now=21,733,945` (from 23,678,779 and 21,732,777 in the
previous sample). The FoF task remains `RUNNING`, while the output sizes and
log are unchanged. This confirms slow forward progress inside the large
component; the catalogue is not yet complete.

## 00977 old-binary input/header cross-check and continued progress (2026-10-21 00:42 KST)

The running pre-header-fix binary did emit the expected `RamsesType` size
warning, but the startup log reports box size 44.1713, H0 67.4, Omega_m 0.315,
Omega_Lambda 0.685, and a=0.821128, consistent with the ASCII `.info` header.
The partial FoF catalogue begins with the seven binary floats that NewGalFinder
reads. This warning alone does not establish bad input or output; structural
validation remains required after completion.

Using the exact loaded executable, rank 8 remains in `pnew_fof_link()` at halo
index 5. Its counters advanced from `ncount=23,777,680`, `now=21,733,945` to
`ncount=24,124,720`, `now=21,743,929`. Slurm still reports the FoF job running
and NewGalFinder dependency-pending; the FoF files remain partial.

## 00977 giant-component progress sample (2026-10-21 00:43 KST)

Slurm reports FoF array element `430418_1` still `RUNNING` (elapsed 3:42:48)
and NewGalFinder elements `430423_1` and `430423_2` still dependency-pending.
Using GDB with the exact executable loaded by rank 8, it remains in
`pnew_fof_link()` at halo index 5; counters advanced to `ncount=24,519,333` and
`now=21,745,492`. The FoF catalogue, member-particle file, and log have not
grown since 00:20:16 KST, so the partial files remain unusable downstream.

## 00977 continued rank-8 sample (2026-10-21 00:44 KST)

The FoF task remains `RUNNING` (3:44:09 elapsed); its dependent NewGalFinder
array elements remain pending. The exact loaded executable still places rank 8
in `pnew_fof_link()` at halo index 5. Counters advanced from the prior sample
to `ncount=24,597,236` and `now=21,747,431`. Output files and log remain at
their 00:20:16 sizes/timestamps; no downstream use is valid yet.

## 00977 continued rank-8 sample (2026-10-21 00:46 KST)

At 3:46:10 elapsed, FoF remains `RUNNING`; NewGalFinder array elements
`430423_1` and `430423_2` remain dependency-pending. Rank 8 is still in
`pnew_fof_link()` at halo index 5, with counters advanced to
`ncount=25,135,474`, `now=21,764,672` (from 24,597,236 and 21,747,431).
The FoF output and log timestamps remain 00:20:16 KST, so downstream processing
is still correctly waiting on the incomplete output.

## 00977 rank-8 sample (2026-10-21 00:47 KST)

FoF remains `RUNNING` at 3:47:16 elapsed; NewGalFinder elements remain
`PENDING (Dependency)`. Rank 8 is still in `pnew_fof_link()` at halo index 5.
Since the 00:46 sample, `ncount` increased from 25,135,474 to 25,157,826 and
`now` from 21,764,672 to 21,766,985. The counter is advancing more slowly in
this interval, but is not stationary. FoF output and log remain unchanged at
their 00:20:16 timestamps; the generated files are still incomplete.

A direct process check on `syn04` at 00:48 KST confirms rank 8 is alive in
running state and consuming 99.4% CPU (elapsed 3:48:23), consistent with the
counter increments rather than a sleeping or terminated process.

## 00977 rank-8 sample (2026-10-21 00:49 KST)

At 3:49:20 elapsed, Slurm still reports FoF `RUNNING`; both reported NewGalFinder
elements remain dependency-pending. Rank 8 remains in `pnew_fof_link()` at
halo index 5, currently traversing the neighbor-distance check (source line 850),
with counters `ncount=25,869,750`, `now=21,779,093`. This is an increase of
711,924 from the previous sample (`ncount=25,157,826`, `now=21,766,985`); process
CPU remains 99.4%. The FoF outputs/log still have their 00:20:16 timestamps.

## 00977 rank-8 sample (2026-10-21 00:50 KST)

At 3:50:51 elapsed, Slurm still reports FoF `RUNNING`; NewGalFinder elements
`430423_1` and `430423_2` remain dependency-pending. Rank 8 is alive at 99.4%
CPU and remains in `pnew_fof_link()` at halo index 5. Counters advanced from
`ncount=25,869,750`, `now=21,779,093` to `ncount=26,021,985`, `now=21,872,209`.
The FoF catalogue, member file, and log remain unchanged at their 00:20:16 KST
timestamps, so the output is still partial and cannot yet feed NewGalFinder.

## 00977 FoF and four-output dependency check (2026-10-21 00:52 KST)

The live queue confirms the 00977 FoF task `430418_1` is `RUNNING`; its
NewGalFinder task `430423_1`, matching extraction `430427_1`, and analysis
`430431_1` remain dependency-pending. The other output chains are also present
and waiting on their matching FoF replacements; no duplicate jobs were added.
Rank 8 remains in `pnew_fof_link()` at halo index 5, with `ncount=26,524,709`
and `now=21,875,757`; its process is still consuming 99.4% CPU. FoF files and
log remain unchanged at 00:20:16 KST and are still partial.

## Methods drafting while catalogues are pending (2026-10-21 00:53 KST)

Created `methods_draft.md` from the implemented extraction and analysis
conventions. It documents central/satellite assignment, the particle-ID audit,
RAMSES birth-metallicity proxy, global and half-mass-radius internal stellar
dispersions, the four-output host-fixed-effects models, host bootstrap, and
limits of the FoF/subhalo mass definitions. It contains no sample counts or
scientific results; those remain reserved for validated NewGalFinder outputs.

## 00977 FoF progress after methods draft (2026-10-21 00:54 KST)

The dependency chain is unchanged: FoF `430418_1` is `RUNNING` at 3:54:29;
NewGalFinder and downstream extraction/analysis tasks remain pending on their
matching prerequisites. Rank 8 remains at halo index 5 in `pnew_fof_link()`;
its counters advanced to `ncount=27,719,180`, `now=21,882,004` from
26,524,709 and 21,875,757. It is alive at 99.4% CPU. The catalogue, member
file, and log still have their 00:20:16 timestamps, so no catalogues can yet
be extracted from this output.

## 00977 FoF sample (2026-10-21 00:55 KST)

At 3:55:36 elapsed, Slurm still reports `430418_1` `RUNNING`; the matching
NewGalFinder, extraction, and analysis tasks remain dependency-pending. Rank 8
continues in `pnew_fof_link()` at halo index 5, with `ncount=28,291,899` and
`now=21,889,209` (previously 27,719,180 and 21,882,004). Its CPU use remains
99.4%. FoF catalogue/member outputs and the log retain their 00:20:16 KST
timestamps and are not complete inputs yet.

## 00977 FoF sample (2026-10-21 00:56 KST)

The queue still shows 00977 FoF `430418_1` `RUNNING` (3:56:44 elapsed), with
its dependent NewGalFinder/extraction/analysis tasks waiting; the other output
chains have not started. Rank 8 remains in the halo-5 `pnew_fof_link()` path,
with counters advanced to `ncount=28,759,569`, `now=21,901,822` from
28,291,899 and 21,889,209. Process state remains running at 99.4% CPU. The
FoF files and log have not changed since 00:20:16 KST and remain partial.

## 00977 FoF progress sample (2026-10-21 00:57 KST)

Slurm reports `430418_1` `RUNNING` at 3:57:51; all dependent galaxy-finder,
extraction, and analysis tasks remain pending. Rank 8 is alive at 99.4% CPU in
`pnew_fof_link()` at halo index 5. Counters advanced from `ncount=28,759,569`,
`now=21,901,822` to `ncount=29,458,807`, `now=22,463,536`. FoF catalogue,
member file, and log timestamps remain 00:20:16 KST; output is still partial.

## 00977 FoF progress and downstream path check (2026-10-21 00:59 KST)

The live queue still shows 00977 FoF `430418_1` `RUNNING` at 3:59:04, with its
NewGalFinder/extraction/analysis chain waiting on dependencies. Rank 8's counters
advanced from `ncount=29,458,807`, `now=22,463,536` to `ncount=30,168,537`,
`now=23,235,491`; it remains at halo index 5 and 99.4% CPU. FoF outputs remain
unchanged at their 00:20:16 KST timestamps. I checked the four expected primary
satellite catalogue/relation output paths (`00972`, `00977`, `00980`, `00984`);
none exists, so the queued extraction and analysis jobs will not be blocked by
stale output files that their scripts refuse to overwrite.

## 00977 rank-8 caller-frame sample (2026-10-21 01:00 KST)

FoF `430418_1` remains `RUNNING` at 4:00:07; its downstream jobs remain
`PENDING (Dependency)`. A GDB sample initially landed in nested `pfof_open()`;
its caller frame confirmed rank 8 is still in `pnew_fof_link()` at halo index 5.
The caller counters are now `ncount=30,819,585`, `now=27,144,394`, advanced
from `ncount=30,168,537`, `now=23,235,491`. The rank remains CPU-bound. The
latest output-file check still shows only the old partial files from 00:20:16.

## 00977 live-process check (2026-10-21 01:01 KST)

Slurm still reports FoF `430418_1` `RUNNING` at about 4:01 elapsed; downstream
jobs remain dependency-pending. GDB samples caught rank 8 in the `pnew_fof_link()`
loop (including nested `pfof_open()` calls), and `ps` confirms PID 386215 is
still running at 99.4% CPU. These samples landed outside the counter variables'
visible source scope, so I do not infer a new counter value from them. The FoF
catalogue/member output and log remain at their 00:20:16 timestamps and are
still partial.

## 00977 rank-8 counter sample (2026-10-21 01:02 KST)

The FoF job remains `RUNNING` at 4:02:27, and the downstream galaxy-finder,
extraction, and analysis jobs remain dependency-pending. A targeted GDB read in
rank 8's `pnew_fof_link()` frame gives `ncount=34,264,861`, `now=34,208,153`,
up from 30,819,585 and 27,144,394 at 01:00. The process remains at 99.4% CPU.
The output files still have their 00:20:16 KST timestamps and are incomplete.

## 00977 rank-8 caller counter sample (2026-10-21 01:03 KST)

FoF `430418_1` remains `RUNNING` at 4:03:33; downstream NewGalFinder,
extraction, and analysis tasks remain pending. The first GDB stop landed in
nested `pfof_open()`; reading its `pnew_fof_link()` caller gave
`ncount=37,208,466`, `now=37,195,622`, up from 34,264,861 and 34,208,153 at
01:02. PID 386215 is running at 99.4% CPU. The FoF files/log remain unchanged
at 00:20:16 KST and partial.

## 00977 FoF sample (2026-10-21 01:04 KST)

Slurm reports FoF `430418_1` `RUNNING` at 4:04:37. The corresponding
NewGalFinder/extraction/analysis tasks remain dependency-pending. Rank 8 is in
`pnew_fof_link()` with nested `pfof_open()` calls; caller counters advanced to
`ncount=39,053,116`, `now=38,578,788` from `37,208,466`, `37,195,622` at 01:03.
PID 386215 remains at 99.4% CPU. The FoF files and log remain unchanged at
00:20:16 KST; the outputs are still incomplete.

## 00977 FoF aggregation wait diagnosis (2026-10-21 01:04 KST)

The unchanged log ends after P0 announces that it is preparing to write
6,101,967 isolated halos. Ranks 1--7 and 9 have sent their data and logged
completion; rank 8 has not. GDB on rank 0 (PID 1762664, `syn01`) confirms it is
inside `WriteIsolatedHalo` at `PMPI_Probe(source=MPI_ANY_SOURCE, tag=987)`,
waiting for the remaining rank-8 message. Rank 0 is busy-spinning at 99.5% CPU;
rank 8 is still alive at 99.4% CPU in `pnew_fof_link()` and its counters
continue increasing. This explains why the large output files remain
unchanged: P0 has not received rank 8's isolated-halo data, so it cannot finish
writing them. There is no evidence of a dead rank or completed FoF output.

## 00977 rank-8 halo-index advance (2026-10-21 01:06 KST)

Slurm still reports FoF `430418_1` `RUNNING` at 4:06:45. A GDB sample now
shows rank 8 at `nhalo=8` in `pnew_fof_link()` (line 849), whereas prior
samples were at `nhalo=5`; its current-call counters are `ncount=1,220,366`,
`now=379,788`. These counters reset between halo calls, so they are not directly
comparable with the earlier `nhalo=5` values; the changed halo index is the
clear progress marker. Rank 0 remains in `PMPI_Probe` for tag 987 inside
`WriteIsolatedHalo`, while the log still lacks rank 8's send/completion lines.
The large FoF outputs remain unchanged at 00:20:16 KST.

## 00977 continuing halo-8 work (2026-10-21 01:07 KST)

Slurm still reports `430418_1` `RUNNING` at 4:07:49; dependent jobs have not
started. Rank 8 remains at `nhalo=8`, and the current-call counters increased
from `ncount=1,220,366`, `now=379,788` to `ncount=2,338,756`, `now=2,295,075`.
Rank 8 remains at 99.4% CPU. A fresh rank-0 backtrace again ends at
`PMPI_Probe(source=MPI_ANY_SOURCE, tag=987)`, so P0 is still waiting for rank
8; the FoF log/output timestamps remain unchanged at 00:20:16 KST.

## 00977 rank-8 traversal advance (2026-10-21 01:08 KST)

Slurm still shows FoF `430418_1` `RUNNING` at 4:08:54; its dependent tasks
remain pending. Rank 8 has advanced from `nhalo=8` to `nhalo=462226` in
`pnew_fof_link()`. The per-call counters at this halo are `ncount=1`, `now=1`,
so they reset between calls and do not indicate total work. Rank 0's fresh
backtrace still ends in `PMPI_Probe` for tag 987 inside the isolated-halo
writer, waiting for rank 8. The log and FoF output timestamps remain unchanged
at 00:20:16 KST; no completion/send marker is present for rank 8.

## 00977 rank-8 seed-loop advance (2026-10-21 01:08 KST)

Using the source corresponding to the loaded 22:22 FoF binary, I confirmed that
`nhalo` is the per-rank outer-loop count of previously unassigned seed
particles, incremented after each `pnew_fof_link()` call; the `ncount`/`now`
variables reset for each such call. Rank 8 has advanced from `nhalo=8` at the
previous check to `nhalo=462,226` and then `nhalo=992,206` in the current
samples. This shows that the extremely large early component has been passed
and the rank is progressing through later seeds. The total remaining seed
count is not yet known, so I make no completion-time estimate. P0 still waits
for rank 8's tag-987 isolated-halo message; outputs remain partial and
unchanged.

## 00977 rank-8 seed count continued (2026-10-21 01:11 KST)

The FoF task is still `RUNNING` at 4:11:11; downstream tasks remain
`PENDING (Dependency)`. Rank 8's outer-loop `nhalo` count has reached 1,665,681,
up from 992,206 at 01:08, confirming continued movement past the early large
component. At the current sample `ncount=1,244,056`, `now=864,326` for the seed
being linked. P0 remains inside `WriteIsolatedHalo` waiting for tag 987; its
log still lists rank 8 as the only worker not yet reported complete. FoF files
remain unchanged at 00:20:16 KST. No total seed bound is available, so I do not
estimate completion time.

## 00977 continued rank-8 large-seed processing (2026-10-21 01:12 KST)

FoF `430418_1` remains `RUNNING` at 4:12:13; all downstream tasks remain
pending on dependencies. Rank 8 is still on outer-loop seed `nhalo=1,665,681`,
with per-seed counters increased from `ncount=1,244,056`, `now=864,326` to
`ncount=2,564,462`, `now=1,450,867`. Rank 8 remains near 99.4% CPU. P0's
stack is still in `PMPI_Probe` for tag 987 while assembling the isolated-halo
catalogue. The Slurm stderr contains only the startup RAMSES-parameter
precision warning and expected missing gas/sink component notices; no MPI
error is reported there. Output file timestamps remain 00:20:16 KST.

## 00977 continued rank-8 seed (2026-10-21 01:13 KST)

The Slurm task remains `RUNNING` at 4:13:18, with all dependent tasks still
pending. Rank 8 remains on seed `nhalo=1,665,681`; its counters increased from
`ncount=2,564,462`, `now=1,450,867` to `ncount=3,087,543`, `now=1,624,256`.
Rank 8 is at 99.4% CPU. Rank 0's stack still shows `PMPI_Probe` for tag 987
inside `WriteIsolatedHalo`, while waiting for rank 8. The FoF output sizes and
log timestamp are unchanged at 00:20:16 KST.

## 00977 rank-8 seed progress (2026-10-21 01:14 KST)

The queue still shows FoF `430418_1` `RUNNING` at 4:14:22 and its dependent
jobs waiting. Rank 8 remains on `nhalo=1,665,681`; its per-seed counters
advanced to `ncount=3,654,628`, `now=1,780,749` from 3,087,543 and 1,624,256.
It remains at 99.3% CPU. Rank 0 is still probing for tag 987 in the isolated-
halo writer. The FoF output files remain unchanged at their 00:20:16 KST
partial-output timestamps.

## 00977 rank-8 seed counter update (2026-10-21 01:15 KST)

FoF `430418_1` remains `RUNNING` at 4:15:28; downstream tasks remain
`PENDING (Dependency)`. Rank 8 is still on seed `nhalo=1,665,681`, with its
per-seed counters advanced from `ncount=3,654,628`, `now=1,780,749` to
`ncount=4,418,517`, `now=2,175,106`; CPU use remains 99.3%. Rank 0's stack is
still in `PMPI_Probe` under `WriteIsolatedHalo`, waiting for tag 987. The FoF
catalogue, member file, and log remain unchanged at 00:20:16 KST.

## 00977 rank-8 FoF activity check (2026-10-21 01:19 KST)

FoF array task `430418_1` remains `RUNNING` at 4:17:47 elapsed; extraction and
analysis tasks remain dependency-pending. The worker process (PID 386215 on
`syn04`) is alive at about 99.3% CPU. A GDB sample still places it in the
`pnew_fof_link()` call for outer seed `nhalo=1,665,681`, but the per-seed
counters have increased to `ncount=7,049,995`, `now=6,928,438` from
`4,418,517`, `2,175,106` at 01:15. A separate sample caught it inside
`pfof_open()`, consistent with continued tree traversal. Rank 0 (PID 1762664
on `syn01`) was still polling for the rank-8 isolated-halo message at tag 987
in the 01:17 sample. The FoF log and output files remain unchanged at their
00:20:16 partial-write timestamps; the catalogue is not complete, and no
downstream task has started. No completion-time estimate is supported by the
available counters.

## 00977 FoF progress and aperture-sample correction (2026-10-21 01:22 KST)

Slurm still reports `430418_1` `RUNNING` at 4:22:04; the dependent NewGalFinder,
extraction, and relation-analysis tasks remain pending, while the replacement
FoF array tasks for snapshots 00972, 00980, and 00984 are pending scheduler
priority. A fresh GDB sample on `syn04` shows rank 8 alive at 99.3% CPU and
advancing its outer seed counter from `nhalo=1,665,681` at 01:19 to
`nhalo=2,313,784` at 01:22. The 00977 FoF log and output files remain at their
00:20:16 partial-write timestamps, so no catalogue is complete yet.

Before the dependent relation-analysis jobs start, I corrected
`analyse_satellite_relations.py`: the minimum count of valid stars inside the
three-dimensional half-mass radius now filters only the within-radius
 dispersion. It no longer removes galaxies from the global-dispersion or
number-weighted global-dispersion samples. This aligns the implementation
with `methods_draft.md`; the script was not executed in this session.

## 00977 rank-8 FoF progress check (2026-10-21 01:23 KST)

Slurm still reports `430418_1` `RUNNING` at 4:23:48; its dependent jobs remain
pending, and replacement FoF tasks for the other three snapshots remain pending
priority. Rank 8 (PID 386215 on `syn04`) is alive at 99.3% CPU. GDB now places
it in `pfof_open()` during `pnew_fof_link()` for outer seed `nhalo=2,967,736`,
up from 2,313,784 at 01:22. The FoF log and output files have not changed from
their 00:20:16 timestamps; no completion marker is present and the catalogue
is still partial.

## 00977 rank-8 FoF continues advancing (2026-10-21 01:25 KST)

At 4:25:47 elapsed, `430418_1` remains `RUNNING`; the matching downstream
jobs are still dependency-pending, and the replacement FoF tasks for 00972,
00980, and 00984 remain pending priority. Rank 8 on `syn04` is alive at 99.3%
CPU and has advanced to outer seed `nhalo=5,094,830`, from 2,967,736 at 01:23.
The FoF log and partial catalogue/member files still have 00:20:16 timestamps,
so the isolated-halo aggregation has not completed. No completion-time estimate
is supported.

## Host-clustered inference correction and 00977 progress (2026-10-21 01:27 KST)

I corrected the CR1 finite-sample factor in `analyse_satellite_relations.py` to
include host-specific intercepts in the residual degrees of freedom for the
host-fixed-effects models. The fit now reports coefficients/standard errors as
unmeasurable when there are not enough residual degrees of freedom or fewer
than two hosts. The secondary stellar-mass-bin fit now also accounts for its
host-by-mass-bin intercepts and does not report a clustered error from a single
host. These edits are in place before the dependent analysis jobs start; I did
not execute the script in this session.

FoF `430418_1` remains `RUNNING` at 4:27:23. Rank 8 on `syn04` is alive at
99.3% CPU and has advanced to `nhalo=8,597,490`, from 5,094,830 at 01:25.
The 00977 log and partial catalogue/member files remain unchanged at 00:20:16;
NewGalFinder and later stages are still dependency-pending. No completion-time
estimate is supported.

## Simulation metadata added to methods draft (2026-10-21 01:28 KST)

Read the four RAMSES `info_*.txt` files used by the extraction code and added
verified output metadata to `methods_draft.md`: snapshots 00972/00977/00980/00984
have redshifts 0.2226/0.2178/0.2132/0.2087, a span of 0.0139, and the same
comoving box size of 44.171261 h^-1 Mpc. The files report H0=67.4 km s^-1
Mpc^-1, Omega_m=0.315, Omega_Lambda=0.685, and Omega_b=0.045. Particle-mass
resolution and validated satellite counts remain to be supplied from the
catalogues; none are inferred from these metadata.

## 00977 FoF moved into tree construction (2026-10-21 01:29 KST)

The FoF log advanced at 01:28:29 and the output files grew: the halo catalogue
is now 16,240,708 bytes and the member file is 19,691,946,792 bytes. The log
shows ranks P0-P6 and P9 have passed tree building; P7 and P8 have not yet
reported that marker. A fresh GDB sample of rank 8 found it in recursive
`FoF_divide_node()` tree construction, rather than the earlier long
`pnew_fof_link()` traversal. Slurm still reports `430418_1` `RUNNING` at
4:29:32, and all downstream jobs remain dependency-pending. The files are still
not verified complete and must not yet be consumed as final catalogues.

## Literature framing and 00977 next FoF pass (2026-10-21 01:31 KST)

Updated `paper_outline.md` to acknowledge the direct metallicity--velocity
literature more precisely. The MNRAS article by Fattahi et al. (2018) explicitly
examines metallicity versus $V_{1/2}$ in Local Group dwarfs and discusses how
tidal stripping flattens the low-velocity relation; Spindler & Wake (2017)
compare central and satellite metallicities and core dispersions in SDSS. The
outline now limits any novelty claim to the conditional, internal-stellar-
dispersion test and its mass/host controls. Corrected the first author's name
in the Spindler & Wake BibTeX entry from Allison to Ashley.

The 00977 FoF log now shows all ten MPI ranks have passed tree building for
their current slab. Rank 8 is back in `pnew_fof_link()` with its per-slab halo
counter reset to zero, consistent with the next local linking pass. Slurm
still reports `430418_1` `RUNNING` at 4:31:21. The last logged output sizes are
from 01:28:29; the final FoF completion marker and downstream starts are still
absent. Literature sources checked: Spindler & Wake,
https://doi.org/10.1093/mnras/stx427; Fattahi et al.,
https://doi.org/10.1093/mnras/sty408.

## 00977 rank-local FoF linking in progress (2026-10-21 01:32 KST)

At 4:32:35 elapsed, `430418_1` remains `RUNNING`; NewGalFinder and the later
stages are still dependency-pending. The log now shows all ten ranks completed
tree construction for their current slabs. Rank 2 has finished its local link
pass for 14,624,800 particles, finding 5,367,904 local groups, and has reached
the halo-boundary check. Rank 8 remains alive at 99.3% CPU inside
`pnew_fof_link()` for its local pass. The last FoF output-file timestamps are
still 01:28:29, and no final catalogue marker has appeared.

## 00977 large-slab linking counters (2026-10-21 01:34 KST)

The FoF job was last confirmed `RUNNING` at 4:33:54 elapsed; its consumers
remain dependency-pending. Rank 8 is alive at 99.3% CPU in `pnew_fof_link()`
for its current slab (`nhalo=0` in this per-slab pass). A GDB sample reports
`ncount=2,718,389`, `now=1,131,583`; because no prior counter sample was taken
for this same call, this is a point-in-time value, not a measured rate. The
stdout remains at 01:32:31 and both FoF files remain at 01:28:29, with P2 at
its halo-boundary check and no final completion marker.

## 00977 link-pass counters advancing (2026-10-21 01:35 KST)

Slurm confirms `430418_1` `RUNNING` at 4:35:18; its extraction/analysis
consumers remain dependency-pending and the other FoF snapshots remain pending
priority. Rank 8 is still in the same per-slab call (`nhalo=0`) at 99.3% CPU,
with `ncount=4,557,224`, `now=4,553,106`, up from 2,718,389 and 1,131,583 at
01:34. The log has advanced: P3 also completed its local link pass (20,104,367
particles; 6,083,819 local groups) and reached the halo-boundary check. The
FoF output files remain at their 01:28:29 timestamps; the final catalogue is
not yet complete.

## EAGLE environmental-metallicity context and 00977 FoF check (2026-10-21 01:36 KST)

Added Bahé et al. (2017) to the literature framing in `paper_outline.md`: their
EAGLE analysis reports stellar-metallicity enhancement in satellites relative
to equal-stellar-mass field galaxies and examines environmental enrichment and
gas-removal processes. The outline now explicitly separates this known
environmental metallicity effect from the present test of a residual relation
with internal stellar dispersion. DOI: https://doi.org/10.1093/mnras/stw2329.

At 4:36:37 elapsed, FoF `430418_1` remained `RUNNING`; its downstream jobs
remain dependency-pending. Rank 8 is alive at 99.3% CPU. A GDB sample places it
back in `pnew_fof_link()` at local counter `nhalo=1`, up from `nhalo=0` at the
previous per-slab sample, indicating it returned from one local linking call.
The log and partial FoF files have not advanced since 01:35:05 and 01:28:29,
respectively. Ranks 2 and 3 have reached their halo-boundary checks; no final
FoF completion marker is present.

## Fattahi kinematic definition checked; 00977 linking continues (2026-10-21 01:38 KST)

Checked the Fattahi et al. (2018) MNRAS article's definition directly: their
$V_{1/2}$ is the circular velocity at the stellar 3D half-light radius
estimated from the luminosity-weighted line-of-sight dispersion, with
$V_{1/2}=\sqrt{3}\,\sigma_{\rm los}$. The distinction already stated in
`paper_outline.md` is supported by the primary article (Section 2.1,
https://academic.oup.com/mnras/article/476/3/3816/4875945).

FoF `430418_1` remains `RUNNING` at 4:37:52; all dependent stages remain
pending. Rank 8 is alive at 99.3% CPU in `pnew_fof_link()` at local counter
`nhalo=1`, with `ncount=5,293,194`, `now=4,127,886`. The stdout and partial FoF
files still have 01:35:05, 01:28:29 timestamps, respectively; completion has
not been reached.

## Working introduction drafted (2026-10-21 01:38 KST)

Created `introduction_draft.md` using the verified literature framing and the
implemented measurement definitions. The draft distinguishes the prior SDSS
core-dispersion and Local Group $V_{1/2}$ studies from the simulation's
mass-weighted internal stellar dispersions, states the conditional test and
metallicity/subhalo-mass limitations, and makes no claim about uncomputed
results or causal mechanisms. It is a working draft pending the final audited
catalogues and results.

## Introduction draft and 00977 status (2026-10-21 01:39 KST)

`introduction_draft.md` now contains four working paragraphs: environmental
metallicity and tidal effects, the direct Spindler/Fattahi precedents and their
kinematic definitions, this study's controlled question, and the limits on
causal interpretation. It contains no unmeasured sample counts or result
claims and remains subject to revision after the audited results are known.

Slurm still reports `430418_1` `RUNNING` at 4:39:07; its NewGalFinder,
extraction, and analysis consumers remain dependency-pending, while the other
FoF outputs wait on scheduler priority. Rank 8 (PID 386215 on `syn04`) is alive
at 99.3% CPU; GDB sampled it inside `pfof_open()` in the local linking pass.
The log and partial FoF files have not advanced since 01:35:05 and 01:28:29,
respectively. There is still no completion marker.

## 00977 FoF remains active (2026-10-21 01:40 KST)

Rechecked the live rank directly: PID 386215 on `syn04` is still running
`opfof.exe 977 100` at 99.3% CPU, with 4:40:16 elapsed. Its recent stdout
still ends after ranks P2 and P3 reached halo-boundary checks; no completion
marker or final catalogue is present. Slurm confirms the job uses 60 CPU cores
on three nodes and requests no GPU. Downstream NewGalFinder, extraction, and
analysis remain dependency-pending, so no satellite catalogue can yet be
validated for snapshot 00977.

## Four-snapshot catalogue outputs still pending (2026-10-21 01:41 KST)

Rechecked the four requested snapshots. None yet has a completed
`GALFIND.DATA` file or the corresponding extracted satellite-properties or
relation catalogue. The 00977 FoF process remains alive on `syn04` at 99.3%
CPU (4:41 elapsed), with no new output beyond the halo-boundary checks for
ranks P2 and P3. Slurm still shows NewGalFinder, extraction, and analysis
waiting on dependencies; replacement FoF jobs for 00972, 00980, and 00984
remain pending scheduler priority.

## 00977 active FoF rank sampled in particle-link search (2026-10-21 01:43 KST)

Attached GDB with the running executable's symbols and sampled PID 386215 on
`syn04`. The rank is executing `new_fof_link()` at the particle-pair distance
test in `Treewalk.fof.ordered.c` (line 919), consistent with a long local
linking pass. Optimized local counter values from this sample are not reliable
and are not treated as progress measurements. No FoF completion marker,
final halo catalogue, or NewGalFinder output is available yet.

## 00977 FoF rank boundary checks advancing (2026-10-21 01:44 KST)

The 00977 FoF process remains `RUNNING` at 4:44 elapsed, and PID 386215 on
`syn04` remains active at 99.3% CPU. Unlike the previous check, the job log
has advanced: ranks P5 and P1 have now also completed their local linking
passes and halo-boundary checks, joining P2 and P3. The 19-GiB member file and
halo catalogue still have their 01:28 timestamps; these are partial outputs,
not a completed FoF catalogue. NewGalFinder and downstream extraction and
analysis remain dependency-pending. FoF jobs for the other requested outputs
are still pending scheduler priority.

## All 00977 FoF MPI ranks remain CPU-active (2026-10-21 01:45 KST)

Checked the process list directly on all three allocated nodes. All ten
`opfof.exe 977 100` MPI ranks are alive and each is using about 99% CPU at
4:45 elapsed. Thus, the lack of additional per-rank completion messages does
not indicate that the whole job is stalled; the ranks are still doing CPU
work. The most recent log advance remains 01:43:55, when P1 completed its
boundary check. Final files and downstream galaxy catalogues are still
pending.

## 00977 FoF additional rank completed (2026-10-21 01:46 KST)

The Slurm job remains `RUNNING` at 4:46 elapsed. Its log advanced at 01:46:03:
rank P6 completed its local linking pass and halo-boundary check. Ranks P2,
P3, P5, P1, and P6 have now reported boundary checks; the remaining ranks
have not yet reported completion. The existing FoF files are still the
partial 01:28 outputs. None of the four snapshots yet has a NewGalFinder
`GALFIND.DATA` file or extracted satellite-properties or relation catalogue.

## 00977 FoF master rank and isolated-halo writes advanced (2026-10-21 01:48 KST)

Slurm still reports `430418_1` running at 4:48 elapsed; its dependent
NewGalFinder and later jobs remain pending. The log has advanced substantially:
rank P0 completed its pass for 31,638,988 particles, finding 5,754,494 local
groups, and reached the halo-boundary check. P0 then began writing isolated
halos; P5, P2, P3, P1, and P6 reported sending their isolated-halo data to P0
and completing their writes. The partial FoF files grew to 19.9 MB and 26.8
GB. These are still intermediate outputs: no final FoF completion marker,
`GALFIND.DATA`, or satellite-property catalogue exists yet.

## 00977 FoF master rank is writing; storage headroom is ample (2026-10-21 01:49 KST)

The job remains `RUNNING` at 4:49 elapsed. The log has not advanced past P0's
isolated-halo write phase since 01:48:37, and the FoF catalogue files remain
at 19.9 MB and 26.8 GB. However, the live P0 process on `syn01` is still at
99.5% CPU and its `/proc` I/O counters show 35.8 GB written through write
system calls, consistent with active output work. The `/gpfs` filesystem has
237 TB available; the user's scratch quota reports no block or file limit.
This rules out a current space/quota shortage but does not establish FoF
completion. Downstream NewGalFinder remains dependency-pending.

## 00977 FoF active CPU, output counters unchanged in latest interval (2026-10-21 01:50 KST)

Slurm still reports `430418_1` `RUNNING` at 4:50 elapsed; P0 on `syn01` is
alive at 99.5% CPU. Since the 01:49:15 process sample, its `wchar` counter has
remained at 35.8 GB, and the log and two FoF file sizes have not changed since
01:48:37. The current evidence is CPU activity without recent file-output
progress; it is not enough to conclude completion or failure. The 00977
NewGalFinder and satellite-property jobs remain dependency-pending.

## 00977 slow-rank bottleneck localized (2026-10-21 01:51 KST)

A fresh symbolized GDB sample shows rank 8 (PID 386215 on `syn04`) still in
`new_fof_link()` at the particle-pair distance check (source line 919). Rank 0
(PID 1762664 on `syn01`) is CPU-active, holds the FoF output files open at the
last written offsets, and its program counter maps into UCX's `libuct.so`.
In `WriteIsolatedHalo`, rank 0 receives worker messages with `MPI_Probe` after
its local halo pass. Together with the log showing that P4, P7, P8, and P9 have
not reported their isolated-halo sends, this is consistent with rank 0
waiting for the remaining workers while rank 8 continues local linking. The
stack did not resolve to a named MPI routine, so this remains a source-guided
inference. Disk space is ample; no output or completion marker has changed.

## 00977 remains in the same large local linking phase (2026-10-21 01:53 KST)

At 4:53 elapsed, Slurm still reports the FoF task `RUNNING`. Rank 8 remains
alive at 99.3% CPU and a fresh symbolized sample still places it inside
`new_fof_link()` at the particle-pair distance test. Rank 0 remains alive at
99.5% CPU in UCX communication progress; its output file offsets and write
counters have not changed since the 01:49 sample. The FoF log and file sizes
are unchanged since 01:48:37. This is consistent with the master waiting for
slow worker ranks, but no named MPI routine or completion has been observed.
All downstream catalogue jobs remain dependency-pending.

## 00977 FoF rank 4 completed; master output resumed (2026-10-21 01:56 KST)

At 4:56 elapsed, Slurm still reports `430418_1` `RUNNING`; all satellite
catalogue jobs remain dependency-pending. The log advanced at 01:55:27: rank P4
completed its pass for 42,336,870 particles (6,831,223 local groups), passed
the halo-boundary check, and sent its isolated-halo data to P0. The partial
FoF halo and member files grew to 20.8 MB and 28.4 GB. P0's live `wchar`
counter increased to 37.4 GB, confirming renewed writes. Ranks P7, P8, and P9
have not yet reported their isolated-halo sends; the snapshot is not complete.

## 00977 FoF file-shard progression clarified (2026-10-21 01:56 KST)

Clarification to the preceding rank updates: `opfof.c` loops over 100 input
files, assigning ten consecutive file shards to each MPI rank. A message such
as “P4 has checked halo boundary” / “P4 have sent all isolated halo data”
marks completion of that rank's current file shard, not completion of the
rank's full ten-file range. The most recent visible pass is file suffix `3`
(e.g. P4 reads `SN.00977.00043`): P0–P6 have completed and sent this shard;
P7, P8, and P9 have not yet reported completion for suffix `3`. P8's current
shard contains 136,713,212 particles, P7's 87,653,295, and P9's 68,881,593.
Thus several input shards remain for every rank, including the heaviest
current shard. The FoF catalogue remains partial and NewGalFinder has not
started. Earlier references to a rank being “complete” should be read as
completion of that file shard only.

## 00977 FoF still active without catalogue writes (2026-10-21 02:04 KST)

At 5:04 elapsed, Slurm still reports `430418_1` as `RUNNING`; the dependent
NewGalFinder, extraction, and analysis tasks remain pending. A direct process
check on the allocated nodes found all opFoF ranks alive. Rank 8 remains at
99.3% CPU, while rank 0 remains at 99.5% CPU. However, the log and both FoF
output files have not changed since 01:55:27, and rank 0's cumulative
`wchar` remains 37.35 GB, so no further catalogue output has been written in
that interval. The run is still computing, but there is no evidence yet of
completion or advancement to the downstream catalogue stages.

## 00977 rank-8 computation localized to particle-linking traversal (2026-10-21 02:06 KST)

A targeted stack sample of rank 8 (PID 386215 on `syn04`) resolves to
`new_fof_link()` at `Treewalk.fof.ordered.c:919`, the particle-coordinate
separation calculation during FoF linking. Its most recent input in the log is
shard `SN.00977.00083` (136,713,212 particles), and it has not logged that
shard's completion. The sample therefore indicates active work in the heavy
shard, rather than an MPI wait; it does not establish how much of the shard
remains. Rank 0's output counter and all output-file sizes remain unchanged,
and downstream jobs remain dependency-pending.

## 00977 remains in heavy linking; other outputs queued (2026-10-21 02:09 KST)

A fresh scheduler check still shows FoF task `430418_1` (`z=0.2178`, output
00977) `RUNNING` at 5:08 elapsed. Its dependent galaxy-finder, extraction, and
analysis tasks remain pending on dependencies. Rank 8 remains at 99.3% CPU in
the same shard/linking phase; its `/proc` character-read counter is unchanged
from 02:06, while the process CPU time continues to advance. No FoF log or
catalogue file has advanced since 01:55. The other current FoF array tasks for
outputs 00972, 00980, and 00984 are pending for priority. The earlier 00972
attempt `430417_0` was cancelled and left no FoF catalogue, so the queued task
for that output is still required.

## Correction and 00977 FoF shard progress (2026-10-21 02:10 KST)

Correction to the 02:09 note: the log and FoF partial files did advance at
02:09:10. Rank P9 completed and sent its isolated-halo data for input shard
`SN.00977.00093` (68,881,593 particles). The halo catalogue grew to 21,966,916
bytes and the member-particle file to 31,466,513,192 bytes. This is one input
shard only; the snapshot-level FoF job remains `RUNNING`, and no GALFIND or
satellite-properties catalogue exists yet. Rank 8's process remains CPU-active
in the previously identified `.00083` linking phase; ranks P7 and P8 have not
yet logged completion of their current shards. The 01:55-to-02:09 interval had
no output growth, but it is no longer accurate to say that the files have not
advanced since 01:55.

## 00977 heavy rank still compute-bound (2026-10-21 02:16 KST)

At 5:15 elapsed, Slurm still reports `430418_1` as `RUNNING`; its dependent
catalogue jobs remain pending. Ranks 7, 8, and 9 on `syn04` are all alive and
using about 99% CPU each. Rank 8's `/proc` read counter is unchanged at
34.19 GB from the 02:06 and 02:14 samples, while its CPU time advances; the
last logged shard for this rank remains `.00083`. The FoF log and output-file
sizes have not changed since 02:09:10. This confirms a prolonged CPU-bound
linking calculation but gives no reliable fraction-complete estimate.

## Rank-8 stack changed to FoF tree construction (2026-10-21 02:17 KST)

A new stack sample for rank 8 resolves to `FoF_divide_node()` at
`Treewalk.fof.ordered.c:964`, whereas the 02:06 sample was in the particle
separation test inside `new_fof_link()`. The process has therefore moved to a
different FoF stage; the stale stdout does not establish which input shard it
is currently handling, since the reader/linker may have buffered messages. Its
read counter remains unchanged and CPU time advances, consistent with
CPU-bound work on already-read particles. The log's latest visible P8 shard
remains `.00083`; this is not proof that `.00083` is still the active shard.
No new snapshot-level output or completion marker is present.

## Rank-8 sample inside tree-node particle aggregation (2026-10-21 02:18 KST)

A symbolized sample placed rank 8 in `FoF_divide_node()` at
`Treewalk.fof.ordered.c:998`, immediately after a particle-list pass that
counts particles and accumulates their coordinates for a tree node. This
confirms ongoing CPU work in tree construction. The sample does not identify
the active input shard or a completion fraction; the log and FoF output sizes
remain unchanged, and the snapshot job is still `RUNNING`.

## Rank 8 returned from tree construction to FoF linking (2026-10-21 02:21 KST)

A later symbolized sample places rank 8 back in `new_fof_link()` at
`Treewalk.fof.ordered.c:919`, with `nhalo=1`; the 02:18 sample was inside
`FoF_divide_node()` tree construction. This confirms that the live process has
moved between the tree-building and particle-linking stages despite unchanged
stdout and catalogue sizes. The exact shard and amount of remaining work are
still unknown; the snapshot task remains `RUNNING` and downstream jobs remain
dependency-pending.

## Rank-8 remains in the same FoF-link call across samples (2026-10-21 02:25 KST)

Symbolized samples at 02:21 and 02:24 both placed rank 8 in
`new_fof_link()` at the particle-separation calculation, with the call argument
`nhalo=1`. The rank remains near 99% CPU and its read counter has not changed;
there is no log or catalogue growth after 02:09. This suggests one prolonged
linking call, but the optimized local counters are not reliable enough to
measure its progress or conclude that it is stuck. Slurm still reports the
snapshot task `RUNNING`; downstream catalogue jobs remain dependency-pending.

## Same rank-8 FoF-link invocation persists (2026-10-21 02:28 KST)

The 02:28 stack sample again places rank 8 in `new_fof_link()` at line 919,
with the same `nhalo=1` argument and the same tree, particle-list, and linked-
array addresses as the 02:24 sample. This is consistent with the same
long-running link call; it still does not expose a trustworthy within-call
counter. The rank remains CPU-active, but no FoF log or output growth has
occurred since 02:09. At 5:28 elapsed, Slurm still reports `430418_1`
`RUNNING`, with its downstream stages dependency-pending.

## Rank-8 stack returned to FoF tree construction (2026-10-21 02:32 KST)

The 02:31 symbolized sample places rank 8 in `FoF_divide_node()` at
`Treewalk.fof.ordered.c:976`, after the 02:28 sample in `new_fof_link()`.
This is another observed transition between tree construction and linking,
consistent with continued FoF work. The current shard remains unknown because
stdout and output files are stale. Slurm reports `430418_1` `RUNNING` at 5:31
elapsed; the downstream catalogue jobs remain dependency-pending.

## Correction: repeated rank-8 pointers do not identify one persistent shard

A source check of `opfof.c` shows that each rank loops over multiple `nowfile`
shards and resets `nhalo=0` for each shard before FoF linking. The rank reuses
its allocated buffers across that loop. Therefore the repeated rank-8 values
(`nhalo=1` and pointer addresses) in the 02:24/02:28 samples do not prove that
the same invocation or input shard persisted. Stack samples show only that the
rank was in linking at those instants; the later tree-construction sample shows
a phase transition. The current shard cannot be identified from the available
logs or these pointer values.

## 00977 FoF and dependent pipeline status (2026-10-21 02:35 KST)

A fresh Slurm query reports FoF array task `430418_1` (`SNAP=977`) still
`RUNNING` after 5:35:12, on `syn[01,03-04]`. Its stdout and the two partial
FoF catalogue files are unchanged since 02:09:10; the log still ends with P9
completing shard `.00093`, which is not snapshot completion. A one-time check
on `syn04` found ranks 7–9 alive and each using about 99% CPU, with no physical
file reads or writes reported by `/proc/<pid>/io` at the sample. This supports
continued CPU computation but gives no reliable remaining-time estimate or
active-shard identity. No snapshot-level completion marker or GALFIND output
exists for 00977. NewGalFinder tasks `430423_1/_2`, extraction tasks
`430427_1/_2`, and relation-analysis tasks `430431_1/_2` remain pending on
dependencies. FoF replacements for 00972, 00980, and 00984 remain pending for
priority. Do not interpret these partial 00977 files as valid input to the
next stage.

## Global stellar-ID audit added to the extractor (2026-10-21 02:42 KST)

Static review found that the original duplicate-ID check covered galaxies only
within the same FoF host. `extract_satellite_properties.py` now retains the
unique stellar IDs from each host and checks them against the full output
before publishing the catalogue. Within-host duplicate galaxies remain
flagged for exclusion by the relation analysis; a repeated ID across different
FoF hosts now makes extraction fail and removes the temporary CSV, preventing
an unverified catalogue from being treated as complete. This code path has
not yet run on a NewGalFinder output. At 02:42 KST, the 00977 FoF task remained
`RUNNING`, its log and partial catalogue files were unchanged since 02:09, and
the dependent extraction and analysis tasks remained pending.

## Four-output catalogue chain audit (2026-10-21 02:48 KST)

A combined `squeue`/`sacct` check confirms that all four outputs still have a
scheduled path through the analysis stage. Snapshot 00977 FoF task `430418_1`
is `RUNNING`; its NewGalFinder, extraction, and analysis tasks are dependency-
pending. The replacement FoF jobs for 00972 (`430610_0`), 00980 (`430611_2`),
and 00984 (`430622_3`) are pending priority, with their downstream
NewGalFinder, extraction, and analysis tasks pending on dependencies. The
cancelled 00972 extraction/analysis array entries (`430426_0`, `430430_0`) are
superseded by the explicit replacement chain `430491_0` -> `430492_0` ->
`430493_0`; the 00984 chain is `430476_3` -> `430477_3` -> `430479_3`. No
snapshot-level satellite catalogue exists yet, and the 00977 FoF log and
partial output files remain unchanged since 02:09.

## 00977 FoF rank diagnostic (2026-10-21 02:51 KST)

The last Slurm poll at 02:49 KST reported `430418_1` running at 5:49 elapsed.
A targeted process sample immediately afterward on `syn04` found ranks 7--9
in runnable state at about 99 per cent CPU;
their character-read/write counters were unchanged from earlier samples. A
single symbolized GDB backtrace for rank 8 placed it at entry to
`FoF_divide_node()` in `Treewalk.fof.ordered.c:964`, a tree-construction
routine. This confirms the rank is spending CPU in FoF tree work, but does not
identify its input shard, tree size, or remaining work and is not proof of
forward progress. The 00977 log and partial catalogue files remain unchanged
since 02:09; no downstream input is ready.

## Current queue has no usable start estimate (2026-10-21 02:54 KST)

A fresh `squeue --start` query returns `N/A` for the priority-pending FoF
replacement tasks (`430610_0`, `430611_2`, `430622_3`) and for their downstream
dependency tasks. There is therefore no current scheduler ETA to report; the
older projected start time recorded on 2026-10-20 is stale. At the preceding
02:53 KST status poll, 00977 FoF `430418_1` remained running and its log and
partial catalogue files were still unchanged since 02:09.

## 00977 FoF rank stage sample (2026-10-21 02:55 KST)

Slurm reports `430418_1` running at 5:55:54 elapsed. A symbolized sample of
rank 8 on `syn04` placed it in `new_fof_link()` at the particle-distance
comparison (`Treewalk.fof.ordered.c:919`), with local `nhalo=1`. Source review
shows `nhalo` is reset for each `nowfile` shard, so this value does not
identify the shard or quantify total snapshot progress. The changed stack
location documents the sampled stage only; it is not proof of forward
progress. The stdout and partial catalogue files remain unchanged since 02:09.
The latest `squeue --start` query gave no ETA for priority-pending replacement
FoF jobs or their dependent stages.

## Repeated 00977 rank-8 linking sample (2026-10-21 02:57 KST)

A second symbolized sample about two minutes after the 02:55 sample again
placed rank 8 at `new_fof_link()` line 919 with the same visible tree, particle,
and linked-buffer addresses and `nhalo=1`; its CPU time advanced by about 2.5
minutes and it remains near 99 per cent CPU. The local `ncount`/`now` values
reported by GDB are not trustworthy in this optimized/corrupt-unwind frame.
Because `opfof.c` reuses buffers and resets `nhalo` for every `nowfile`, these
samples cannot prove one invocation or shard remained active. They also do not
provide a reliable progress fraction. The Slurm job remains `RUNNING`; stdout
and partial output files have not changed since 02:09.

## 00977 FoF status and shard-loop audit (2026-10-21 03:02 KST)

A fresh queue check reports FoF task `430418_1` still `RUNNING` at 6:00:45;
replacement FoF tasks for 00972, 00980, and 00984 remain pending with scheduler
reason `Priority`, and all downstream tasks remain dependency-pending. The
00977 stdout and partial catalogue files are still unchanged since 02:09:10;
no `satellite_properties_*.csv` or relation result for the new NewGalFinder
chain exists. Reading `opfof.c` confirms that the job loops through each rank's
assigned shard range (`nfile/nid`); with 100 files and 10 MPI ranks, each rank
handles ten shards. The log ending with rank 9 completing shard `.00093`
therefore does not mark snapshot completion. Local FoF halo files remain
partial and must not be passed downstream. No `--gres=gpu` request appears in
the submitted FoF script, and Slurm currently labels the three replacement
FoF jobs as `Priority` pending; this check does not implicate a GPU quota.

## NewGalFinder binary-layout and unit audit (2026-10-21)

Read-only inspection of the exact executable's source tree at
`/gpfs/kjhan/Others/Jaehyun/2026.Apr/NewGalFinder` confirms that the extractor's
168-byte `StarType` layout matches `ramses.h` and the active Makefile flags:
`dptype=double`, default 32-bit `idtype`, no `OUTPUT_PARTICLE_POTENTIAL`, no
`NBODY`, and `NCHEM=9`. The source documents particle positions in
`cMpc/h`, velocities in `km/s`, current and initial stellar masses in
`Msun/h`, and `tp,zp` in RAMSES units. This validates the extractor's field
offsets and the analysis conversion from `Msun/h` to `Msun`; it does not yet
validate the parsed output values, which still requires a completed
NewGalFinder catalogue and its numerical consistency checks.

## Match the stellar-particle cut to the dispersion weights (2026-10-21)

Static review found that the 100-star selection counted any finite velocity,
while the primary velocity dispersion uses only stars with both finite
velocities and positive current masses. The extractor now records this
mass-valid velocity count separately as `nstar_velocity_mass_valid`, and the
analysis applies the 100-particle cut to that field. The existing
`nstar_velocity_valid` still records the broader finite-velocity count. The
Methods draft now states the actual primary and half-mass-aperture sample
criteria. This change has not yet been exercised on a NewGalFinder output; no
tests were run.

## 00977 FoF CPU-activity diagnostic (2026-10-21 03:08 KST)

`430418_1` remained `RUNNING` at 6:07 elapsed. A single targeted process
sample on `syn04` found three `opfof.exe 977 100` MPI ranks, each runnable and
using 99.3--99.5 per cent CPU; their process CPU times were about 6:05--6:06.
This rules out an idle or dependency wait on those ranks, but does not quantify
FoF progress. The log still ends after the tree/linking stage for the records
read from shards `.00003`, `.00013`, ..., `.00093`, and the two FoF files still
have their 02:09 timestamps. No complete-output marker or downstream catalogue
is present. Keep the job running; no cancellation or resubmission is warranted
by this sample alone.

## Introduction revised while catalogue jobs run (2026-10-21)

Condensed `introduction_draft.md` into four paragraphs: physical motivation,
prior measurements, the two conditional models, and scope/limitations. Kept
explicit that the analysis uses internal stellar-particle dispersion, that
present-day finder-assigned dark-matter mass is not verified bound mass, and
that attenuation after this mass control cannot exclude tidal effects. No
simulation result is asserted before the new catalogues are validated.

## Correct opFoF CPU request for future submissions (2026-10-21)

The active FoF executable is built from MPI-only C sources: a targeted source
search found no OpenMP directives or runtime calls, and `Rules.make` does not
enable OpenMP. Yet `opfof_array.sbatch` requested six CPUs per MPI rank and set
OpenMP environment variables. Updated the reusable script to request one CPU
per rank and removed the unused OpenMP settings. This affects future
submissions only; Slurm has already captured the resource requests of the
currently running 00977 job and the queued jobs. No active job was cancelled
or resubmitted.

## Clarify uncertainty estimation in Methods (2026-10-21)

Reworded the Methods description to distinguish the host-clustered standard
errors from the percentile confidence intervals obtained by resampling whole
hosts 2,000 times. This matches the two separate uncertainty calculations in
`analyse_satellite_relations.py`.

## Align planned diagnostics with available NewGalFinder fields (2026-10-21)

The primary NewGalFinder adapter explicitly sets recent SFR and stellar age to
missing values. Removed activity-class comparisons from the planned results
and figure description unless those values are later derived and validated;
the Methods draft now records that limitation. Stellar size remains available
as a secondary diagnostic. This prevents the planned paper from promising
measurements that the current primary catalogue does not contain.

## Clarify unavailable diagnostics as not yet derived (2026-10-21)

Clarified the Methods and outline wording: the NewGalFinder stellar records
contain birth-time information, but the current analysis has not derived
recent SFRs or stellar ages from it. Activity and age comparisons remain out
of the planned primary results until their time-unit conversion and derived
measurements are validated.

## Birth-time field provenance check (2026-10-21)

Read-only source inspection confirms NewDD reads each stellar `tp` field from
RAMSES `part` files and carries it unchanged into the split particle records;
it scales mass, position, velocity, and `mass0`, but not `tp`. The RAMSES
particle documentation identifies `tp` as the birth epoch, and the yt RAMSES
documentation says cosmological formation times are stored as conformal times
([particle fields](https://ramses-organisation.readthedocs.io/en/latest/dev_docs/lecture/particles.html),
[yt RAMSES fields](https://yt-project.org/doc/examining/loading_data.html)).
This makes a future age or recent-SFR diagnostic possible in principle, but
requires the cosmological conformal-to-proper-time conversion to be validated
before use; no such derived quantity is included in the current catalogue.

## RAMSES birth-time conversion source check (2026-10-21)

The official RAMSES `part2birth` utility clarifies the required conversion:
for a cosmological output it builds conformal-time and proper-time tables from
the output cosmology, interpolates each stellar birth field `tp` in conformal
time, and reports the elapsed proper age between birth and the output in Gyr
([official RAMSES source](https://raw.githubusercontent.com/ramses-organisation/ramses/dev/utils/f90/part2birth.f90)).
The RAMSES documentation separately notes that non-RHD cosmological outputs
store stellar birth times in conformal units
([RAMSES RHD parameters](https://ramses-organisation.readthedocs.io/en/latest/wiki/RHD_params.html)).
This confirms the conversion method and rules out treating `tp` as a proper
age or multiplying it directly by a time unit. We have not yet applied the
conversion to the simulation particles or compared the result with
`part2birth`/yt on a validation sample; stellar ages and recent SFRs therefore
remain absent from the primary catalogue and analysis.

## 00977 FoF resource and progress check (2026-10-21 03:22 KST)

Slurm still reports `430418_1` as `RUNNING` at 6:23:49. The active step's
cumulative CPU time increased to 21:12:47 and its maximum resident memory was
about 89.7 GB, so the computation remains CPU-active despite no log or FoF
catalogue growth since 02:09. The allocation requests 60 CPUs across three
nodes and has no GPU resource in `ReqTRES` or `AllocTRES`; `a10` is the
partition name, not evidence of a GPU allocation. Do not interrupt or
resubmit this task based on the unchanged output files alone.

The FoF source shows that each slab writes isolated halo data, synchronizes
all MPI ranks, then traverses boundary-linked particles to append its
`BottomFaceContactHalo` file before stacking particles for the next slab.
Those auxiliary files are already 1.4--8.9 GB each. Their latest modification
times are 01:28, while the main FoF outputs last changed at 02:09 and the log
ends after rank 9 saved isolated halo data. This identifies a potentially
expensive stage in the slab pipeline, but does not prove which rank or slab is
currently active or provide a completion fraction.

## 00977 FoF and GPU-request audit (2026-10-21 03:38 KST)

The same FoF task `430418_1` remains `RUNNING` (elapsed 6:37:35). Slurm maps
the array task to job `430419`; its active MPI step reports cumulative CPU
time 21:59:21, maximum RSS about 86.8 GiB, and cumulative disk read/write of
about 71.5/44.5 GiB. The log and both main FoF products remain unchanged since
02:09, and the final `GALFIND.DATA.00977` catalogue has not been created. The
CPU and I/O counters show continued work, but do not establish a completion
fraction or identify its current MPI stage. Leave the running task untouched.

Resource records for this FoF job show `ReqTRES=cpu=60,mem=768G,node=3` and
`AllocTRES=cpu=60,mem=768G,node=3`, with no GPU/GRES field. The queued galaxy-
finder jobs also request CPU and memory only; their current pending reasons
are dependencies, while the replacement FoF jobs are pending for priority.
Their use of partitions named `a10` or `a100` does not constitute a GPU
request. Keep these expensive stages under Slurm rather than running them on
the shared login node.

A follow-up sample at 03:40 KST showed cumulative MPI-step CPU time rising to
22:07:43 while the log and FoF products remained unchanged. This confirms
that the job is still doing work, although its output has not yet advanced to
the next visible checkpoint.

At 03:41 KST, the output advanced: rank P7 completed another slab and saved
its isolated halos. The halo catalogue grew from 21,966,916 to 22,769,732
bytes and the member-particle file from 31,466,513,192 to 32,654,885,672
bytes; both were updated at 03:40:55. FoF is therefore progressing through
slabs. `GALFIND.DATA.00977` is still absent, so the dependent galaxy-finding
and analysis stages remain correctly queued.

At 03:44 KST, Slurm still reported `430418_1` as `RUNNING`. The MPI step's
cumulative CPU time had risen to 22:19:18, but the FoF catalogue and member
file remained at their 03:40:55 sizes and the auxiliary boundary files had
not changed. The log still ends with P7 saving its isolated halos. The job is
active, but there is not yet a newer visible slab-output checkpoint.

The 03:46 KST sample still showed the task `RUNNING` at 6:46:21 and cumulative
MPI-step CPU time 22:27:34. No FoF output or log lines had advanced since
03:40:55; the last logged state therefore remains P8 past tree construction
without a completed halo-count message.

At 03:47 KST, the task was still `RUNNING` and cumulative MPI-step CPU time
had reached 22:30:05. A targeted check of all four requested outputs found
only the partial 00977 FoF files; 00972, 00980, and 00984 have no main FoF
files yet, and none of the four has a `GALFIND.DATA` or extracted satellite
catalogue. All galaxy-finding and analysis jobs remain dependency-pending.

At 03:52 KST, `430418_1` remained `RUNNING` (elapsed 6:52:48). Its cumulative
MPI-step CPU time rose to 22:48:58, while the log and FoF file sizes were
unchanged from the 03:40:55 checkpoint. The most recent log line remains P7
after saving isolated halos; P8 has not yet reported its halo count for this
slab. The rising CPU counter confirms that the task remains active, but no
new output checkpoint is available yet.

Inspection of the current slab's ordered log and `opfof.c` narrows the
unfinished rank: P8 logged `passed Tree building` for 136,713,212 particles,
but has not yet logged its halo count. The code next traverses the particle
array with the FoF-linking routine before printing that count; all other ranks
have already reached the isolated-halo synchronization. This points to P8's
large FoF-linking pass as the current bottleneck, rather than a stalled output
write. The inference is from the ordered log and source sequence; no live
per-rank stack trace was taken.

## Recent tidal mass-loss context added to Introduction (2026-10-21)

Verified the 2026 MNRAS Auriga Streams III paper (Riley et al.,
doi:10.1093/mnras/stag029). It finds that a tight observed Local Group
satellite mass--metallicity relation can coexist with substantial tidal
stellar-mass loss. Added this limited context to `introduction_draft.md`,
`paper_outline.md`, and `references.bib`; it does not claim that this paper
tests the same conditional metallicity--internal-dispersion relation.

## Primary-results figure pipeline (2026-10-21)

Added `make_satellite_relation_figure.py` to plot the percent change in
internal dispersion per 0.3-dex metallicity-proxy difference, with whole-host
bootstrap intervals, comparing the host+stellar-mass model with the added
finder-assigned dark-matter-mass model. The four outputs are categorical
replication checks, not a redshift track. Slurm job `431161` is queued with
`afterok` dependencies on all four relation-analysis array tasks; it requests
one CPU and no GPU. It will run after upstream FoF, galaxy-finding, extraction,
and analysis stages complete successfully.

Clarified the pending figure's labels in `make_satellite_relation_figure.py`:
the all-member-star dispersion is no longer labeled only as “global,” the
y-axis names the internal stellar dispersion and metallicity-proxy interval,
and the added control is explicitly called the finder-assigned subhalo dark-
matter mass.

## Galaxy-finder CPU allocation audit (2026-10-21)

`gfind_array.sbatch` requests five MPI ranks, 20 CPUs per rank, and 600 GB of
memory per node across two nodes. `gfind_wrapper.sh` runs one thread on rank 0
and 20 OpenMP threads on each of ranks 1--4, for 81 configured CPU threads in
total; it pins each worker to CPU cores. Neither script requests a GPU/GRES.
The queued 00977 galaxy-finder job is pending on its FoF dependency. Its use
of the `a100` CPU partition does not by itself mean a GPU was allocated.

## Four-output catalogue pipeline status (2026-10-21 04:00 KST)

The live FoF task `430418_1` for snapshot 00977 remains `RUNNING` at 7:00:18;
its MPI step's cumulative CPU time has risen to 23:11:17. The main FoF files
and log have not advanced beyond the 03:40:55 checkpoint, whose last entries
show P7 completing its isolated-halo save. This confirms continued CPU work,
but there is no completed 00977 FoF catalogue yet.

An exact-path audit found no `GALFIND.DATA`, `GALCATALOG.LIST`, extracted
`satellite_properties_*.csv`, or relation-analysis `satellite_relations_*.csv`
for any of the four target outputs (00972, 00977, 00980, 00984). The 00972
`satellite_properties_ascii_v2_00972.csv` and its associated preliminary
tables are older legacy products; their schema lacks the finder-assigned
subhalo dark-matter mass and the new extraction's particle-ID audit, so they
cannot substitute for the queued validated pipeline.

The replacement FoF tasks for 00972, 00980, and 00984 remain pending with
reason `Priority`; the 00977 galaxy-finding job has not started, and all other
galaxy-finding, extraction, relation-analysis, and figure jobs remain pending
on dependencies. At this check all cluster nodes were allocated. The Slurm
scripts request CPU and memory but no GPU/GRES, so the current wait is due to
node availability/priority and upstream dependencies, not a GPU allocation
request. No job was restarted or duplicated.

## Running FoF allocation versus current script (2026-10-21 04:02 KST)

`scontrol show job 430418_1` reports the active 00977 allocation as 10 MPI tasks
with 60 CPUs across three nodes (`ReqTRES=cpu=60,mem=768G,node=3`), a 48-hour
time limit, and no GPU/GRES request. The current `opfof_array.sbatch` requests
one CPU per MPI task because opFoF is MPI-only; the already-running allocation
retains its earlier six-CPU-per-task request. The job has used about seven
hours of its two-day limit and remains active. Do not restart it to change the
allocation; replacement tasks will use the current one-CPU-per-rank script.

## Pending replacement FoF job resource requests (2026-10-21 04:03 KST)

`scontrol` confirms each replacement task (00972 `430610_0`, 00980 `430611_2`,
00984 `430622_3`) is pending with `Partition=a100`, two nodes, ten MPI tasks,
one CPU per task, 448 GB total memory, and a 24-hour time limit. Their
`ReqTRES` contains CPU, memory, and nodes only; it contains no GPU/GRES
request. The partition name `a100` identifies the node class and does not
indicate that a GPU was allocated. These are the actual scheduler records;
they are authoritative even though the current script's default directives
now specify a different partition, node count, memory, and time limit. No
pending job was modified or resubmitted.

## Extractor and analysis contract review (2026-10-21 04:04 KST)

Static review of the current extractor confirmed that NewGalFinder's `idtype=int`
star IDs are read as signed 32-bit integers, and the build's documented C
layout uses the same four-byte ID field. The extractor computes metallicity
using current stellar particle masses, computes the primary internal
one-dimensional dispersion from the same positive-mass stars with finite
velocities, and separately reports a half-mass-radius aperture dispersion.
The analysis reader requires the particle-ID duplicate audit fields before it
will use the NewGalFinder schema, converts the stored `M_sun/h` masses using
that output's `h`, and excludes galaxies with repeated or multiply assigned
stellar IDs. The implementation matches the current methods draft on these
points; no code change or test run was needed. Runtime validation still depends
on the queued extraction completing successfully.

## Live 00977 FoF process sample (2026-10-21 04:10 KST)

The scheduler still reports array task `430418_1` as `RUNNING` at 7:10 elapsed.
The MPI-step counters continue to rise, while the FoF output files remain at
03:40:55. A read-only `srun --overlap` process sample on `syn04` found three
`opfof.exe` processes in runnable state, each using about 99 per cent of one
CPU. This confirms that the current FoF stage is actively computing even
though its log and output checkpoint have not advanced. The array element maps
to underlying Slurm job `430419`; no process was interrupted or resubmitted.

## Host-mass boundary sensitivity pipeline added (2026-10-21 04:14 KST)

The manuscript plan specified nearby host-mass boundary checks, but the primary
analysis array used only the interval `1e13 <= Mhost < 1e14 Msun`. I added a
separate wrapper and CPU-only Slurm array to rerun the same host-fixed models
on four prespecified alternatives, shifting one boundary at a time by 0.2 dex:
`[10^12.8, 10^14)`, `[10^13.2, 10^14)`, `[10^13, 10^13.8)`, and
`[10^13, 10^14.2) Msun`. The fits retain the primary analysis's stellar-mass
thresholds, particle-count cut, mass-bin and matching settings, and 2,000
whole-host bootstrap draws; each output row records its exact mass limits and
variant label.

The new job is `431211`, submitted with `afterok` dependencies on all four
satellite extractions (`430492_0`, `430427_1`, `430427_2`, and `430477_3`). It
requests one CPU and 16 GB of memory, with no GPU/GRES request. Its results
will be written separately as `satellite_host_mass_sensitivity_NNNNN.csv`, so
the primary analysis output is unchanged. The scripts were statically
reviewed; no tests or calculations were run because the required catalogues
are not yet available.

## Current catalogue-pipeline state (2026-10-21 04:21 KST)

The active 00977 FoF array element `430418_1` (underlying job `430419`) is
still `RUNNING`; `scontrol` reports partition `a10`, three nodes, ten MPI
tasks, six CPUs per task, 768 GB, and no GPU/GRES allocation. `sstat` reports
about 24 hours 18 minutes of cumulative CPU time and a maximum RSS near 90.9
GB at 7 hours 20 minutes elapsed. The FoF log currently ends after P9 and P7
saved their isolated halos. The expected catalogue outputs are under
`FoF_Data/FoF.00977/` (not the raw snapshot directory `output_00977/`):
`FoF_halo_cat.00977` (22,769,732 bytes) and
`FoF_member_particle.00977` (32,654,885,672 bytes), both last modified at
03:40:55 KST. They have not yet advanced to a newer checkpoint. The Galaxy
Finder task for 00977 waits on `aftercorr:430418_*`; extraction waits on the
corresponding Galaxy Finder task. Replacement FoF tasks for 00972, 00980, and
00984 remain pending for priority, while the host-mass sensitivity array
waits on all four extractions. No running or queued task was modified or
restarted.

## NewGalFinder source and executable synchronized (2026-10-21 04:26 KST)

Before any NewGalFinder array task became runnable, I checked the binary against
the exact build target. `gfind.exe` dated 2026-09-26, while
`src/find_galaxies.c` and `src/basin_members.c` were newer (2026-09-28 and
2026-09-29); `make -n gfind.exe` showed both objects would be rebuilt. No
galaxy-finder task was running. I preserved the old executable as
`/gpfs/kjhan/Others/Jaehyun/2026.Apr/NewGalFinder/gfind.exe.before_source_sync_20261021`,
then rebuilt `gfind.exe` with the configured Intel compiler/MPI modules. The
build completed, `make -q gfind.exe` now reports it up to date, and `ldd`
reports no missing libraries. This is a compile check only; no galaxy-finder
calculation or catalogue validation has run yet.

Source review of `find_galaxies.c`, `basin_members.c`, and `binding.c` confirms
that satellite candidates receive an iterative tidal-radius cut and a
self-binding energy filter; the remaining stellar members then undergo the
finder's stellar-linking step. I updated `methods_draft.md` and the working
manuscript to describe these finder-defined memberships and to state that we
do not independently rerun the membership procedure. The active 00977 FoF job
was still running at 7:25 elapsed, and the downstream galaxy-finder tasks were
still dependency-pending when checked. Its two FoF catalogue files remained
unchanged since 03:40:55 KST. No Slurm task was modified or restarted.

## 00977 FoF process remains active (2026-10-21 04:27 KST)

A fresh scheduler check reports `430418_1` still `RUNNING` at 7:26 elapsed.
The FoF catalogue files remain at their 03:40:55 modification time, but a
read-only process sample on `syn04` found three `opfof.exe` ranks in runnable
state, each using about 99.5 per cent of one CPU. This confirms ongoing
computation rather than an idle process. No signal was sent and the task was
not restarted.

## 00977 FoF output advanced (2026-10-21 04:29 KST)

The next check found `430418_1` still running at 7:28 elapsed. Both expected
FoF catalogue files advanced at 04:27:57 KST, from 22.8 to 24.4 MB for the
halo catalogue and from 32.7 to 37.2 GB for the member-particle file. The
latest standard output shows further NewDD slabs being read and tree-building
ranks passing their checkpoints. `GALFIND.DATA.00977` and the extracted
satellite CSV do not exist yet; the galaxy-finder task remains dependency-
pending. The FoF run is progressing, so no intervention was made.

## Literature framing checked against journal articles (2026-10-21 04:30 KST)

I checked the cited papers against their MNRAS articles. Spindler & Wake
(2017, DOI `10.1093/mnras/stx427`) directly study stellar metallicity versus
core velocity dispersion for SDSS central and satellite galaxies, reporting
little central--satellite metallicity difference after splitting by star
formation activity. Fattahi et al. (2018, DOI `10.1093/mnras/sty408`) directly
show a Local Group `[Fe/H]`--`V_{1/2}` relation and attribute its low-velocity
flattening to stripped satellites; their inferred-progenitor correction
assumes that stripping leaves mean stellar metallicity unchanged. Neither
study tests the present paper's narrower conditional question: the residual
metallicity--internal-dispersion association among satellites after controlling
for both stellar mass and host in a cosmological simulation. I revised
`introduction_draft.md` and the integrated manuscript to state this distinction
explicitly and to describe the earlier results more precisely. The current
00977 FoF job remained active at the last scheduler check; its files had
advanced at 04:27:57 KST, while the galaxy-finder stage was still pending.

## Host-mass literature added to introduction (2026-10-21 04:34 KST)

A further MNRAS study directly motivates the host control: Gallazzi et al.
(2021, DOI `10.1093/mnras/stab265`) find a mild increase in satellite
stellar metallicity with host halo mass at fixed stellar mass below
`10^10.5 Msun`. I added this result and a verified BibTeX entry to both
`introduction_draft.md` and the integrated manuscript. The paper studies
stellar populations, not internal stellar dispersion, so it supports the
need for host control without replacing the present paper's distinct
kinematic question.

At the latest check, FoF task `430418_1` remained `RUNNING` at 7:34 elapsed.
Its log had advanced to P2 finding 5.56 million halos from 19.7 million
particles and checking the halo boundary. The two final FoF catalogue files
had not changed since 04:27:57 KST, and the 00977 Galaxy Finder still waited
on the FoF dependency. The process was left running.

## Stellar birth-time interpretation checked (2026-10-21 05:24 KST)

The canonical extractor reads each star particle's `tp`, but does not yet
derive stellar age or recent star-formation rate. I checked the official
RAMSES documentation before adding such quantities: it identifies `tp` as a
birth epoch, and states that non-radiation-hydrodynamics runs store birth
times in conformal units; conversion to proper time is needed when
post-processing those runs ([particle fields](https://ramses-organisation.readthedocs.io/en/latest/dev_docs/lecture/particles.html),
[RAMSES guide, `convert_birth_times`](https://ramses-organisation.readthedocs.io/en/latest/Ramses.pdf)).
The local NewGalFinder header only labels `tp` as being in RAMSES units, while
the separate `README4GADGET.md` describes the GADGET input mapping and does
not establish the convention of these native RAMSES snapshots. Therefore,
`tp` must not be treated as a scale factor or directly converted to age/SFR
until the exact output time convention and its conversion are verified. This
prevents adding misleading age or activity controls to the canonical
metallicity--dispersion analysis.

At 05:24 KST, FoF task `430418_1` was still `RUNNING` at 8:24 elapsed, but
its output catalogue files and log had not advanced since 05:21:39 KST. The
downstream Galaxy Finder, extraction, and analysis tasks remained pending on
dependencies. No job was stopped or restarted.

At 05:25 KST, `scontrol` still reported `430418_1` as `RUNNING`. A read-only
process sample across its three allocated nodes found all ten `opfof.exe`
MPI ranks running at about 99.4--99.5 per cent CPU, so the unchanged output
files reflect ongoing computation rather than an idle MPI job. Slurm reports
60 CPUs and 768 GB of memory allocated, with no GPU in `AllocTRES`; the `a10`
partition name did not allocate a GPU. Output files still had their 05:21:39
timestamps. The Galaxy Finder dependency remains correctly pending until FoF
finishes.

## 00977 FoF output stage identified (2026-10-21 05:32 KST)

The latest log ends while rank P0 is preparing to write isolated halos from a
slab containing about 5.9 million halo candidates; ranks P1--P9 have sent
their isolated-halo data. Source inspection of `opfof.c` and
`Treewalk.fof.ordered.c` shows that, in the active default writer branch, P0
first walks its local halo memberships and writes them, then receives and
writes the other ranks' records. Only after that function returns do ranks
enter the barrier, write bottom-face contact records, and advance to the next
slab. A read-only check found P0 still running with the FoF output files open;
their file positions and sizes were unchanged over a one-minute interval.
This is a slow active output-construction stage, not evidence that FoF has
finished or that it is safe to launch NewGalFinder on partial outputs. The
chain was left untouched.

At 05:33 KST, the latest `has read` entries show ranks P0--P9 on shards
`00004`, `00014`, ..., `00094`, respectively: the fifth input shard in each
rank's ten-shard block. This is halfway through the 100 shard files by count,
not necessarily by bytes or runtime because shard particle counts differ. The
FoF catalogue file sizes and modification times were unchanged at this check;
the current slabs were still being linked and written.

At 05:35 KST, a fresh check shows FoF still `RUNNING` and the outputs have
advanced: the halo catalogue grew from 29,626,404 to 30,560,292 bytes, while
the member-particle file grew from 49,185,323,176 to 50,042,427,560 bytes.
The latest log reports P4 completing its isolated-halo write for a slab with
6.36 million halo candidates and 71.8 million particles. The downstream
NewGalFinder and extraction jobs remain dependency-pending until the whole FoF
catalogue is complete.

At 05:38 KST, the FoF job remained `RUNNING`. Since the 05:34:35 output
update, the file sizes, output offsets, and log tail had not advanced. A
read-only sample still found rank P0 running at 99.5 per cent CPU; its
accumulated CPU time rose from 08:34:43 to 08:36:09 across the observation
interval. The log remains in the isolated-halo output stage, so the process is
still doing CPU work on this slab even though no new output was visible in
that interval. No task was interrupted.

At 05:42 KST, all ten FoF processes still used about 99.4--99.5 per cent CPU,
but the files and log remained unchanged since 05:34:35. The latest worker
messages show ranks P1 and P4 have sent their isolated-halo records and
returned from that write call, while P0 has not reported the end of the
isolated-halo output. The default `WriteIsolatedHalo` source path makes P0
scan its local halo memberships and then serially receive/write records from
the other ranks before the barrier and bottom-face step. This is consistent
with the current delay being concentrated in the master-writer stage; active
CPU use alone does not establish progress through its loop. The task remains
live and was not interrupted.

At 05:45 KST on 2026-10-21, a targeted `srun --overlap` inspection of the
live allocation found all ten `opfof.exe 977 100` MPI processes still running
at about 99.4--99.5 per cent CPU. The FoF catalogue files and log remained at
their 05:34:35 sizes/timestamps. Source inspection confirms that each rank
still scans its local halo list; rank P0 then receives and writes worker
records serially. The log shows P1, P4, and P9 have completed their sends,
while the remaining ranks have not yet reported completion. Thus the unchanged
files do not mean the job is idle or complete. The live FoF allocation requests
only CPU, memory, and nodes; NewGalFinder also has no GPU GRES request. Its
Slurm `PENDING` reason is the FoF `aftercorr` dependency, and it requests two
A100-partition nodes with 600 GiB per node for the large RAM working set.
No running task was interrupted and no partial FoF output was passed downstream.

A static format/unit audit of the canonical NewGalFinder source was completed
on 2026-10-21. `NewGalFinder/ramses.h` defines particle positions in cMpc/h,
velocities in km/s, and current particle masses in Msun/h; its `StarType`
contains the stored `tp`, `zp`, elemental abundances, and initial mass. The
NewGalFinder README documents catalogue component masses in Msun/h and bulk
velocities in km/s. These conventions match the extractor's unit handling:
it converts catalogue stellar masses by the snapshot h parameter, converts
its measured half-mass radius to physical kpc, and uses the stored particle
velocities directly in km/s. The `zp` values are retained as the simulation's
total birth-metallicity mass fraction and are not interpreted as [Fe/H]. This
source audit supports the measurement definitions but does not replace the
pending particle-ID, particle-mass, and catalogue-integrity checks on completed
outputs.

At 05:48 KST, Slurm still reported FoF task 430418_1 as `RUNNING`; its
catalogue files had the same sizes as at 05:34:35. A targeted overlapping
process sample at 05:48:23 again found all ten MPI ranks at 99.4--99.5 per cent
CPU, with their accumulated CPU times increased from the 05:45:55 sample. This
confirms continued computation but not completed output. Downstream jobs remain
held by the FoF dependency, and no partial catalogue was consumed.

The NewGalFinder build configuration was also compared with the extractor's
binary record layout. Its Makefile uses `NCHEM=9`, `NDUST=4`, and does not define
`DOUBLE_PRECISION`, `OUTPUT_PARTICLE_POTENTIAL`, or `NBODY`; `ramses.h` therefore
uses the 32-bit particle ID and double-precision particle fields assumed by the
extractor. The executable `NewGalFinder/gfind.exe` was built at 04:25 KST on
2026-10-21, before the pending galaxy-finder job. This is a source/build
configuration audit, not a substitute for post-run record-length and catalogue
integrity checks.

## 00977 late-rank stage advance (2026-10-21 05:52 KST)

The earlier 05:49 GDB sample placed rank 8 in `FoF_divide_node()` while it
built a particle tree. A second read-only GDB sample at 05:51:xx placed the
same confirmed MPI rank (`PMI_RANK=8`, PID 386215 on syn04) in
`new_fof_link()` at `Treewalk.fof.ordered.c:919`, linking particles against
the tree. This is a later stage than the prior sample and shows the rank
advanced through tree construction. Slurm reported the FoF job `RUNNING` at
8:51 elapsed at 05:51; the latest log and 00977 output files were still
unchanged from 05:34:35 at 05:52. The pending NewGalFinder step must continue
to wait for a complete FoF catalogue.

## Correction to the 05:52 rank-8 interpretation (2026-10-21 05:53 KST)

The 05:51 GDB sample located rank 8 in `new_fof_link()`, but a further sample
at 05:52 located it in `FoF_divide_node()`. Since these are point-in-time
samples and optimized local variables/backtraces are unreliable, they do not
establish monotonic stage progress between samples. The prior 05:52 note's
claim that rank 8 had definitively advanced through tree construction is
withdrawn. Slurm still reports FoF `430418_1` as `RUNNING` at 8:53 elapsed;
the log and output files remain unchanged from 05:34:35. The calculation is
live, but recent file-level progress is not confirmed.

## Multi-output Slurm dependency audit (2026-10-21 05:55:01 KST)

The pending chain was inspected directly in Slurm. NewGalFinder task 430423_1
(00977) is waiting on `aftercorr:430418_*`; task 430423_2 (00980) is waiting
on `afterok:430611_*`. Their extraction and relation-analysis jobs remain
chained to those galaxy-finder tasks. The FoF jobs for 00972, 00980, and 00984
(array tasks 430610_0, 430611_2, and 430622_3) request CPU, memory, and nodes
only; their `ReqTRES` contains no `gres/gpu`, and Slurm reports `Priority`, not
a GPU quota, as their pending reason. They each request 10 CPUs and 448 GB
across two A100-partition nodes. The 00977 job remains live on A10 nodes with
no GPU allocation. No pending task was canceled or resubmitted; the dependency
chain is valid and the current wait is attributable to FoF completion and
scheduler priority.

## Rank-8 hot-loop sample (2026-10-21 05:55 KST)

A new read-only GDB sample of confirmed MPI rank 8 (PID 386215) again located
it in `new_fof_link()` at line 919 with argument `nhalo=0`, matching the 05:51
sample. It is therefore still processing the first seed call visible at those
samples; optimized local loop counters were nonsensical and are not used as
progress evidence. Slurm still showed the parent FoF task `RUNNING` at 8:55
elapsed, while the output and log timestamps remained at 05:34:35. No single
sample proves a hang, but this call has not yielded a new log marker or output
write during the interval observed. Downstream catalogue processing remains
correctly dependency-gated.

## Rank-8 CPU/I/O interval sample (2026-10-21 ~05:59 KST)

Two process-counter samples roughly one minute apart on rank 8 (PID 386215 on
syn04) show the process still `Rl` at 99.4 per cent CPU; accumulated CPU time
rose from 08:54:16 to 08:55:36. Its `/proc/386215/io` counters did not change
(`rchar=42080966616`, `wchar=9130533123`, `syscr=19004`, `syscw=109293` at
both samples), while `read_bytes` and `write_bytes` remained zero. The FoF
catalogue files and job log also remained unchanged. Combined with the recent
GDB sample in `new_fof_link()`, this supports an in-memory, CPU-bound linking
calculation during the interval, rather than a filesystem stall. It does not
show how much work remains or justify interrupting the live run.

## Repeated rank-8 linking-call sample (2026-10-21 ~06:00 KST)

A further GDB attach at 05:59 again placed rank 8 (PID 386215) at
`new_fof_link()` line 919 with the same argument addresses and `nhalo=0` as the
05:55 sample. Combined with the unchanged process I/O counters and increasing
CPU time, this call has remained active for at least four minutes without a
visible file or log update. This identifies a long-running in-memory linking
call, but optimized counters do not reveal its iteration count or remaining
work; it is not sufficient evidence of a deadlock. The live job was left
untouched.

## 00977 live-resource check (2026-10-21 06:03 KST)

Slurm still reports FoF task 430418_1 `RUNNING` at 9:03 elapsed. `sstat`
reports `MaxRSS=90908028K` (about 87 GiB) against 768 GiB allocated, so the
current delay is not explained by exhausting the requested memory. The `/gpfs`
filesystem reports about 236 TiB available and the user's GPFS quota has no
configured block limit, so storage capacity is not the current constraint.
The log and FoF output modification times remain at 05:34:35. Together with
rank 8's CPU-active, unchanged-I/O samples, the evidence still points to a
CPU-bound linking phase; completion time remains unknown.

At 06:03 KST, the latest FoF log still places each rank on its fifth input
shard of ten: P0--P9 are on shards `00004`, `00014`, ..., `00094`. This is
about halfway by shard count only; shard sizes and processing time vary. Rank 8
has loaded shard `00084` (142,836,000 reported particles) and passed tree
construction, while the most recent GDB sample places it inside linking. The
log and output files have not advanced since 05:34:35, so the shard count must
not be treated as a runtime-completion percentage.

## Rank-8 GDB location updated (2026-10-21 06:04 KST)

A 06:04 GDB sample of the same live rank 8 (PID 386215) places its instruction
pointer in `FoF_divide_node()` at line 964, whereas the 05:59 sample placed it
in `new_fof_link()` at line 919. This means the earlier linking call cannot be
claimed to remain active now, nor do these isolated samples establish a
monotonic transition between stages; optimized locals and backtraces are not
reliable enough to infer shard progress. Slurm still reports the job
`RUNNING` (9:04 elapsed), but the output files and log remain unchanged since
05:34:35. No estimate of remaining work is justified.

The 00977 FoF stderr was searched at about 06:06 KST for fatal errors,
segmentation faults, MPI aborts, overflows, and NaN diagnostics. None were
found. Its tail contains only the expected notices that absent GAS and SINK
component files are treated as empty. This supports treating the live CPU-bound
job as still computing rather than as a failed run; completion is not inferred.

## FoF input-read marker interpretation (2026-10-21 06:07 KST)

The RAMSES reader prints each `P%d has read ...` component marker and flushes
stdout immediately after it. Therefore, the latest logged component-read markers
(P0--P9 on shards `00004` through `00094`) reliably mean no later component
read has completed and been logged as of this check. This narrows what the stale
log tells us, but it does not identify the rank's current inner computation:
other FoF progress messages can remain buffered. The 00977 array task is still
`RUNNING` at 9:07 elapsed; no failure or completion is inferred from the lack of
new read markers.

## Live FoF and downstream gate (2026-10-21 06:08 KST)

The scheduler still reports array task `430418_1` (underlying job `430419`)
`RUNNING`, with 60 CPUs and 768 GiB allocated and no GPU GRES requested. At
09:08 elapsed, `sstat` reports 90,960,444 KiB maximum RSS and accumulated CPU
and disk I/O, consistent with an active compute step. Its FoF log and both
`FoF.00977` products remain timestamped 05:34:35; the directory contains the
FoF halo and member-particle files but no NewGalFinder catalogue. Jobs
`430423_1` and `430427_1` remain dependency-pending, so extraction and relation
analysis cannot yet run on the canonical 00977 output.

A static source review confirms that the prepared extractor computes the
current-mass-weighted `zp`, global and half-mass-radius stellar dispersions,
particle-mass consistency, and within/across-galaxy star-ID audits. The
analysis code fits equal-weight host-fixed-effects models with stellar-mass
control and host bootstrap intervals, with separate subhalo-mass and aperture
sensitivity models. No tests or new calculations were run during this review;
canonical numerical conclusions remain unavailable until the dependent jobs
produce and validate the catalogues.

## Additional snapshot FoF tasks started (2026-10-21 06:10 KST)

A new scheduler check shows FoF array task `430610_0` (snapshot 00972) and
`430611_2` (snapshot 00980) have started, each on two nodes with 10 CPUs and
448 GiB, with no GPU GRES requested. The 00977 task `430418_1` remains active
at 09:10 elapsed. The 00984 task `430622_3` is still pending for resources.
NewGalFinder and extraction remain dependency-pending for 00977; canonical
satellite catalogues and relation results are not yet available.

## 00972/00980 live FoF input progress (2026-10-21 06:11 KST)

The live logs for `430610_0` (00972) and `430611_2` (00980) have advanced:
rank P7 has read the `.00070` GAS shard (zero particles) and about 54 million
particles in the corresponding shard, while several ranks report tree building.
Their 28-byte halo-catalogue placeholders and zero-byte member-file
placeholders are not completed FoF products; both jobs remain `RUNNING` and
must finish successfully before downstream NewGalFinder can use them. Stderr
shows only the expected missing GAS/SINK component notices treated as empty.

## FoF queue and fresh progress check (2026-10-21 06:12 KST)

At 06:12, 00977 (`430418_1`) remains `RUNNING` at 09:12 elapsed; 00972
(`430610_0`) and 00980 (`430611_2`) remain `RUNNING` at 02:36. Their latest
logs now show tree-build completion on several ranks and P1 has checked/saved
isolated halos for a shard, but the 28-byte/zero-byte files are still only
initial placeholders. 00984 (`430622_3`) remains pending for resources. All
NewGalFinder, extraction, primary-analysis, and host-sensitivity jobs shown in
`squeue` remain dependency-pending. No canonical galaxy catalogue or relation
CSV exists for these outputs yet.

## Fresh FoF progress (2026-10-21 06:14 KST)

At 06:14, 00972 (`430610_0`) and 00980 (`430611_2`) are still running
(about 4 minutes elapsed). Their log timestamps have advanced to 06:13:53 and
06:13:48, respectively; both logs now show P6 sending and saving isolated
halos for a shard. The output files remain 28-byte halo headers and empty
member-file placeholders, so final catalogue generation is still in progress.
00977 remains `RUNNING` at 09:14 elapsed, with its log and 50-GB member file
unchanged since 05:34:35. NewGalFinder/extraction/analysis tasks are still
waiting on dependencies; 00984 FoF is pending for resources.

## FoF activity revalidated (2026-10-21 06:15 KST)

At 06:15 the scheduler still reports 00977, 00972, and 00980 FoF tasks as
`RUNNING`; downstream catalogue and analysis tasks remain dependency-pending,
and 00984 remains pending for resources. The 00972/00980 logs advanced again
(to 06:14:54 and 06:15:01) and show P2 completing a shard save/check. The
00977 log and products remain unchanged since 05:34:35, but its `sstat` CPU
counter rose from 1d06:24:51 to 1d06:41:52 since the prior measurement, with
MaxRSS 93,811,260 KiB. This is evidence of continuing CPU activity, not a
completion-time estimate. The 00972/00980 catalogue files remain empty/header
placeholders; no NewGalFinder or canonical relation products are present yet.

## FoF check (2026-10-21 06:16 KST)

At 06:16 all three launched FoF tasks (00972, 00977, 00980) are still
`RUNNING`; the 00984 task remains pending for resources. Logs for 00972 and
00980 advanced to 06:16:40 and 06:16:18, respectively, and show another rank
checking halo boundaries. Their catalogue files remain at the initial 28-byte
and zero-byte sizes. The 00977 log/products remain unchanged since 05:34:35;
its prior `sstat` CPU increase remains the evidence that it is computing. All
downstream jobs are still dependency-pending, and no canonical satellite list
or relation result has been produced.

## FoF jobs rechecked (2026-10-21 06:18 KST)

The 00972 and 00980 FoF tasks remain active at about 8 minutes; their logs
advanced to 06:17:57 and 06:17:30 and show P4 checking halo boundaries. Their
halo/member files are still 28-byte/zero-byte initial placeholders. The 00977
task remains active at 09:17 elapsed, with its log and existing large FoF
files unchanged since 05:34:35. The 00984 task is still resource-pending.
NewGalFinder, extraction, primary relation analysis, and host-mass sensitivity
jobs remain dependency-pending; no canonical galaxy catalogues or analysis
outputs are present.

## Analysis-label clarification (2026-10-21)

Static review found that the analyzer's docstring incorrectly said the
secondary host orbital-dispersion comparison did not match stellar mass. The
code actually pairs low- and high-metallicity satellites within each host
within a stellar-mass caliper. Updated the docstring to state this procedure
and distinguish the galaxies' orbital velocities from the primary internal
stellar dispersion. Renamed the CSV field to
`median_log10_orbital_sigma_highZ_over_lowZ` to make that distinction explicit.
No tests or calculations were run.

## CPU activity during buffered FoF phases (2026-10-21 06:21 KST)

Although the 00972/00980 log files have not advanced since 06:17, Slurm `sstat`
confirms active CPU and I/O for both live steps: task 430610 has accumulated
57:24 CPU time with 27.96 GB disk reads; task 430611 has 57:40 CPU time with
28.09 GB disk reads, each after about 11:41 wall time. MaxRSS is about 45 GB
for both. The 00977 step has also accumulated CPU time (1d07:03:01) and has
MaxRSS 94,016,060 KiB. These measurements support ongoing computation despite
buffered/stale stdout; they do not identify completion time. No final FoF or
NewGalFinder catalogue has appeared.

## FoF live-step recheck (2026-10-21 06:23 KST)

All three started FoF tasks remain `RUNNING`; 00984 remains pending for
resources, and the galaxy-finder/extraction/analysis dependency chain remains
pending. The 00980 log advanced to 06:22:34 and reports P8 checking the halo
boundary. The 00972 log remains buffered since 06:17:57, but its CPU counter
rose to 01:05:21 and read volume remains about 27.96 GB. The 00980 CPU counter
rose to 01:05:40 with about 28.09 GB read. The 00977 CPU counter rose to
1d07:08:21, with MaxRSS 94,149,180 KiB. These counters confirm continued work;
no final products or completion estimate can yet be claimed.

## FoF resource check (2026-10-21 06:24 KST)

The 00977 task remains `RUNNING` at 09:24, with 1d07:12:57 accumulated CPU
and MaxRSS 94,288,444 KiB. Tasks 00972 and 00980 remain `RUNNING` at 14:40;
their CPU totals rose to 01:12:12 and 01:12:33, while each has read about
28.0 GB and reached about 45 GB MaxRSS. Their log and placeholder-output times
have not advanced since the previous check. No diagnostic error was found in
the stderr tails; only absent GAS/SINK components treated as empty. The 00984
FoF task remains resource-pending, and downstream tasks are dependency-pending.

## FoF activity (2026-10-21 06:26 KST)

At 06:26, 00972/00977/00980 remain `RUNNING`; 00984 remains pending for
resources and downstream jobs are dependency-pending. Logs for 00972 and 00980
advanced to 06:24:13 and 06:25:26 and show P8/P9 checking halo boundaries.
Their member files are still empty and halo files remain 28 bytes. The 00977
log and product timestamps remain at 05:34:35; its `sstat` CPU increased to
1d07:17:13, with MaxRSS 94,870,076 KiB. The 00972/00980 CPU totals increased
to 01:18:35 and 01:18:58, with disk reads still about 28 GB each. FoF work is
active, but no final products or downstream catalogues are available.

## Manuscript contribution statement updated (2026-10-21)

Replaced the acknowledgements placeholder with separate attribution: K.J.H.
provided the MPI+OpenMP halo-finding and galaxy-finding codes; AI-assisted
tools helped develop catalogue extraction/validation, analysis, and figure
scripts and assisted with manuscript drafting/editing. Kept data availability
as a placeholder pending confirmation of simulation-data access conditions.
No tests or calculations were run.

## FoF and manuscript check (2026-10-21 06:27 KST)

The 00972/00977/00980 FoF tasks remain `RUNNING`; 00984 remains resource-
pending and all downstream jobs remain dependency-pending. FoF logs/products
have not advanced since the previous check: 00972/00980 still have only
28-byte/zero-byte placeholders, while 00977 retains its pre-existing products.
No NewGalFinder, extraction, or relation result is present. Separately, the
working manuscript now distinguishes the user's supplied MPI+OpenMP halo and
galaxy finders from AI-assisted catalogue, analysis, figure-script, and
writing contributions; data-access terms remain unresolved.

## AI-use disclosure aligned with MNRAS guidance (2026-10-21)

Checked the current MNRAS author instructions and OUP AI-use policy. MNRAS
requires AI use to be disclosed in the cover letter and in Methods or
Acknowledgements; OUP asks the manuscript disclosure to identify tool/version,
use, and the authors' role in directing/reviewing the work. Updated the draft
acknowledgement to name OpenAI Codex (GPT-6), Claude Opus 5.5, and GPT-6 Astra,
describe their work, and state K.J.H.'s direction and author responsibility.
The final cover-letter disclosure remains due at submission. See the MNRAS
instructions and OUP policy for the current requirements.

## Active FoF confirmed (2026-10-21 06:30 KST)

Tasks for 00972, 00977, and 00980 are still `RUNNING`; 00984 remains pending
for resources. `sstat` counters increased for every active step: CPU totals are
1d07:30:04 for 00977, 01:37:47 for 00972, and 01:38:15 for 00980. Disk reads
are about 94.3 GB for 00977 and 28.0 GB for each newer task. No NewGalFinder,
extracted-satellite, or relation CSVs are present. Existing pending jobs remain
blocked on their declared dependencies; do not bypass those dependencies while
FoF products are still being written.

## FoF progress rechecked (2026-10-21 06:31 KST)

At 06:31 all three launched FoF tasks remain `RUNNING`; 00984 is pending for
resources and downstream jobs remain dependency-pending. Logs for 00972/00980
advanced to 06:29:18/06:28:57 and show P5 checking and saving isolated halos
for another shard. Their live steps' CPU totals increased to 01:44:28 and
01:44:58 with roughly 28 GB read each. The 00977 step remains active at
09:31 elapsed; its CPU total rose to 1d07:34:33 and MaxRSS is 95,300,156 KiB,
though its log/products remain unchanged since 05:34. The output catalogue
placeholders for 00972/00980 remain incomplete.

## FoF activity (2026-10-21 06:32 KST)

Slurm still reports FoF tasks for 00972, 00977, and 00980 as `RUNNING`; 00984
is resource-pending, and all downstream jobs are dependency-pending. `sstat`
CPU totals rose to 1d07:38:42 (00977), 01:50:40 (00972), and 01:51:12
(00980), with disk reads unchanged at about 94.3 GB and 28.0 GB per newer
snapshot. No satellite catalogue files have appeared. `scontrol listpids` is
not available from this login node (no local slurmd domain socket); this
failed process-list query does not change the scheduler or accounting evidence
that the jobs are running.

## Introduction literature claims checked (2026-10-21)

Checked the cited SDSS and simulation studies against their primary article
records. The current draft's summaries are supported: Gallazzi et al. report a
mild host-halo-mass dependence in stellar metallicity for low-mass satellites;
Spindler & Wake find little central--satellite metallicity offset after
separating star-forming and quiescent galaxies; NIHAO XV reports overlapping
satellite/field trends with lower satellite dispersions at fixed stellar mass;
and Fattahi et al. discuss a Local Group [Fe/H]--$V_{1/2}$ relation altered by
stripping. These use differing metallicity and kinematic definitions, so the
draft correctly limits the comparison to context rather than claiming prior
work tested the same host- and stellar-mass-controlled internal-dispersion
relation. No text change was required.

At the same check (06:33 KST), the 00972/00977/00980 FoF tasks were still
`RUNNING`, while 00984 waited for resources and downstream tasks waited on
dependencies. Logs last advanced at 06:29:18 (00972), 06:28:57 (00980), and
05:34:35 (00977). No NewGalFinder, satellite-property, or relation outputs
were present.

## FoF CPU check (2026-10-21 06:35 KST)

All three FoF tasks remain `RUNNING`; 00984 remains pending for resources,
and NewGalFinder/extraction/analysis tasks remain dependency-pending. Since
06:33, `sstat` CPU counters rose to 1d07:47:22 (00977), 02:03:36 (00972), and
02:04:11 (00980). The 00972/00980 logs have not advanced since 06:29 and
00977's since 05:34; their output files remain unchanged. This is continuing
CPU work with delayed log/output progress, not evidence of a hang or a runtime
estimate. No canonical satellite catalogue is available yet.

## FoF status revalidated (2026-10-21 06:36 KST)

The 00972, 00977, and 00980 FoF tasks remain `RUNNING`; 00984 remains pending
for resources. CPU totals increased to 02:10:30 (00972), 02:11:07 (00980), and
1d07:51:59 (00977). Disk reads remain about 28.0 GB, 28.1 GB, and 94.3 GB,
respectively. Logs have not advanced since 06:29 (00972/00980) or 05:34
(00977), and 00972/00980 outputs are still 28-byte/zero-byte placeholders.
00977 still has its existing 30.6-MB/50.0-GB files, but the job has not exited,
so they are not yet accepted as final outputs. Galaxy finding and all
satellite-analysis tasks remain dependency-pending; no canonical galaxy or
relation files exist.

## FoF status revalidated (2026-10-21 06:38 KST)

The live Slurm check still shows FoF 00972, 00977, and 00980 as `RUNNING`, with
00984 pending for resources. CPU counters increased since 06:36 to 02:19:24
(00972), 02:20:04 (00980), and 1d07:57:57 (00977), confirming ongoing CPU
execution. The 00972 and 00980 log files still end at the P5 isolated-halo
save message (06:29); the 00977 log remains at 05:34. Exact FoF output-file
stats show 00972 and 00980 still have 28-byte halo-catalogue placeholders and
empty member files (mtime 06:10); 00977 remains at 30.6 MB and 50.0 GB (mtime
05:34) while its job is live, so none of these are accepted as completed
catalogues. FoF output files for 00984 do not yet exist. No restart or
termination action was taken.

## FoF pipeline check (2026-10-21 06:39 KST)

A new filtered Slurm query confirms 00972, 00977, and 00980 are still live;
their CPU counters rose again to 02:24:54, 1d08:01:38, and 02:25:36,
respectively. The downstream galaxy-finder, extraction, relation-analysis,
host-mass-sensitivity, and figure jobs checked remain dependency-pending.
Exact stats show no FoF output size or timestamp has advanced since the 06:38
check; the logs also remain unchanged at the previous save messages. Therefore
this is a verified wait on active CPU jobs, with no completed catalogue to
analyze yet. 00984 is still queued for resources.

## FoF stage diagnosis (2026-10-21 06:42 KST)

The filtered queue and `sstat` still show 00972, 00977, and 00980 running.
Their current CPU totals are 02:35:46 (00972), 02:36:30 (00980), and
1d08:08:55 (00977); these rose since the prior poll. The 00972 and 00980 logs
show all ten ranks have built their local trees, eight have completed halo
boundary checks, and three have finished writing/sending isolated-halo data.
Inspection of the active `opfof.c` call sequence shows that ranks then wait at
an MPI barrier before writing bottom-face contact particles; the expected
contact files are not yet present for 00972/00980. Thus the snapshots are
still inside FoF's per-slab halo assembly/communication stage, not ready for
NewGalFinder. Their disk counters have not changed from the earlier poll while
CPU continues to accrue. The expected missing GAS/SINK warnings remain the
only messages in the inspected error logs.

## Downstream resource check (2026-10-21 06:44 KST)

The queued NewGalFinder, extraction, analysis, sensitivity, and figure jobs
still show `Dependency` as their pending reason. The scheduler reports no GRES
request for these jobs, and `gfind_array.sbatch` contains no `--gres` option;
they are CPU-only MPI/OpenMP jobs despite using the A100 partition. The
NewGalFinder array requests 600 GB per node and, for the main configuration,
100 CPUs across two nodes. That memory request is large and may affect when it
starts after its dependency clears, but it is not the current pending reason.
No scheduler or resource changes were made.

## FoF output progress (2026-10-21 06:45 KST)

The latest Slurm poll still shows 00972, 00977, and 00980 running. Unlike the
previous check, exact file stats now show 00977's FoF halo catalogue and member
file have grown to 31.4 MB and 52.1 GB, respectively, with a new mtime of
06:42:34; its log has also advanced to P7 saving isolated halos. These are
ongoing outputs, not yet complete while the job remains live. The 00972 and
00980 output files and logs have not advanced; their CPU counters continue to
rise. Downstream NewGalFinder and analysis tasks remain dependency-pending.

## FoF output progress (2026-10-21 06:49 KST)

The latest queue and accounting checks still show 00972, 00977, and 00980
running, with CPU totals increasing. New output progress is visible for 00972:
at 06:49:23 its FoF halo file had grown from the 28-byte placeholder to
172,060 bytes and its member file to 10,862,592 bytes; its log now shows
additional ranks completing isolated-halo writes. These files are still
partial because the FoF job remains active. 00977 remains at 31,404,068 and
52,130,289,832 bytes (last update 06:42:34), while 00980 remains at its
placeholder sizes. Downstream galaxy-finder and analysis jobs remain
dependency-pending; no partial FoF output was consumed.

## FoF output progress (2026-10-21 06:52 KST)

All three active FoF jobs remain `RUNNING` and their CPU counters continue to
increase. The 00980 snapshot has now also begun writing FoF output: at 06:52:48
its halo file was 172,060 bytes and its member file was 10,838,016 bytes. These
are partial files while the job remains live. Snapshot 00972 remains at
172,060/10,862,592 bytes (last update 06:49:23), and 00977 remains at
31,404,068/52,130,289,832 bytes (last update 06:42:34). No `GALFIND.DATA` file
is ready for any of these snapshots; downstream galaxy-finding and analysis
jobs remain dependency-pending.

## FoF writer-stage check (2026-10-21 06:55 KST)

The filtered queue and `sacct` still show all three FoF jobs `RUNNING`, with
their CPU totals increasing; downstream jobs remain dependency-pending. The
00980 log identifies its current expensive stage: rank 0 has begun writing
isolated-halo data for 5,162,663 halos containing 41,364,410 particles, after
other ranks sent their isolated-halo data. The 00980 partial output remains
172,060/10,838,016 bytes at 06:52:48, and 00972 remains
172,060/10,862,592 bytes at 06:49:23. This matches the serial master-writer
section in `WriteIsolatedHalo`; CPU continues to accrue while that loop runs.
00977 remains at 31,404,068/52,130,289,832 bytes from 06:42:34. No output has
passed its completion checks and no partial data were analyzed.

## FoF rank-0 write check (2026-10-21 06:56 KST)

The next live queue and accounting poll again shows all three FoF tasks running;
CPU totals rose to 1d08:58:02 (00977), 03:49:06 (00972), and 03:50:12
(00980). In both newer logs, nonzero ranks have sent and saved their isolated
haloes, while rank 0 remains in `WriteIsolatedHalo`: it reports 5,188,551
haloes for 00972 and 5,162,663 for 00980 before its serial write loop. The
partial FoF output sizes have not advanced since their first writes. This
pinpoints the active stage more precisely; it is not an observed job
termination or catalogue completion.

## FoF and dependent-pipeline check (2026-10-21 07:01 KST)

The project directory was confirmed as
`/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020`. Slurm
still reports FoF tasks 430418_1 (00977), 430610_0 (00972), and 430611_2
(00980) as `RUNNING`; 430622_[3] (00984) is `PENDING (Resources)`. The two
10-CPU tasks' accumulated CPU times increased from 03:49:06/03:50:12 at 06:56
to 04:13:00/04:14:12 at 07:01, confirming continued CPU work. The FoF logs and
partial output files have not advanced: 00972 remains at 172,060/10,862,592
bytes (06:49:23), 00977 at 31,404,068/52,130,289,832 bytes (06:42:34), and
00980 at 172,060/10,838,016 bytes (06:52:48). Thus none of these catalogues
has passed completion checks. The preceding queue check showed downstream
NewGalFinder and analysis arrays dependency-gated; they were not re-queried in
this update. No partial FoF output was used and no job was restarted or
cancelled.

## FoF straggler diagnosis (2026-10-21 07:03 KST)

The live queue still shows 00972, 00977, and 00980 FoF tasks running. For
00972 and 00980, ranks 1, 2, 3, 4, 5, 6, 8, and 9 report that they sent their
isolated-halo data to rank 0; rank 7 has not reported a halo count, boundary
check, or send completion. Its last log message in each job is `passed Tree
building` (about 54 million particles assigned), whereas rank 0 and the other
workers have advanced farther. This points to rank 7's particle-linking walk
as the likely straggler, although the log has no progress counters to prove
its exact current instruction. The active `WriteIsolatedHalo` path also has
no progress messages inside its long loops. CPU accounting continues to rise,
but output sizes and logs have not changed since the previous check. For
00977, rank 7 has completed the isolated-halo write, but the job has not
reached a completion marker. These live jobs remain untouched; partial FoF
files are still excluded from analysis.

## FoF and galaxy-finder dependency check (2026-10-21 07:06 KST)

Slurm still reports FoF 00977 (430418_1), 00972 (430610_0), and 00980
(430611_2) as running; 00984 (430622_[3]) is pending for resources. Galaxy
finder tasks 430423_1 (00977) and 430423_2 (00980), along with the queued
extraction tasks, remain dependency-pending. The 00972 and 00980 CPU totals
rose from 04:22:32/04:23:47 at 07:03 to 04:38:05/04:39:25 at 07:06, while
all three FoF catalogue/member files retain their previous sizes and mtimes.
Thus the jobs remain live but no snapshot has yet cleared the FoF completion
gate; no downstream task can safely consume the partial files.

## Analysis-sample mask correction (2026-10-21 07:08 KST)

While the FoF jobs run, I reviewed the canonical relation-analysis code. The
previous common validity mask required finite galaxy-centre velocities, even
for the primary relation between stellar metallicity and internal stellar
velocity dispersion. Those centre velocities are only needed for the separate
orbital-dispersion comparison. The code now builds the primary sample from
metallicity and stellar mass validity, and applies a finite-centre-velocity
mask only to the orbital comparison. It also reports the orbital-check sample
count explicitly. This prevents missing orbital velocities from needlessly
selecting the primary internal-dispersion sample. The canonical catalogues
have not been produced yet, so this adjustment has not generated results; no
tests were run.

## FoF dependency recheck (2026-10-21 07:13 KST)

Slurm still reports 00977 (430418_1), 00972 (430610_0), and 00980
(430611_2) as `RUNNING`; 00984 (430622_[3]) remains `PENDING (Resources)`.
NewGalFinder and extraction tasks remain `PENDING (Dependency)`. CPU totals for
00972 and 00980 rose from 05:06:09/05:07:36 at 07:11 to 05:12:23/05:13:51
at 07:13. However, all FoF catalogue and member-file sizes and log timestamps
remain unchanged from 06:42:34 (00977), 06:49:23 (00972), and 06:52:48
(00980). Exact checks also confirm that no `GALFIND.DATA`, satellite-property,
or relation-result file exists for snapshots 00972, 00977, 00980, or 00984.
The active outputs remain partial and were not consumed.

## Rank-7 process sample (2026-10-21 07:15--07:16 KST)

I sampled the two exact rank-7 PIDs from their startup logs on `syn102` once,
then again 65 seconds later. PID 1037709 (00972) stayed in `Rl` state at
96 per cent CPU and its process CPU time advanced from 01:02:25 to 01:03:27;
PID 1037703 (00980) stayed in `Rl` at 99 per cent CPU and advanced from
01:04:35 to 01:05:40. Their last application log line is still `passed Tree
building`. In `opfof.c`, the next stage scans each local particle and calls
`new_fof_link`/`pnew_fof_link`; these ranks each hold about 54 million
particles. This confirms active CPU computation in the likely expensive
particle-linking stage, rather than an idle process. Slurm still shows the
FoF jobs running and their dependent galaxy-finder tasks pending, so I left
the jobs untouched and did not use partial catalogues.

## Per-rank CPU check (2026-10-21 07:17 KST)

The 00977 FoF task remains `RUNNING`; a read-only sample of all ten exact MPI
PIDs from its startup log across `syn01`, `syn03`, and `syn04` showed every
rank in `Rl` state at 99.4--99.5 per cent CPU. For the 00972 and 00980 tasks,
rank 7 on `syn102` was sampled at 07:15 and 07:16: its CPU times advanced by
about one second per wall-clock second while it stayed at 96--99 per cent
CPU. These samples confirm that all three FoF tasks are doing CPU work even
though their application logs and catalogue files have not advanced. The
galaxy-finder and extraction jobs remain dependency-pending; no partial FoF
outputs were consumed or altered.

## FoF worker recheck (2026-10-21 07:21 KST)

The fresh Slurm query still shows all three FoF jobs running and every linked
NewGalFinder/extraction task dependency-pending. A second exact-PID sample
shows all ten 00977 MPI workers still in `Rl` at 99.4--99.5 per cent CPU; the
two tracked rank-7 processes for 00972 and 00980 remain in `Rl` at 96 and 99
per cent CPU. Their process CPU clocks increased nearly one-for-one with wall
time since the prior sample. The FoF catalogues/member files and application
logs remain unchanged. The 00977 bottom-face files also have no mtimes later
than 04:27:58, so the inferred current contact scan has not yet produced a
new visible write. No process was interrupted; these outputs are still
unvalidated partials.

## FoF output progress observed (2026-10-21 07:22 KST)

The 07:22 Slurm check still showed 00972, 00977, and 00980 running, with
downstream galaxy-finder/extraction tasks dependency-pending. Unlike the
previous checks, snapshot 00980 has now advanced: its halo/member files grew
to 198,028/12,065,760 bytes at 07:21:25, and its log was updated through
07:22:33. The new log lines show rank 7 reading input slab 00071 and then
multiple ranks passing tree construction. This is real output progress, but
the catalogue remains partial while the FoF task runs. Snapshot 00972 remains
at 172,060/10,862,592 bytes from 06:49:23, and 00977 remains at
31,404,068/52,130,289,832 bytes from 06:42:34. No partial catalogue has been
used.

The last 00977 log block shows every MPI rank has reported `saved isolated
halos`, but none has yet reported the following `saved Bottom Faced Halo`
message. In `opfof.c`, the next operation is `WriteBottomFaceContact`, which
scans each rank's halo/member lists before printing that completion message.
Together with all ten ranks being CPU-runnable at the 07:17 sample, this makes
the bottom-face contact scan the likely current 00977 stage. This is inferred
from the sequential log and source; the executable has no progress counter
inside that scan. The 00972/00980 rank-7 tasks remain earlier, in the
per-particle linking pass after tree construction.

## FoF output progress recheck (2026-10-21 07:24 KST)

A fresh Slurm query at 07:24 shows FoF tasks 430418_1 (00977), 430610_0
(00972), and 430611_2 (00980) still running; 00984 remains pending for
resources and the galaxy-finder/extraction tasks remain dependency-pending.
The 00972 FoF halo and member files have advanced to 197,548 and 12,046,608
bytes, with filesystem timestamps at 07:24:13. Its log now shows all ranks
through tree construction and several ranks completing isolated-halo writes.
Snapshot 00980 remains at 198,028 and 12,065,760 bytes (07:21:25), with its
latest log showing ranks reaching isolated-halo writes. Snapshot 00977 remains
at 31,404,068 and 52,130,289,832 bytes (06:42:34). None has a FoF completion
marker, and no `GALFIND.DATA` exists for these outputs; all catalogue files
remain partial and are excluded from analysis.


## FoF progress recheck (2026-10-21 07:26 KST)

Slurm still reports FoF tasks 00972, 00977, and 00980 as `RUNNING`; 00984
remains pending for resources, and galaxy-finder/extraction jobs remain
pending on dependencies. `sstat` shows accumulated CPU times of 06:17:53 for
00972 and 06:19:41 for 00980, with maximum RSS values of about 58 and 61 GiB,
respectively. The 00972 log file timestamp advanced to 07:25:40, but its size
and last lines are unchanged from the 07:24 sample; all FoF catalogue/member
file sizes are unchanged from the preceding check. The 00980 and 00977 logs
and outputs also show no new writes. No completion marker or `GALFIND.DATA`
exists for these outputs, so the current files remain partial and are not used
for satellite measurements.

## FoF rank and output check (2026-10-21 07:27 KST)

A fresh queue check still shows FoF 00972, 00977, and 00980 running; all
NewGalFinder and extraction jobs remain dependency-pending. The two exact
rank-7 processes on `syn102` are runnable and actively using CPU: snapshot
00972 process 1037709 used 1:14:32 CPU during 1:17:34 elapsed, and snapshot
00980 process 1037703 used 1:17:07 CPU during 1:17:34 elapsed. The 00980 log
advanced at 07:27:36 with ranks 1 and 4 completing isolated-halo writes, but
its catalogue files remain at their 07:21 sizes. Snapshot 00972's latest log
and FoF output remain at their earlier sizes; rank 7 is still computing. The
00977 files and log also remain unchanged. No FoF completion marker or
NewGalFinder catalogue is present, so no downstream catalogue extraction or
science analysis can yet run.

## Satellite-extraction integrity gate review (2026-10-21)

A static review of `extract_satellite_properties.py` confirms the extraction
stage checks repeated stellar IDs within each galaxy, assignments of the same
ID to multiple galaxies within a FoF host, and reuse of an ID across different
FoF hosts. Cross-host ID reuse aborts extraction and removes the temporary
CSV; the final CSV is published only after the complete file and global ID
check succeed. The batch script requires a non-empty `GALFIND.DATA` input and
refuses to overwrite an existing satellite CSV. These safeguards are ready
for the Darwin-1 NewGalFinder outputs; no catalogue was available to exercise
them in this review, and no tests were run.

## FoF compute check (2026-10-21 07:29 KST)

The queue still shows FoF tasks 00972, 00977, and 00980 running and all
NewGalFinder/extraction tasks pending on dependencies. Catalogue and member
file sizes have not changed since the previous recorded write. The exact
rank-7 processes on `syn102` remain runnable: 00972 (PID 1037709) has used
1:16:07 CPU in 1:19:12 elapsed, and 00980 (PID 1037703) has used 1:18:44
CPU in 1:19:12 elapsed. Their logs have no lines beyond the 07:25:40 and
07:27:36 updates, respectively; 00977's files and log remain at the 06:42:34
write. These observations indicate continuing CPU work for the two sampled
FoF ranks, with no new catalogue output yet. No partial catalogue was passed
to NewGalFinder or the science analysis.

A file-only follow-up later in the 07:29 check found one additional 00980 log
update at 07:29:21: rank 2 completed its halo-boundary check. The FoF output
file sizes had not changed by that check. This is application-stage progress,
not a completed catalogue; the prior Slurm query was less than 60 seconds old,
so I did not query the scheduler again.

## FoF boundary-stage progress (2026-10-21 07:30 KST)

Slurm continues to show FoF tasks 00972, 00977, and 00980 as running; 00984 is
pending for resources and NewGalFinder/extraction remain dependency-pending.
The latest logs show new per-rank work: 00972 rank 1 completed its halo
boundary check at 07:30:29, and 00980 rank 6 did so at 07:30:00. Their
catalogue/member file sizes have not advanced since the previous recorded
writes. Snapshot 00977's log and outputs are still unchanged from 06:42:34.
The observed rank messages indicate the two newer tasks are advancing through
FoF halo assembly, but none has reached the wrapper's completion marker; their
outputs remain partial and are not used for galaxy finding or analysis.

## FoF rank progress recheck (2026-10-21 07:31 KST)

The 07:31 Slurm query still shows 00972, 00977, and 00980 FoF tasks running;
00984 remains resource-pending and downstream jobs remain dependency-pending.
Snapshot 00980's log advanced to 07:31:26: rank 3 completed its halo-boundary
check after ranks 1, 2, 4, and 6. Its halo/member files are still unchanged
from 07:21:25. The 00972 log and files have not advanced since 07:30:29 and
07:24:13, respectively; 00977 remains unchanged since 06:42:34. No task has
printed `opfof complete`, so no FoF output has passed the gate for NewGalFinder.

## 00977 exact-rank activity and stage check (2026-10-21 07:32 KST)

I re-read the startup log to obtain all ten exact MPI PIDs for task 430418_1,
then sampled only those PIDs on their allocated nodes. All ten remain runnable
at 99.4--99.5 per cent CPU, with 10:29--10:30 CPU time over 10:33 elapsed.
The latest application log block has every rank saved its isolated halos for
the current slab, with no subsequent `saved Bottom Faced Halo` line. The
source calls `WriteBottomFaceContact` immediately after that point. Its ten
contact files have no mtimes later than 04:27:58, and the main FoF files and
log remain at their 06:42:34 sizes/timestamps. This is consistent with an
active, expensive contact-processing pass; the source has no progress counter
inside the function. Since every rank is still consuming CPU, I left the task
running. It has not reached FoF completion and its partial outputs remain
unusable.

## Storage headroom for active FoF outputs (2026-10-21 07:34 KST)

The cluster reports about 235 TiB available on the shared `/gpfs` filesystem;
`mmlsquota` shows no finite per-user block or file limit for `scratch`. The
current multi-gigabyte FoF outputs therefore are not blocked by the user's
quota or current filesystem free space. This check does not estimate the final
catalogue volume. The 07:34 Slurm query still showed all three FoF jobs active,
and a targeted sample of all ten 00977 MPI ranks showed each using about
99.4--99.5 per cent CPU. The latest 00977 log is still after isolated-halo
writes and before the subsequent bottom-face-contact completion messages.

## 00977 compute-versus-I/O diagnostic (2026-10-21 07:36 KST)

The 07:35:44 Slurm query still showed 00977, 00972, and 00980 FoF tasks
running, with all galaxy-finder/extraction tasks waiting on dependencies. I
sampled `/proc/<pid>/io` for the ten exact 00977 rank PIDs, waited 60 seconds,
and sampled them again. Their write-syscall counters did not increase, while
all ten process CPU times advanced by about one minute and remained at
99.4--99.5 per cent CPU. This points to computation rather than a storage
write wait during that interval; it does not identify a specific instruction.
The log still ends after isolated-halo saves and before `WriteBottomFaceContact`
reports completion. The 00977 contact/catalogue files did not change. By
contrast, the 00972 and 00980 logs advanced to 07:34:53 and 07:35:27,
respectively, while their FoF catalogue/member files remain at their earlier
sizes. No FoF completion marker or NewGalFinder catalogue is present.

## 00977 FoF per-rank slab position (2026-10-21 07:38 KST)

A source-and-log cross-check clarifies the 00977 progress position. The wrapper
passes `nfile=100`; `opfof.c` divides this into ten sequential files per MPI
rank (`nfile/mid`, with ten ranks). The latest read records are file indices
00004, 00014, ..., 00094, so each rank is processing its fifth assigned file
(out of ten), rather than the fifth file out of all 100. The latest tree,
halo-save, and contact-save messages correspond to the preceding iteration;
for the current file, all ranks have saved isolated haloes but have not logged
completion of bottom-face-contact handling. This establishes that the job is
midway through its per-rank input ranges, but the very different slab sizes and
unreported function-level progress prevent a reliable completion-time
estimate. The task remains active and its outputs remain partial.

## 00977 repeated CPU and syscall sample (2026-10-21 07:39 KST)

Between the 07:36 and 07:39 targeted samples, all ten 00977 MPI ranks stayed
in runnable state at 99.4--99.5 per cent CPU, and each process CPU clock
advanced by about three minutes. Their `/proc/<pid>/io` write-syscall counts
were identical in both samples, and no FoF/contact output file changed. The
07:38 Slurm query still showed the job running. Thus the sampled interval was
CPU-active without observable file writes; whether the ranks are inside the
contact scan, list compaction, or a busy MPI wait is not established by these
counters. No interruption or restart was made.

## 00977 boundary-check log update (2026-10-21 07:45 KST)

Although the FoF catalogue and member file remain at their 06:42:34 sizes, the
00977 application log has advanced: at 07:44:46 rank 8 completed the halo
boundary check for its current input file. The other recent rank messages show
this slab's tree and isolated-halo steps; no `saved Bottom Faced Halo` or
`opfof complete` marker follows yet. The job was still `RUNNING` in the
07:45:02 Slurm query. The new log line confirms continued per-rank progress,
while the FoF output remains partial and unusable for NewGalFinder.

## 00977 FoF output growth (2026-10-21 07:46 KST)

The 07:46 Slurm query still showed FoF tasks 00972, 00977, and 00980 running;
00984 remained resource-pending and NewGalFinder/extraction remained
pendency-pending. Snapshot 00977 has now made visible output progress: its
halo catalogue/member file grew from 31,404,068/52,130,289,832 bytes at 06:42
to 32,901,988/62,901,606,600 bytes at 07:46:09. The application log also
advanced through 07:46:31; ranks 7 and 9 are reading files 00075 and 00095,
the sixth inputs in their respective ten-file ranges, and rank 2 has passed
tree construction for its current input. This confirms continued multi-rank
FoF progress and output growth, but the job has no completion marker yet, so
the catalogue remains partial. Snapshot 00972 and 00980 FoF files have not
changed since their prior recorded writes.

## Darwin-1 FoF status recheck (2026-10-21 07:47 KST)

The 07:47:34 Slurm query still showed FoF snapshots 00972, 00977, and 00980
`RUNNING`; snapshot 00984 remained pending for resources, and the NewGalFinder
array remained dependency-pending. The 00977 catalogue and member file were
unchanged from 07:46:09 (32,901,988 and 62,901,606,600 bytes), and its log tail
still ended with ranks completing tree building, without an `opfof complete`
marker. The 00972 and 00980 log tails show later per-rank halo and boundary
checks, but neither task has published a completion marker. All current FoF
outputs therefore remain provisional; none was passed to galaxy finding.
The simulation is Darwin-1; `Cocoa_analysis_20261020` is only the working
folder name.

## Darwin-1 FoF log advancement (2026-10-21 07:48 KST)

At the 07:48:51 check, Slurm still reported FoF tasks 00972, 00977, and 00980
as `RUNNING`; 00984 remained resource-pending, and NewGalFinder remained
blocked on its declared job dependencies. The 00977 application log was
updated at 07:47:50 and now records rank 6 reading input file 00065, followed
by all ten ranks passing tree building. Its FoF files had not grown since
07:46:09 and no completion marker was present. The 00980 log was updated at
07:47:33 with further halo-boundary checks. The 00972 log and outputs had no
new updates since the earlier check. These records show continuing activity
for 00977/00980 but do not establish completion; no partial catalogue was
used downstream.

## Darwin-1 FoF status poll (2026-10-21 07:50 KST)

The 07:50:14 Slurm query still showed 00972, 00977, and 00980 as `RUNNING`;
00984 remained pending for resources, and the NewGalFinder tasks remained
pending on dependencies. No `opfof complete`, MPI abort, or error marker was
found in the three active FoF logs. The three halo-catalogue/member-file pairs
were unchanged from the prior size check. Latest log writes were 07:37:47 for
00972, 07:47:50 for 00977, and 07:49:08 for 00980. Because the products are
still partial and no FoF task has completed, none was supplied to NewGalFinder.

## FoF input-header precision audit (2026-10-21 07:52 KST)

The 00977 FoF stderr contains a `RamsesType` size-mismatch warning. The Slurm
log shows that this task started at 21:00:10 on 2026-10-20, whereas the current
`opfof.exe` and object files were rebuilt at 22:22 that night. The running task
therefore cannot be audited solely from the current binary/source pair. The
current `opfof.c` documents that NewDD `.info` files contain only an ASCII
header and parses that header; it is not sufficient evidence for the older
executable already in memory.

I compared the 00977 NewDD header values with the original RAMSES
`output_00977/info_00977.txt`: NewDD writes `aexp=0.821128` versus
`0.821128359613652`, `H0=67.400002` versus `67.4000015258789`, and rounds
`omega_m`/`omega_l` from `0.314999997615814`/`0.685000002384186` to
`0.315000`/`0.685000`. These differences are small, but this comparison alone
does not establish how the older FoF executable consumed the header. I have
not interrupted the live task; its 00977 catalogue must receive an explicit
post-run consistency check before being accepted into the scientific sample.

## Darwin-1 FoF activity recheck (2026-10-21 07:53 KST)

At the 07:53:35 scheduler check, FoF tasks 00972, 00977, and 00980 remained
`RUNNING`; 00984 was still resource-pending, and NewGalFinder remained
pendency-pending. All three live application logs have advanced since the
07:50 check: 00977 rank 2 completed a halo-boundary check, 00972 rank 5
completed one, and 00980 rank 9 completed one. Their FoF catalogue/member
files have not changed since the prior check, and no completion marker is
present. This is application-log progress, not a completed catalogue. The header warning
is traced to the pre-fix ASCII-header reader; its measured rounding is
negligible and the inspected source diff leaves the FoF algorithm unchanged.
The 00977 output still awaits successful completion and structural validation
before acceptance.

### 00977 warning follow-up: source and header path (2026-10-21)

The pre-fix executable backup at `opFoF/opfof.exe.pre_ascii_header_fix_20261020`
was timestamped 20:29 on 2026-10-20, before task 430418_1 started at 21:00,
and contains the exact `RamsesType` warning strings found in that task's stderr.
A source diff between `opfof.c.bak` and the current `opfof.c` shows changes to
command-line arguments, input filenames, and the conditional read of a binary
`RamsesType` block; it shows no change to the FoF tree-building or linking
algorithm. The NewDD `.info` file ends after its ASCII header and has no
`SIZEOF(RamsesType)` field or binary struct payload. `read_head()` scans the
header's numeric fields into `double` values; the old mismatch branch skips the
absent binary block. Its warning text about degrading to float is therefore
not an accurate description of the ASCII fallback. The six-decimal header
rounding is confirmed by comparison with the original RAMSES metadata (for
00977, relative differences are about $4.4\times10^{-7}$ in scale factor and
$7\times10^{-9}$ in $H_0$; density parameters differ by less than
$8\times10^{-9}$ fractionally). This resolves the warning as a small header
precision/reader-version issue rather than evidence of changed FoF membership
logic. Final acceptance still requires successful job completion and output
validation.

## Darwin-1 FoF check after header audit (2026-10-21 07:56 KST)

Slurm at 07:56:33 still reported FoF tasks 00972, 00977, and 00980 as
`RUNNING`; 00984 remained resource-pending and NewGalFinder tasks remained
dependency-pending. The 00972 log advanced at 07:56:17 with rank 8 completing
another halo-boundary check. The 00977 and 00980 logs had not advanced since
07:52:55 and 07:53:07, respectively, but Slurm still reported both tasks as
running. All three output pairs remained unchanged and incomplete; no partial
catalogue was used.

## Darwin-1 FoF quiet-log recheck (2026-10-21 07:57 KST)

At 07:57:56, Slurm continued to report FoF 00972, 00977, and 00980 as
`RUNNING`; 00984 remained pending for resources, and all NewGalFinder tasks
remained dependency-pending. The 00977 and 00980 logs had not advanced since
07:52:55 and 07:53:07; 00972's last log write was 07:56:17. The three FoF
output pairs were unchanged, and no `GALCATALOG.LIST` or `GALFIND.DATA` file
existed for these snapshots. I will check Slurm accounting next to distinguish
long-running computation from a quiet or stalled task; current `RUNNING` state
alone is not proof of active CPU progress.

## Darwin-1 FoF log advancement (2026-10-21 08:01 KST)

The 08:01 Slurm query still showed FoF 00972, 00977, and 00980 as `RUNNING`;
00984 remained pending for resources and NewGalFinder remained dependency-
pending. The 00977 log advanced at 08:00:48 with another rank completing halo
boundary checks, and 00972 advanced at 08:00:15 with rank 9 completing one.
The 00980 log has not advanced since 07:53:07. None of the FoF output files
has changed since their earlier partial writes, and the NewGalFinder outputs
remain absent. The corrected `sstat` step query is still pending.

## Darwin-1 FoF live-step accounting (2026-10-21 08:02 KST)

The corrected `sstat` query at 08:02:45 returned the MPI job steps rather than
only their batch wrappers. Reported fields were: 00977 (step 430419.0),
AveCPU=1-12:37:34, MaxRSS=95,300,156 KiB, MaxDiskRead=109,989,125,171 bytes,
and MaxDiskWrite=79,200,291,586 bytes; 00972 (430610.0), AveCPU=09:16:55,
MaxRSS=68,728,524 KiB, MaxDiskRead=47,650,851,731 bytes, and MaxDiskWrite=
5,624,694,333 bytes; 00980 (430611.0), AveCPU=09:19:26, MaxRSS=69,326,308 KiB,
MaxDiskRead=47,924,259,643 bytes, and MaxDiskWrite=5,533,299,841 bytes.
`MaxDiskRead` and `MaxDiskWrite` are per-step maxima, not cumulative I/O totals.
FoF output files remain partial, so no NewGalFinder task has been released.


## Darwin-1 CPU-accounting recheck (2026-10-21 08:05 KST)

A second `sstat` sample at 08:05:17 showed the `AveCPU` counters rising from
the 08:02:45 sample: by 00:08:23 for 00977, 00:12:32 for 00972, and 00:12:34
for 00980. This confirms CPU use despite the quiet 00980 log. The `MaxRSS`,
`MaxDiskRead`, and `MaxDiskWrite` fields were unchanged; the latter two are
per-task maxima and do not measure cumulative I/O. Logs and FoF output files
were unchanged at the 08:05 check, and no completion marker was present.

## Darwin-1 FoF CPU recheck (2026-10-21 08:06 KST)

At 08:06:30, `sstat` reported AveCPU=1-12:50:00 for 00977,
09:35:30 for 00972, and 09:38:04 for 00980. Relative to the 08:05:17
sample, these counters increased by 00:04:03, 00:06:03, and 00:06:04 over
73 seconds. MaxRSS and the per-task maximum disk counters did not change.
The application logs and catalogue/member files also remained unchanged at
this check. The rising CPU counters show ongoing computation, but the output
catalogues remain partial and NewGalFinder has not started.

## Darwin-1 snapshot 00972 FoF output growth (2026-10-21 08:07 KST)

Snapshot 00972 made its first large visible output update at 08:06:43: the FoF
halo file grew from 197,548 to 7,603,116 bytes, and the member-particle file
from 12,046,608 to 5,255,295,248 bytes. Its log shows ranks 1--9 sending
isolated-halo data to rank 0 and reporting that their isolated halos were
saved; rank 0 was preparing to write 6,869,691 halos. Slurm still reported
00972 `RUNNING` at 08:07:05, with no `opfof complete` marker, so the new files
are still partial. The 00977 and 00980 products have not changed, and all
NewGalFinder jobs remain dependency-pending.

## Darwin-1 snapshots 00972/00980 enter output-write phase (2026-10-21 08:09 KST)

At 08:08:00, snapshot 00980's FoF halo file grew from 198,028 to 7,546,252
bytes and its member-particle file from 12,065,760 to 5,114,887,136 bytes.
Snapshot 00972's corresponding files remain at 7,603,116 and 5,255,295,248
bytes after their 08:06:43 increase. Both logs show rank 0 preparing to write
about 6.8 million isolated haloes after the other ranks sent their data and
reported saving their portions. No `opfof complete` marker or NewGalFinder
output exists yet. The 08:09:15 Slurm query still showed both tasks running;
00977 also remained running, 00984 was resource-pending, and NewGalFinder
remained dependency-pending. Treat all three FoF products as partial.

## Darwin-1 FoF follow-up (2026-10-21 08:10 KST)

At 08:10:32, all three FoF tasks were still `RUNNING`. The 00977 log advanced
to rank 6 completing a halo-boundary check; 00972 and 00980 logs now show
ranks 1--9 reporting that they sent and saved isolated-halo data to rank 0.
The 00972 and 00980 catalogues/member files remain at their newly grown sizes
from 08:06:43 and 08:08:00, respectively; 00977's outputs remain unchanged
since 07:46:09. No `opfof complete` or error marker is present, so the apparent
root-write phase is not yet complete and none of the three catalogues is
released to NewGalFinder.

## Darwin-1 FoF root-write wait recheck (2026-10-21 08:12 KST)

At 08:12:02, Slurm still reported all three active FoF tasks as `RUNNING`;
00984 remained resource-pending and NewGalFinder remained dependency-pending.
The 00972 and 00980 logs still end with ranks 1--9 saving their isolated
haloes, while rank 0 is in its output phase; neither log nor output file has
changed since 08:06:43/08:08:00. Snapshot 00977's output and log are likewise
unchanged since 07:46:09/08:07:38. No completion or error marker is present.
At 08:12, `sstat` AveCPU had risen since 08:06:30 by 00:18:41 (00977),
00:27:55 (00972), and 00:28:02 (00980); MaxRSS and the per-task maximum disk
counters were unchanged. Thus the jobs continue consuming CPU, but there is no
new catalogue output to validate yet.

## FoF rank-0 serialization diagnosis (2026-10-21 08:13 KST)

Inspection of `opFoF/Treewalk.fof.ordered.c::WriteIsolatedHalo` explains the
quiet final phase for 00972 and 00980. Rank 0 first loops over its local halo
array and serially writes each qualifying halo record and its member particles;
it then receives the other ranks' halo/particle payloads and serially writes
those to the same two files. Only after all payloads have been written does
rank 0 print `have written all the isolated halo data`. The other ranks' `has
saved isolated halos` lines follow their sends and do not mean rank 0 has
finished. The current log pattern is consistent with rank 0 still serializing
the large output, and the rising `AveCPU` counters support continued compute.
This source-level diagnosis does not substitute for the completion marker or
validation of the resulting catalogue.

## Darwin-1 snapshot 00977 FoF output resumes (2026-10-21 08:14 KST)

After remaining unchanged since 07:46, snapshot 00977's FoF files grew at
08:13:50: the halo catalogue increased from 32,901,988 to 35,515,236 bytes,
and the member-particle file from 62,901,606,600 to 69,901,768,904 bytes.
The latest log block reports rank 0 checking a slab with 5,735,826 halos and
preparing to write its isolated-halo data; several other ranks have sent their
payloads, while not all ranks have reported completion yet. Slurm still showed
00977 `RUNNING` at 08:14:07. Tasks 00972 and 00980 are likewise `RUNNING` in
their isolated-halo output phase; their files have not changed since 08:06:43
and 08:08:00. No FoF completion marker is present and NewGalFinder remains
dependency-pending.

## Darwin-1 FoF status and CPU recheck (2026-10-21 08:17 KST)

The 08:17:53 Slurm query still showed FoF 00972, 00977, and 00980 as
`RUNNING`, with 00984 resource-pending and NewGalFinder dependency-pending.
Since the 08:15:48 sample, AveCPU rose by 00:06:31 for 00977, 00:09:45 for
00972, and 00:09:47 for 00980. The FoF output files and logs had not changed
since their last recorded writes (00972 at 08:06:43, 00980 at 08:08:00, and
00977 at 08:13:50). The processes continue consuming CPU, but the current
outputs remain partial and no NewGalFinder catalogue is available.

## Darwin-1 FoF activity recheck (2026-10-21 08:20 KST)

At 08:20:26, `sstat` AveCPU had risen since 08:17:53 by 00:08:28 for
00977, 00:12:39 for 00972, and 00:12:43 for 00980. The corresponding
MaxRSS and per-task maximum disk counters were unchanged. FoF catalogue and
member-file sizes, as well as all three application-log timestamps, were also
unchanged from the previous sample. The MPI steps continue using CPU, but the
rank-0 isolated-halo writes have not yet produced a new visible file update or
completion marker; NewGalFinder remains gated on FoF success.

## Darwin-1 FoF CPU recheck (2026-10-21 08:21 KST)

The 08:21:49 queue check still showed 00972, 00977, and 00980 `RUNNING`;
00984 remained resource-pending, with NewGalFinder tasks dependency-pending.
Between `sstat` samples at 08:20:26 and 08:21:49, AveCPU increased by
00:04:57 for 00977, 00:07:25 for 00972, and 00:07:25 for 00980. Their FoF
catalogue/member files remain at the previously recorded sizes; there is still
no completion marker or galaxy-finder output. The tasks continue consuming
CPU in the final isolated-halo output phase, but the scientific catalogues are
not yet ready.

## Darwin-1 FoF CPU-only progress sample (2026-10-21 08:24 KST)

The 08:24:17 `sstat` sample showed AveCPU increasing since 08:21:49 by
00:07:50 for 00977, 00:11:41 for 00972, and 00:11:44 for 00980. MaxRSS and
the per-task maximum disk counters were unchanged. The FoF files and logs also
remain unchanged since their last writes. This confirms ongoing CPU activity
but no visible output progress; the three catalogues remain partial, and
NewGalFinder is still waiting on successful FoF completion.

## Darwin-1 FoF CPU recheck (2026-10-21 08:25 KST)

At 08:25:36, the MPI-step `AveCPU` counters had risen since 08:24:17 by
00:04:21 (00977), 00:06:31 (00972), and 00:06:32 (00980). MaxRSS and the
per-task maximum disk counters were unchanged. The latest FoF file sizes and
log timestamps remain 7,603,116/5,255,295,248 bytes (00972, 08:06:43),
35,515,236/69,901,768,904 bytes (00977, 08:13:50), and
7,546,252/5,114,887,136 bytes (00980, 08:08:00). No completion marker is
present. CPU use continues, but all catalogue files remain partial.

## Darwin-1 identity and snapshot 00984 NewDD validation (2026-10-21 08:27 KST)

The simulation being analysed is **Darwin-1**. `Cocoa_analysis_20261020` is
the analysis-directory name; it is not the simulation name.

NewDD for snapshot 00984 finished at 22:34:36 KST on 2026-10-20. I checked
all 100 DM slabs and all 100 STAR slabs: every file exists, is non-empty, and
has a size divisible by the 168-byte particle-record size. Their totals are
1,073,741,824 DM records and 1,326,175,827 stellar records, matching the sums
from all 6,000 CPU entries in the NewDD log. The log contains the snapshot
completion marker. The corresponding FoF task (array 430622, index 3) is
pending in Slurm with no GPU requested; its reported reason changed from
compute resources to priority by 08:29 KST. Its NewGalFinder successor is
still dependency-pending.

At 08:26:47 KST on 2026-10-21, FoF tasks for 00972, 00977, and 00980 were
still running. Their output files and logs had not advanced since the last
recorded writes, and none had a FoF completion marker. NewGalFinder catalogues
and the new satellite-relation analysis therefore remain unavailable; no
physical conclusion is drawn from these in-progress products.

## Darwin-1 FoF input-chunk progress (2026-10-21 08:28 KST)

The opFoF source confirms that the 100 NewDD input slabs are split across ten
MPI ranks, with ten slabs processed sequentially by each rank. Thus the
repeated `P0 have written all the isolated halo data` messages mark completed
input chunks, not completion of the whole snapshot. At 08:28:08 KST, the
00977 log had five completed rank-0 chunks and was writing the next; the 00972
and 00980 logs each had one completed rank-0 chunk and were processing the
next. Their MPI jobs were still `RUNNING`, and the 00972/00980 `AveCPU`
counters had increased since the previous sample. All three FoF outputs are
still intermediate files and must not be passed to NewGalFinder until the
batch scripts report successful completion.

At 08:30 KST, the three FoF jobs remained `RUNNING`; the completed-chunk
counts and output-file sizes had not changed since 08:28. Their Slurm `AveCPU`
counters nevertheless continued to rise (00977: 1-14:06:52 to 1-14:10:03;
00972: 11:29:47 to 11:35:09; 00980: 11:32:48 to 11:38:11). This confirms
ongoing computation, but the current FoF catalogues are still incomplete.

At 08:33 KST, Slurm still reported all three FoF tasks as `RUNNING`. `AveCPU`
had risen to 1-14:20:38 (00977), 11:50:58 (00972), and 11:54:04 (00980).
However, their log timestamps, rank-0 completed-chunk counts (5, 1, and 1),
and member-file sizes were unchanged from the 08:30 sample. No completion or
error marker appeared; downstream NewGalFinder tasks remain gated.

## opFoF build-performance finding (2026-10-21 08:35 KST)

The current `/gpfs/kjhan/Others/Jaehyun/2026.Apr/opFoF/opfof.exe` was built
with Intel oneAPI DPC++/C++ Compiler 2025.3.0. The embedded compile command in
`opfof.o` has `-g` but no explicit `-O` option; `Rules.make` sets
`OPT = -DINTEL -g`. Intel documents that for this compiler, `-g` without an
explicit optimization level selects `-O0`, despite the compiler's normal
default being `-O2` ([Intel compiler guidance](https://www.intel.com/content/www/us/en/developer/articles/technical/compiler-optimization-debugging-trade-off.html)).
This is a credible contributor to the long CPU-heavy FoF write phase, although
it does not by itself prove the whole delay is caused by compiler settings.
The three running jobs keep using the binary already loaded by their MPI
processes; do not interrupt them or change their output directories. If a
retry becomes necessary, build a separate optimized binary with explicit
`-O2`, preserve this build and all partial products, and validate output
compatibility before adopting it for the catalogue pipeline.

At 08:36 KST, the same three FoF tasks were still `RUNNING`; their member-file
sizes and modification times remained 5,255,295,248 bytes (00972; 08:06:43),
69,901,768,904 bytes (00977; 08:13:50), and 5,114,887,136 bytes (00980;
08:08:00). No completion or error marker was present, and NewGalFinder
remained dependency-pending.

At 08:37 KST, I changed `/gpfs/kjhan/Others/Jaehyun/2026.Apr/opFoF/Rules.make`
to include explicit `-O2` while retaining debug symbols. I did not rebuild or
replace `opfof.exe`; the current jobs and the pending 00984 task therefore
continue to use the existing binary. Any future rebuild will use the optimized
setting, and its FoF output must still pass the same catalogue checks.

At 08:39 KST, a separate optimized executable was successfully built at
`/gpfs/kjhan/Others/Jaehyun/2026.Apr/opFoF/build_o2_20261021/opfof.exe`.
Its embedded compile record confirms `-g -O2`, its MPI library resolves, and
the original `opFoF/opfof.exe` remains unchanged. Compilation emitted existing
format-specifier warnings and assignment-in-condition warnings in the FoF
source. This candidate has not been run or validated against the current
catalogue; do not submit it until its output compatibility has been checked.

## Downstream catalogue handoff review (2026-10-21 08:40 KST)

I rechecked the queued stage scripts. NewGalFinder requires non-empty FoF
catalogue and membership files, then requires non-empty `GALCATALOG.LIST` and
`GALFIND.DATA` outputs before reporting success. The extractor consumes that
same `GALFIND.DATA` path and refuses to overwrite its CSV; the relation job
waits for the extracted CSV and also refuses to overwrite its output. These
filename and completion guards match across the scripts, so no handoff change
is needed. The NewGalFinder Slurm script requests CPUs and memory but contains
no explicit `--gres=gpu` request; it does target the A100 partition.

The optimized-build warning review found that the active `BIG_MPI_Send` and
`BIG_MPI_Recv` call sites pass `MPI_BYTE`; those calls take the explicit
one-byte branch before the assignment-based branches. Thus the compiler's
assignment-in-condition warning is a latent problem for other datatype uses,
but it does not change the byte-count behavior at these current call sites.
The `size_t`/`%d` format warnings remain in snapshot-path and diagnostic
formatting; they were not changed in the candidate build. This static review
still does not establish output equivalence, so the optimized binary has not
been used for catalogue production.

## Darwin-1 FoF output advances (2026-10-21 08:42 KST)

The 08:41--08:42 checks show visible output progress in all three running FoF
jobs. The rank-0 completed-write marker count is now 2 for snapshots 00972 and
00980 (up from 1), and remains 5 for 00977 while that snapshot's next chunk is
in progress. Current member-file sizes are 5,867,741,376 bytes (00972),
70,825,547,976 bytes (00977), and 5,739,794,088 bytes (00980); the 00972 and
00980 files have grown since their earlier samples. Slurm `AveCPU` is
1-14:47:21 (00977), 12:30:52 (00972), and 12:34:05 (00980). All three jobs
were still `RUNNING`, and no completion or error marker had appeared. These
files remain partial and are not yet valid NewGalFinder inputs.

## Darwin-1 FoF output sample (2026-10-21 08:43 KST)

At 08:43:14 KST, Slurm still reported the FoF jobs for snapshots 00972
(430610_0), 00977 (430418_1), and 00980 (430611_2) as `RUNNING`; 00984
(430622_[3]) remained pending with reason `Priority`. A filesystem sample at
08:43:17 found the 00977 member file at 72,205,695,176 bytes, up
1,380,147,200 bytes from the 08:42 sample. Snapshot 00972 remained at
5,867,741,376 bytes and 00980 at 5,739,794,088 bytes. The 00972 and 00980
files therefore showed no growth in this sample; this alone does not establish
that either job is stalled or failed. Snapshot 00977's halo catalogue was
37,301,092 bytes; the 00972 and 00980 halo catalogues were 8,416,348 and
8,351,428 bytes, respectively. The outputs are still being written and must
not be passed to NewGalFinder until the writer and wrapper completion checks
succeed.

## Darwin-1 FoF unchanged-file follow-up (2026-10-21 08:45 KST)

At 08:45:37 KST, Slurm still reported snapshots 00972 (430610_0), 00977
(430418_1), and 00980 (430611_2) as `RUNNING`; snapshot 00984 (430622_[3])
remained pending for `Priority`. The FoF member and halo-catalogue file sizes
for 00972, 00977, and 00980 were unchanged from the 08:43 sample. The
`AveCPU` values for 00972 and 00980 had increased to 12:50:30 and 12:53:49,
with maximum RSS values of 71,716,416 KiB and 75,693,368 KiB, respectively.
Thus, unchanged output-file sizes in this interval do not by themselves show
that those two jobs are idle or failed. The 00977 task did not appear in this
`sstat` sample; Slurm still reported it as running. No downstream catalogue
step should consume these outputs before the FoF wrapper's completion checks
pass.

## Darwin-1 FoF active-rank check (2026-10-21 08:46 KST)

At 08:46:24 KST, the three FoF jobs were still `RUNNING` and 00984 remained
pending for `Priority`. The three output file pairs had not changed size since
08:43. However, the 00972 and 00980 `AveCPU` values rose to 12:54:14 and
12:57:34, with maximum RSS values of 72,083,008 KiB and 75,982,136 KiB. The
latest log tails show ongoing rank-level halo-boundary checks and isolated
halo saves (00977), additional boundary checks (00980), and active tree-building
output (00972). This is evidence of continued computation despite unchanged
files; no completion marker or wrapper success line has yet been observed.

## Darwin-1 FoF worker-log follow-up (2026-10-21 08:47 KST)

At 08:47:10 KST, Slurm continued to report jobs 430610_0 (00972), 430611_2
(00980), and 430418_1 (00977) as `RUNNING`; 430622_[3] (00984) remained
pending for `Priority`. The 00972 and 00980 `.out` logs had advanced to 08:47
and their latest lines report additional halo-boundary checks, although their
member-file sizes remained unchanged. The 00977 log and member file had not
advanced since 08:42:17. Slurm still reported 00977 as running, so this is a
lack of recent observable output rather than proof of job failure; preserve
the job and inspect its accounting again before taking any intervention.

## Darwin-1 FoF accounting check (2026-10-21 08:48 KST)

At 08:47:54 KST, Slurm still listed 430418_1 (00977), 430610_0 (00972), and
430611_2 (00980) as running and 430622_[3] (00984) as pending. `sacct`
reported the three parent tasks and their active steps as running with exit
code 0:0; no terminal job result was recorded. The 00972 and 00980 logs had
recent worker output through 08:47:13 and 08:47:09, respectively, including
new boundary-check records. The 00977 log remained unchanged at 08:42:17 and
its member file remained 72,205,695,176 bytes. No success, fatal-error, or
MPI-abort marker was present in the three jobs' current logs. The 00977 task
therefore remains unresolved but live; continue observing without cancelling
or restarting it.

## Darwin-1 FoF output-writer inspection (2026-10-21 08:49 KST)

The 08:48:45 scheduler check still showed all three active FoF tasks running
and 00984 pending. `sstat` gave current samples for 00972 and 00980 only:
maximum RSS 73,946,688 KiB / 76,537,668 KiB and `AveCPU` 13:05:22 /
13:08:45, respectively. Their member-file sizes remained unchanged. The 00977
step still had no `sstat` row and no log or file update since 08:42:17, while
Slurm listed it as running.

A read-only inspection of `opFoF/Treewalk.fof.ordered.c` found that its
non-`USE_MASTER` `WriteIsolatedHalo` path has rank 0 receive and serially write
isolated-halo records sent by the other ranks; rank 0 prints its completion
line only after it has received all senders' data and closed both output files.
The build rules inspected do not define `USE_MASTER`, so this path is a
plausible explanation for a long interval without a new rank-0 completion line
while a very large member file is being assembled. This is a source-code-based
explanation, not proof that it is the cause of the 00977 delay. Keep the job
running and use subsequent Slurm accounting and log/output timestamps to decide
whether it advances or reaches a terminal state.

## Darwin-1 FoF status follow-up (2026-10-21 08:50 KST)

At 08:49:57 KST, Slurm still showed 00972 (430610_0), 00977 (430418_1), and
00980 (430611_2) as running, with 00984 (430622_[3]) pending for `Priority`.
The 00972/00980 `AveCPU` samples rose to 13:11:19 and 13:14:44, and their
maximum RSS values were 75,189,824 KiB and 76,738,372 KiB. Their logs and
member-file sizes did not change since the 08:47 sample. The 00977 step still
had no `sstat` row, and its log and member file remained unchanged since
08:42:17. No terminal state or error marker appeared. These are continued
live-job observations, not catalogue-completion evidence; the rank-0 serial
writer remains one possible explanation for the output delay, not a verified
cause.

## Darwin-1 snapshot 00977 Slurm allocation audit (2026-10-21 08:50 KST)

`scontrol show job 430418_1` resolves array element 1 to concrete Slurm JobId
430419. At 08:50:43 KST it was `RUNNING`, with `Requeue=1`, `Restarts=0`, an
elapsed time of 11:50:35, and end time 2026-10-22 21:00:08. The allocation is
60 CPUs and 768 GiB across three A10-partition nodes (`syn[01,03-04]`); its
`ReqTRES` and `AllocTRES` list CPU, memory, and nodes, with no GPU GRES. This
confirms that the job is live and has not restarted, despite no new log or
member-file timestamp since 08:42:17. Use the concrete JobId (430419) when
querying job-step accounting in the next sample; the array-form identifier did
not expose its active step in the previous `sstat` query.

## Darwin-1 snapshot 00977 concrete-step accounting (2026-10-21 08:51 KST)

The concrete Slurm JobId 430419 resolves the accounting gap seen when querying
array-form ID 430418_1. At 08:51:22 KST, step 430418_1.0 was `RUNNING` across
three nodes with three Hydra bootstrap proxy tasks. `sstat -j 430419` returned
MaxRSS 95,300,156 KiB, AveRSS 56,258,440 KiB, and accumulated AveCPU
1-15:18:46. The parent job remained `RUNNING`, with 60 CPUs and 768 GiB
allocated and no restart. Thus the absence of recent application log output
is not evidence of an idle job; its step remains active and resource-accounted.
The 00977 member file and application log were still unchanged at
72,205,695,176 bytes and 08:42:17, respectively. Meanwhile, the 00972 and
00980 logs had advanced again to 08:50:47 and 08:50:17 with additional halo
boundary-check lines, though their member-file sizes remained unchanged. The
FoF wrapper has not yet reported completion for these snapshots.

## Darwin-1 FoF resource-accounting follow-up (2026-10-21 08:52 KST)

At 08:51:56 KST, Slurm still listed 00972 (430610_0), 00977 (array element
430418_1 / concrete job 430419), and 00980 (430611_2) as `RUNNING`; 00984
(430622_[3]) remained pending for `Priority`. Querying concrete job 430419
returned its active step with MaxRSS 95,300,156 KiB and accumulated AveCPU
1-15:20:39, up from 1-15:18:46 at 08:51:22. The 00972/00980 steps also showed
increased AveCPU, to 13:20:34 and 13:24:01. Their application logs and all
three member-file sizes were unchanged from the prior sample. This confirms
continued CPU work but not completion; rank-0 output closure and the wrapper's
success checks remain outstanding.

## Darwin-1 FoF live-job poll (2026-10-21 08:52 KST)

At 08:52:23 KST, Slurm still reported the 00972, 00977, and 00980 FoF tasks
as running and 00984 as pending for `Priority`. Querying concrete JobId 430419
again returned active step 430419.0, with MaxRSS 95,300,156 KiB and AveCPU
1-15:22:11 (up from 1-15:20:39 at 08:51:56). The 00972/00980 `AveCPU`
samples were 13:22:51 and 13:26:18, also higher than at the previous sample.
The member-file sizes remained 5,867,741,376 bytes (00972), 72,205,695,176
bytes (00977), and 5,739,794,088 bytes (00980). The jobs are still consuming
CPU, but the current FoF outputs remain incomplete until the application and
wrapper finish writing and validate them.

## Darwin-1 snapshot 00977 FoF output resumes (2026-10-21 08:54 KST)

The 00977 FoF log and member file advanced at 08:54:12 KST. Rank 4 completed
its halo-boundary check for 74,410,415 particles, sent its isolated-halo data
to rank 0, and logged that it had saved the isolated halos. The member file
increased from 72,205,695,176 bytes (08:52 sample) to 78,231,427,272 bytes,
a growth of 6,025,732,096 bytes. This directly confirms that the long quiet
interval included continued rank-by-rank processing and that rank 0 is writing
new data. In the current log tail, ranks 7, 8, and 9 have not yet reported
completion of this output round; the job remains `RUNNING`, so the file is
still partial and cannot yet be consumed downstream.

## Darwin-1 FoF post-write interval (2026-10-21 08:55 KST)

At 08:55:30 KST, all three FoF jobs (00972, 00977, 00980) remained `RUNNING`
and 00984 remained pending for `Priority`. The 00977 concrete step 430419.0
still showed MaxRSS 95,300,156 KiB; AveCPU increased from 1-15:28:16 at
08:54:13 to 1-15:32:30. The 00972 and 00980 AveCPU values also increased to
13:38:15 and 13:41:48. Since rank 4's 08:54:12 transfer, 00977's log and
output files had not changed size; the member file remained 78,231,427,272
bytes. No success or error marker was present. This is a further live-resource
sample, not completion evidence; do not launch the dependent galaxy finder
from these partial outputs.

## Darwin-1 FoF quiet-interval check (2026-10-21 08:56 KST)

At 08:56:03 KST, Slurm still reported 00972, 00977, and 00980 as `RUNNING`
and 00984 as pending for `Priority`. The 00977 active step's AveCPU increased
to 1-15:34:21 (from 1-15:32:30 at 08:55:30), while MaxRSS remained
95,300,156 KiB. Its log, halo catalogue, and member file had not advanced
since 08:54:12; the latest application messages still show rank 4 as the last
rank to finish this output round. The 00972/00980 AveCPU values also continued
to rise. This is another active-computation sample, but there is no new FoF
output completion evidence.

## Darwin-1 snapshot 00977 slow-rank phase diagnosis (2026-10-21 08:57 KST)

At 08:57:24 KST, Slurm still showed all three active FoF tasks running and
00984 pending for `Priority`. Job 430419.0 AveCPU had risen to 1-15:38:50;
its MaxRSS remained 95,300,156 KiB. The 00977 log and output files had not
advanced since 08:54:12, and the member file remained 78,231,427,272 bytes.
The log's current processing round (lines 1088--1137) shows all ranks passing
tree construction, but ranks 7, 8, and 9 have not yet printed their halo
counts or boundary-check results. In `opfof.c`, after the `passed Tree
building` message, each rank loops over its local particles and calls
`pnew_fof_link` before printing its halo count (lines 656--665). The current
input records contain 97,786,946 particles for rank 7
(`SN.00977.00075`), 88,726,609 for rank 8 (`SN.00977.00085`), and 99,824,420
for rank 9 (`SN.00977.00095`). The unchanged output and logs are therefore consistent with those
ranks still doing the CPU-intensive linking pass; this is supported by the
code path and rising CPU accounting, though no per-rank stack trace was
collected. The active job continues; do not restart it based on quiet output.

## Darwin-1 FoF slow-rank accounting (2026-10-21 08:59 KST)

At 08:58:46 KST, Slurm still listed the three FoF tasks 00972, 00977, and
00980 as `RUNNING`, with 00984 pending for `Priority`. The 00977 step's AveCPU
rose to 1-15:43:22 while MaxRSS remained 95,300,156 KiB. The 00972 and 00980
AveCPU samples were 13:54:27 and 13:58:05. The 00977 log and its output files
had not changed since rank 4's transfer at 08:54:12; the last visible stage
still has ranks 5, 2, 6, 1, 3, and 4 reporting boundary checks/transfers,
while ranks 7, 8, and 9 have not reported their halo counts for this round.
Given the source order (`pnew_fof_link` loop before halo-count logging) and
continued job CPU accounting, the likely state remains CPU-intensive linking
on the last large chunks; this remains an inference without rank-level traces.

## Darwin-1 snapshot 00977 rank-9 output progress (2026-10-21 08:59 KST)

At 08:59:39 KST, Slurm still listed the three FoF jobs as `RUNNING`; 00984
remained pending for `Priority`. The 00977 log advanced at 08:59:32: rank 9
reported 9,329,380 halos from 99,824,420 particles, passed its halo-boundary
check, and sent its isolated-halo data to rank 0. The halo-catalogue file grew
to 39,381,860 bytes and the member file to 83,729,905,864 bytes, up
5,498,478,592 bytes from the 08:54:12 sample. No success, MPI-abort, or fatal
error marker was found. In the current log, ranks 7 and 8 still have no
completion message for this round. The 00977 output remains partial while the
job is active.

## Darwin-1 FoF status at 09:00 KST (2026-10-21)

At 09:00:22 KST, Slurm still showed the 00972, 00977, and 00980 FoF jobs as
`RUNNING`; 00984 remained pending for `Priority`. Concrete step 430419.0 for
00977 had AveCPU 1-15:48:39 and MaxRSS 95,300,156 KiB. The 00972 and 00980
AveCPU samples were 14:02:21 and 14:06:01. The 00977 log and both output
files were unchanged since 08:59:32, when rank 9 completed its transfer; no
new completion or error marker appeared. The active jobs remain in progress,
with no snapshot yet ready for the dependent galaxy-finder chain.

## Darwin-1 FoF live-resource follow-up (2026-10-21 09:02 KST)

At 09:01:46 KST, jobs 430610_0 (00972), 430419 (array element 430418_1,
00977), and 430611_2 (00980) remained `RUNNING`; 430622_[3] (00984) remained
pending for `Priority`. The 00977 step's AveCPU rose to 1-15:53:17 from
1-15:48:39 at 09:00:22, with MaxRSS still 95,300,156 KiB. The other two
steps' AveCPU values rose to 14:09:16 and 14:12:58. The 00977 application log
and output files remained unchanged since 08:59:32; ranks 7 and 8 still had no
current-round completion lines in the log. No job completion or error was
reported, so the catalogue output remains partial.

## Darwin-1 FoF live poll (2026-10-21 09:02 KST)

At 09:02:34 KST, 00972 (430610_0), 00977 (array 430418_1 / JobId 430419),
and 00980 (430611_2) remained `RUNNING`; 00984 (430622_[3]) remained pending
for `Priority`. The 00977 step's AveCPU was 1-15:55:57 (up from 1-15:53:17 at
09:01:46), with MaxRSS 95,300,156 KiB. AveCPU for 00972 and 00980 rose to
14:13:15 and 14:16:58. The 00977 log and outputs remained unchanged since
08:59:32; rank 9 is the last current-round completion message, with no current
round completion from ranks 7 and 8. No completion or error marker was found.

## Darwin-1 FoF CPU-progress sample (2026-10-21 09:03 KST)

At 09:03:10 KST, the active FoF jobs for 00972, 00977, and 00980 remained
`RUNNING`; 00984 remained pending for `Priority`. The 00977 active step's
AveCPU reached 1-15:57:56, up from 1-15:55:57 at 09:02:34; MaxRSS remained
95,300,156 KiB. The 00972/00980 AveCPU values reached 14:16:12 and 14:19:56.
The 00977 application log and output files remained unchanged since 08:59:32,
with rank 9 still the last recorded transfer. This confirms continued CPU
accounting without new output or completion evidence.

## Darwin-1 FoF wait verification (2026-10-21 09:04 KST)

At 09:03:47 KST, the three FoF jobs remained `RUNNING`, and snapshot 00984
remained pending for `Priority`. The 00977 step's AveCPU increased to
1-15:59:58 (from 1-15:57:56 at 09:03:10); MaxRSS remained 95,300,156 KiB.
The 00972/00980 AveCPU samples also increased to 14:19:15 and 14:22:59. The
00977 log and output files remained unchanged since 08:59:32, with rank 9
still the last current-round completion message. This confirms the process is
still consuming CPU but produced no new output in this interval; no terminal
or error state was observed.

## Darwin-1 FoF per-snapshot log check (2026-10-21 09:04 KST)

At 09:04:25 KST the three FoF jobs were still `RUNNING`, while 00984 remained
pending for `Priority`. No output file had reached a newer timestamp. The
00972 log advanced through 09:02:06 and its last lines show ranks 4 and 5
reporting halo counts and boundary checks; the 00980 log advanced through
08:57:09 and likewise shows ranks 4 and 5 in boundary checks. Their member
files remain at 5,867,741,376 bytes (00972) and 5,739,794,088 bytes (00980).
For 00977, concrete step 430419.0 remains CPU-active (AveCPU 1-16:02:05,
MaxRSS 95,300,156 KiB) but its last log/output update remains 08:59:32, with
rank 9 the latest current-round sender. None of these snapshots has a
completion marker, so downstream galaxy finding remains gated.

## Darwin-1 FoF status and downstream gate (2026-10-21 09:05 KST)

At 09:05:06 KST, FoF tasks for 00972, 00977, and 00980 were still
`RUNNING`; 00984 remained pending for `Priority`. Their AveCPU samples rose
to 14:25:48 (00972), 1-16:04:22 (00977, concrete step 430419.0), and
14:29:35 (00980). The 00972 and 00980 logs reached rank-5 halo-boundary
checks (09:02:06 and 08:57:09), but their member files remained unchanged at
5,867,741,376 and 5,739,794,088 bytes. The 00977 log/output remained at the
08:59:32 sample (rank 9 sender), with
member file 83,729,905,864 bytes. No FoF completion marker exists, so no
snapshot yet satisfies the input gate for NewGalFinder.

## Darwin-1 FoF active-job check (2026-10-21 09:06 KST)

At 09:06:02 KST, FoF tasks for 00972, 00977, and 00980 remained `RUNNING`;
00984 remained pending for `Priority`. AveCPU had risen to 14:30:24 (00972),
1-16:07:26 (00977, JobId 430419.0), and 14:34:11 (00980). The 00977 log and
both files remained at the 08:59:32 state (39,381,860-byte halo catalogue;
83,729,905,864-byte member file). The 00972 and 00980 logs showed rank-5
boundary checks, with outputs unchanged from their earlier timestamps. No
success, fatal-error, or MPI-abort marker was present. No downstream catalogue
step can safely use these partial FoF outputs yet.

## Darwin-1 FoF status recheck (2026-10-21 09:07 KST)

At 09:06:35 KST, the FoF jobs for 00972, 00977, and 00980 remained `RUNNING`;
00984 remained pending for `Priority`. AveCPU samples rose to 14:33:06,
1-16:09:15 (JobId 430419.0), and 14:36:55, respectively. Application-log
mtimes and member-file sizes were unchanged: 00977 remained at 83,729,905,864
bytes with rank 9 the last current-round sender; 00972 and 00980 remained at
5,867,741,376 and 5,739,794,088 bytes, with rank 5 last in their current log
tails. No completion or error marker was observed. The outputs remain partial.

## Darwin-1 FoF verified wait (2026-10-21 09:07 KST)

At 09:07:03 KST, 00972, 00977, and 00980 were still `RUNNING`, with 00984
pending for `Priority`. AveCPU rose to 14:35:27, 1-16:10:49 (concrete job
430419.0), and 14:39:16, respectively. The 00977 log and member-file sizes
remained unchanged since 08:59:32; 00972 and 00980 also showed no member-file
growth since their earlier timestamps. There was no completion/error marker.
The output remains partial; this check verifies ongoing CPU activity only.

## Darwin-1 FoF verified wait (2026-10-21 09:08 KST)

At 09:08:17 KST, the 00972, 00977, and 00980 FoF tasks were still `RUNNING`;
00984 remained pending for `Priority`. AveCPU rose to 14:41:33 (00972),
1-16:14:55 (00977, concrete step 430419.0), and 14:45:24 (00980). The 00977
log and member file remained unchanged since 08:59:32 (83,729,905,864 bytes),
and no completion/error marker appeared. This remains an active calculation
with no newly written FoF data in the observed interval.

## Darwin-1 FoF live-job check (2026-10-21 09:09 KST)

At 09:09:35 KST, the FoF jobs for 00972, 00977, and 00980 remained
`RUNNING`, with 00984 pending for `Priority`. AveCPU values continued to rise:
14:47:59 (00972), 1-16:19:14 (00977, JobId 430419.0), and 14:51:52
(00980). The 00977 log and member output remained unchanged since 08:59:32,
with rank 9 still the latest current-round transfer; no completion or error
marker appeared. Output is still partial and downstream NewGalFinder remains
gated on a successful FoF wrapper exit.

## Darwin-1 FoF scheduler check (2026-10-21 09:10:55 KST)

Jobs 430610_0 (snapshot 00972), 430418_1 (00977), and 430611_2 (00980)
remained `RUNNING`; 430622_[3] (00984) remained `PENDING` with Slurm reason
`Priority`. The concrete 00977 step (430419.0) continued accumulating CPU
time (AveCPU 1-16:23:38; MaxRSS 95,300,156 KiB). No FoF output file had a
newer timestamp: the 00977 member file remained 83,729,905,864 bytes (last
write 08:59:32), 00972 remained 5,867,741,376 bytes (last write 08:41:45),
and 00980 remained 5,739,794,088 bytes (last write 08:38:04). The three job
logs last advanced at 08:59:32, 09:02:06, and 08:57:09, respectively. No
complete-output marker was present, so these files remain gated from
NewGalFinder. Slurm reported the pending 00984 allocation on partition
`a100` with a predicted start at 2026-10-22 06:10:13; the reason shown was
priority. The checked-in opFoF batch script currently names `a10` and requests
no GPU GRES, so the submitted job's partition and the current script differ
and should be reconciled before any future submission; active jobs were left
untouched.

## Darwin-1 downstream pipeline queued (2026-10-21 09:12 KST)

For each output, I submitted a separate dependency chain so the catalogue
work starts as soon as that output's FoF job succeeds: NewGalFinder, satellite
property extraction, then relation analysis. The job IDs are:

| Snapshot | FoF prerequisite | NewGalFinder | Extraction | Relation analysis |
|---|---:|---:|---:|---:|
| 00972 | 430610_0 | 431520 | 431524 | 431528 |
| 00977 | 430418_1 | 431521 | 431525 | 431529 |
| 00980 | 430611_2 | 431522 | 431526 | 431530 |
| 00984 | 430622_3 | 431523 | 431527 | 431531 |

Each job depends on successful completion of the preceding job for the same
snapshot. Before submission, I confirmed that the expected galaxy-finder,
property-catalogue, and relation-catalogue outputs do not already exist. The
batch scripts request CPU and memory resources but no GPU GRES. The NewGalFinder
script still specifies the `a100` partition; this reserves CPU/memory on those
nodes but does not request a GPU. No FoF jobs were modified or restarted.

At 09:13:40 KST, Slurm showed all four NewGalFinder jobs pending on their
respective `afterok` FoF dependencies, with extraction and relation-analysis
jobs likewise dependency-pending. The 00972, 00977, and 00980 FoF jobs were
still running; 00984 was still pending for `Priority`. For NewGalFinder job
431520, `ReqTRES` was `cpu=100,mem=1200G,node=2,billing=100` and contained no
GPU request, confirming that the CPU-only job is queued on the A100 partition
without reserving a GPU.

## Darwin-1 pipeline dependency verification (2026-10-21 09:14:24 KST)

Slurm still showed the FoF jobs for 00972, 00977, and 00980 as `RUNNING`; the
00984 FoF job remained `PENDING (Priority)`. Their measured CPU times continued
to increase: 15:11:53, 1-16:35:14, and 15:15:52, respectively. All four
NewGalFinder jobs remained pending on their individual FoF `afterok`
dependencies; the extraction and relation-analysis jobs remained dependency-
pending behind NewGalFinder and extraction. Thus the queued chains have not
started against incomplete FoF files.

## Darwin-1 FoF progress check (2026-10-21 09:16 KST)

All three active FoF jobs continued accumulating CPU time: 15:20:50 for 00972,
1-16:41:14 for 00977, and 15:24:52 for 00980. Their output files had not
changed since the 09:14:59 sample. The 00972 log shows ranks 2, 3, 6, 1, 4,
and 5 reporting halo counts and boundary checks, followed by rank 9; no wrapper
completion is present. In the latest 00977 round, ranks 8 and 9 have sent and
saved their isolated halos, while rank 7 remains without a halo-count result
after its tree-building message. In 00980, rank 9 has sent and saved its
isolated halos, while no wrapper completion is present. These remain active
calculations; the per-snapshot downstream jobs are still gated on FoF success.

## Darwin-1 scheduler recheck (2026-10-21 09:18:16 KST)

The three running FoF jobs remained active, with AveCPU increasing to 15:30:59
(00972), 1-16:48:01 (00977), and 15:35:03 (00980). The 00972 log had advanced
to 09:16:12, but the FoF member-file sizes for all three snapshots were
unchanged from the prior check; 00977 and 00980 logs also remained unchanged.
The 00984 task's pending reason changed from `Priority` to `Resources`. All
four NewGalFinder jobs still waited on their own FoF success dependencies,
with extraction and relation analysis queued behind them.

## Darwin-1 verified wait (2026-10-21 09:19:32 KST)

All three active FoF jobs remained `RUNNING`; the 00984 task remained pending
for `Resources`. AveCPU rose to 15:37:17 (00972), 1-16:52:15 (00977), and
15:41:23 (00980), while the FoF logs and member-file sizes had not changed
since 09:16:12, 09:11:32, and 09:12:16, respectively. The dependent
NewGalFinder, extraction, and analysis jobs remained queued on their success
gates. This confirms ongoing CPU activity without a new catalogue or a
completed FoF output.

The earlier Darwin-1 NewGalFinder run for 00967 provides a runtime reference:
its 333-GB `GALFIND.DATA.00967` was produced from a 395-GB FoF member file in
about 10 h 20 min, using five worker ranks with 60 OpenMP threads each. The
queued four-worker configuration uses 20 threads per worker and has a 48-h
limit. A linear inverse-thread-count extrapolation would be about 39 h for a
similar workload; this is only a rough planning estimate, since I/O, snapshot
size, and parallel scaling can differ. We should watch the first queued
NewGalFinder run's elapsed time before assuming the same limit is sufficient
for the remaining outputs.

## Darwin-1 verified wait (2026-10-21 09:22:41 KST)

Slurm continued to report 00972, 00977, and 00980 as `RUNNING`, with AveCPU
15:52:54, 1-17:02:43, and 15:57:05, respectively. Snapshot 00984 remained
pending for `Resources`; the per-snapshot downstream chains remained gated.
The three FoF member files had not changed from their 09:19:51 sizes. The
existing 00967 run indicates that NewGalFinder itself may take many hours, so
its 48-h limit will be reviewed against the first actual queued run before
changing any later submissions.

## Darwin-1 verified wait (2026-10-21 09:23:57 KST)

The 00972, 00977, and 00980 FoF jobs remained `RUNNING`, with AveCPU
15:59:09, 1-17:06:55, and 16:03:22. The 00984 job remained pending for
`Resources`. Their output files and logs still had the 09:19:51-or-earlier
timestamps recorded above, and all satellite-catalogue pipeline jobs remained
gated on successful FoF completion. CPU activity continues, but there is no
new FoF output or scientific result to report yet.

## Darwin-1 live-job check (2026-10-21 09:20:44 KST)

The three FoF jobs were still `RUNNING`; 00984 remained pending for
`Resources`. AveCPU continued to rise to 15:43:15 (00972), 1-16:56:15
(00977), and 15:47:23 (00980). At 09:20:44, all three application logs and
member files still had the same modification times and sizes recorded at
09:19:51. All downstream jobs remained dependency-pending. The FoF jobs are
active, but have not produced a new output sample or exited successfully.

## Darwin-1 FoF output growth (2026-10-21 09:14:59 KST)

The 00977 and 00980 FoF jobs wrote additional data since the previous sample.
For 00977, the halo catalogue grew to 40,397,668 bytes and the member file to
90,641,668,296 bytes (both last written 09:11:32); its log now shows rank 8
finished sending and saving its isolated halos. For 00980, the halo catalogue
grew to 13,921,988 bytes and the member file to 17,425,829,544 bytes (both
last written 09:12:16); rank 9 reported its isolated halos saved. Slurm still
reported both jobs as `RUNNING`, with no wrapper completion marker, so their
outputs remain incomplete. Snapshot 00972's files were unchanged from the
08:41:45 sample despite the job continuing to accrue CPU time. The dependency-
gated NewGalFinder jobs therefore remain the correct next step only after each
FoF wrapper exits successfully.

In the latest visible 00977 output round, rank 7 last reported `passed Tree
building` (log line 1096) and has not yet reported its halo count or sent its
isolated halos, while ranks 8 and 9 have completed those steps. This is
consistent with rank 7's larger 97.8-million-particle input taking longer in
the linking phase; aggregate CPU time is still rising, so the evidence does
not indicate a stopped job.

## Darwin-1 FoF wait verification (2026-10-21 09:25:22 KST)

The 00972, 00977, and 00980 FoF jobs remained `RUNNING`, with AveCPU
16:06:09, 1-17:11:36, and 16:10:25. Snapshot 00984 remained pending for
`Resources`. The four NewGalFinder arrays and all eight extraction/analysis
jobs were still dependency-pending. No current snapshot had yet passed the FoF
success gate, so no satellite catalogue or relation result was available.

## Darwin-1 pipeline and storage check (2026-10-21 09:26:45 KST)

The three active FoF jobs remained `RUNNING`; AveCPU reached 16:12:57 (00972),
1-17:16:10 (00977), and 16:17:15 (00980). Snapshot 00984 was again pending
for `Priority`; all downstream jobs remained dependency-pending. No FoF output
file had advanced past the timestamps listed above. The GPFS quota report
showed no user quota limit, and `/gpfs` reported 234 TB available, so storage
capacity is not currently blocking the planned catalogues. This is a storage
availability check only; it does not estimate the final four output sizes.

## Darwin-1 FoF status (2026-10-21 09:28:13 KST)

The 00972, 00977, and 00980 jobs remained `RUNNING`, with AveCPU 16:20:16,
1-17:21:04, and 16:24:36. Snapshot 00984 returned to pending reason
`Priority`. All per-snapshot NewGalFinder, extraction, and relation-analysis
jobs remained dependency-pending. Logs and member files were unchanged from
the previous stat check, so the active FoF jobs had not yet exposed another
completed output block.

## Darwin-1 FoF active-job check (2026-10-21 09:29:27 KST)

The 00972, 00977, and 00980 jobs remained `RUNNING`; their AveCPU values rose
to 16:26:21, 1-17:25:08, and 16:30:42. Snapshot 00984 was pending for
`Priority`. All four NewGalFinder jobs and their extraction/analysis chains
remained dependency-pending. No new FoF catalogue or member-file write was
available at this check.

## Darwin-1 verified wait (2026-10-21 09:30:52 KST)

Slurm again showed the 00972, 00977, and 00980 FoF jobs `RUNNING`, with AveCPU
16:33:22, 1-17:29:50, and 16:37:45. The 00984 FoF job remained pending for
`Priority`. All per-snapshot catalogue jobs remained behind their FoF
dependencies. No completion or downstream start was observed.

## Darwin-1 active FoF verification (2026-10-21 09:32:52 KST)

The 00972, 00977, and 00980 FoF jobs remained `RUNNING`, with AveCPU
16:43:17, 1-17:36:29, and 16:47:44. Snapshot 00984 remained pending for
`Priority`. The application logs and output files were unchanged from the
09:31:22 sample, and no completion or error marker appeared. All four
dependency-gated NewGalFinder jobs remained pending, with their extraction
and relation-analysis jobs queued behind them.

## Darwin-1 supplementary 00967 pipeline queued (2026-10-21 09:34:33 KST)

To make progress while the four primary FoF outputs are running, I queued a
separate supplementary extraction and relation analysis using the completed
00967 NewGalFinder file (`GALFIND.DATA.00967`, 332,983,471,392 bytes). It writes
distinct `satellite_properties_00967.csv` and `satellite_relations_00967.csv`
files and will not overwrite existing results. Job 431549 requests one CPU,
32 GB RAM, and no GPU GRES; Slurm shows it pending for `Priority`, with a
predicted start at 13:52:24 KST. This output is ancillary and will not be used
to change the four-snapshot primary sample unless its provenance and the
result are validated.

## Darwin-1 scheduler check (2026-10-21 09:37:32 KST)

The three primary FoF jobs remained `RUNNING`; AveCPU reached 17:06:16 (00972),
1-17:51:58 (00977), and 17:10:56 (00980). The 00980 MaxRSS sample was
78,150,940 KiB. The 00984 FoF job and supplementary 00967 job 431549 were
pending for `Priority`; the four primary NewGalFinder jobs and their
downstream chains remained dependency-pending. The 00972 files and log had
not advanced since 09:35:13, and job 431549 had not started (so its log files
were not yet created).

## Darwin-1 FoF output growth for 00972 (2026-10-21 09:35:12 KST)

Snapshot 00972 entered another output-writing round. Its FoF halo catalogue
grew to 12,880,988 bytes and the member file to 15,498,526,912 bytes. The log
records a new rank-0 halo count of 6,121,044 and shows ranks 5, 6, 2, 9, 3, 1,
and 4 sending and saving their isolated halo data. The job had no completion
marker at this sample, so the 00972 NewGalFinder dependency remains gated.

At 09:36:26 KST, the 00972 halo catalogue had grown further to 14,003,292
bytes and its member file to 17,685,520,576 bytes. The latest log block shows
rank 0's 6,121,044-halo count and ranks 5, 6, 2, 9, 3, 1, and 4 saving their
isolated halos. Ranks 7 and 8 had passed tree building but had not yet reported
halo counts in this block. The job had not reached its wrapper completion
marker, so NewGalFinder remained gated.

## Darwin-1 pipeline status (2026-10-21 09:36:05 KST)

Slurm still showed FoF jobs 00972, 00977, and 00980 as `RUNNING`, with AveCPU
16:59:07, 1-17:47:08, and 17:03:42. The 00980 MaxRSS sample rose to
77,878,556 KiB. Snapshot 00984 and supplementary 00967 job 431549 remained
pending for `Priority`; all NewGalFinder/extraction/analysis jobs for the
primary snapshots remained dependency-pending. The new 00972 output block had
not yet passed the FoF wrapper's success gate.

## Darwin-1 active jobs (2026-10-21 09:38:48 KST)

Slurm continued to show 00972, 00977, and 00980 as `RUNNING`, with AveCPU
17:12:31, 1-17:56:10, and 17:17:14. The 00980 MaxRSS sample was
78,470,428 KiB. Snapshot 00984 and supplementary job 431549 remained pending
for `Priority`, while all primary downstream jobs remained dependency-pending.
No newer FoF output appeared after the 00972 writes at 09:35:13.

## Darwin-1 extractor dtype correction (2026-10-21 09:42 KST)

The supplementary 00967 job 431549 failed immediately before reading the input:
`STAR_DTYPE` described the six-byte C alignment gap as a fifteenth data field,
which shifted the later format entries and made NumPy require 232 bytes while
the declared record size was 168 bytes. I removed that padding pseudo-field;
the explicit offsets retain the alignment, and the intended 15 fields now map
to the 168-byte record. This fixes the shared extractor used by the primary
snapshots as well as the supplementary run. I resubmitted 00967 as job 431591.
Its output is still ancillary pending input provenance and successful
extraction checks. No downstream science results have been produced by this
retry yet.

## Darwin-1 00967 schema provenance correction (2026-10-21 09:44 KST)

Job 431591 also failed before producing a catalogue. The 00967 `GALFIND.DATA`
file was created on 2026-09-25 with the earlier NewGalFinder build. Its existing
stellar reader identifies that legacy product as a 120-byte, NCHEM=3 record
with unreliable 32-bit particle IDs. The active NewGalFinder build was rebuilt
on 2026-10-21 from the current source and uses `NCHEM=9`, `NDUST=4`; its
particle and gas records have different sizes. The current extractor correctly
rejects the old 00967 stream when its halo sequence becomes inconsistent. I
removed only the 1.9-kB partial CSV left by the failed attempt. I will exclude
00967 rather than mix its incompatible layout or untrusted IDs into the
primary Darwin-1 sample. The four requested primary snapshots use the current
NewGalFinder build and remain the analysis target.

At 09:44 KST, FoF jobs for 00972, 00977, and 00980 were still `RUNNING`; 00984
was pending for priority. The existing outputs had not grown since their
previous samples (00972 at 09:35:13; 00977 at 09:11:32; 00980 at 09:12:16).
All four NewGalFinder jobs, extractions, and relation analyses remained
queued behind their successful FoF dependencies. No satellite results are
available from this run yet.

## Darwin-1 primary FoF check and legacy-job guard (2026-10-21 09:48 KST)

Slurm still reports the 00972, 00977, and 00980 FoF tasks as `RUNNING`; 00984
is pending for priority. All four NewGalFinder, extraction, and relation-analysis
jobs remain dependency-pending. The three current FoF output pairs and logs are
unchanged from 09:35:13 (00972), 09:11:32 (00977), and 09:12:16 (00980), so
there is no newer file-level progress to report. A `sstat` query of the batch
shell steps returned no useful CPU data (and no running batch step for 00977);
I make no process-activity inference from it and did not restart any task.

I added an explicit refusal at the start of `supplement_00967.sbatch` so this
legacy, ID-unreliable 00967 product cannot be resubmitted accidentally. It is
retained only as a record of the failed supplementary attempt; the four current
Darwin-1 snapshots remain the analysis sample.

## Darwin-1 FoF live-step resource sample (2026-10-21 09:49 KST)

A targeted `sstat --allsteps` sample confirms the three FoF MPI steps remain
active and accumulating CPU time despite unchanged output files: 00972 step
430610.0 reports 18:06:05 AveCPU and 75,522,524 KiB MaxRSS; 00977 concrete
step 430419.0 reports 1-18:32:15 AveCPU and 95,300,156 KiB MaxRSS; 00980
step 430611.0 reports 18:11:18 AveCPU and 84,132,636 KiB MaxRSS. Slurm still
lists all three array tasks as `RUNNING`; 00984 remains pending for priority.
This is evidence of live computation, not a completion estimate. The dependent
catalogue and science-analysis jobs remain held until FoF exits successfully.

## 00977 node-level FoF process check (2026-10-21 09:54 KST)

A targeted read-only `srun --overlap` sample on `syn04` found all three local
`opfof.exe 977 100` MPI ranks in runnable (`Rl`) state at 99.4--99.5% CPU,
with elapsed process time 12:54:05. Together with the rising Slurm CPU counters
for all three active snapshots, this confirms that the lack of new catalogue
writes does not reflect idle MPI ranks. No task was restarted or interrupted;
the FoF completion markers are still required before downstream processing.

## Darwin-1 FoF output advances for 00980 (2026-10-21 09:55 KST)

Snapshot 00980 made measurable output progress: its FoF halo catalogue grew
from 13,921,988 bytes at 09:12:16 to 14,692,036 bytes at 09:54:00, and its
member-particle file grew from 17,425,829,544 to 18,472,816,296 bytes. The
latest log shows rank P7 completed its boundary check and saved its isolated
halos after rank P9. Slurm still reports 00980 `RUNNING`; its final output
marker has not appeared. The 00972 and 00977 files remain unchanged at the
latest sample, but their Slurm MPI-step CPU counters continue to rise. All
catalogue consumers remain dependency-gated.

## Darwin-1 FoF progress and literature wording audit (2026-10-21 10:02 KST)

At the 10:02 KST scheduler check, FoF jobs `430418_1` (00977), `430610_0`
(00972), and `430611_2` (00980) were still `RUNNING`; 00984 remained
`PENDING (Priority)`. The NewGalFinder, extraction, and relation-analysis
chains remain dependency-pending. The NewGalFinder dependency for 00972 is
explicitly `afterok:430610_0`, so it will not consume the in-progress FoF
files. The 00977 log records all ten ranks repeatedly completing tree
building and saving isolated halos across successive passes; it is still
advancing and has no `opFoF complete` marker. The 00977 halo and member files
advanced to 41,270,668 and 92,222,079,936 bytes (10:00:25 KST). The 00972
and 00980 pairs remain at 14,003,292 / 17,685,520,576 bytes and 14,692,036 /
18,472,816,296 bytes, respectively. These are still in-progress products and
are not accepted as inputs to the satellite extraction.

The introduction was revised after checking the primary papers: Gallazzi et
al. (2021) report an increase in satellite metallicity with host mass at
fixed stellar mass below $10^{10.5} M_\odot$, with much of the central--
satellite difference tied to the higher passive fraction. Buck et al. (2019)
report nearly coincident satellite/field stellar-mass--metallicity and
stellar-mass--dispersion relations, alongside lower overall dispersions for
satellites and nearby dwarfs. The draft no longer claims a lower dispersion
at fixed stellar mass for the latter comparison.

## Darwin-1 00977 FoF live-process verification (2026-10-21 10:05 KST)

Slurm still reports 00977 FoF as `RUNNING` (array task `430418_1`, concrete
allocation `430419`). An overlapping read-only process sample on its three
allocated nodes found all ten `opfof.exe 977 100` MPI ranks in runnable
`Rl` state at 99.4--99.5 per cent CPU. The member and halo files have not
grown since 10:00:25 KST, and the log last advanced at 10:02:27 KST, but the
live rank sample confirms the job is actively computing. It was left running;
no restart or cancellation was made. Jobs 00972 and 00980 remain scheduler-
active, 00984 remains priority-pending, and all downstream NewGalFinder,
extraction, and analysis tasks remain dependency-gated.

## Darwin-1 FoF scheduler recheck (2026-10-21 10:06 KST)

The 10:06 KST scheduler poll still showed 00972 (`430610_0`), 00977
(`430418_1`), and 00980 (`430611_2`) running. The four catalogue-generation
and analysis chains remain dependency-pending; 00984 FoF remains
priority-pending. The 00972, 00977, and 00980 FoF file sizes are unchanged
from the preceding check. The 00977 compute ranks were directly confirmed
runnable and CPU-active at 10:05 KST, so no process was restarted or
cancelled.

## Darwin-1 00972/00980 FoF live-process verification (2026-10-21 10:07 KST)

A read-only process sample through the active allocation found five
`opfof.exe 972 100` ranks and five `opfof.exe 980 100` ranks in runnable
`Rl` state, each using 99.2--99.5 per cent CPU. This explains why their
output files and logs have not changed recently: both jobs are still in
compute phases. The 10:07 scheduler check also showed 00977 running, with its
log advancing at 10:06:37 KST; all three FoF tasks remain active and their
catalogues remain unvalidated until the completion markers and final checks.

## Darwin-1 FoF live status (2026-10-21 10:08 KST)

At 10:08 KST, Slurm still showed 00972 (`430610_0`), 00977 (`430418_1`),
and 00980 (`430611_2`) as `RUNNING`; 00984 remained priority-pending and
all downstream NewGalFinder/extraction/analysis tasks remained pending on
dependencies. A fresh process sample again found five ranks each for 00972
and 00980 in runnable state at 96.1--99.4 per cent CPU. Their output files
were unchanged. The 00977 log had last advanced at 10:06:37 KST, with no
completion marker; its output files also remained unchanged. None of the
partial FoF products was passed to the extractor.

## Darwin-1 00977 FoF active-compute recheck (2026-10-21 10:09 KST)

A fresh overlapping process sample through concrete Slurm allocation `430419`
again found all ten 00977 FoF ranks in runnable `Rl` state at 99.4--99.5
per cent CPU. Slurm showed 00972, 00977, and 00980 still running; the
catalogue-generation arrays remained dependency-pending and 00984 remained
priority-pending. The three output pairs were unchanged at this check, so the
FoF products remain provisional and no downstream job was started manually.

## Darwin-1 FoF recheck (2026-10-21 10:10 KST)

The 10:10 KST scheduler poll still showed FoF snapshots 00972, 00977, and
00980 running; 00984 remained priority-pending. NewGalFinder, extraction,
and relation-analysis jobs were still waiting on dependencies. All three
FoF log and output-file timestamps/sizes were unchanged since the 10:09
check. The most recent direct process sample (10:09 KST) found the 00977
MPI ranks actively using CPU. No FoF output is treated as complete yet.

## Darwin-1 host-sensitivity dependency repaired (2026-10-21 10:13 KST)

The pending host-mass sensitivity array `431211` was still attached to four
superseded extraction jobs (`430492_0`, `430427_1`, `430427_2`, `430477_3`),
while the active per-output chain uses extraction array `431524` (tasks 0--3).
I changed the pending sensitivity array dependency to `aftercorr:431524`, so
each sensitivity task now waits for the successful extraction of the
corresponding snapshot. `scontrol show job` confirms tasks 0 and 1 now report
`aftercorr:431524_* (unfulfilled)`. This changes only the pending sensitivity
job dependency; no running FoF task was altered.

## Darwin-1 obsolete duplicate analysis chains removed (2026-10-21 10:15 KST)

A dependency audit found that four earlier downstream chains were still pending
alongside the current chains, and several were tied directly to the same active
FoF jobs. They could therefore launch duplicate NewGalFinder/extraction runs
and race for the same snapshot output filenames. Before cleanup, `squeue`
showed all 12 old tasks pending on dependencies and all canonical galaxy and
satellite catalogue files absent. I cancelled only those pending obsolete
tasks: NewGalFinder `430423_1`, `430423_2`, `430491_0`, `430476_3`; extraction
`430427_1`, `430427_2`, `430492_0`, `430477_3`; and relation analysis
`430431_1`, `430431_2`, `430493_0`, `430479_3`. They are no longer in
`squeue`. No running FoF task was cancelled or changed. The current four
chains (`431520`--`431531`) remain pending on their corresponding FoF,
NewGalFinder, and extraction jobs. The host-mass sensitivity array `431211`
now uses `aftercorr:431524`, pairing each task with its current extraction.

## Darwin-1 00972 FoF output advancement (2026-10-21 10:16 KST)

The 10:15 KST Slurm check still showed FoF jobs 00972, 00977, and 00980
running; 00984 remained priority-pending. Snapshot 00972 advanced at 10:12:36:
its FoF halo catalogue grew to 14,781,532 bytes and member-particle file to
18,738,323,648 bytes. The log records P0 preparing to write 6,121,044
isolated halos and the ranks sending/saving their isolated-halo data. There is
still no `opFoF complete` marker and Slurm reports the job running, so this
is a progress checkpoint, not an accepted catalogue. The 00977 and 00980
FoF output pairs were unchanged; all current downstream jobs remain pending.

## Darwin-1 FoF live compute check (2026-10-21 10:17 KST)

The 10:16 KST scheduler poll showed 00972, 00977, and 00980 FoF jobs running;
00984 remained priority-pending. Direct overlapping process samples at 10:17
found all ten 00977 ranks, all five 00972 ranks, and all five 00980 ranks in
runnable state at roughly 99 per cent CPU. Snapshot 00972's latest FoF output
growth remains timestamped 10:12:36; the 00977/00980 files have not advanced.
No job has the `opFoF complete` marker, and the current downstream catalogue
chains remain dependency-pending.

## Darwin-1 current chains and host-sensitivity dependency recheck (2026-10-21 10:20 KST)

Slurm recheck at 10:19 KST showed FoF tasks 00972, 00977, and 00980 still
`RUNNING`; 00984 remained `PENDING (Priority)`. Read-only `srun` samples at
10:20 found all ten 00977 ranks and all five ranks each for 00972 and 00980
runnable at approximately 96--99.5 per cent CPU. Their recent FoF products
have not advanced since the 10:12:36 checkpoint for 00972; the other two
remain unchanged. No `GALFIND.DATA` files exist for any output. A direct
`scontrol show job 431211_0` check confirms the host-mass sensitivity task
still waits on `aftercorr:431524_*`, as intended; all current catalogue jobs
remain dependency-pending.

## Darwin-1 snapshot 00977 FoF catalogue advances (2026-10-21 10:22 KST)

The active FoF task for snapshot 00977 (`430418_1`) wrote a new output
checkpoint at 10:22:12 KST. Its halo catalogue is now 42,524,044 bytes and
its member-particle file is 94,162,813,888 bytes; the log reports ranks saving
isolated-halo data. The outputs for 00972 and 00980 have not advanced since
10:12:36 and 09:54:00, respectively. All three FoF jobs remain `RUNNING`, and
none of their logs contains the completion marker, so none of these partial
catalogues is ready for NewGalFinder or scientific analysis. The simulation
and all four intended snapshots are Darwin-1; the `Cocoa_analysis_20261020`
string is only the analysis-directory name.

## Darwin-1 final relation-figure dependency repaired (2026-10-21 10:23 KST)

The pending plot job `431161` was permanently blocked by `afterok` dependencies
on the obsolete analysis tasks `430431_1` and `430431_2`, which had been
cancelled during duplicate-chain cleanup. Its script requires all four current
per-output relation CSVs before plotting. I updated the pending job dependency
to `afterok:431528`, the current four-task satellite-analysis array, so the
figure can start only after that complete array succeeds. The PDF does not yet
exist; this change does not launch analysis early or alter any running job.

## Darwin-1 00977 FoF output continues; plotting dependency confirmed (2026-10-21 10:24 KST)

The 10:24 KST scheduler check found FoF tasks 00972, 00977, and 00980 still
running, with 00984 priority-pending. Snapshot 00977's FoF products advanced
again at 10:24:08 to 43,228,556 bytes (halo catalogue) and 98,093,728,704
bytes (member-particle file). The 00972 and 00980 products remain at their
10:12:36 and 09:54:00 checkpoints. No `GALFIND.DATA` exists for any of the
four outputs, and all current NewGalFinder, extraction, and relation-analysis
tasks are still dependency-pending. The plot job now reports
`afterok:431528_* (unfulfilled)`, correctly waiting for the current four-output
analysis array; no plot has run.

## Darwin-1 live FoF process check (2026-10-21 10:27 KST)

At 10:26 KST, Slurm still showed FoF jobs for 00972, 00977, and 00980
`RUNNING`; 00984 remained priority-pending, with all downstream steps waiting
on dependencies. A read-only overlapping process sample on the allocated nodes
found all five opFoF ranks for 00972 and all five ranks for 00980 actively
running at about 96--99.5 per cent CPU. Their output files have not advanced
since 10:12:36 and 09:54:00, respectively, but the ranks are computing rather
than idle; neither task should be restarted or treated as failed. Snapshot
00977's latest output checkpoint remains 10:24:08. No NewGalFinder catalogue
is present yet.

## Darwin-1 FoF rank-output checkpoint for 00977 (2026-10-21 10:28 KST)

The 10:28 KST Slurm check still showed 00972, 00977, and 00980 FoF running;
00984 remained priority-pending. Snapshot 00977 wrote another FoF checkpoint at
10:27:30: the halo-catalogue file is 44,105,100 bytes and the member-particle
file is 100,206,355,392 bytes. The log records rank P5 completing its
isolated-halo write. Snapshot 00972 remains at its 10:12:36 checkpoint and
00980 at 09:54:00. The current NewGalFinder, extraction, analysis, and figure
jobs remain dependency-pending; no `GALFIND.DATA` catalogue is available yet.

## Darwin-1 FoF compute remains active (2026-10-21 10:29 KST)

A fresh Slurm poll at 10:29 KST still showed 00972, 00977, and 00980 running;
00984 remained priority-pending. The 00972 and 00980 FoF files and logs had not
changed, so I sampled their allocated nodes again with an overlapping,
read-only process check. All five ranks for each snapshot were still runnable
at roughly 96--99.5 per cent CPU. This confirms continued computation during
the quiet output interval; the live jobs remain untouched and downstream
catalogue tasks remain dependency-gated.

## Darwin-1 snapshot 00984 FoF input gate rechecked (2026-10-21 10:33 KST)

Snapshot 00984's prepared NewDD input still passes the `opfof_array.sbatch`
preflight: all 100 required DM shard files are non-empty, all 100 stellar
shard paths exist (empty stellar shards are allowed), and the snapshot info
file is present. Its FoF task remains priority-pending, so no FoF or galaxy
catalogue output has been produced for this snapshot yet.

## Darwin-1 snapshot 00977 live-rank recheck (2026-10-21 10:41 KST)

Snapshot 00977's FoF log and output files have not advanced since 10:27:30,
but a direct read-only process check on allocated node `syn01` found five
`opfof.exe 977 100` ranks actively running at 99.4--99.5 per cent CPU. The
10:41 Slurm poll still reports 00977 `RUNNING`; 00972 and 00980 are also
`RUNNING`, while 00984 remains priority-pending. The unchanged files therefore
reflect ongoing compute, not a completed or failed task. No NewGalFinder output
is available yet.

## Darwin-1 FoF long-output interval explained by shard workload (2026-10-21 10:47 KST)

The FoF logs are quiet, but their last shard records show why this can take a
long time: each MPI rank processes its assigned NewDD shards sequentially,
then builds a local tree, links particles, checks halo boundaries, and writes
members/contact particles. The active binary's Makefile links `opfof.c` and
`Treewalk.fof.ordered.c`, matching the inspected loop. In the last logged
shards, ranks handled roughly 108--116 million particles each (for example,
00972 rank P8: 115,829,557 particles; 00980 rank P8: 115,305,736; 00977 rank
P9: 111,862,141). Those large per-shard tree/FoF workloads follow the last
output lines. This source/log evidence plus live CPU use supports leaving the
running jobs intact; it does not make their partial catalogues complete.

## Darwin-1 snapshot 00977 FoF checkpoint advances (2026-10-21 10:49 KST)

Snapshot 00977 advanced after the long quiet interval. At 10:49:02 the FoF
halo catalogue reached 44,809,612 bytes and the member-particle file reached
105,586,262,976 bytes, up from 44,105,100 and 100,206,355,392 bytes at
10:27:30. The new log records rank P4 finding 5,803,041 halos from
58,774,664 particles, checking their boundaries, and saving isolated-halo
data. The Slurm task remains `RUNNING`; it has no completion marker yet, and
its NewGalFinder catalogue has not been created.

## Darwin-1 FoF outputs advance for 00972 and 00980 (2026-10-21 10:53 KST)

A targeted Slurm poll at 10:53 KST shows FoF jobs 00972, 00977, and 00980
still `RUNNING`; 00984 remains priority-pending. The 00972 and 00980 halo
catalogues and member-particle files were updated at 10:52 KST and currently
occupy about 16 MB/20 GB and 16 MB/19 GB, respectively. Their latest logs show
ranks progressing through the assigned shards after local tree construction.
Snapshot 00977's products remain at the 10:49 checkpoint (43 MB/99 GB). No
FoF completion marker or `GALFIND.DATA` was observed, so dependent satellite
catalogue and relation-analysis jobs remain correctly pending. These are
Darwin-1 outputs; `Cocoa_analysis_20261020` is only the analysis-directory
name.

## Legacy 00972 relation files are not accepted as final inputs (2026-10-21 10:55 KST)

Inspection of the existing 00972 exploratory relation CSVs found
`duplicate_star_id_audit_available=False`; those files predate the current
catalogue extraction and are dated 2026-10-20. The current extraction writes
per-galaxy stellar-ID duplicate counts, audits IDs shared between galaxies,
and the current relation-analysis path requires those audit fields for its
catalogue schema. Therefore the earlier 00972 slopes remain exploratory and
are not evidence for the final physical claim. The working manuscript already
excludes legacy-catalogue values pending validation of all four Darwin-1
outputs.

## Darwin-1 FoF processes remain CPU-active during quiet output interval (2026-10-21 10:57 KST)

The 10:57 KST Slurm poll still reports 00972, 00977, and 00980 `RUNNING`,
with 00984 priority-pending and all dependent arrays waiting. FoF output files
have not advanced since 10:52 (00972/00980) and 10:49 (00977), but a single
read-only process sample on allocated nodes found ten 00972/00980 `opfof.exe`
processes and four sampled 00977 processes each using about 99 per cent CPU.
The unchanged files therefore do not indicate idle tasks; the active FoF jobs
were left untouched. No `GALFIND.DATA` or final satellite relation catalogue
is available yet.

## Darwin-1 FoF logs advance into halo-boundary checks (2026-10-21 11:02 KST)

The 00972 and 00980 FoF logs advanced after the 10:57 checkpoint, with new
entries at 10:59:50 and 11:00:03 KST. Ranks P2 and P3 report completing local
tree construction, identifying several million candidate haloes in their
current particle batches, and checking halo boundaries. Their FoF output files
still have the earlier 10:52 timestamps, and no `GALFIND.DATA` or completion
marker is present. These log updates show work advancing inside the active FoF
stage, but they do not yet make either snapshot ready for galaxy finding.

## Darwin-1 snapshot 00977 FoF output advances again (2026-10-21 11:03 KST)

At 11:03:45 KST, Slurm still showed 00972, 00977, and 00980 FoF jobs
`RUNNING`, with 00984 priority-pending and downstream arrays dependency-pending.
Snapshot 00977's FoF products advanced at 11:03:26: the halo catalogue grew
to 45,522,316 bytes and the member-particle file to 110,176,781,248 bytes,
about 4.6 GB above its 10:49 checkpoint. The log records ranks P4 and P1
sending and saving isolated-halo data. Snapshot 00972 and 00980 output files
remain at their 10:52 checkpoints. No FoF completion marker or `GALFIND.DATA`
exists yet, so the catalogue and relation analysis are still pending.

## Darwin-1 FoF stays CPU-active while outputs are quiet (2026-10-21 11:04 KST)

The 11:04:58 KST Slurm poll still showed 00972, 00977, and 00980 `RUNNING`;
00984 was priority-pending and downstream arrays were dependency-pending. Since
the 11:03:26 checkpoint, no FoF output file has advanced. A read-only process
sample nevertheless found ten 00972/00980 ranks and four sampled 00977 ranks
running at about 99 per cent CPU; their CPU elapsed times had advanced by about
4.8 minutes since the previous process sample. The active computations were
left untouched. No `GALFIND.DATA` or relation-analysis output is available.

## Darwin-1 snapshot 00977 FoF output advances at 11:06 KST

The 11:06:17 KST Slurm check still showed FoF 00972, 00977, and 00980
`RUNNING`; 00984 remained priority-pending, and satellite-finder and analysis
arrays remained dependency-pending. Snapshot 00977's output files advanced at
11:06:08: the halo catalogue reached 46,431,628 bytes and the member-particle
file reached 116,740,670,400 bytes, an increase of about 6.6 GB since 11:03.
The log records rank P8 checking boundaries and saving isolated-halo data.
The 00972 and 00980 files remain at their 10:52 checkpoints. No completion
marker or `GALFIND.DATA` exists yet.

## Darwin-1 downstream catalogue arrays now advance per snapshot (2026-10-21 11:10 KST)

The original dependencies made all four property-extraction tasks wait for the
entire NewGalFinder array, and all four relation analyses wait for the entire
extraction array. Since each array task reads only its matching snapshot, I
changed these two dependencies to Slurm `aftercorr`: extraction task *i* can
start after NewGalFinder task *i*, and analysis task *i* can start after
extraction task *i*. `scontrol show job` confirms both new dependencies and
shows them still correctly unfulfilled. The final figure job retains
`afterok:431528_*`, so it waits for all four relation analyses. No running job
was changed or restarted.

## Darwin-1 snapshot 00980 FoF rank logs progress (2026-10-21 11:11 KST)

A targeted 11:11:32 KST Slurm poll still showed FoF 00972, 00977, and 00980
`RUNNING`, with 00984 priority-pending. Snapshot 00980's log advanced at
11:11:14: ranks P1 and P5 reported halo identification and boundary checks.
Its FoF files remain at the 10:52 checkpoint, so the snapshot is not yet ready
for NewGalFinder. After the `aftercorr` change, the matching downstream array
tasks are represented individually and remain pending on their own upstream
snapshots; no extraction or analysis has started.

## Darwin-1 FoF remains active and per-snapshot dependencies are intact (2026-10-21 11:12 KST)

At 11:12:49 KST, the three active FoF jobs were still `RUNNING`; the sampled
ranks on `syn101` and `syn01` were using about 99 per cent CPU. Output files
remain at their latest checkpoints (00977 at 11:06; 00972/00980 at 10:52), and
no `GALFIND.DATA` exists. `scontrol show job` reconfirmed the per-snapshot
`aftercorr` dependencies for extraction and relation analysis, while the final
figure continues to wait for all four relation tasks to succeed.

## Darwin-1 snapshot 00972 FoF log advances (2026-10-21 11:14 KST)

At 11:14:08 KST, Slurm still showed FoF 00972, 00977, and 00980 `RUNNING`,
with 00984 priority-pending; all downstream tasks were still dependency-pending.
The 00972 log advanced at 11:13:20: ranks P5 and P1 identified haloes in new
particle batches and checked their boundaries. Its FoF output files remain at
the 10:52 checkpoint, and no `GALFIND.DATA` exists for any of the four
snapshots. The log progress is not a completion signal.

## Darwin-1 snapshot 00972 FoF output advances at 11:15 KST

At 11:15:28 KST, the FoF jobs for 00972, 00977, and 00980 were still
`RUNNING`; 00984 remained priority-pending, and downstream work remained
pending. Snapshot 00972's products advanced at 11:15:01: the halo catalogue
reached 19,295,052 bytes and the member-particle file reached 26,241,271,800
bytes, about 5.8 GB above its 10:52 checkpoint. The log says rank P0 began
writing isolated-halo data for 5,763,660 haloes after several ranks reported
saving their data. This is output progress, not job completion; no
`GALFIND.DATA` has been produced.

## Darwin-1 FoF output and rank progress continue (2026-10-21 11:16 KST)

At 11:16:21 KST, Slurm still reported 00972, 00977, and 00980 `RUNNING`;
00984 remained priority-pending and all downstream tasks were pending. Snapshot
00972's halo catalogue and member-particle file advanced at 11:15:54 to
19,999,564 bytes and 27,480,623,096 bytes, respectively; rank P6 reported
saving isolated-halo data. Snapshot 00980's log advanced at 11:15:36 with rank
P6 checking halo boundaries, although its output files remain at the 10:52
checkpoint. Snapshot 00977 remains at its 11:06 checkpoint. No FoF completion
marker or `GALFIND.DATA` exists.

## Darwin-1 snapshot 00980 FoF output advances at 11:18 KST

At 11:18:57 KST, the FoF jobs for 00972, 00977, and 00980 remained
`RUNNING`; 00984 remained priority-pending and all downstream arrays were
pending. Snapshot 00980's files advanced at 11:18:18: the halo catalogue
reached 19,907,220 bytes and the member-particle file reached 26,850,819,040
bytes, about 7.0 GB above its 10:52 checkpoint. The log shows rank P0 preparing
to write data for 5,745,168 isolated haloes after several other ranks saved
their data. This is not job completion; `GALFIND.DATA` is still absent.

## Darwin-1 snapshot 00977 FoF output advances at 11:21 KST

At 11:21:35 KST, the 00972, 00977, and 00980 FoF jobs remained `RUNNING`,
with 00984 priority-pending and all dependent jobs pending. Snapshot 00977's
files advanced at 11:21:28: the halo catalogue reached 47,947,148 bytes and
the member-particle file reached 123,267,900,352 bytes, about 6.5 GB above
its 11:06 checkpoint. Rank P9 reported completing boundary checks and saving
isolated-halo data. `GALFIND.DATA` is still absent for all four snapshots.

## Darwin-1 snapshot 00980 FoF output advances at 11:22 KST

At 11:22:53 KST, FoF 00972, 00977, and 00980 were still `RUNNING`; 00984
remained priority-pending and all downstream tasks remained pending. Snapshot
00980's files advanced at 11:22:18: its halo catalogue reached 20,800,148 bytes
and member-particle file reached 28,435,045,344 bytes, about 1.6 GB above the
11:18 checkpoint. Rank P4 reported checking boundaries and saving isolated
haloes. No `GALFIND.DATA` exists yet.

## Darwin-1 snapshot 00972 FoF output advances at 11:23 KST

At 11:24:13 KST, FoF 00972, 00977, and 00980 were still `RUNNING`, while
00984 remained priority-pending and the downstream arrays remained pending.
Snapshot 00972's output files advanced at 11:23:10: the halo catalogue reached
20,892,492 bytes and the member-particle file reached 29,052,790,776 bytes,
about 1.6 GB above its 11:15 checkpoint. Rank P4 reported saving
isolated-halo data. `GALFIND.DATA` remains absent for all four outputs.

## Darwin-1 ID audit memory bounded before catalogue extraction (2026-10-21 11:31 KST)

A targeted 11:31:17 KST Slurm check showed FoF snapshots 00972, 00977, and
00980 still running; 00984 remained priority-pending. NewGalFinder and
extraction tasks remained pending on their per-snapshot dependencies, and no
catalogue extraction had started. Before those jobs become eligible,
`extract_satellite_properties.py` was changed to replace the snapshot-wide list
and final concatenation of unique stellar-particle IDs with a fixed 512 MiB
bitset keyed by the NewGalFinder 32-bit particle ID. The per-host audit still
checks duplicates across galaxies and the bitset rejects IDs repeated across
FoF hosts. The patch was reviewed by inspection only; it has not yet been run
against a generated Darwin-1 catalogue.

## Darwin-1 FoF jobs confirmed active by Slurm step accounting (2026-10-21 11:33 KST)

The targeted Slurm poll still showed FoF snapshots 00972, 00977, and 00980
running, with 00984 priority-pending; no `GALFIND.DATA` or extracted satellite
catalogue exists yet. `sstat --allsteps` reported the three active FoF task
steps with 1.0--2.0 accumulated CPU-days and about 87--91 GiB maximum RSS.
This confirms the jobs are actively doing substantial work. At that
accounting check, the FoF output files and rank logs had not advanced since
11:23 KST; the subsequent file check recorded below found newer output from
00972 and 00980. No live job was modified or restarted.

## Darwin-1 FoF 00972 and 00980 output growth; 00977 remains active (2026-10-21 11:34 KST)

A targeted 11:33:40 KST check showed all three FoF tasks still `RUNNING` and
00984 still priority-pending. Snapshot 00972's files advanced at 11:33:06 to
22,006,604 bytes (halo catalogue) and 31,546,124,280 bytes (member particles),
up from 20,892,492 and 29,052,790,776 bytes at 11:23. Snapshot 00980 advanced
at 11:31:58 to 21,914,260 and 30,840,994,784 bytes, up from 20,800,148 and
28,435,045,344 bytes at 11:22. Snapshot 00977's files and rank log last
advanced at 11:21:28; its Slurm task remains active and its accounting shows
substantial CPU and memory use. No NewGalFinder catalogue has been created,
so the dependent extraction and science analyses have not started.

## Introduction framed as a complementary conditional test (2026-10-21 11:38 KST)

After checking the cited primary studies, the manuscript now acknowledges their
existing metallicity--kinematics and environmental results, then poses this
study's narrower question directly: whether satellite metallicity is associated
with internal stellar dispersion within the same cosmological host, controlling
for stellar mass and present-day finder-assigned subhalo dark-matter mass. This
avoids implying that a metallicity--kinematics relation itself is new. The
revision is in `satellite_metallicity_kinematics_manuscript.md`; it does not
change or anticipate the still-pending Darwin-1 measurements.

## MNRAS submission requirements added to the manuscript plan (2026-10-21 11:40 KST)

The official MNRAS author instructions were checked while the Darwin-1 FoF jobs
continued running. `paper_outline.md` now records the final-format requirements
relevant to this study: a single-paragraph abstract of at most 250 words,
1--6 approved keywords, numbered sections and final section order, sequential
figure/table citations, vector quantitative plots with labelled axes, figure
alt text, and machine-readable catalogue supplements. The current status poll
still showed 00972/00977/00980 FoF running, 00984 priority-pending, and all
NewGalFinder/extraction/analysis arrays dependency-pending. No satellite
catalogue exists yet.

## MNRAS keywords provisioned; Darwin-1 FoF remains CPU-active (2026-10-21 11:43 KST)

The official MNRAS keyword list was checked and five provisional terms were
added to `paper_outline.md`: `galaxies: abundances`, `galaxies: clusters:
general`, `galaxies: kinematics and dynamics`, `galaxies: evolution`, and
`methods: numerical`. The final set remains subject to the validated results.
At 11:43:37 KST, Slurm still showed FoF 00972, 00977, and 00980 running and
00984 priority-pending; downstream NewGalFinder, extraction, and analysis jobs
were dependency-pending. The FoF products had not advanced since 11:33 (00972),
11:31 (00980), and 11:21 (00977), but `sstat` task CPU counters continued to
rise. No `GALFIND.DATA` or satellite-property catalogue exists yet.

## Darwin-1 FoF task I/O continues while catalogues remain incomplete (2026-10-21 11:45 KST)

At 11:45:26 KST, Slurm still showed FoF 00972, 00977, and 00980 running;
00984 remained priority-pending, with NewGalFinder and downstream analysis
arrays still dependency-pending. `sstat` task CPU counters had advanced since
the 11:43 sample. The 00977 FoF task's cumulative disk-write counter increased
by about 1.03 GiB, although its FoF output files had not yet received a newer
mtime. No `GALFIND.DATA` or satellite catalogue exists yet; do not treat the
step-level I/O counter as a completed FoF product.

## Darwin-1 snapshot 00977 FoF output advances again (2026-10-21 11:45 KST)

At 11:45:44 KST, Slurm still showed FoF snapshots 00972, 00977, and 00980
running; 00984 remained priority-pending. Snapshot 00977's FoF output and rank
log advanced at 11:44:03: its halo catalogue grew from 47,947,148 to 49,151,372
bytes, and its member-particle file grew from 123,267,900,352 to
126,578,590,656 bytes (about 3.08 GiB). The log records rank P3 processing
6,503,004 haloes from 107,787,892 particles, checking their boundaries, and
saving isolated-halo data. This is a completed shard checkpoint, not FoF job
completion; no NewGalFinder catalogue exists yet.

## Darwin-1 simulation identity and catalogue pipeline check (2026-10-21 11:47 KST)

The simulation identity is Darwin-1. The analysis directory and historical
RAMSES output path contain `Cocoa` in their storage names, but the job scripts,
methods draft, manuscript, and existing project notes identify the simulation
as Darwin-1. A targeted check found no completed `GALFIND.DATA`, extracted
satellite-property catalogue, or relation-analysis CSV for snapshots 00972,
00977, 00980, or 00984. Slurm still showed FoF 00972, 00977, and 00980 running,
00984 pending priority, and the dependent catalogue/analysis stage pending.
The latest visible 00977 FoF log/output update remains 11:44:03 KST; this is a
checkpoint only, so no satellite-galaxy measurements are available to report.

## Darwin-1 FoF jobs remain active; downstream catalogues have not started (2026-10-21 11:50 KST)

Slurm accounting at 11:50 KST confirms FoF tasks 430418_1 (snapshot 00977),
430610_0 (00972), and 430611_2 (00980) are still running within their
48-hour wall-time requests. Snapshot 00984 remains pending on its dependency.
The four NewGalFinder array tasks and downstream extraction/analysis tasks
remain pending. Targeted checks found no `GALFIND.DATA`, extracted satellite
catalogues, or relation tables. The FoF logs and output files have not advanced
since the 00977 checkpoint at 11:44:03 KST. No job was restarted or modified.

## Darwin-1 FoF activity and stellar-velocity units verified (2026-10-21 11:52 KST)

At 11:52:04 KST, `squeue` still showed FoF 00972, 00977, and 00980 running;
00984 was priority-pending, and all NewGalFinder, extraction, and relation
analysis tasks were dependency-pending. `sstat --allsteps` showed active MPI
steps for the three FoF runs, with cumulative CPU use of 2 d 1 h (00977),
1 d 4 h (00972), and 1 d 4 h (00980), plus tens to hundreds of GB of recorded
disk writes. These counters confirm active work but are not completion markers;
the latest FoF log/output update remains 11:44:03 KST and no galaxy catalogue
has appeared.

Before interpreting the extracted dispersion fields as km/s, their unit path
was checked in source: `NewDD/rd_info.c` defines `kmscale_v` in km/s,
`NewDD/rd_part.c` applies it to stellar-particle velocities, and
`NewGalFinder/gfind.c` copies those velocities into the galaxy records without
another conversion. `extract_satellite_properties.py` now documents that it
preserves the km/s values, avoiding an unsupported second conversion.

## Signed particle-ID duplicate audit aligned with NewDD tracking convention (2026-10-21 11:57 KST)

Static source review found that `NewDD/rd_part.c` preserves signed particle
IDs but documents `abs(ID)` as the downstream tracking key. The extractor had
compared signed values, which could miss a repeated tracking key with opposite
signs. It now widens IDs to 64-bit integers before taking their absolute value
and uses that key for per-galaxy, within-host, and cross-host duplicate checks;
the bitmap also covers the largest possible absolute signed 32-bit ID without
overflow. The methods draft and working manuscript now state this convention.
The extraction jobs are still pending, so they will read the updated script at
runtime. This was a static code/documentation change; no test or calculation
was run.

## Darwin-1 FoF work continues while output files remain unchanged (2026-10-21 11:59 KST)

The 11:59 KST Slurm check still reports FoF 00972, 00977, and 00980 as
`RUNNING`; 00984 is priority-pending, and all NewGalFinder, extraction, and
analysis tasks remain dependency-pending. Since the 11:56 check, the MPI-step
CPU counters rose by about 10.5 minutes for 00977 and 15.8 minutes each for
00972 and 00980; 00980's recorded maximum RSS also rose by about 1 GiB. The
three FoF logs and output files have not advanced since their 11:31--11:44
checkpoints, and no galaxy catalogue exists. Increased CPU and memory use with
unchanged output supports continued in-memory processing; there is no evidence
to restart these jobs.

## Darwin-1 FoF still active; shared-filesystem capacity checked (2026-10-21 12:01 KST)

At 12:01 KST, Slurm continued to report FoF 00972, 00977, and 00980 as
`RUNNING`, with their MPI-step CPU counters increasing; 00984 remained
priority-pending and dependent galaxy finding, extraction, and analysis stayed
queued. The FoF logs and output files were still at the previously recorded
checkpoints. The GPFS quota report showed no per-user block or file limit, and
`df` reported 231 TB available on the shared filesystem, so current storage
capacity is not the limiting resource. No job or output was changed.

## Darwin-1 FoF tasks remain CPU-active (2026-10-21 12:03 KST)

At 12:03 KST, the three running FoF MPI steps had accumulated about 7 minutes
more CPU for 00977 and about 10.5 minutes more each for 00972 and 00980 since
the 12:01 check. Their recorded peak resident memory also rose, by about 0.3
GiB for 00972 and 1.8 GiB for 00980. The jobs remain `RUNNING`, but FoF logs
and output files still show their earlier 11:31--11:44 checkpoints. Snapshot
00984 and the dependent galaxy-catalogue chain remain pending. This confirms
continued computation without indicating completion or a remaining-time
estimate; no job was restarted.

The current FoF stderr logs contain warnings that the SINK and GAS component
files are missing. These are expected for this pipeline: `newdd_array.sbatch`
sets `NEWDD_SKIP_AMR=1`, intentionally omitting AMR, gas, and sinks while
retaining dark-matter and stellar particles. The warnings do not indicate a
failed FoF input check.

## Darwin-1 FoF resource counters continue to rise (2026-10-21 12:06 KST)

The 12:06 KST scheduler/accounting check still reports FoF 00972, 00977, and
00980 `RUNNING`; the MPI-step CPU counters increased since the prior sample,
and recorded peak RSS rose to about 91 GiB, 94 GiB, and 96 GiB, respectively.
Snapshot 00984 remains priority-pending; all NewGalFinder and later stages
remain dependency-pending. The FoF logs and products are still at the earlier
11:31--11:44 checkpoints, so none of the dependent catalogue outputs is ready.
The three active jobs were left untouched.

## Darwin-1 FoF 00984 waits for resources; other snapshots continue (2026-10-21 12:11 KST)

The 12:11 KST Slurm check still shows FoF 00972, 00977, and 00980 `RUNNING`,
with their MPI-step CPU counters and recorded maximum RSS increasing. FoF
00984 (`430622_3`) remains pending; its scheduler reason changed from
`Priority` to `Resources`. NewGalFinder, extraction, and relation analysis
remain dependency-pending. No new FoF or satellite-catalogue files appeared;
the three active jobs were left untouched.

## Darwin-1 FoF 00977 output advances into tree building (2026-10-21 12:13 KST)

At 12:13 KST, snapshot 00977's halo catalogue grew from 49,151,372 to
49,913,668 bytes and its member-particle file grew from 126,578,590,656 to
130,435,463,088 bytes (about 3.59 GiB). The FoF log advanced at 12:13:28 and
reports new particle slabs and successful tree-building messages from every
rank P0--P9. The 12:14:22 Slurm poll still showed its FoF task running; the
files therefore remain intermediate outputs, not a completed FoF catalogue.
Snapshots 00972 and 00980 had no new output, and no NewGalFinder catalogue or
satellite-property table exists yet.

## Darwin-1 snapshot 00984 scheduler reservation (2026-10-21 12:16 KST)

`scontrol show job 430622_3` reports state `PENDING (Resources)` and a
projected start of 2026-10-22 06:10 KST, coinciding with the 24-hour end time
reserved for the currently running 00972/00980 jobs on A100 nodes syn101--102.
The queued 00984 FoF task requests 2 nodes, 10 CPUs, and 448 GiB of memory
(224 GiB per node); `ReqTRES` contains no GPU/GRES request. The A100 partition
currently has no idle whole nodes. This start time is a scheduler forecast,
not a guarantee; the independent 00972/00977/00980 chains can proceed before
00984. No pending job was changed.

## Snapshot 00984 memory alternative reviewed (2026-10-21 12:19 KST)

The free A100-PCIe node syn103 has 512 GiB of physical memory but is a
single-node partition; the 00984 FoF request is for two nodes. An `sstat`
sample of the analogous running 00972/00980 FoF MPI steps reported an AveRSS
near 69 GiB and MaxRSS of 97--98 GiB on task 1 of syn102. These are step-level
task summaries, not separate per-node memory totals; they do not establish
that the workload would fit safely on syn103. The queued 00984 job was left
unchanged. The existing multi-node job is forecast to start after the A100
nodes' current reservations release.

## Darwin-1 identity clarification and 00977 FoF progress (2026-10-21 12:23 KST)

The user clarified that the simulation is **Darwin-1**, not Cocoa. `Cocoa_analysis_20261020` is only the working-directory name; it does not identify the simulation. At 12:23:14 KST, Slurm still showed FoF jobs for snapshots 00972 (`430610_0`), 00977 (`430418_1`), and 00980 (`430611_2`) running; 00984 (`430622_3`) remained pending for resources, and downstream NewGalFinder tasks remained dependency-pending. The 00977 log had advanced beyond P7's boundary check: P6 also reported 6,595,238 halos from 27,238,900 particles and completed its halo-boundary check. This is rank-level progress only; it does not establish that the full FoF catalogue is complete. No jobs were changed.

## Darwin-1 downstream gate and FoF product check (2026-10-21 12:24 KST)

A static review of `gfind_array.sbatch`, `extract_array.sbatch`, and `analyse_array.sbatch` confirms the staged guards are in place: NewGalFinder requires both FoF products and checks for its catalogue outputs; extraction requires `GALFIND.DATA`; relation analysis requires the extracted satellite table. A targeted product check found no `GALCATALOG.LIST`, `GALFIND.DATA`, satellite-property CSV, or relation CSV for snapshots 00972, 00977, 00980, or 00984. FoF product sizes/timestamps were unchanged from the 12:22 checkpoint for 00972/00980 and 12:12 checkpoint for 00977; 00984 still has no FoF products. The Slurm dependencies therefore remain correctly closed, and no downstream job was released manually.

## Darwin-1 00977 FoF remains actively computing (2026-10-21 12:26 KST)

A follow-up scheduler check at 12:26:41 KST still showed FoF tasks 00972, 00977, and 00980 `RUNNING`, with 00984 pending for resources and the downstream arrays waiting on dependencies. Snapshot 00977's FoF files are unchanged since 12:12:41, and its log is unchanged since P6 completed the halo-boundary check. Although the `sstat` query under the array alias returned no row, resolving the actual array element showed its live step (`430419.0`, `hydra_bstrap_proxy`) at 15:27 elapsed, 2-03:16:15 aggregate CPU, 95.3 GB maximum RSS, and about 140 GB disk read / 147 GB disk write. Thus the lack of new log/file timestamps is not evidence that the job has stopped; it remains resource-active, so it was left running.

## Darwin-1 00977 FoF log-stage audit (2026-10-21 12:28 KST)

The current 00977 FoF log contains 1,485 lines. In its latest phase, all ten ranks report `Tree building` complete; ranks P7 and P6 have subsequently reported halo counts and completed their boundary checks for their current particle slabs (slab 00077). Earlier output shows the same boundary-check messages for other ranks during a prior processing phase, so these messages are per-stage/per-slab, not whole-job completion markers. At 12:28:03, Slurm still showed 00977 running and the downstream arrays dependency-pending. The three available FoF output files remain unchanged from their recorded timestamps; no galaxy catalogue is ready.

## Darwin-1 00977 FoF advances to P2 boundary check (2026-10-21 12:29 KST)

At 12:29:52 KST, the 00977 FoF job remained `RUNNING`; its live MPI step's aggregate CPU counter increased from 2-03:16:15 at 12:26 to 2-03:23:44, and the log advanced to rank P2's current halo-boundary check. This confirms ongoing computation despite the unchanged final FoF-file timestamps. Snapshots 00972 and 00980 also remain `RUNNING`; 00984 is still resource-pending, while all dependent NewGalFinder/extraction/analysis tasks remain pending. No catalogue or satellite-property CSV exists yet, and no job was modified.

## Darwin-1 three FoF jobs remain resource-active (2026-10-21 12:33 KST)

The 12:33 scheduler check still showed 00972, 00977, and 00980 running, while 00984 remained pending for resources and the follow-on arrays remained dependency-pending. Snapshot 00972 and 00980 logs and FoF member files are unchanged since 11:33 and 11:31, respectively; their live MPI steps nevertheless report about 1-07:31 and 1-07:44 aggregate CPU, maximum RSS near 99 GB, and roughly 100 GB read / 51 GB written. Snapshot 00977's live step reports 2-03:34 aggregate CPU and about 95 GB maximum RSS, with its log currently at P2's boundary check. These are active jobs, not completed catalogues; no NewGalFinder or satellite-property output exists.

## Darwin-1 00977 FoF root writer starts final halo consolidation (2026-10-21 12:34 KST)

The 12:34:37 scheduler poll still showed FoF 00972, 00977, and 00980 running; 00984 remained pending for resources and the NewGalFinder/extraction/analysis arrays remained dependency-pending. Snapshot 00977's live MPI-step CPU counter rose to 2-03:39:27 since 12:31:42, and its log advanced: P2 sent its isolated-halo data to P0 and saved its local output; P0 reports preparing to write the consolidated data for 5,585,046 halos. The FoF halo catalogue and member file grew to 52,543,300 bytes and 136,562,309,040 bytes, respectively, with 12:33:00 timestamps. Since the Slurm task is still running and the final job marker has not appeared, these remain in-progress products; downstream jobs were not started manually.

## Darwin-1 00977 consolidated FoF writer remains active (2026-10-21 12:36 KST)

At 12:36:11 KST, Slurm still showed 00977 FoF `RUNNING`; the matching NewGalFinder task and later stages remained dependency-pending. The FoF files stayed at the 12:33 sizes while the live MPI-step aggregate CPU counter rose from 2-03:39:27 at 12:34:37 to 2-03:45:06. The log remains at P0's consolidated-write stage for 5,585,046 halos, with P7, P2, and P6 having sent and saved their isolated-halo data. No final `opFoF complete` marker or GALFIND output is present. The compute step is active, so it was left untouched.

## Darwin-1 00977 FoF output grows again after P4 completes (2026-10-21 12:37 KST)

The 12:37:52 check shows the 00977 FoF task still running, with its dependent NewGalFinder/extraction/analysis tasks still pending. Since 12:36, the halo and member files have grown to 53,206,852 bytes and 138,008,041,392 bytes (12:36:38 timestamps); the live step's aggregate CPU counter rose to 2-03:50:15 and its disk-write counter increased by about 1.45 GB. The log now records P4's halo-boundary check and transfer/save of its isolated-halo data after the P7, P2, and P6 transfers. The full FoF job has not yet completed; no galaxy-finder output is ready.

## Darwin-1 00977 FoF consolidation continues (2026-10-21 12:39 KST)

At 12:39:13 KST, FoF 00977 remained `RUNNING`, and the corresponding NewGalFinder task was still dependency-pending. Its live MPI-step aggregate CPU counter rose from 2-03:50:15 at 12:37:52 to 2-03:54:44. The consolidated FoF files and log have not changed since 12:36:38; the latest log still shows P4, P6, P7, and P2 transferring/saving isolated-halo outputs. The job remains active, with no final completion marker or galaxy catalogue yet.

## Darwin-1 00977 FoF output advances through rank P1 (2026-10-21 12:39 KST)

At 12:39:38 KST, all three active FoF jobs remained `RUNNING`; 00984 remained resource-pending. Snapshot 00977's halo and member files advanced again to 53,976,900 bytes and 140,091,209,648 bytes (12:39:21 timestamps). Its log records P1 completing its boundary check and transferring/saving its isolated-halo data after P4. The live MPI-step CPU counter is 2-03:56:06; the job has not yet emitted its final completion marker. Snapshot 00972/00980 output files remain at their earlier timestamps, although their live-step CPU counters are still increasing. No GALFIND or extracted satellite-property products are present.

## Darwin-1 FoF tasks remain live while 00977 consolidation waits (2026-10-21 12:41 KST)

The 12:41:04 KST poll still shows FoF 00972/00977/00980 running and 00984 pending for resources; all downstream arrays remain dependency-pending. Snapshot 00977's output and log have not changed since P1 saved its isolated-halo data at 12:39:21, but the live MPI-step aggregate CPU counter rose from 2-03:56:06 to 2-04:00:53. The 00972 and 00980 MPI-step CPU counters also continued increasing. No completion marker or satellite catalogue exists yet.

## Darwin-1 00972/00980 FoF checkpoint audit (2026-10-21 12:42 KST)

A targeted log search shows snapshots 00972 and 00980 have also reached their root-writer consolidation phases: P0's last logged message says it is preparing isolated-halo data for 5,763,660 and 5,745,168 halos, respectively, with several ranks reporting their local data saved. Those log files and FoF products have not been updated since 11:33 (00972) and 11:31 (00980), so these are stale stage markers, not evidence of current output completion. However, at 12:41:37 their live-step aggregate CPU counters were 1-08:13 and 1-08:25, up by about ten minutes from the 12:39:38 sample; both Slurm jobs remain running. No intervention was made.

## Darwin-1 FoF consolidation CPU progress without file refresh (2026-10-21 12:44 KST)

At 12:44:40 KST all three active FoF tasks were still running and 00984 was still pending for resources; NewGalFinder and later stages remained dependency-pending. Snapshot 00977's files and log remain unchanged since 12:39:21, with P1's transfer/save as the last marker, but its live MPI-step CPU counter increased from 2-04:08:34 to 2-04:12:49 between 12:41:04 and 12:44:40. The 00972 and 00980 step CPU counters also rose. This is continued computation without a new output checkpoint; no job was altered.

## Darwin-1 00984 FoF resource queue rechecked (2026-10-21 12:46 KST)

`scontrol show job 430622_3` still reports `PENDING (Resources)`, forecast start 2026-10-22 06:10 KST, partition `a100`, two nodes, 10 CPUs, and 448 GiB memory. `ReqTRES` contains no GPU/GRES request, so the wait is not due to GPU allocation quota. `sinfo -s` shows both A100 nodes allocated and one idle A100-PCIe node; the latter is a single-node partition and does not match the pending two-node request. The other A10, A40, H100, and H200 nodes are allocated. Since the single-node alternative has not been shown to safely satisfy the workload's memory/runtime needs, the pending job and its dependency chain were left unchanged.

## Darwin-1 FoF live-step revalidation (2026-10-21 12:46 KST)

At 12:46:05 KST, the three FoF tasks remained `RUNNING`; the 00984 job remained resource-pending, and its downstream arrays remained dependency-pending. The 00977 FoF files/log and the 00972/00980 FoF files/logs were unchanged since their last timestamps, but live MPI-step aggregate CPU counters continued to rise: 00972 1-08:35, 00980 1-08:47, and 00977 2-04:17. Their latest visible logs are still the root-consolidation/local-save phase. No final FoF completion marker or GALFIND catalogue exists.

## Darwin-1 00977 FoF output advances through rank P8 (2026-10-21 12:47 KST)

At 12:47:46 KST, all three active FoF tasks remained running, with 00984 still resource-pending and downstream catalogue/analysis jobs pending on dependencies. Snapshot 00977's FoF files grew again at 12:47:27 to 54,853,444 bytes (halo catalogue) and 141,569,775,536 bytes (member-particle file). The log records P8 completing its boundary check and sending/saving its isolated-halo data. Its live MPI-step CPU counter reached 2-04:23:07. The full task is still running and no GALFIND output is ready.

## Darwin-1 FoF live recheck (2026-10-21 12:49 KST)

The 12:49:04 poll still shows 00972/00977/00980 FoF jobs running, 00984 pending for resources, and downstream jobs waiting on their dependencies. The 00977 log and files have not changed since P8's save at 12:47:27; its live step's aggregate CPU has risen to 2-04:27:27. No satellite catalogue is ready.

## Darwin-1 FoF output-path capacity check (2026-10-21 12:50 KST)

A targeted filesystem check found about 230 TiB free on the `/gpfs` scratch mount; the user's GPFS quota has no configured block or file limit. Thus shared-filesystem capacity/quota is not the current reason the FoF output files have paused between updates. At 12:49:27, all three FoF jobs were still running and their live-step CPU counters continued to rise; no satellite catalogue was present.

## Darwin-1 00977 FoF active-step check (2026-10-21 12:51 KST)

At 12:50:55 KST, Slurm still reports the three FoF jobs running and 00984 pending for resources. Snapshot 00977's file/log timestamps remain at 12:47:27, but its live MPI-step aggregate CPU counter increased from 2-04:28:41 at 12:49:27 to 2-04:33:35. Its latest log still shows P1/P8 local halo saves; no final completion marker or downstream catalogue has appeared. `/gpfs` has ample space, so the job was left to continue.

## Darwin-1 downstream dependency and NewGalFinder memory gate verified (2026-10-21 12:53 KST)

`scontrol` checks confirm the four NewGalFinder jobs depend on the matching FoF array tasks via `afterok`; each extraction job depends on its matching NewGalFinder task, and each relation-analysis job depends on its matching extraction task. The NewGalFinder requests two A100 nodes and 600 GiB per node (1.2 TiB total); current A100 nodes provide 1,028,000 MiB each, so the per-node memory request is feasible on that partition. Both A100 nodes are currently allocated. The idle A100-PCIe node has 512,000 MiB and cannot meet that same per-node request or two-node layout; the other idle/active partition alternatives therefore do not establish a safe replacement. The existing dependency chain remains valid; no job was changed.

## Darwin-1 FoF remains active after dependency audit (2026-10-21 12:54 KST)

The 12:54:31 KST poll still shows FoF 00972, 00977, and 00980 running; 00984 remains resource-pending. Snapshot 00977's files/log remain unchanged since 12:47:27, but its live MPI-step CPU counter increased from 2-04:39:07 at 12:52:36 to 2-04:45:29. The downstream dependencies are intact and will release the matching catalogue/extraction/analysis jobs only after each upstream success. No NewGalFinder catalogue has appeared.

## Darwin-1 FoF job time-limit check (2026-10-21 12:58 KST)

At 12:57:50 KST, the 00972 and 00980 FoF jobs were still running with 6:47 elapsed and a 24-hour limit (scheduled end 2026-10-22 06:10 KST). The 00977 FoF job was running with 15:57 elapsed and a 48-hour limit (scheduled end 2026-10-22 21:00 KST). The live-step CPU counters continued increasing, but no output/log update followed the 12:47:27 00977 checkpoint. All NewGalFinder jobs remain pending on the matching FoF successes, so the pipeline is not at risk of processing partial files. No satellite catalogue exists yet.

## Darwin-1 manuscript comparison literature rechecked (2026-10-21 12:59 KST)

A targeted review of the primary MNRAS papers already cited confirms the draft's stated distinction: Spindler & Wake (2017) compare central and satellite metallicity/kinematics in SDSS, with little central-satellite metallicity difference after star-formation separation; Fattahi et al. (2018) discuss the Local Group dwarf [Fe/H]--$V_{1/2}$ relation and interpret its low-velocity flattening through tidal stripping, assuming the mean metallicity is unchanged by tides. These are useful context but differ from the planned satellite-only, within-host conditional test using total birth-metallicity mass fraction and three-dimensional stellar dispersion. Existing manuscript wording already makes that distinction, so no premature results-dependent edits were made.

## Darwin-1 FoF status and downstream gate revalidated (2026-10-21 13:01 KST)

At 13:00:49 KST, the 00972, 00977, and 00980 FoF jobs remained `RUNNING`; 00984 remained pending for resources, with its NewGalFinder task also waiting. The live MPI-step CPU counters continued to increase (00977: 2-05:06:24), while the 00977 FoF files and log were unchanged since 12:47:27. No catalogue/property/relation output exists yet. The documented after-success dependency chain remains intact.

## Darwin-1 FoF root-consolidation stage diagnosed from source (2026-10-21 13:02 KST)

The job's executable path is `/gpfs/kjhan/Others/Jaehyun/2026.Apr/opFoF/opfof.exe`. Its Makefile selects `Treewalk.fof.ordered.c`; the corresponding object and executable were built after that source file. In the P0-root `WriteIsolatedHalo` branch, P0 first walks and writes its local isolated halos, then waits for the other MPI ranks, receives their halo/particle payloads, writes them to the shared files, and prints the completion message only after all ranks finish sending. The log sequence (P0 preparing 5,585,046 halos; ranks P1/P2/P4/P6/P7/P8 saving at different times) matches this long consolidation stage. This explains why output/log timestamps can pause while rank-side CPU counters continue increasing; it does not yet prove completion. No job was modified.

## Darwin-1 00977 second consolidation phase waiting on remaining ranks (2026-10-21 13:02 KST)

The latest 00977 log phase begins at line 1494 with P0 preparing the consolidated data for 5,585,046 halos. In this phase, P1, P2, P4, P6, P7, and P8 have sent their isolated-halo data; P3, P5, and P9 have not yet logged their send, and P0 has not logged `have written all the isolated halo data`. The source routine explains that P0 waits for and writes each remaining rank's halo/particle payload before emitting that marker. At 13:02:23 the Slurm job remained `RUNNING` and its aggregate CPU counter was 2-05:11:38, while the output files remained at their 12:47:27 sizes. The next completion signal to watch is those remaining rank transfers followed by P0's write-complete marker and Slurm success.

## Darwin-1 00977 remaining-rank workload identified (2026-10-21 13:04 KST)

In the current 00977 consolidation phase, P3, P5, and P9 are the three ranks not yet recorded sending isolated-halo data to P0. Their latest log messages show they have passed tree building but have not yet logged halo counts/boundary completion. Their current slabs contain about 101.8 million particles (P3, slab 00037), 85.1 million (P5, slab 00057), and 108.9 million (P9, slab 00097), substantially larger than several already completed ranks. The root writer is waiting for their boundary checks and payloads; this is consistent with the continuing MPI-step CPU usage and explains the output pause. The live job remains running; it was not interrupted.

## Darwin-1 00977 boundary-check cost confirmed from source (2026-10-21 13:06 KST)

The active `CheckHaloBound` routine is linear in each rank's particle count plus its halo count: it initializes per-halo boundary records, walks every local particle once to build halo membership and z limits, then scans the halo list. The current-phase ranks P3/P5/P9 each passed tree building but have not yet printed their halo counts, which are emitted only after this routine returns. Their 85--109 million-particle slabs can therefore account for the current quiet-log interval while the other ranks and P0 wait. This code-path evidence supports continued computation rather than a dependency or output-file error; no process was interrupted.

## Darwin-1 00977 remaining boundary checks continue (2026-10-21 13:07 KST)

At 13:07:13 KST, the FoF tasks remained running, and their downstream jobs remained dependency-pending. Snapshot 00977's last log/file timestamp is still 12:47:27; P3, P5, and P9 have not yet emitted their current-phase halo-count/boundary markers. The live MPI-step aggregate CPU counter rose from 2-05:22:51 at 13:05:46 to 2-05:27:40. The source shows this stage scans each rank's large particle slab to assign halo membership and boundary flags; no sign of a failed job or disk-space issue has appeared.

## Darwin-1 00977 large-rank boundary interval persists (2026-10-21 13:09 KST)

At 13:09:20 KST, 00977 remained running; its output files and log remain unchanged since 12:47:27. The live MPI-step aggregate CPU counter rose from 2-05:27:40 at 13:07:13 to 2-05:34:39. The current log still shows six non-root ranks sent their isolated halos; P3/P5/P9 have not yet reached their halo-count/boundary output, matching their much larger particle shards. NewGalFinder/extraction/analysis remain pending on upstream success; no catalogues exist yet.

## Darwin-1 FoF stderr warnings classified (2026-10-21 13:10 KST)

The targeted stderr review shows only warnings that missing GAS/SINK component files are being treated as empty. This matches the NewDD extraction setup, which deliberately writes only dark-matter and stellar particles, and the FoF submission guard checks both required DM and STAR chunks. These messages are expected for this run, not fatal errors; they also reinforce the documented limitation that FoF host masses include DM+stars but omit gas and sinks. No corrective change is needed.

## Darwin-1 00977 long boundary-check interval revalidated (2026-10-21 13:11 KST)

At 13:11:16 KST, all three active FoF tasks remained running; 00984 remained resource-pending and downstream tasks remained dependency-pending. The 00977 log and FoF files are unchanged since 12:47:27, while its live MPI-step CPU counter increased from 2-05:36:44 at 13:09:58 to 2-05:41:04. The log still indicates P3/P5/P9 have not completed their current boundary-check stage. Stderr contains only the expected missing-GAS/SINK-as-empty warnings; no fatal errors or storage limit were found.

## Darwin-1 00984 FoF queue reason updated (2026-10-21 13:13 KST)

At 13:12:42 KST, `squeue` changed 00984 (`430622_3`) from `PENDING (Resources)` to `PENDING (Priority)`. `scontrol show job` now reports `StartTime=Unknown` and retains the two-A100-node, 448-GiB, 10-CPU request with no GPU/GRES. The scheduler no longer provides the earlier 06:10 forecast; the job and dependency chain were left unchanged. The other three FoF jobs remain running, and 00977's log/output have not changed since 12:47:27 despite increasing CPU use.
