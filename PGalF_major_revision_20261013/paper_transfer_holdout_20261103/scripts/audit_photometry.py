"""Recover archived flux normalization by reproducing every SSP grid node."""
import os
# Set SPS_HOME to the FSPS stellar-library directory before invoking this script.
import json
import hashlib
import sys
from pathlib import Path
import numpy as np
import fsps
import fsps._fsps as fsps_extension
from astropy.cosmology import FlatLambdaCDM

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT/'supplement'
sp = fsps.StellarPopulation(zcontinuous=1, imf_type=1, sfh=0, dust_type=0,
                          dust2=0., add_neb_emission=False)
libraries = [value.decode().strip() if isinstance(value, bytes)
             else str(value).strip() for value in sp.libraries]
build = dict(libraries=libraries, extension_path=fsps_extension.__file__,
             extension_sha256=hashlib.sha256(Path(fsps_extension.__file__).read_bytes()).hexdigest(),
             python_version=sys.version,
             population_initial_parameters={key: sp.params[key] for key in sp.params.all_params},
             parameter_note="Initial population parameters. The grid loop sets logzsol and supplies tage and redshift to get_mags. zcontinuous=1 and compute_vega_mags=False.")
(INPUT/'fsps_build.json').write_text(json.dumps(build, indent=2)+'\n')
sim = FlatLambdaCDM(H0=67.4, Om0=.315)
# Reproduce FSPS src/get_lumdist.f90 and the interpolation table in sps_setup.
scale = np.linspace(1./1001., 1., 500)
ztable = 1./scale - 1.
clight = 2.9979e18  # Angstrom/s, exact constant in sps_vars.f90
dl = []
for z in ztable:
    zz = np.arange(1, 10001)/10000.*z
    dl.append(np.trapezoid(1./np.sqrt(.27*(1+zz)**3+.73), zz)
              *(1+z)*clight/1e13/72.*1e6)
records = []
for name in ['fsps_rest_z0.npz'] + [f'fsps_z{z}.npz' for z in
                                 ['0.190','0.227','0.500','1.002','1.990','3.002']]:
    path = INPUT/name
    with np.load(path) as grid:
        z = float(grid['z'])
        restored = np.empty_like(grid['flux'])
        for j, lz in enumerate(grid['logz']):
            sp.params['logzsol'] = float(lz)
            for i, age in enumerate(grid['age']):
                mag = sp.get_mags(tage=float(age), redshift=z, bands=['v','sdss_g'])
                restored[i,j] = 10.**(-.4*mag)
        delta = -2.5*np.log10(restored/grid['flux'])
    maxerr = float(np.max(np.abs(delta)))
    record = dict(file=name, z=z, max_abs_magnitude_difference=maxerr,
                  sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    if z > 0:
        fsps_pc = float(np.interp(z, ztable[::-1], np.asarray(dl)[::-1]))
        sim_pc = float(sim.luminosity_distance(z).value)*1e6
        record.update(fsps_distance_pc=fsps_pc, simulation_distance_pc=sim_pc,
                      flux_factor=(fsps_pc/sim_pc)**2,
                      magnitude_shift=5*np.log10(sim_pc/fsps_pc))
    print(record, flush=True)
    records.append(record)
report = dict(python_fsps_version=fsps.__version__,
              compiled_libraries=libraries,
              recovered_by_grid_reproduction=True,
              grid_reproduction_passed=all(r['max_abs_magnitude_difference']<1e-6 for r in records),
              records=records)
out = ROOT/'supplement/photometry_audit.json'
out.write_text(json.dumps(report, indent=2)+'\n')
if not report['grid_reproduction_passed']:
    raise SystemExit('Archived grids differ from the installed FSPS. Do not assume a distance correction.')
