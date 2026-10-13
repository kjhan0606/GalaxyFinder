"""Record the bounded producing configuration without catalogues or credentials."""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
import sys,json,hashlib,platform,socket
from pathlib import Path
import numpy,scipy
ROOT=Path(__file__).resolve().parents[1]
def main():
    files=[ROOT/'scripts/measure_halo.py',ROOT/'scripts/pilot_properties.py',ROOT/'scripts/source/sb_vs_b.py',
           ROOT/'scripts/source/sb_revise.py',ROOT.parent/'opFoF/variant/libopfof_sb.so',
           ROOT/'inputs/supplement/photometry_audit.json',ROOT/'inputs/supplement/halo_metadata.tsv']
    config=dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,platform=platform.platform(),host=socket.gethostname(),
        base_manuscript_commit=(ROOT/'base_commit.txt').read_text().strip(),threads=1,
        grid=[.02,.03,.05,.06,.07,.08,.09,.10,.11,.12,.13,.14,.15,.18,.20],
        sampling=dict(factors=[2,4,8],mass_flux_reweighting='inverse probability',nested=True,seed_sequence=['20261013','integer output number','catalogue total star count']),
        boundary_width_arcsec_per_side=[1,2],controls=['primary','threshold_bright','threshold_faint','wing_off','light_peak'],
        distance='simulation cosmology for both angular and luminosity distances',
        source_sha256={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else 'original linking shared library':hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (ROOT/'results/opus_calculation_config.json').write_text(json.dumps(config,indent=2)+'\n')
    print('Configuration recorded')
if __name__=='__main__':main()
