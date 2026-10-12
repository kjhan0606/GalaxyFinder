#!/usr/bin/env python3
"""Summarize per-halo read-only SSP audits generated on Lageunha.

Usage: python scripts/aggregate_clipping.py /path/to/audit/results
The per-halo input JSON files require the original particle catalogues; the
released TSV/JSON files can be read directly without rerunning the audit.
"""
import csv,json,sys
from pathlib import Path
import numpy as np
SUPP=Path(__file__).resolve().parents[1]/'supplement'
def main():
 source=Path(sys.argv[1])
 meta=list(csv.DictReader([s for s in (SUPP/'halo_metadata.tsv').read_text().splitlines() if s and not s.startswith('#')],delimiter='\t'))
 rows=[]
 for r in meta:
  d=json.loads((source/f"clip_{r['step']}_{r['nstar']}.json").read_text());out={k:d[k] for k in ['step','nstar','z','n_input']}
  for kind in ['young','old','metal_low','metal_high']:
   out[kind+'_count']=d['counts'][kind]
   for j,band in enumerate(['V','g']):out[kind+'_flux_fraction_'+band]=d['flux_fractions'][kind][j]
  rows.append(out)
 with (SUPP/'ssp_clipping.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 summary=dict(n_halos=len(rows),n_input_stars=sum(r['n_input'] for r in rows),counts={k:sum(r[k+'_count'] for r in rows) for k in ['young','old','metal_low','metal_high']},percentiles={k:np.quantile([r[k] for r in rows],[.16,.5,.84]).tolist() for k in rows[0] if 'fraction' in k},maximum={k:max(r[k] for r in rows) for k in rows[0] if 'fraction' in k},per_output={})
 for step in sorted(set(r['step'] for r in rows)):
  rr=[r for r in rows if r['step']==step];summary['per_output'][step]={k:float(np.median([r[k] for r in rr])) for k in rr[0] if 'fraction' in k}
 (SUPP/'ssp_clipping_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print('Summarized',len(rows),'particle audits.')
if __name__=='__main__':main()
