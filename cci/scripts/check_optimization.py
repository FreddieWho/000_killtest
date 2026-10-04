"""Conclusion-affecting regression checks against frozen data and matched fits."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
root=Path('results/optimization_20260930')
checks=[]
for name in ['CellChat','LIANA_CellPhoneDB']:
 with np.load(root/'full'/f'{name}.npz') as a, np.load(root/'stream4096'/f'{name}.npz') as b:
  for k in a.files:
   np.testing.assert_allclose(a[k],b[k],rtol=1e-12,atol=1e-12)
 checks.append(name+': all scores, support, tested feature set, beta, SE, p identical')
with np.load('results/kang_run_03/sufficient_statistics.npz') as a, np.load(root/'kang/sufficient_statistics.npz') as b:
 ix=pd.Index(a['genes']).get_indexer(b['genes']);assert (ix>=0).all()
 for k in ['mean','variance','detection','raw_sum']:
  np.testing.assert_allclose(a[k][:,ix],b[k],rtol=1e-12,atol=1e-12)
 for k in ['library_size','n_cells']:
  np.testing.assert_allclose(a[k],b[k],rtol=1e-12,atol=1e-12)
checks.append('All retained gene moments/raw sums and full-library totals preserved')
assert json.loads((root/'kang/receipt.json').read_text())['status']=='DONE'
for n in ['programs.tsv','null_replicates.tsv','paired_fit.tsv','factorial_fit.tsv','synthetic_score_power.tsv']:
 a=pd.read_csv(Path('results/kang_run_03')/n,sep='\t');b=pd.read_csv(root/'kang'/n,sep='\t')
 # Minimum covers may choose different equally optimal targets: compare counts, not arbitrary representatives.
 if 'targets' in a:a=a.drop(columns='targets');b=b.drop(columns='targets')
 pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=1e-12,atol=1e-12)
checks.append('Full streamed Kang: IC programs, null, paired fit, factorial and synthetic power unchanged')
(root/'checks.json').write_text(json.dumps({'status':'PASS','checks':checks},indent=2)+'\n')
print(json.dumps(checks))
