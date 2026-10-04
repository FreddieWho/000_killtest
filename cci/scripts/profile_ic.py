from pathlib import Path
import numpy as np,pandas as pd
from cci_core import exact_cover,exact_hitting_set
root=Path('results/kang_run_03');programs=pd.read_csv(root/'programs.tsv',sep='\t');rows=[]
for name,p in programs.groupby('database'):
 z=np.load(root/f'{name}_communication.npz');r=pd.read_csv(root/f'{name}_mapped_resource.tsv',sep='\t');f=pd.read_csv(root/f'{name}_features.tsv',sep='\t');labels=f.sender.unique().tolist();samples=z['samples'].tolist();k=len(r);simple=~(r.ligand_subunits.str.contains('|',regex=False)|r.receptor_subunits.str.contains('|',regex=False));timing={}
 for row in p.itertuples():
  i=(labels.index(row.sender)*len(labels)+labels.index(row.receiver))*k;rr=r[z['support'][samples.index(row.sample),i:i+k]&simple.to_numpy()];n,_=exact_cover(rr.ligand_subunits,rr.receptor_subunits,timing)
  assert n==row.ic100
 for stage,sec in timing.items():rows.append(dict(database=name,algorithm='exact_MVC',stage=stage,seconds=sec,n_programs=len(p)))
 if name=='CellChat':
  timing={};h=pd.read_csv(root/'complex_hitting_sets.tsv',sep='\t')
  for row in h.itertuples():
   i=(labels.index(row.sender)*len(labels)+labels.index(row.receiver))*k;rr=r[z['support'][samples.index(row.sample),i:i+k]]
   sets=[str(v.ligand_subunits).split('|')+str(v.receptor_subunits).split('|') for v in rr.itertuples()];answer=exact_hitting_set(sets,profile=timing)
   assert answer['status']=='OPTIMAL' and answer['value']==row.value
  for stage,sec in timing.items():rows.append(dict(database=name,algorithm='exact_hitting_set',stage=stage,seconds=sec,n_programs=len(h)))
pd.DataFrame(rows).to_csv('results/ic_phase_profile.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
