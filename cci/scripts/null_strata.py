import itertools,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
root=Path(sys.argv[1]);rows=[]
for name in ['CellChat','LIANA_CellPhoneDB']:
 z=np.load(root/f'{name}_communication.npz');ss=z['samples'].tolist();ds=sorted({s.split('__')[0] for s in ss});ctrl=[ss.index(d+'__ctrl') for d in ds];stim=[ss.index(d+'__stim') for d in ds]
 y=z['scores'];keep=np.isfinite(y).all(axis=0)&(np.nanstd(y,axis=0)>1e-12);dif=y[stim][:,keep]-y[ctrl][:,keep]
 signs=np.array(list(itertools.product([-1.,1.],repeat=8)));mu=signs@dif/8;v=np.maximum(((dif*dif).sum(axis=0)-8*mu*mu)/7,0)
 with np.errstate(divide='ignore',invalid='ignore'):t=mu/np.sqrt(v/8)
 p=np.nan_to_num(2*stats.t.sf(abs(t),7),nan=1);ep=stats.rankdata(-np.round(abs(mu),12),axis=0,method='max')/256
 n=(abs(dif)>1e-12).sum(axis=0)
 for count in range(9):
  ix=n==count
  if ix.any():rows.append(dict(resource=name,nonzero_paired_donors=count,n_features=int(ix.sum()),sample_t_FPR=float((p[:,ix]<.05).mean()),exact_signflip_FPR=float((ep[:,ix]<.05).mean())))
pd.DataFrame(rows).to_csv(root/'null_information_strata.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
