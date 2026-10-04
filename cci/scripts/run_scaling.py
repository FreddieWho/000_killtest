import time,json,itertools,resource,sys
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import stats
from statsmodels.stats.multitest import multipletests
from cci_core import aggregate,mapped_resource
OUT=Path(sys.argv[1]);OUT.mkdir(exist_ok=False)
a=ad.read_h5ad('data/raw/kang_2018.h5ad');obs=pd.DataFrame({'celltype':a.obs.cell_type.astype(str),'patient':a.obs.replicate.astype(str),'condition':a.obs.label.astype(str)});obs['sample']=obs.patient+'__'+obs.condition
r=pd.read_csv('data/resources/cellchat_human.tsv',sep='\t');r=r[r.annotation.isin(['Secreted Signaling','Cell-Cell Contact'])];r=r[~r.ligand_subunits.str.contains('|',regex=False)&~r.receptor_subunits.str.contains('|',regex=False)].drop_duplicates(['ligand_subunits','receptor_subunits']).reset_index(drop=True)
r,li,ri=mapped_resource(r,a.var_names);li=np.array(li).ravel();ri=np.array(ri).ravel()
labels=sorted(obs.celltype.value_counts().head(6).index);donors=sorted(obs.patient.unique());samples=[d+'__'+c for d in donors for c in ['ctrl','stim']]
signs=np.array(list(itertools.product([-1.,1.],repeat=len(donors))));rng=np.random.default_rng(20260928);outer=rng.choice(256,100,replace=False);signs=signs[outer]
w=np.zeros((100,16));w[:,::2]=(signs<0);w[:,1::2]=(signs>0)
# Ordered cell lists allow nested stratified subsampling and preserve whole-cell covariance.
groups=[rng.permutation(v) for v in obs.groupby(['sample','celltype'],sort=True).indices.values()]
cache=[];profiles=[]
sizes=[len(obs)] if '--full-only' in sys.argv else [2000,5000,10000,len(obs)]
for size in sizes:
 frac=np.array([len(g)*size/len(obs) for g in groups]);take=np.floor(frac).astype(int)
 for j in np.argsort(-(frac-take))[:size-take.sum()]:take[j]+=1
 ids=np.sort(np.concatenate([v[:n] for v,n in zip(groups,take)]));np.save(OUT/f'cells_{size}.npy',ids)
 t=time.perf_counter();ag=aggregate(a.X[ids],obs.iloc[ids]);elapsed=time.perf_counter()-t
 gi={tuple(row):i for i,row in enumerate(ag['groups'].itertuples(index=False,name=None))}
 means=np.zeros((16,len(labels),a.n_vars));var=np.zeros_like(means);nn=np.zeros((16,len(labels)))
 for i,s in enumerate(samples):
  for j,c in enumerate(labels):
   index=gi.get((s,c))
   if index is not None:means[i,j]=ag['mean'][index];var[i,j]=ag['variance'][index];nn[i,j]=ag['n_cells'][index]
 cache.append((size,ids,ag,means,var,nn));profiles.append(dict(stage='aggregation',n_cells=size,seconds=elapsed,peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
valid_types=np.all(np.stack([c[5] for c in cache])>=5,axis=(0,1))
features=[(ia,ib,k) for ia in range(len(labels)) for ib in range(len(labels)) for k in range(len(r)) if ia!=ib and valid_types[ia] and valid_types[ib]]
ii=np.array(features);aa=ii[:,0];bb=ii[:,1];ll=li[ii[:,2]];rr=ri[ii[:,2]]
full=cache[-1][3];keep=np.std(full[:,aa,ll]*full[:,bb,rr],axis=0)>1e-12
aa=aa[keep];bb=bb[keep];ll=ll[keep];rr=rr[keep];ii=ii[keep]
pd.DataFrame({'sender':[labels[i] for i in aa],'receiver':[labels[i] for i in bb],'ligand':a.var_names[ll],'receptor':a.var_names[rr]}).to_csv(OUT/'common_features.tsv',sep='\t',index=False)
allrows=[]
def add(ps,method,size):
 for q,v in enumerate(ps):
  reject=multipletests(v,method='fdr_bh')[0]
  allrows.append(dict(method=method,n_cells=size,replicate=int(outer[q]),n_features=len(v),FPR=float((v<.05).mean()),FDP=float(reject.any()),n_BH_discoveries=int(reject.sum())))
 probs=np.linspace(.001,.999,999);pd.DataFrame({'expected':probs,'observed':np.quantile(ps,probs)}).to_csv(OUT/f'QQ_{method}_{size}.tsv',sep='\t',index=False)
for size,ids,ag,means,var,nn in cache:
 t=time.perf_counter();c=means[:,aa,ll]*means[:,bb,rr];dif=c[1::2]-c[::2];mu=signs@dif/8;vv=((dif*dif).sum(axis=0)-8*mu*mu)/7
 with np.errstate(divide='ignore',invalid='ignore'):tt=mu/np.sqrt(np.maximum(vv,0)/8)
 pv=np.nan_to_num(2*stats.t.sf(abs(tt),7),nan=1);add(pv,'paired_sample_t',size)
 profiles.append(dict(stage='sample_fit_100_null_after_aggregation',n_cells=size,seconds=time.perf_counter()-t,n_features=len(aa)))
 # Delta-method pooled means, explicitly ignoring donor clustering.
 def pool(gene,types,W):
  n=nn[:,types];m=means[:,types,gene];v=var[:,types,gene]
  total=W@n;sm=W@(n*m);ss=W@((n-1)*v+n*m*m);mu=sm/total
  vm=np.maximum((ss-sm*mu)/(total-1),0)/total
  return mu,vm
 l1,vl1=pool(ll,aa,w);r1,vr1=pool(rr,bb,w);l0,vl0=pool(ll,aa,1-w);r0,vr0=pool(rr,bb,1-w)
 observed=l1*r1-l0*r0
 se=np.sqrt(r1*r1*vl1+l1*l1*vr1+r0*r0*vl0+l0*l0*vr0)
 with np.errstate(divide='ignore',invalid='ignore'):z=observed/se
 add(np.nan_to_num(2*stats.norm.sf(abs(z)),nan=1),'pooled_cell_delta',size)
 # Actual condition-label permutations within each cell type (NOT a CPDB specificity test).
 t=time.perf_counter();genes=np.unique(np.r_[ll,rr]);lookup={g:i for i,g in enumerate(genes)};lig=np.array([lookup[g] for g in ll]);rec=np.array([lookup[g] for g in rr])
 xx=ag['normalized'][:,genes].toarray();types=obs.iloc[ids].celltype.to_numpy();blocks=[xx[types==c] for c in labels];n1=(w@nn).astype(int);n0=(1-w)@nn
 sums=[b.sum(axis=0) for b in blocks];exceed=np.zeros_like(observed,dtype=int)
 for perm in range(199):
  m1=np.zeros((100,len(labels),len(genes)));m0=np.zeros_like(m1)
  for j,b in enumerate(blocks):
   if not valid_types[j]:continue
   cumul=np.cumsum(b[rng.permutation(len(b))],axis=0);total=cumul[n1[:,j]-1]
   m1[:,j]=total/n1[:,j,None];m0[:,j]=(sums[j]-total)/n0[:,j,None]
  delta=m1[:,aa,lig]*m1[:,bb,rec]-m0[:,aa,lig]*m0[:,bb,rec]
  exceed+=abs(delta)>=abs(observed)-1e-12
 pp=(exceed+1)/200;add(pp,'pooled_cell_condition_permutation',size)
 profiles.append(dict(stage='pooled_condition_permutation_100x199',n_cells=size,seconds=time.perf_counter()-t,n_features=len(aa),peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
 pd.DataFrame(profiles).to_csv(OUT/'profile.tsv',sep='\t',index=False);pd.DataFrame(allrows).to_csv(OUT/'null_replicates.tsv',sep='\t',index=False)
 print('DONE size',size,'features',len(aa),flush=True)
(OUT/'receipt.json').write_text(json.dumps(dict(status='DONE',outer_null_assignments=100,inner_cell_permutations=199,seed=20260928,common_features=len(aa),cell_types_used=[c for c,v in zip(labels,valid_types) if v],excluded_types_due_smallest_scale=[c for c,v in zip(labels,valid_types) if not v],scope='Six largest original cell types; common all-scale complete-case simple inter-type LR universe. Real whole cells subsampled; no synthetic duplicates. Conditional paired label randomization, not independent null cohorts.'),indent=2))
