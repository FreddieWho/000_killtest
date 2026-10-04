"""Task-specific numerical evidence, not a broad package test suite."""
import json,time,resource,sys
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from patsy import dmatrix
from cci_core import aggregate,communication,fit_linear
root=Path(sys.argv[1] if len(sys.argv)>1 else 'results/kang_run_03');t=time.perf_counter();a=ad.read_h5ad('data/raw/kang_2018.h5ad');obs=pd.DataFrame({'patient':a.obs.replicate.astype(str),'condition':a.obs.label.astype(str),'celltype':a.obs.cell_type.astype(str)});obs['sample']=obs.patient+'__'+obs.condition
profile=[dict(stage='load',seconds=time.perf_counter()-t)]
t=time.perf_counter();ag=aggregate(a.X,obs);profile.append(dict(stage='aggregation',seconds=time.perf_counter()-t));checks=[]
checks.append(dict(check='raw_count_conservation',pass_=bool(np.allclose(ag['raw_sum'].sum(axis=0),np.asarray(a.X.sum(axis=0)).ravel())),max_abs_error=float(abs(ag['raw_sum'].sum(axis=0)-np.asarray(a.X.sum(axis=0)).ravel()).max())))
for name in ['CellChat','LIANA_CellPhoneDB']:
 r=pd.read_csv(root/f'{name}_mapped_resource.tsv',sep='\t');t=time.perf_counter();out=communication(ag,a.var_names,r,obs);profile.append(dict(stage='communication_'+name,seconds=time.perf_counter()-t,n_features=len(out['features'])))
 prior=np.load(root/f'{name}_communication.npz');v=out['scores'];w=prior['scores'];finite=np.isfinite(v)&np.isfinite(w);err=float(abs(v[finite]-w[finite]).max())
 checks.append(dict(check='optimized_equals_frozen_'+name,pass_=bool(np.allclose(v,w,rtol=1e-12,atol=1e-12,equal_nan=True) and np.array_equal(out['support'],prior['support'])),max_abs_error=err))
 samples=out['samples'];good=np.isfinite(v).all(axis=0)&(np.nanstd(v,axis=0)>1e-12);meta=pd.DataFrame({'patient':[s.split('__')[0] for s in samples],'condition':[s.split('__')[1] for s in samples]});x=dmatrix('~ condition + patient',meta,return_type='dataframe')
 t=time.perf_counter();b,se,p,df=fit_linear(v[:,good],x,x.columns.get_loc('condition[T.stim]'));profile.append(dict(stage='fit_'+name,seconds=time.perf_counter()-t,n_features=int(good.sum())))
 checks.append(dict(check='paired_df_'+name,pass_=bool(df==7),df=int(df)))
 if name=='CellChat':
  ds=sorted(meta.patient.unique());c=[samples.index(d+'__ctrl') for d in ds];s=[samples.index(d+'__stim') for d in ds];diff=v[s][:,good]-v[c][:,good];rng=np.random.default_rng(20260928);target=rng.choice(diff.shape[1],max(1,diff.shape[1]//10),replace=False);sd=np.std(diff,axis=0,ddof=1)
  truth={'n_injected':len(target),'nonzero_injected':int((sd[target]>1e-12).sum()),'all_nonzero':bool((sd[target]>1e-12).all())};(root/'power_truth_audit.json').write_text(json.dumps(truth,indent=2))
  checks.append(dict(check='prototype_power_nonzero_effects',pass_=truth['all_nonzero'],**truth))
for row in profile:row['peak_rss_mb_process_final']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
pd.DataFrame(profile).to_csv('results/prototype_current_profile.tsv',sep='\t',index=False)
Path('results/algebra_checks.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
if not all(c['pass_'] for c in checks):raise RuntimeError('A required algebra or power truth check failed')
