import os,sys,time,json,resource
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
method=sys.argv[1];out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
inputs=Path('data/baseline_inputs');samples=sorted(p.stem for p in inputs.glob('patient*.h5ad'));rows=[];timings=[]
if method=='CellPhoneDB':from cellphonedb.src.core.methods import cpdb_statistical_analysis_method as cpdb
elif method=='FastCCC':from fastccc.core import statistical_analysis_method as fastccc
elif method=='LIANA':import liana as li
for sample in samples:
 d=out/sample;d.mkdir();t=time.perf_counter()
 if method=='CellPhoneDB':
  r=cpdb.call(cpdb_file_path='data/resources/cellphonedb_v5.zip',meta_file_path=str(inputs/f'{sample}.tsv'),counts_file_path=str(inputs/f'{sample}.h5ad'),counts_data='hgnc_symbol',output_path=str(d),iterations=1000,threads=1,debug_seed=20260928,result_precision=8)
  v=r['means'];pair_cols=[c for c in v if '|' in c]
  z=v[['id_cp_interaction']+pair_cols].melt(id_vars='id_cp_interaction',var_name='pair',value_name='score')
  z['feature']=z.id_cp_interaction.astype(str)+'::'+z.pair
 elif method=='FastCCC':
  strength,pvals,percents=fastccc(database_file_path='data/resources/cpdb_v5',celltype_file_path=str(inputs/f'{sample}.tsv'),counts_file_path=str(inputs/f'{sample}.h5ad'),single_unit_summary='Mean',complex_aggregation='Minimum',LR_combination='Arithmetic',min_percentile=.1,save_path=str(d))
  strength.to_csv(d/'strength.tsv',sep='\t');pvals.to_csv(d/'pvalues.tsv',sep='\t')
  z=strength.stack(dropna=False).reset_index();z.columns=['pair','interaction','score'];z['feature']=z.interaction.astype(str)+'::'+z.pair.astype(str)
 elif method=='LIANA':
  a=ad.read_h5ad(inputs/f'{sample}.h5ad')
  li.mt.cellphonedb(a,groupby='celltype',resource_name='cellphonedb',use_raw=False,n_perms=1000,n_jobs=1,seed=20260928,return_all_lrs=True,min_cells=5,expr_prop=.1,inplace=True)
  z=a.uns['liana_res'];z.to_csv(d/'liana_native.tsv',sep='\t',index=False)
  z['feature']=z.source.astype(str)+'|'+z.target.astype(str)+'::'+z.ligand_complex.astype(str)+'|'+z.receptor_complex.astype(str);z['score']=z.lr_means
 z=z[['feature','score']].copy();z['sample']=sample;rows.append(z)
 elapsed=time.perf_counter()-t;timings.append(dict(method=method,sample=sample,seconds=elapsed,peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,n_scores=len(z)))
 pd.DataFrame(timings).to_csv(out/'profile.tsv',sep='\t',index=False)
 pd.concat(rows).to_csv(out/'scores.tsv.gz',sep='\t',index=False)
 print(method,sample,elapsed,flush=True)
(out/'receipt.json').write_text(json.dumps({'status':'DONE','method':method,'samples':len(samples),'iterations':0 if method=='FastCCC' else 1000,'score_null':'sample-level downstream only; native specificity p values are not condition-difference p values'},indent=2))
