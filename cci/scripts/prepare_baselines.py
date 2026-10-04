import json,time
from pathlib import Path
import anndata as ad,pandas as pd,numpy as np
from scipy import io
from cci_core import normalize_counts
out=Path('data/baseline_inputs');out.mkdir(exist_ok=True)
a=ad.read_h5ad('data/raw/kang_2018.h5ad');a.obs['celltype']=a.obs.cell_type.astype(str);a.obs['sample']=a.obs.replicate.astype(str)+'__'+a.obs.label.astype(str)
t=time.perf_counter();a.X,_=normalize_counts(a.X)
cc=pd.read_csv('data/resources/cellchat_human.tsv',sep='\t');gg=pd.read_csv('data/resources/cpdb_v5/gene_table.csv')
keep=set(gg.hgnc_symbol.dropna().astype(str))
for col in ['ligand_subunits','receptor_subunits']:
 for x in cc[col]:keep.update(str(x).split('|'))
a=a[:,a.var_names.isin(keep)].copy();a.write_h5ad(out/'all_samples.h5ad')
a.obs[['sample','celltype']].to_csv(out/'metadata.tsv',sep='\t')
io.mmwrite(out/'expression.mtx',a.X.T);pd.Series(a.var_names).to_csv(out/'genes.txt',index=False,header=False)
for sample in sorted(a.obs['sample'].unique()):
 b=a[a.obs['sample']==sample].copy();b.write_h5ad(out/f'{sample}.h5ad')
 pd.DataFrame({'Cell':b.obs_names,'cell_type':b.obs.celltype.to_numpy()}).to_csv(out/f'{sample}.tsv',sep='\t',index=False)
Path('results/baseline_input_receipt.json').write_text(json.dumps({'cells':a.n_obs,'LR_genes':a.n_vars,'samples':16,'seconds':time.perf_counter()-t,'normalization':'log1p(CP10k) computed on all 15706 measured genes BEFORE complete LR-gene subset; same normalized inputs for native baselines'},indent=2))
print(a.shape)
