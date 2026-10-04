from pathlib import Path
import anndata as ad
import pandas as pd
from scipy import io
from cci_core import normalize_counts
out=Path('data/baseline_inputs_cellchat');out.mkdir(exist_ok=True)
a=ad.read_h5ad('data/raw/kang_2018.h5ad');a.X,_=normalize_counts(a.X)
b=ad.read_h5ad('data/baseline_inputs/all_samples.h5ad',backed='r')
keep=set(b.var_names)|set(Path('data/resources/cellchat_cofactor_genes.txt').read_text().splitlines())
a=a[:,a.var_names.isin(keep)].copy()
io.mmwrite(out/'expression.mtx',a.X.T);pd.Series(a.var_names).to_csv(out/'genes.txt',index=False,header=False)
print(a.shape,'includes required measured cofactor genes')
