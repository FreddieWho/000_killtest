from pathlib import Path
import h5py,numpy as np,pandas as pd,sys
from scipy import sparse,io
from cci_core import normalize_counts
p=Path('data/spatial_baseline_inputs');p.mkdir(exist_ok=True)
with h5py.File('data/raw/visium_filtered_feature_bc_matrix.h5') as f:
 m=f['matrix'];x=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr();genes=np.array([s.decode() for s in m['features/name'][:]]);barcodes=[s.decode() for s in m['barcodes'][:]]
u,inv=np.unique(genes,return_inverse=True);x=x@sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),inv)),shape=(len(genes),len(u)));genes=u
x,_=normalize_counts(x);r=pd.read_csv('data/resources/cellchat_mouse.tsv',sep='\t');keep=set(Path('data/resources/cellchat_mouse_cofactor_genes.txt').read_text().splitlines())
for c in ['ligand_subunits','receptor_subunits']:
 for v in r[c]:keep.update(str(v).split('|'))
ix=np.isin(genes,list(keep));io.mmwrite(p/'expression.mtx',x[:,ix].T);pd.Series(genes[ix]).to_csv(p/'genes.txt',index=False,header=False)
meta=pd.read_csv(sys.argv[1] if len(sys.argv)>1 else 'results/spatial_run_04/spots.tsv',sep='\t',index_col=0).loc[barcodes];meta.to_csv(p/'metadata.tsv',sep='\t')
print(len(barcodes),'spots',int(ix.sum()),'LR/cofactor genes')
