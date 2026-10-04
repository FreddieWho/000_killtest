import time,json,resource,sys
from pathlib import Path
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.spatial import distance,cKDTree
from cci_core import aggregate,communication,opportunity_graph,components,mapped_resource,exact_cover,aggregate_opportunity
OUT=Path(sys.argv[1] if len(sys.argv)>1 else 'results/spatial_run_01');OUT.mkdir(exist_ok=False)
profiles=[]
def rec(stage,start,**kw):
 profiles.append(dict(stage=stage,seconds=time.perf_counter()-start,peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,**kw))
 pd.DataFrame(profiles).to_csv(OUT/'profile.tsv',sep='\t',index=False)
 print(stage,profiles[-1]['seconds'],flush=True)
t=time.perf_counter()
with h5py.File('data/raw/visium_filtered_feature_bc_matrix.h5') as f:
 m=f['matrix'];x=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr()
 genes=np.array([v.decode() for v in m['features/name'][:]]);barcodes=[v.decode() for v in m['barcodes'][:]]
# Collapse identical symbols deterministically; no re-clustering or deconvolution.
u,inv=np.unique(genes,return_inverse=True)
if len(u)!=len(genes):x=x@sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),inv)),shape=(len(genes),len(u)));genes=u
pos=pd.read_csv('data/resources/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','row','col']).set_index('barcode').loc[barcodes]
anno=pd.read_csv('data/resources/visium_annotation.tsv',sep='\t').set_index('barcode')
labels=anno.celltype.reindex(barcodes).fillna('UNANNOTATED')
obs=pd.DataFrame({'sample':'visium_anterior1','celltype':labels.to_numpy()},index=barcodes)
scales=json.load(open('data/resources/scalefactors_json.json'))
pixels=pos[['row','col']].to_numpy()
pixel_pitch=float(np.median(cKDTree(pixels).query(pixels,k=2)[0][:,1]))
coords=pixels*(100./pixel_pitch)
nearest=cKDTree(coords).query(coords,k=2)[0][:,1]
obs.assign(x_um=coords[:,0],y_um=coords[:,1]).to_csv(OUT/'spots.tsv',sep='\t')
rec('read_load',t,n_spots=len(obs))
input_info=dict(spots=len(obs),genes=len(genes),annotated=int((labels!='UNANNOTATED').sum()),unannotated=int((labels=='UNANNOTATED').sum()),annotation_file_spots=len(anno),annotation_barcode_not_in_counts=int((~anno.index.isin(barcodes)).sum()),median_nearest_distance_um=float(np.median(nearest)),spot_diameter_um=55,spot_diameter_fullres_pixels=scales['spot_diameter_fullres'],source='https://www.10xgenomics.com/support/cn/software/space-ranger/3.1/getting-started/space-ranger-glossary',pixel_pitch=pixel_pitch,coordinate_scale_um_per_pixel=100./pixel_pitch,nearest_um_using_55_over_json_diameter=pixel_pitch*55/scales['spot_diameter_fullres'],nearest_um_using_tutorial_65=pixel_pitch*65/scales['spot_diameter_fullres'],tutorial_difference='Official nominal v1 spot diameter is 55um and pitch100um; JSON diameter calibration does not reproduce pitch with55um. Use observed nearest-neighbor pixel pitch and official100um center spacing; retain both alternative diameter estimates.',sample_units=1)
(OUT/'input_audit.json').write_text(json.dumps(input_info,indent=2))
t=time.perf_counter();ag=aggregate(x,obs);rec('aggregation',t)
res=pd.read_csv('data/resources/cellchat_mouse.tsv',sep='\t');res=res[res.annotation.isin(['Secreted Signaling','Cell-Cell Contact'])].drop_duplicates(['ligand_subunits','receptor_subunits']).reset_index(drop=True)
outs={};programs=[]
for mode,rad,kernel in [('mixing',0,'none'),('contact',110,'contact'),('secreted',250,'exponential')]:
 t=time.perf_counter();o=communication(ag,genes,res,obs,mode='mixing' if mode=='mixing' else 'spatial',coords=coords,radius=rad,kernel=kernel);outs[mode]=o
 rec('communication',t,mode=mode,n_lr=len(o['resource']))
 np.savez_compressed(OUT/f'{mode}_communication.npz',scores=o['scores'],support=o['support'])
 rr=o['resource'];simple=~(rr.ligand_subunits.str.contains('|',regex=False)|rr.receptor_subunits.str.contains('|',regex=False));k=len(rr)
 for ai,a in enumerate(o['labels']):
  for bi,b in enumerate(o['labels']):
   sl=slice((ai*len(o['labels'])+bi)*k,(ai*len(o['labels'])+bi+1)*k)
   active=o['support'][0,sl]&simple.to_numpy();z=rr[active];ic,targets=exact_cover(z.ligand_subunits,z.receptor_subunits)
   programs.append(dict(mode=mode,sample='visium_anterior1',sender=a,receiver=b,count=len(z),strength=float(o['scores'][0,sl][active].sum()),ic100=ic,targets=';'.join(targets)))
 o['features'].to_csv(OUT/'features.tsv',sep='\t',index=False)
 rr.to_csv(OUT/'mapped_resource.tsv',sep='\t',index=False)
pd.DataFrame(programs).to_csv(OUT/'programs.tsv',sep='\t',index=False)
# Dense distance is comparator ONLY, never used by the formal sparse communication path.
rr,li,ri=mapped_resource(res,genes);L=components(ag['normalized'],li);R=components(ag['normalized'],ri)
labels_sorted=outs['mixing']['labels'];type_array=obs.celltype.to_numpy();g=sparse.csr_matrix((np.ones(len(obs)),([labels_sorted.index(c) for c in type_array],np.arange(len(obs)))),shape=(len(labels_sorted),len(obs)))
counts=np.asarray(g.sum(axis=1)).ravel()
def dense_aggregation(w):
 blocks=[]
 for b in labels_sorted:
  mask=(type_array==b);values=g@(L*(w[:,mask]@R[mask]))
  blocks.append(values/(counts[:,None]*mask.sum()))
 return np.stack(blocks,axis=1).reshape(-1)
checks=[]
for kernel,rad in [('contact',110),('exponential',250)]:
 t=time.perf_counter();w=opportunity_graph(coords,rad,kernel);rec('sparse_graph',t,mode=kernel,n_spots=len(coords),edges=w.nnz,storage_bytes=w.data.nbytes+w.indices.nbytes+w.indptr.nbytes)
 t=time.perf_counter();d=distance.cdist(coords,coords);rec('dense_distance_comparator',t,mode=kernel,n_spots=len(coords),storage_bytes=d.nbytes)
 wd=((d>0)&(d<=rad)).astype(float)
 if kernel!='contact':wd*=np.exp(-d/100.)
 t=time.perf_counter();dense=aggregate_opportunity(wd,L,R,type_array,labels_sorted).reshape(-1);rec('dense_communication_comparator',t,mode=kernel,n_spots=len(coords),n_lr=len(rr))
 t=time.perf_counter();sp=aggregate_opportunity(w,L,R,type_array,labels_sorted).reshape(-1);rec('sparse_communication',t,mode=kernel,n_spots=len(coords),n_lr=len(rr))
 sparse_values=outs['contact' if kernel=='contact' else 'secreted']['scores'][0]
 checks.append(dict(kernel=kernel,max_absolute_difference=float(np.max(abs(dense-sparse_values))),all_LR_and_all_group_pairs=True,edges=w.nnz,dense_pairs=len(coords)**2))
 del d,wd
# Constant all-to-all opportunity must recover mixing, including complex algebra.
full=dense_aggregation(np.ones((len(coords),len(coords))))
checks.append(dict(kernel='all_to_all_equals_mixing',max_absolute_difference=float(np.max(abs(full-outs['mixing']['scores'][0]))),all_LR_and_all_group_pairs=True))
(OUT/'numerical_equivalence.json').write_text(json.dumps(checks,indent=2))
# Geometry scaling: disconnected copies are computational inputs, never new replicates.
for copies in [1,2,4]:
 xy=np.vstack([coords+[i*(coords[:,0].max()-coords[:,0].min()+1000),0] for i in range(copies)])
 for kernel,rad in [('contact',110),('exponential',250)]:
  t=time.perf_counter();w=opportunity_graph(xy,rad,kernel);rec('sparse_graph_scaling',t,mode=kernel,n_spots=len(xy),edges=w.nnz,storage_bytes=w.data.nbytes+w.indices.nbytes+w.indptr.nbytes,copies=copies)
(OUT/'receipt.json').write_text(json.dumps({'status':'DONE','biological_sample_units':1,'statistical_replication':'NOT_AVAILABLE','inputs':input_info,'same_interface_resource_schema':True,'unannotated_spots_retained':True},indent=2))
print('DONE',OUT,flush=True)
