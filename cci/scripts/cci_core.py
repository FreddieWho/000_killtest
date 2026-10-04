"""Minimal CCC algebra. Cells are measurements; sample rows are inference units."""
from dataclasses import dataclass
import time
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree
from scipy.optimize import milp, Bounds, LinearConstraint
import networkx as nx


def normalize_counts(x, library_size=None):
    x = sparse.csr_matrix(x, dtype=np.float64)
    lib = (np.asarray(x.sum(axis=1)).ravel() if library_size is None
           else np.asarray(library_size, dtype=float))
    if lib.shape != (x.shape[0],) or not np.isfinite(lib).all() or (lib < 0).any():
        raise ValueError('library_size must contain one finite nonnegative total per cell')
    y = sparse.diags(10000 / np.maximum(lib, 1)) @ x
    y = y.tocsr(); y.data = np.log1p(y.data)
    return y, lib


def aggregate(x, obs, library_size=None):
    """All measured genes; raw sums/library size and log1p(CP10k) moments."""
    y, lib = normalize_counts(x, library_size)
    keys = pd.MultiIndex.from_frame(obs[['sample','celltype']].astype(str))
    codes, levels = pd.factorize(keys, sort=True)
    g = sparse.csr_matrix((np.ones(len(codes)), (codes,np.arange(len(codes)))),
                          shape=(len(levels), len(codes)))
    n = np.bincount(codes).astype(float)
    sums = (g @ y).toarray(); means = sums/n[:,None]
    ss = (g @ y.multiply(y)).toarray()
    var = np.maximum(ss-sums*means,0)/np.maximum(n[:,None]-1,1)
    binary = y.copy(); binary.data[:] = 1
    return dict(groups=pd.DataFrame(levels.tolist(),columns=['sample','celltype']),
                mean=means, variance=var, detection=(g@binary).toarray()/n[:,None],
                raw_sum=(g@sparse.csr_matrix(x)).toarray(),library_size=np.asarray(g@lib),
                n_cells=n, normalized=y, group_operator=g)


def mapped_resource(resource, genes):
    index = {str(g):i for i,g in enumerate(genes)}
    keep=[];lig=[];rec=[]
    for i,r in resource.iterrows():
        l=str(r.ligand_subunits).split('|'); t=str(r.receptor_subunits).split('|')
        if all(g in index for g in l+t):
            keep.append(i);lig.append([index[g] for g in l]);rec.append([index[g] for g in t])
    return resource.loc[keep].reset_index(drop=True),lig,rec


def components(matrix, lists):
    if sparse.issparse(matrix):
        matrix=matrix.tocsc()
        cache={}
        for j in lists:
            key=tuple(j)
            if key not in cache:cache[key]=matrix[:,j].min(axis=1).toarray().ravel()
        return np.column_stack([cache[tuple(j)] for j in lists])
    return np.stack([matrix[:,j].min(axis=1) for j in lists],axis=1)


def opportunity_graph(coords, radius, kernel='contact', scale=100.):
    """KD-tree sparse graph, no self loops or dense N by N intermediates."""
    coords=np.asarray(coords,float)
    pairs=cKDTree(coords).query_pairs(radius,output_type='ndarray')
    if len(pairs)==0:return sparse.csr_matrix((len(coords),len(coords)))
    d=np.linalg.norm(coords[pairs[:,0]]-coords[pairs[:,1]],axis=1)
    w=np.ones(len(d)) if kernel=='contact' else np.exp(-d/scale)
    i=np.r_[pairs[:,0],pairs[:,1]];j=np.r_[pairs[:,1],pairs[:,0]]
    return sparse.csr_matrix((np.r_[w,w],(i,j)),shape=(len(coords),len(coords)))


def aggregate_opportunity(graph, L, R, types, labels):
    """One receiver pass; group sender sums with a sparse indicator matrix."""
    codes=np.array([labels.index(t) for t in types]);n=np.bincount(codes,minlength=len(labels))
    group=sparse.csr_matrix((np.ones(len(types)),(codes,np.arange(len(types)))),shape=(len(labels),len(types)))
    result=np.full((len(labels),len(labels),L.shape[1]),np.nan)
    for b in range(len(labels)):
        mask=codes==b
        if not mask.any():continue
        for k in range(0,L.shape[1],128):
            values=group@(L[:,k:k+128]*(graph[:,mask]@R[mask,k:k+128]))
            result[:,b,k:k+128]=values/np.maximum(n[:,None]*n[b],1)
    return result


def communication(agg, genes, resource, obs, mode='mixing', coords=None,
                  radius=110., kernel='contact', cutoff=.1, min_cells=5):
    """Identical sample x sender x receiver x LR schema for both modalities.

    Mixing Omega=1/(n_sender*n_receiver). Spatial Omega=K/(n_sender*n_receiver).
    Graph operation is O(K_LR*|E|), not O(|E|) independently of resource size.
    """
    res, li,ri=mapped_resource(resource,genes)
    labels=sorted(obs.celltype.astype(str).unique());samples=sorted(obs['sample'].astype(str).unique())
    lm=components(agg['mean'],li);rm=components(agg['mean'],ri)
    ld=components(agg['detection'],li);rd=components(agg['detection'],ri)
    # Complexes use cell-level minimum BEFORE sample aggregation in BOTH modes.
    # min(mean(subunits)) would not equal mixing of spatial cell-level components.
    cell_matrix=agg['normalized'].tocsc();complex_cache={}
    for indices,means,detect in [(li,lm,ld),(ri,rm,rd)]:
        for k,jj in enumerate(indices):
            if len(jj)>1:
                key=tuple(jj)
                if key not in complex_cache:
                    values=cell_matrix[:,jj].min(axis=1).toarray().ravel()
                    complex_cache[key]=((agg['group_operator']@values)/agg['n_cells'],(agg['group_operator']@(values>0).astype(float))/agg['n_cells'])
                means[:,k],detect[:,k]=complex_cache[key]
    gi={tuple(r):i for i,r in enumerate(agg['groups'].itertuples(index=False,name=None))}
    rows=[(a,b,k) for a in labels for b in labels for k in range(len(res))]
    features=pd.DataFrame(rows,columns=['sender','receiver','lr_index'])
    scores=np.full((len(samples),len(rows)),np.nan);support=np.zeros_like(scores,dtype=bool)
    graph_info=[]
    for si,s in enumerate(samples):
        ids=np.flatnonzero(obs['sample'].astype(str).to_numpy()==s)
        types=obs.celltype.astype(str).to_numpy()[ids]
        if mode=='spatial':
            if coords is None:raise ValueError('spatial requires observed coordinates')
            graph=opportunity_graph(np.asarray(coords)[ids],radius,kernel)
            yy=cell_matrix[ids]
            cell_l=components(yy,li);cell_r=components(yy,ri)
            spatial_scores=aggregate_opportunity(graph,cell_l,cell_r,types,labels)
            graph_info.append(dict(sample=s,n=len(ids),edges=graph.nnz,graph_bytes=graph.data.nbytes+graph.indices.nbytes+graph.indptr.nbytes))
        for ai,a in enumerate(labels):
            for bi,b in enumerate(labels):
                off=(ai*len(labels)+bi)*len(res);sl=slice(off,off+len(res))
                ia=gi.get((s,a));ib=gi.get((s,b))
                if ia is None or ib is None or min(agg['n_cells'][ia],agg['n_cells'][ib])<min_cells:continue
                if mode=='mixing':score=lm[ia]*rm[ib]
                elif mode=='spatial':
                    score=spatial_scores[ai,bi]
                else:raise ValueError(mode)
                scores[si,sl]=score
                support[si,sl]=(ld[ia]>=cutoff)&(rd[ib]>=cutoff)&(score>0)
    return dict(samples=samples,labels=labels,resource=res,features=features,scores=scores,
                support=support,graphs=graph_info,mode=mode)


def exact_cover(ligands,receptors,profile=None):
    start=time.perf_counter()
    g=nx.Graph();left={('L',str(x)) for x in ligands}
    g.add_edges_from((('L',str(l)),('R',str(r))) for l,r in zip(ligands,receptors))
    if profile is not None:profile['graph_seconds']=profile.get('graph_seconds',0.)+time.perf_counter()-start
    if not len(g):return 0,[]
    start=time.perf_counter()
    matching=nx.algorithms.bipartite.maximum_matching(g,top_nodes=left)
    cover=nx.algorithms.bipartite.to_vertex_cover(g,matching,top_nodes=left)
    if not all(('L',str(l)) in cover or ('R',str(r)) in cover for l,r in zip(ligands,receptors)):
        raise RuntimeError('uncovered LR edge')
    if len(cover)!=len(matching)//2:raise RuntimeError('Konig mismatch')
    if profile is not None:profile['solve_seconds']=profile.get('solve_seconds',0.)+time.perf_counter()-start
    return len(cover),sorted(f'{a}:{b}' for a,b in cover)


def exact_hitting_set(target_sets,time_limit=30,profile=None):
    start=time.perf_counter()
    sets=[set(s) for s in target_sets]
    if not sets:return dict(value=0,targets=[],status='OPTIMAL',gap=0.)
    genes=sorted(set.union(*sets));ix={g:i for i,g in enumerate(genes)}
    rows=[];cols=[]
    for i,s in enumerate(sets):
        for g in s:rows.append(i);cols.append(ix[g])
    a=sparse.csc_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(sets),len(genes)))
    if profile is not None:profile['graph_seconds']=profile.get('graph_seconds',0.)+time.perf_counter()-start
    start=time.perf_counter()
    r=milp(np.ones(len(genes)),integrality=np.ones(len(genes)),bounds=Bounds(0,1),
           constraints=LinearConstraint(a,1,np.inf),options={'time_limit':time_limit,'mip_rel_gap':0})
    if profile is not None:profile['solve_seconds']=profile.get('solve_seconds',0.)+time.perf_counter()-start
    if not r.success:return dict(value=None,targets=[],status=str(r.message),gap=getattr(r,'mip_gap',None))
    targets=[g for g,v in zip(genes,r.x) if v>.5]
    if not all(set(targets)&s for s in sets):raise RuntimeError('hitting set feasibility')
    return dict(value=len(targets),targets=targets,status='OPTIMAL',gap=float(r.mip_gap))


def fit_linear(y,design,coefficient):
    """OLS t inference across independent sample rows; paired fixed effects allowed."""
    x=np.asarray(design,float);y=np.asarray(y,float)
    rank=np.linalg.matrix_rank(x);df=x.shape[0]-rank
    if rank<x.shape[1] or df<=0:raise ValueError('nonidentifiable design or no residual df')
    b=np.linalg.lstsq(x,y,rcond=None)[0];e=y-x@b
    scale=(e*e).sum(axis=0)/df
    se=np.sqrt(np.linalg.inv(x.T@x)[coefficient,coefficient]*scale)
    with np.errstate(divide='ignore',invalid='ignore'):t=b[coefficient]/se
    p=2*stats.t.sf(abs(t),df);p=np.where(np.isfinite(p),p,1.)
    return b[coefficient],se,p,df
