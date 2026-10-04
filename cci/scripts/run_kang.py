"""Real-data exploratory nulls, general designs, and IC; see report for estimands."""
import os,time,json,itertools,resource,sys
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import stats,sparse
from statsmodels.stats.multitest import multipletests
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from patsy import dmatrix
from cci_core import aggregate,communication,exact_cover,exact_hitting_set,fit_linear
from streaming_input import read_lr_counts

OUT=Path(sys.argv[1] if len(sys.argv)>1 else 'results/kang_run_01')
OUT.mkdir(parents=True,exist_ok=False)
rng=np.random.default_rng(20260928);timings=[]
def record(stage,start,**kw):
 timings.append(dict(stage=stage,seconds=time.perf_counter()-start,peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,**kw))
 pd.DataFrame(timings).to_csv(OUT/'profile.tsv',sep='\t',index=False)
 print(stage,timings[-1]['seconds'],flush=True)
def dump(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,default=str))

cc=pd.read_csv('data/resources/cellchat_human.tsv',sep='\t')
cc=cc[cc.annotation.isin(['Secreted Signaling','Cell-Cell Contact'])].copy()
cc=cc.drop_duplicates(['ligand_subunits','receptor_subunits']).reset_index(drop=True)
om=pd.read_csv('data/resources/liana_cellphonedb.tsv',sep='\t')
resources={'CellChat':cc,'LIANA_CellPhoneDB':om}
t=time.perf_counter()
if '--stream' in sys.argv:
 counts,raw_obs,genes,lib=read_lr_counts('data/raw/kang_2018.h5ad',resources.values())
else:
 a=ad.read_h5ad('data/raw/kang_2018.h5ad')
 counts,raw_obs,genes,lib=a.X,a.obs,a.var_names,None
obs=pd.DataFrame({'patient':raw_obs.replicate.astype(str),'condition':raw_obs.label.astype(str),'celltype':raw_obs.cell_type.astype(str)})
obs['sample']=obs.patient+'__'+obs.condition
obs.to_csv(OUT/'cell_metadata.tsv',sep='\t')
record('load',t,n_cells=len(obs))
t=time.perf_counter();agg=aggregate(counts,obs,library_size=lib);record('aggregation',t,n_cells=len(obs))
agg['groups'].assign(n_cells=agg['n_cells'],library_size=agg['library_size']).to_csv(OUT/'groups.tsv',sep='\t',index=False)
np.savez_compressed(OUT/'sufficient_statistics.npz',**{k:agg[k] for k in ['mean','variance','detection','raw_sum','library_size','n_cells']},genes=genes.to_numpy(dtype=str))
outputs={};programs=[];complex_candidates=[]
for name,res in resources.items():
 t=time.perf_counter();o=communication(agg,genes,res,obs);outputs[name]=o
 record('LR_lookup_and_mixing',t,db=name,n_features=len(o['features']))
 o['resource'].to_csv(OUT/f'{name}_mapped_resource.tsv',sep='\t',index=False)
 o['features'].to_csv(OUT/f'{name}_features.tsv',sep='\t',index=False)
 np.savez_compressed(OUT/f'{name}_communication.npz',scores=o['scores'],support=o['support'],samples=np.array(o['samples']))
 t=time.perf_counter();rr=o['resource'];simple=~(rr.ligand_subunits.str.contains('|',regex=False)|rr.receptor_subunits.str.contains('|',regex=False))
 k=len(rr)
 for si,s in enumerate(o['samples']):
  for ai,sender in enumerate(o['labels']):
   for bi,receiver in enumerate(o['labels']):
    sl=slice((ai*len(o['labels'])+bi)*k,(ai*len(o['labels'])+bi+1)*k)
    score=o['scores'][si,sl];active=o['support'][si,sl]
    if not np.isfinite(score).all():continue
    z=rr[active&simple.to_numpy()];ic,targets=exact_cover(z.ligand_subunits,z.receptor_subunits)
    prog=dict(database=name,sample=s,patient=s.split('__')[0],sender=sender,receiver=receiver,
              count=len(z),strength=float(score[active&simple.to_numpy()].sum()),ic100=ic,
              targets=';'.join(targets),role_overlap=len(set(z.ligand_subunits)&set(z.receptor_subunits)))
    programs.append(prog)
    if name=='CellChat' and (active&~simple.to_numpy()).any():
     complex_candidates.append((int(active.sum()),s,sender,receiver,rr[active]))
 record('IC_graph_and_exact_MVC',t,db=name)
p=pd.DataFrame(programs);p.to_csv(OUT/'programs.tsv',sep='\t',index=False)
# Exact gene-deduplicated hitting set, all required ligand/receptor subunits are targetable.
t=time.perf_counter();hits=[]
for count,s,sender,receiver,z in sorted(complex_candidates,key=lambda x:(-x[0],x[1],x[2],x[3]))[:100]:
 sets=[str(r.ligand_subunits).split('|')+str(r.receptor_subunits).split('|') for r in z.itertuples()]
 h=exact_hitting_set(sets);hits.append(dict(sample=s,sender=sender,receiver=receiver,channels=count,**h))
pd.DataFrame(hits).to_csv(OUT/'complex_hitting_sets.tsv',sep='\t',index=False);record('exact_hitting_set_top100',t,n_programs=len(hits))
# Leave-one-donor-out prediction, rather than evaluating the fitted curve on its own samples.
summ=[];quadrants=[];similar=[]
for name,g in p.groupby('database'):
 g=g[g['count']>0].copy()
 for feature in ['count','strength']:
  x=np.log1p(g[feature].to_numpy());y=g.ic100.to_numpy();pred=np.zeros(len(g));lin=np.zeros(len(g))
  for donor in g.patient.unique():
   test=(g.patient==donor).to_numpy();train=~test
   pred[test]=IsotonicRegression(out_of_bounds='clip').fit(x[train],y[train]).predict(x[test])
   lin[test]=LinearRegression().fit(x[train,None],y[train]).predict(x[test,None])
  summ.append(dict(database=name,predictor=feature,n_programs=len(g),pearson=float(stats.pearsonr(g[feature],y).statistic),spearman=float(stats.spearmanr(g[feature],y).statistic),heldout_donor_isotonic_R2=r2_score(y,pred),heldout_donor_loglinear_R2=r2_score(y,lin)))
  # Count similar / strength similar <=5%; IC difference >=2 and >=1.5 fold (descriptive, not a gate).
  z=g.reset_index(drop=True);vals=z[feature].to_numpy();ic=z.ic100.to_numpy();cases=[];involved=set()
  for i in range(len(z)):
   matches=np.flatnonzero((np.arange(len(z))>i)&(abs(vals-vals[i])<=.05*np.maximum(vals,vals[i]))&(abs(ic-ic[i])>=2)&(np.maximum(ic,ic[i])>=1.5*np.maximum(1,np.minimum(ic,ic[i]))))
   for j in matches:
    involved.update([i,int(j)])
    if len(cases)<30:cases.append(dict(database=name,matched_on=feature,program1='|'.join(str(z.iloc[i][c]) for c in ['sample','sender','receiver']),program2='|'.join(str(z.iloc[j][c]) for c in ['sample','sender','receiver']),value1=vals[i],value2=vals[j],ic1=ic[i],ic2=ic[j]))
  similar.extend(cases);summ[-1]['similar_value_different_IC_program_fraction']=len(involved)/len(z)
 g['quadrant']=np.char.add(np.char.add(np.where(g.strength>=g.strength.median(),'high_strength','low_strength'),' / '),np.where(g.ic100>=g.ic100.median(),'high_IC','low_IC'))
 quadrants.append(g.groupby('quadrant',group_keys=False).head(3))
pd.DataFrame(summ).to_csv(OUT/'ic_independence.tsv',sep='\t',index=False)
pd.DataFrame(similar).to_csv(OUT/'similar_value_cases.tsv',sep='\t',index=False)
pd.concat(quadrants).to_csv(OUT/'quadrants.tsv',sep='\t',index=False)
m=p[p.database=='CellChat'].merge(p[p.database=='LIANA_CellPhoneDB'],on=['sample','patient','sender','receiver'],suffixes=('_a','_b'))
valid=(m['count_a']>0)&(m['count_b']>0);m=m[valid]
dump('database_sensitivity.json',{'n_matched_programs':len(m),'IC_rank_spearman':float(stats.spearmanr(m.ic100_a,m.ic100_b).statistic),'count_rank_spearman':float(stats.spearmanr(m.count_a,m.count_b).statistic),'high_low_agreement':float(((m.ic100_a>=m.ic100_a.median())==(m.ic100_b>=m.ic100_b.median())).mean()),'note':'Native catalogues mapped to same measured genes; catalogue coverage differs. Shared LR intersection would be identical by construction.'})
m.to_csv(OUT/'database_matched_programs.tsv',sep='\t',index=False)
# Full enumeration of all 2^8 within-donor label swaps preserves original cells and all within-sample structure.
o=outputs['CellChat'];samples=o['samples'];donors=sorted(obs.patient.unique());signs=np.array(list(itertools.product([-1.,1.],repeat=len(donors))))
np.save(OUT/'null_signs.npy',signs)
ci=[samples.index(d+'__ctrl') for d in donors];si=[samples.index(d+'__stim') for d in donors]
mask=np.isfinite(o['scores']).all(axis=0)&(np.nanstd(o['scores'],axis=0)>1e-12)
y=o['scores'][:,mask];diff=y[si]-y[ci]
features=o['features'][mask].copy();features.to_csv(OUT/'inference_features.tsv',sep='\t',index=False)
meta=pd.DataFrame({'sample':samples,'patient':[s.split('__')[0] for s in samples],'condition':[s.split('__')[1] for s in samples]})
x=dmatrix('~ condition + patient',meta,return_type='dataframe');coef=x.columns.get_loc('condition[T.stim]')
t=time.perf_counter();beta,se,pvals,df=fit_linear(y,x,coef)
pd.DataFrame({'beta':beta,'SE':se,'pvalue':pvals,'FDR':multipletests(pvals,method='fdr_bh')[1]}).to_csv(OUT/'paired_fit.tsv',sep='\t',index=False)
record('paired_OLS',t,n_features=y.shape[1],n_replicates=len(donors))
def null_metrics(dif,name):
 n=dif.shape[0];means=signs@dif/n;sq=(dif*dif).sum(axis=0)
 var=np.maximum((sq[None,:]-n*means*means)/(n-1),0)
 with np.errstate(divide='ignore',invalid='ignore'):tt=means/np.sqrt(var/n)
 pv=2*stats.t.sf(abs(tt),n-1);pv=np.where(np.isfinite(pv),pv,1.)
 exact=stats.rankdata(-np.round(abs(means),12),axis=0,method='max')/len(signs)
 rows=[]
 for method,ps in [(name,pv),(name+'_exact_signflip',exact)]:
  for i,v in enumerate(ps):
   rej=multipletests(v,method='fdr_bh')[0]
   rows.append(dict(method=method,replicate=i,n_features=len(v),FPR=float((v<.05).mean()),FDP=float(rej.any()),n_BH_discoveries=int(rej.sum())))
  q=np.linspace(.001,.999,999);pd.DataFrame({'expected':q,'observed':np.quantile(ps,q)}).to_csv(OUT/f'QQ_{method}.tsv',sep='\t',index=False)
 return rows
rows=[];t=time.perf_counter()
for name,o2 in outputs.items():
 good=np.isfinite(o2['scores']).all(axis=0)&(np.nanstd(o2['scores'],axis=0)>1e-12)
 rows.extend(null_metrics(o2['scores'][si][:,good]-o2['scores'][ci][:,good],name+'_paired_product'))
pd.DataFrame(rows).to_csv(OUT/'null_replicates.tsv',sep='\t',index=False)
record('null_256_label_swaps',t,n_null=256)
# Score-level alternatives: explicitly synthetic, not independent perturbation validation.
t=time.perf_counter();alt=rng.choice(diff.shape[1],max(1,diff.shape[1]//10),replace=False);is_alt=np.zeros(diff.shape[1],bool);is_alt[alt]=True
base_sd=np.std(diff,axis=0,ddof=1);power=[]
for effect in [.5,1.,2.]:
 for i,sg in enumerate(signs):
  d=diff*sg[:,None];d[:,alt]+=effect*base_sd[alt]
  _,pv=stats.ttest_1samp(d,0,axis=0);pv=np.nan_to_num(pv,nan=1);reject=multipletests(pv,method='fdr_bh')[0]
  power.append(dict(effect_in_donor_SD=effect,replicate=i,n_features=len(pv),n_injected=int(is_alt.sum()),power=float(reject[is_alt].mean()),FDP=float(reject[~is_alt].sum()/max(1,reject.sum()))))
pd.DataFrame(power).to_csv(OUT/'synthetic_score_power.tsv',sep='\t',index=False);record('score_level_power',t)
# General factorial design with batch, independent sample rows.
metadata=pd.DataFrame(list(itertools.product([0,1],[0,1,2],['A','B'],range(4))),columns=['treatment','time','batch','rep'])
x=dmatrix('~ treatment * time + batch',metadata,return_type='dataframe');truth=np.zeros((x.shape[1],200));j=x.columns.get_loc('treatment:time');truth[j,:100]=1
ys=np.asarray(x)@truth+rng.normal(size=(len(x),200));b,s,pv,df=fit_linear(ys,x,j)
pd.DataFrame({'truth':truth[j],'beta':b,'SE':s,'pvalue':pv,'FDR':multipletests(pv,method='fdr_bh')[1]}).to_csv(OUT/'factorial_fit.tsv',sep='\t',index=False)
metadata.to_csv(OUT/'factorial_metadata.tsv',sep='\t',index=False)
dump('designs.json',{'paired_columns':list(dmatrix('~ condition + patient',meta,return_type='dataframe').columns),'paired_df':len(meta)-len(donors)-1,'factorial_columns':list(x.columns),'factorial_df':df,'factorial_rank':int(np.linalg.matrix_rank(x))})
dump('receipt.json',{'status':'DONE','seed':20260928,'cells':len(obs),'genes':len(genes),'input_mode':'stream' if '--stream' in sys.argv else 'full','donors':len(donors),'samples':len(samples),'null_assignments':256,'unique_two_sided_sign_patterns':128,'null_interpretation':'conditional randomization of observed paired profiles; not sharp-null biological data nor independent cohorts','biology_blockade':'NOT_VALIDATED','remaining':'native baselines, cell scaling, spatial, consolidated reporting'})
print('DONE',OUT,flush=True)
