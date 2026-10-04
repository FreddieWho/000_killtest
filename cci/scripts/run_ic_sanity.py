import json,sys,time
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import r2_score
from cci_core import exact_cover,exact_hitting_set
src=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=False)
p=pd.read_csv(src/'programs.tsv',sep='\t');rows=[]
for name,pp in p.groupby('database'):
 z=np.load(src/f'{name}_communication.npz');samples=z['samples'].tolist();res=pd.read_csv(src/f'{name}_mapped_resource.tsv',sep='\t');features=pd.read_csv(src/f'{name}_features.tsv',sep='\t');labels=features.sender.unique().tolist();k=len(res)
 simple=~(res.ligand_subunits.str.contains('|',regex=False)|res.receptor_subunits.str.contains('|',regex=False))
 for row in pp.itertuples():
  si=samples.index(row.sample);ai=labels.index(row.sender);bi=labels.index(row.receiver);sl=slice((ai*len(labels)+bi)*k,(ai*len(labels)+bi+1)*k)
  act=z['support'][si,sl]&simple.to_numpy();rr=res[act]
  if row.role_overlap:
   h=exact_hitting_set([set([r.ligand_subunits,r.receptor_subunits]) for r in rr.itertuples()])
   val=h['value'];status=h['status'];targets=h['targets']
  else:
   val=row.ic100;status='EXACT_BY_DISJOINT_BIPARTITE_MVC';targets=[t.split(':',1)[1] for t in str(row.targets).split(';') if ':' in t]
  rows.append(dict(database=name,sample=row.sample,patient=row.patient,sender=row.sender,receiver=row.receiver,count=row.count,strength=row.strength,role_ic=row.ic100,molecular_ic=val,role_overlap=row.role_overlap,status=status,targets=';'.join(targets)))
q=pd.DataFrame(rows);q.to_csv(out/'gene_deduplicated_ic.tsv',sep='\t',index=False)
summ=[]
for name,g in q.groupby('database'):
 g=g[(g['count']>0)&g.molecular_ic.notna()]
 for feat in ['count','strength']:
  x=np.log1p(g[feat].to_numpy());y=g.molecular_ic.to_numpy();pred=np.zeros(len(g))
  for donor in g.patient.unique():
   test=(g.patient==donor).to_numpy();model=IsotonicRegression(out_of_bounds='clip').fit(x[~test],y[~test]);pred[test]=model.predict(x[test])
  summ.append(dict(database=name,predictor=feat,n_programs=len(g),molecular_IC_heldout_donor_R2=r2_score(y,pred),spearman=stats.spearmanr(g[feat],y).statistic,programs_with_role_overlap=int((g.role_overlap>0).sum()),programs_with_changed_IC=int((g.role_ic!=g.molecular_ic).sum())))
pd.DataFrame(summ).to_csv(out/'molecular_ic_independence.tsv',sep='\t',index=False)
# Exact synthetic certificates: same count and total strength, contrasting bottlenecks.
examples={}
for name,ls,rs in [('star',[f'L{i}' for i in range(12)],['R0']*12),('independent',[f'L{i}' for i in range(12)],[f'R{i}' for i in range(12)])]:
 v,targets=exact_cover(ls,rs);degrees={t:sum(t==f'L:{l}' or t==f'R:{r}' for l,r in zip(ls,rs)) for t in targets}
 examples[name]=dict(count=12,total_strength=12,IC100=v,certificate=targets,best_single_target_fraction_removed=max(degrees.values())/12)
examples['boundary']='Graph deletion identities only. No downstream response, dose, compensation or experimental blockade has been validated.'
(out/'synthetic_intervention.json').write_text(json.dumps(examples,indent=2))
# Stimulation sanity at the donor unit, separate from endogenous CCC blockade.
z=np.load(src/'sufficient_statistics.npz');groups=pd.read_csv(src/'groups.tsv',sep='\t');genes=z['genes'].tolist();markers=['ISG15','IFIT1','IFIT3','MX1','OAS1'];san=[]
for ct,g in groups.groupby('celltype'):
 for gene in markers:
  if gene not in genes:continue
  values=dict(zip(g['sample'],z['mean'][g.index,genes.index(gene)]));delta=[]
  for donor in sorted({s.split('__')[0] for s in values}):
   c=donor+'__ctrl';s=donor+'__stim'
   if c in values and s in values:
    valid=g[g['sample'].isin([c,s])].n_cells.min()>=5
    if valid:delta.append(values[s]-values[c])
  if len(delta)>=2:
   test=stats.ttest_1samp(delta,0);san.append(dict(celltype=ct,gene=gene,n_donors=len(delta),mean_stim_minus_ctrl=float(np.mean(delta)),donors_positive=int((np.array(delta)>0).sum()),pvalue=float(test.pvalue)))
s=pd.DataFrame(san);s['FDR']=multipletests(s.pvalue,method='fdr_bh')[1];s.to_csv(out/'IFN_response_sanity.tsv',sep='\t',index=False)
(out/'receipt.json').write_text(json.dumps({'status':'DONE','biology_blockade':'NOT_VALIDATED','same_gene_both_roles':'Molecular targets deduplicated by exact ILP when necessary; otherwise exact MVC certificate applies.'},indent=2))
