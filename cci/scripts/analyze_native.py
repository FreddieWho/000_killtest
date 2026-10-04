import sys,json,time,itertools
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from patsy import dmatrix
from cci_core import fit_linear
out=Path(sys.argv[1]);out.mkdir(exist_ok=False)
nulls=[];summaries=[];fits=[];powers=[]
for path in sys.argv[2:]:
 root=Path(path);receipt=json.load(open(root/'receipt.json'))
 if receipt.get('status')!='DONE' or receipt.get('samples')!=16:raise ValueError('Incomplete native method '+path)
 method=receipt['method'];z=pd.read_csv(root/'scores.tsv.gz',sep='\t')
 if z.groupby(['sample','feature']).score.nunique().max()>1:raise ValueError('Inconsistent duplicate feature')
 z=z.drop_duplicates(['sample','feature']);wide=z.pivot(index='sample',columns='feature',values='score').sort_index()
 good=np.isfinite(wide.to_numpy()).all(axis=0)&(wide.std().to_numpy()>1e-12);y=wide.iloc[:,good].to_numpy();cols=wide.columns[good]
 samples=wide.index.tolist();donors=sorted({s.split('__')[0] for s in samples});ci=[samples.index(d+'__ctrl') for d in donors];si=[samples.index(d+'__stim') for d in donors]
 dif=y[si]-y[ci];signs=np.array(list(itertools.product([-1.,1.],repeat=8)))
 meta=pd.DataFrame({'sample':samples,'patient':[s.split('__')[0] for s in samples],'condition':[s.split('__')[1] for s in samples]})
 x=dmatrix('~ condition + patient',meta,return_type='dataframe');t=time.perf_counter();b,se,p,df=fit_linear(y,x,x.columns.get_loc('condition[T.stim]'))
 fits.append(dict(method=method,n_features=y.shape[1],seconds=time.perf_counter()-t))
 pd.DataFrame({'feature':cols,'beta':b,'SE':se,'pvalue':p,'FDR':multipletests(p,method='fdr_bh')[1]}).to_csv(out/f'{method}_paired_fit.tsv',sep='\t',index=False)
 mu=signs@dif/8;v=np.maximum(((dif*dif).sum(axis=0)-8*mu*mu)/7,0)
 with np.errstate(divide='ignore',invalid='ignore'):ts=mu/np.sqrt(v/8)
 ps=np.nan_to_num(2*stats.t.sf(abs(ts),7),nan=1)
 for i,pv in enumerate(ps):
  reject=multipletests(pv,method='fdr_bh')[0]
  nulls.append(dict(method=method,replicate=i,n_features=len(pv),FPR=float((pv<.05).mean()),FDP=float(reject.any()),n_BH_discoveries=int(reject.sum())))
 q=np.linspace(.001,.999,999);pd.DataFrame({'expected':q,'observed':np.quantile(ps,q)}).to_csv(out/f'QQ_{method}.tsv',sep='\t',index=False)
 summaries.append(dict(method=method,all_features=wide.shape[1],complete_nonconstant_features=y.shape[1],fraction_p_equals_1=float((ps==1).mean()),uniform_KS_D=float(stats.kstest(ps.ravel(),'uniform').statistic)))
 # Additional score-level power, no claim of biological intervention ground truth.
 sd=np.std(dif,axis=0,ddof=1);eligible=np.flatnonzero(sd>1e-12);rng=np.random.default_rng(20260928)
 targets=rng.choice(eligible,max(1,len(eligible)//10),replace=False);truth=np.zeros(len(sd),bool);truth[targets]=True
 for effect in [.5,1.,2.]:
  for i,sg in enumerate(signs):
   d=dif*sg[:,None];d[:,targets]+=effect*sd[targets]
   pv=np.nan_to_num(stats.ttest_1samp(d,0,axis=0).pvalue,nan=1);r=multipletests(pv,method='fdr_bh')[0]
   powers.append(dict(method=method,effect_in_donor_SD=effect,replicate=i,power=float(r[truth].mean()),FDP=float(r[~truth].sum()/max(1,r.sum())),n_features=len(sd),n_injected=int(truth.sum())))
 print(method,'analyzed',y.shape,flush=True)
pd.DataFrame(nulls).to_csv(out/'null_replicates.tsv',sep='\t',index=False)
pd.DataFrame(summaries).to_csv(out/'uniformity.tsv',sep='\t',index=False)
pd.DataFrame(fits).to_csv(out/'profile.tsv',sep='\t',index=False)
pd.DataFrame(powers).to_csv(out/'synthetic_score_power.tsv',sep='\t',index=False)
(out/'receipt.json').write_text(json.dumps({'status':'DONE','methods':[s['method'] for s in summaries],'null_replicates_per_method':256,'unit':'8 donors /16 sample units','estimand':'native per-sample score condition differences, not native specificity p-values','power':'score-level standardized injections in nonzero-variance features only'},indent=2))
