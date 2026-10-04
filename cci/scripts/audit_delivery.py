"""Completion evidence for the requested kill-test, separate from scientific PASS."""
from pathlib import Path
import csv,json,hashlib
import numpy as np,pandas as pd
C=json.load(open('results/selected_runs.json'));K=Path(C['kang']);S=Path(C['spatial'])
rows=[]
def add(item,verdict,evidence,note):rows.append(dict(requirement=item,work_status='COMPLETED',scientific_verdict=verdict,evidence=evidence,note=note))
for name in ['00_DECISION.md','01_collision_matrix.md','02_statistical_killtest.md','03_engineering_profile.md','04_intervention_complexity.md','05_spatial_homology.md','results/benchmark.tsv','README.md']:
 assert Path(name).is_file() and Path(name).stat().st_size>0
assert Path('collision_matrix.md').resolve()==Path('01_collision_matrix.md').resolve()
add('required_deliverables','DELIVERED','00–05 reports; results/benchmark.tsv; README.md; scripts/','Section12 numbering is canonical;section2 collision_matrix.md is an actual symlink to01_collision_matrix.md.')
text=Path('01_collision_matrix.md').read_text()
for name in ['CellPhoneDB','CellChat','FastCCC','LIANA','MultiNicheNet','Tensor-cell2cell','dominoSignal','scDiffCom','SpatialDM','COMMOT','scRICH','Renoir','MOSANIC']:assert name in text
add('collision_audit','PARTIAL_COLLISION;NO_FIRST_CLAIM','01_collision_matrix.md; sources/CellChat_spatial_tutorial.Rmd','All named methods reviewed through original papers/official sources; MOSANIC fulltext403 and U fields remain explicit. No overall GO or absence theorem.')
input=json.load(open('results/kang_input_audit.json'));assert input['shape']==[24673,15706] and input['donors']==8 and input['sample_units']==16 and input['max_row_sum_difference_nCount_RNA']==0
with np.load(K/'sufficient_statistics.npz') as z:
 assert set(['mean','variance','detection','raw_sum','library_size','n_cells','genes'])<=set(z.files)
 assert z['raw_sum'].shape==z['mean'].shape==z['variance'].shape==z['detection'].shape
add('Kang_and_sufficient_statistics','VERIFIED','results/kang_input_audit.json; '+str(K/'sufficient_statistics.npz'),'Original cell labels retained;8 independent donors,16 samples.')
sg=np.load(K/'null_signs.npy');assert sg.shape==(256,8) and len(np.unique(sg,axis=0))==256
add('100_to_500_null_replicates','RUN_CONDITIONAL_NULL',str(K/'null_signs.npy')+'; results/pooled_full_run_01/null_replicates.tsv','256 whole-pair swaps;100 external assignments x199 actual inner cell permutations;not independent null cohorts.')
for method,path in C['native'].items():
 p=Path(path);r=json.load(open(p/'receipt.json'));assert r['status']=='DONE' and r['samples']==16
 z=pd.read_csv(p/'scores.tsv.gz',sep='\t');assert z['sample'].nunique()==16 and np.isfinite(z.score).all()
 add('native_'+method,'EXECUTED',str(p/'receipt.json')+'; '+str(p/'scores.tsv.gz'),'Native API; condition null evaluated on cached sample score, not native specificity p.')
seen=[]
for path in C['native_analyses']:
 p=Path(path);r=json.load(open(p/'receipt.json'));assert r['status']=='DONE'
 n=pd.read_csv(p/'null_replicates.tsv',sep='\t');assert n.groupby('method').size().eq(256).all();seen+=list(n.method.unique())
assert set(seen)==set(C['native'])
add('native_condition_FPR_FDR_QQ_power','EXECUTED','results/native_*analysis_01/; results/figures/null_QQ.png','Family-specific universes disclosed; score-level synthetic power only.')
for f in ['paired_fit.tsv','factorial_fit.tsv']:
 z=pd.read_csv(K/f,sep='\t');assert {'beta','SE','pvalue','FDR'}<=set(z.columns) and len(z)>0
add('paired_and_factorial_design','RUN;NO_EXCLUSIVE_ADVANTAGE',str(K/'designs.json'),'C~condition+patient; synthetic treatment*time+batch. Optional mixed model not run.')
for size in [2000,5000,10000,24673]:
 ix=np.load(Path(C['scaling'])/f'cells_{size}.npy');assert len(ix)==size and len(np.unique(ix))==size
n=pd.read_csv(Path(C['pooled_full'])/'null_replicates.tsv',sep='\t');assert n.groupby('method').size().eq(100).all();assert set(n.n_features)=={4951}
add('KT1_cell_number_scaling','POOLED_INFLATION;SAMPLE_CONSERVATIVE',C['scaling']+'; '+C['pooled_full'],'450-feature two-type common-scale universe is not substituted for full six-type4951-feature experiment.')
h=pd.read_csv(K/'complex_hitting_sets.tsv',sep='\t');assert len(h)==100 and h.status.eq('OPTIMAL').all() and h.gap.eq(0).all()
add('IC_exact_MVC_and_complex_hitting_set','EXACT_CONDITIONAL_ON_GRAPH',str(K/'programs.tsv')+'; '+str(K/'complex_hitting_sets.tsv'),'All graph cover certificates checked by execution;100 ILP optima;role/gene distinction audited.')
q=pd.read_csv(K/'quadrants.tsv',sep='\t');assert q.groupby('database').quadrant.nunique().eq(4).all()
ic=pd.read_csv(K/'ic_independence.tsv',sep='\t');v=ic[(ic.database=='LIANA_CellPhoneDB')&(ic.predictor=='count')].iloc[0];assert v.heldout_donor_isotonic_R2>.9 and v.spearman>.95
add('KT4_IC_independence','FAIL_IN_ALTERNATIVE_RESOURCE',str(K/'ic_independence.tsv')+'; '+str(K/'similar_value_cases.tsv')+'; '+str(K/'quadrants.tsv'),'Four quadrants and matched-value examples exist, but do not erase count-rescale failure.')
d=json.load(open(K/'database_sensitivity.json'));assert d['n_matched_programs']==747
add('KT5_database_sensitivity','NOT_RANDOM;NOT_FULLY_STABLE',str(K/'database_matched_programs.tsv'),'Same measured genes;different native catalogue coverage;8 donors rather than747 independent replicates.')
s=json.load(open(S/'input_audit.json'));assert s['spots']==2695 and s['annotated']==1072 and s['unannotated']==1623 and s['annotation_barcode_not_in_counts']==1
for c in json.load(open(S/'numerical_equivalence.json')):assert c['all_LR_and_all_group_pairs'] and c['max_absolute_difference']<1e-10
sp=json.load(open(S/'typed_summary.json'));assert sp['removed_by_geometry']==9188 and sp['added_by_geometry']==0
add('KT6_spatial_homology','ENGINEERING_PASS;BIOLOGICAL_REPLICATION_UNTESTABLE',str(S/'numerical_equivalence.json')+'; '+str(S/'typed_summary.json'),'Full2695 spots;2 kernels;one section;100um pitch calibration;no invented labels or donor replication.')
assert 'cdist' not in Path('scripts/cci_core.py').read_text()
add('sparse_formal_operator','VERIFIED_SOURCE_AND_NUMERIC', 'scripts/cci_core.py; '+str(S/'profile.tsv'),'Dense N² only comparator script;formal KD-tree/CSR;O(K E) aggregation plus graph build.')
r=json.load(open(Path(C['spatial_native'])/'receipt.json'));assert r['status']=='DONE' and r['spots']==2695 and r['nboot']==20
add('spatial_native_reference','PROFILING_ONLY',C['spatial_native'],'CellChat2.2.0 spatial mode;20bootstrap;not v3 experiment or replicate-level calibration.')
for f in ['results/prototype_current_profile.tsv','results/ic_phase_profile.tsv',str(S/'profile.tsv')]:assert len(pd.read_csv(f,sep='\t'))>0
add('KT3_stage_runtime_RAM','PARTIAL_ENGINEERING_SUPPORT','03_engineering_profile.md; results/benchmark.tsv','Stage-matched spatial gain measured;no claim of same native statistical outputs with universally lower end-to-end RAM.')
add('intervention_validity','PERMITTED_SYNTHETIC_FALLBACK;NOT_VALIDATED',C['ic_sanity']+';04_intervention_complexity.md','Public perturbation/stimulation options reviewed;no matched real blockade test. Synthetic proof andKangsanity only, as control section6 allows.')
for row in csv.DictReader(open('infra/bioinf-data-index/BIOINF_DATA_INDEX.tsv'),delimiter='\t'):
 p=Path(row['path']);assert p.stat().st_size==int(row['bytes']) and hashlib.file_digest(p.open('rb'),'sha256').hexdigest()==row['sha256']
add('external_data_index','VERIFIED','infra/bioinf-data-index/BIOINF_DATA_INDEX.tsv','Actual bytes/SHA256 rechecked;derived files not counted as independent replication.')
add('optional_branches','NOT_RUN_AS_OPTIONAL','README.md;TODO.md','Skin,IC90,mixed model,50k/100k duplication andgreedy approximation optional;none substituted for required work.')
add('final_decision','KILL','00_DECISION.md;04_intervention_complexity.md','Kill-test completion is not all hypotheses passing. No formal package/paper/training follows.')
pd.DataFrame(rows).to_csv('results/completion_audit.tsv',sep='\t',index=False)
print('Audited',len(rows),'requirements;scientific failures and limitations retained.')
