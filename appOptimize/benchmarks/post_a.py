"""After real W2: count actual bacterial placements, replay threads, compare jplace, sample functions."""
import csv,json,os,pathlib,subprocess,sys
from wait_receipt import wait_receipt
ROOT=pathlib.Path(__file__).resolve().parents[1];PROFILE=ROOT/os.environ.get('APPOPT_A_PROFILE','profiles/A/W2_8')
OUT=ROOT/os.environ.get('APPOPT_A_SCALING','profiles/A_scaling');OUT.mkdir(exist_ok=True)
wait_receipt(PROFILE/'status.json',int(sys.argv[1]))
subprocess.run([sys.executable,str(ROOT/'benchmarks/analyze_process_trace.py'),str(PROFILE/'process.trace'),'--out',str(PROFILE/'analysis')],check=True)
intervals=json.loads((PROFILE/'analysis/exec_intervals.json').read_text())['completed']
def placements(path):
 doc=json.loads(path.read_text());fields=doc['fields'];result={}
 for p in doc['placements']:
  names=p.get('n',[item[0] for item in p.get('nm',[])])
  values=[dict(zip(fields,row)) for row in p['p']]
  for name in names:result[name]=values
 return result
cases=[];all_names=set();missing=[]
for row in intervals:
 argv=row['argv']
 if pathlib.Path(row['executable']).name!='pplacer' or '-o' not in argv or 'bac120' not in ' '.join(argv):continue
 file=pathlib.Path(argv[argv.index('-o')+1])
 if not file.exists():missing.append(str(file));continue
 placed=placements(file);all_names.update(placed);cases.append((len(placed),row,file))
summary=ROOT/'benchmarks/A_W2_retained_20261002/classify/gtdbtk.bac120.summary.tsv'
old_summary=ROOT/'benchmarks/A_W2_8/classify/gtdbtk.bac120.summary.tsv'
if summary.exists():
 def taxonomy(path):
  with path.open() as f:return {r['user_genome']:r for r in csv.DictReader(f,delimiter='\t')}
 old_tax=taxonomy(old_summary);new_tax=taxonomy(summary)
 (OUT/'taxonomy_comparison.json').write_text(json.dumps({'genomes':len(new_tax),'same_genome_set':set(old_tax)==set(new_tax),'all_summary_fields_exact':old_tax==new_tax,'taxonomy_exact':{k:v['classification'] for k,v in old_tax.items()}=={k:v['classification'] for k,v in new_tax.items()},'scope':'Original classify versus same binary with keep_intermediates; not a modified placement algorithm'},indent=2))
qualification={'distinct_bacterial_placed_queries':len(all_names),'names':sorted(all_names),'pplacer_calls':len(cases),
 'missing_placement_files':missing,'task_requires_at_least_50':len(all_names)>=50,'denominator':'Union of bacterial jplace names, not downloaded genomes'}
(OUT/'W2_placement_denominator.json').write_text(json.dumps(qualification,indent=2))
if missing:raise RuntimeError('Placement files were removed; cannot infer zero placements: '+str(missing))
if len(all_names)<50:raise RuntimeError('W2 has fewer than 50 actual bacterial placements; add real MAGs before qualification')
_,case,baseline=max(cases,key=lambda c:c[1]['wall_s'])
(OUT/'case_selection.json').write_text(json.dumps({'baseline':str(baseline),'wall_s':case['wall_s'],'queries':len(placements(baseline)),'reason':'Replay slowest actual placement call; query count alone does not identify the critical path.','all_call_wall_s':sum(c[1]['wall_s'] for c in cases)},indent=2))
original=list(case['argv']);original[0]=str(ROOT/'.envs/gtdbtk/bin/pplacer')
env=dict(os.environ,PATH=str(ROOT/'.envs/gtdbtk/bin')+':/usr/bin:/bin')
def command(threads,output):
 argv=list(original);argv[argv.index('-j')+1]=str(threads);argv[argv.index('-o')+1]=str(output)
 if '--mmap-file' in argv:argv[argv.index('--mmap-file')+1]=str(output.parent/'scratch.mmap')
 return argv
for threads in [1,8]:
 path=OUT/f'threads_{threads}'
 if path.exists():raise RuntimeError('Refuse existing replay output '+str(path))
 path.mkdir();argv=command(threads,path/'result.jplace')
 subprocess.run([sys.executable,str(ROOT/'benchmarks/run_profile.py'),'--out',str(path),'--cwd',str(path),'--',*argv],env=env,check=True)
one=placements(OUT/'threads_1/result.jplace');eight=placements(OUT/'threads_8/result.jplace');base=placements(baseline)
comparison={'baseline':str(baseline),'queries':len(base),'query_sets_equal':set(base)==set(one)==set(eight),
 'baseline_vs_1_exact':base==one,'baseline_vs_8_exact':base==eight,'one_vs_eight_exact':one==eight,
 'prototype':'NONE; thread count replay of unchanged binary','taxonomy_equivalence':'NOT_RUN; baseline taxonomy exists, no modified full classify workflow'}
details=[]
for name in sorted(set(base)|set(one)|set(eight)):
 b=base.get(name,[]);x=one.get(name,[]);y=eight.get(name,[])
 edge=lambda v:[z.get('edge_num') for z in v]
 def maxdiff(left,right):
  if len(left)!=len(right):return None
  return max([abs(a[k]-b[k]) for a,b in zip(left,right) for k in a if k!='edge_num' and isinstance(a[k],(int,float)) and isinstance(b.get(k),(int,float))] or [0])
 weights=sorted([v.get('like_weight_ratio',0) for v in b],reverse=True)
 details.append({'query':name,'baseline_edges':edge(b),'one_edges':edge(x),'eight_edges':edge(y),
  'baseline_like_weight_ratios':[z.get('like_weight_ratio') for z in b],
  'one_like_weight_ratios':[z.get('like_weight_ratio') for z in x],
  'eight_like_weight_ratios':[z.get('like_weight_ratio') for z in y],
  'max_abs_numeric_diff_1':maxdiff(b,x),'max_abs_numeric_diff_8':maxdiff(b,y),
  'retained_candidate_count':len(b),'top_two_weight_gap':weights[0]-weights[1] if len(weights)>1 else None})
comparison['per_query']=details
(OUT/'jplace_comparison.json').write_text(json.dumps(comparison,indent=2))
perfout=OUT/'perf_8'
spec={'out':str(perfout.relative_to(ROOT)),'argv':command(8,perfout/'result.jplace'),
 'cpus':8,'memory':'512g','path_prefix':str(ROOT/'.envs/gtdbtk/bin'),
 'scope':'Slowest actual bacterial W2 pplacer call; diagnostic function profile, not end-to-end speed'}
specfile=OUT/'perf_spec.json';specfile.write_text(json.dumps(spec,indent=2))
subprocess.run([sys.executable,str(ROOT/'benchmarks/perf_replay.py'),str(specfile)],check=True)
print('A placement count, thread replay, jplace comparison and function diagnostic complete; full route decision still requires analysis.',flush=True)
