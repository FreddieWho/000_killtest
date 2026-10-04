import csv,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
def timing(path):
 if not path.exists():return {}
 d={}
 for line in path.read_text(errors='replace').splitlines():
  if ': ' in line:
   k,v=line.strip().split(': ',1);d[k]=v
 def number(k):
  try:return float(d[k])
  except (KeyError,ValueError):return None
 wall=next((v for k,v in d.items() if k.startswith('Elapsed (wall clock)')),None)
 if wall:
  parts=wall.split(':');wall=sum(float(x)*60**i for i,x in enumerate(reversed(parts)))
 return {'wall_s':wall,'user_s':number('User time (seconds)'),'system_s':number('System time (seconds)'),'peak_rss_kib':number('Maximum resident set size (kbytes)'),'major_faults':number('Major (requiring I/O) page faults'),'minor_faults':number('Minor (reclaiming a frame) page faults'),'fs_input_bytes':None if number('File system inputs') is None else number('File system inputs')*512,'fs_output_bytes':None if number('File system outputs') is None else number('File system outputs')*512,'exit_code':number('Exit status')}
rows=[]
def add(route,workload,path,n,unit):
 p=ROOT/path;r={'route':route,'workload':workload,'n':n,'unit':unit,'run_state':'NOT_RUN','profile_path':str(p.relative_to(ROOT)),'storage':'local_rotational_sda','measurement_scope':'instrumented_baseline','equivalence':'NOT_RUN'}
 r.update(timing(p/'time.txt'))
 if (p/'invocation.json').exists():r['run_state']='RUNNING'
 if (p/'status.json').exists():r['run_state']=json.loads((p/'status.json').read_text())['status']
 if r.get('exit_code') is not None:r['run_state']='EXIT_0' if r['exit_code']==0 else 'FAILED'
 if (p/'tracer_detach.json').exists():
  r['measurement_scope']='diagnostic_strace_overhead_not_clean_baseline'
  if r.get('exit_code') is None:r['run_state']='RUNNING_DIAGNOSTIC'
 if (p/'bounded_stop.json').exists():
  r['run_state']='BOUNDED_3_ROUNDS'
  r['measurement_scope']='rounds1_3_R1_trace_contaminated_R4_partial_excluded_see_round_metrics'
  for key in ['wall_s','user_s','system_s','peak_rss_kib','fs_input_bytes','fs_output_bytes','major_faults','minor_faults']:r[key]=None
 rows.append(r)
for w,n in [('W1',1),('W2',100)]:add('A',w+'_pplacer8','profiles/A/'+w+'_8',n,'input_genomes_not_placement_denominator')
add('B','rice_16','profiles/B_rice',374422835,'genome_bp')
r={'route':'B','workload':'dmel_16_interrupted','n':143726002,'unit':'genome_bp','run_state':'INTERRUPTED_ROUND4_INCOMPLETE','profile_path':'profiles/B_dmel.time.txt','storage':'local_rotational_sda','measurement_scope':'interrupted_original_not_full_workflow','equivalence':'NOT_RUN'};r.update(timing(ROOT/'profiles/B_dmel.time.txt'))
if r.get('exit_code') is not None:r['run_state']='EXIT_0' if r['exit_code']==0 else 'FAILED'
rows.append(r)
add('B','dmel_16_fresh_full','profiles/B_dmel_full_20260928',143726002,'genome_bp')
for kind in ['WGS','WES']:add('C',kind+'_NA12878_t8','profiles/C/'+kind+'_typing_8',1,'sample_assay_same_donor')
rows.append({'route':'C','workload':'RNA_HLA_HD','n':1,'unit':'sample_assay_same_donor','run_state':'BLOCKED_SOFTWARE','equivalence':'NOT_RUN'})
for n in [50,200]:
 for mode in ['native_import','gatk_java','gatk_native_reader']:add('D',f'N{n}_{mode}',f'profiles/D_20260928/N{n}_{mode}',n,'distinct_samples_chr22')
 # Native import above excludes separately recorded initialization; do not silently compare it as total ETL.
 rows[-3]['measurement_scope']='native_import_only_add_init_for_ETL'
fields=['route','workload','n','unit','run_state','wall_s','user_s','system_s','peak_rss_kib','fs_input_bytes','fs_output_bytes','major_faults','minor_faults','exit_code','storage','measurement_scope','equivalence','profile_path']
with (ROOT/'benchmarks/summary.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for r in rows:w.writerow({k:'NA' if r.get(k) is None else r.get(k,'NA') for k in fields})
print('summary rows',len(rows))
