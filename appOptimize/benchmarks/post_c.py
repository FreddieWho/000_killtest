"""After complete BAM baselines, replay their actual extracted reads and compare calls."""
import csv,json,os,pathlib,re,subprocess,sys
from collections import Counter
from wait_receipt import wait_receipt
ROOT=pathlib.Path(__file__).resolve().parents[1];OUT=ROOT/os.environ.get('APPOPT_C_SCALING','profiles/C_scaling');OUT.mkdir(exist_ok=True)
wait_receipt(ROOT/'profiles/C/WES_typing_8/status.json',int(sys.argv[1]))
cases=[]
for assay in ['WGS','WES']:
 baseline=ROOT/'profiles/C'/f'{assay}_typing_8'
 subprocess.run([sys.executable,str(ROOT/'benchmarks/analyze_process_trace.py'),str(baseline/'process.trace'),'--out',str(baseline/'analysis')],check=True)
 rows=json.loads((baseline/'analysis/exec_intervals.json').read_text())['completed']
 matches=[r for r in rows if pathlib.Path(r['executable']).name=='HLA-LA' and '--action' in r['argv'] and r['argv'][r['argv'].index('--action')+1]=='HLA']
 if len(matches)!=1:raise RuntimeError('Expected one actual '+assay+' HLA backend invocation, found '+str(len(matches)))
 cases.append({'assay':assay,**matches[0]})
chosen=max(cases,key=lambda r:r['wall_s']);assay=chosen['assay']
(OUT/'case_selection.json').write_text(json.dumps({'selected_assay':assay,'cases':cases,
 'reason':'Profile the slower of the two actual backends; report both full baselines separately, not as a population-average claim.'},indent=2))
original=chosen['argv'];original[0]=str(ROOT/'.envs/hla/opt/hla-la/bin/HLA-LA')
def command(threads,out):
 argv=list(original);argv[argv.index('--maxThreads')+1]=str(threads)
 argv[argv.index('--outputDirectory')+1]=str(out)
 argv[argv.index('--sampleID')+1]=assay+'_kernel_'+str(threads)
 return argv
env=dict(os.environ,PATH=str(ROOT/'.envs/hla/bin')+':/usr/bin:/bin')
for threads in [1,8]:
 out=OUT/f'threads_{threads}'
 if out.exists():raise RuntimeError('Refuse existing replay directory '+str(out))
 out.mkdir();argv=command(threads,out/'typing')
 subprocess.run([sys.executable,str(ROOT/'benchmarks/run_profile.py'),'--out',str(out),'--cwd',str(ROOT/'.envs/hla/opt/hla-la/src'),'--',*argv],env=env,check=True)
with (ROOT/'data/C_hla/truth_20181129.tsv').open() as f:
 truth=next(r for r in csv.DictReader(f,delimiter='\t') if r['Sample ID']=='NA12878')
def calls(path):
 with path.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 by_locus={}
 for r in rows:by_locus.setdefault(r['Locus'].removeprefix('HLA-'),[]).append(r['Allele'])
 return {locus:sorted(alleles) for locus,alleles in by_locus.items()}
paths={'WGS_baseline':ROOT/'benchmarks/C_20260928/typing/WGS_NA12878_t8/hla/R1_bestguess.txt',
       'WES_baseline':ROOT/'benchmarks/C_20260928/typing/WES_NA12878_t8/hla/R1_bestguess.txt',
       assay+'_kernel_1':OUT/'threads_1/typing/hla/R1_bestguess.txt',
       assay+'_kernel_8':OUT/'threads_8/typing/hla/R1_bestguess.txt'}
all_calls={name:calls(path) for name,path in paths.items()}
all_G_calls={name:calls(path.with_name('R1_bestguess_G.txt')) for name,path in paths.items()}
details=[]
for name,typed in all_calls.items():
 for locus in ['A','B','C','DQB1','DRB1']:
  expected=sorted([truth[f'HLA-{locus} 1'],truth[f'HLA-{locus} 2']])
  alternatives=[sorted(set(re.findall(r'(?<![\d:])(\d+:\d+)(?::\d+)*(?:[A-Z])?',x))) for x in typed.get(locus,[])]
  unambiguous=len(alternatives)==2 and all(len(x)==1 for x in alternatives)
  observed=sorted(x[0] for x in alternatives) if unambiguous else None
  compatible=(len(alternatives)==2 and ((expected[0] in alternatives[0] and expected[1] in alternatives[1]) or
                                      (expected[1] in alternatives[0] and expected[0] in alternatives[1])))
  confirmed=[x[0] for x in alternatives if len(x)==1] if len(alternatives)<=2 else []
  confirmed_count=sum((Counter(expected)&Counter(confirmed)).values())
  padded=alternatives+[[]]*(2-len(alternatives))
  compatible_count=max(sum(e[i] in padded[i] for i in [0,1]) for e in [expected,expected[::-1]]) if len(padded)==2 else 0
  details.append({'run':name,'locus':locus,'reported':typed.get(locus,[]),'two_field_alternatives':alternatives,
   'truth':expected,'unambiguous_two_field_match':observed==expected,'compatible_with_truth':compatible,
   'confirmed_two_field_alleles':confirmed_count,'compatible_two_field_alleles':compatible_count,'allele_denominator':2,
   'expected_zygosity':'homozygous' if expected[0]==expected[1] else 'heterozygous',
   'observed_zygosity':('homozygous' if observed[0]==observed[1] else 'heterozygous') if observed else 'AMBIGUOUS_OR_MISSING'})
summaries={name:{'locus_denominator':5,'allele_denominator':10,
 'unambiguous_matching_genotypes':sum(d['unambiguous_two_field_match'] for d in details if d['run']==name),
 'confirmed_two_field_alleles':sum(d['confirmed_two_field_alleles'] for d in details if d['run']==name),
 'compatible_two_field_alleles':sum(d['compatible_two_field_alleles'] for d in details if d['run']==name)} for name in all_calls}
report={'selected_assay':assay,'calls':all_calls,'G_calls':all_G_calls,'truth_comparison':details,'truth_summary':summaries,
 'baseline_kernel1_exact_allele_sets':all_calls[assay+'_baseline']==all_calls[assay+'_kernel_1'],
 'baseline_kernel8_exact_allele_sets':all_calls[assay+'_baseline']==all_calls[assay+'_kernel_8'],
 'baseline_kernel1_exact_G_calls':all_G_calls[assay+'_baseline']==all_G_calls[assay+'_kernel_1'],
 'baseline_kernel8_exact_G_calls':all_G_calls[assay+'_baseline']==all_G_calls[assay+'_kernel_8'],
 'four_field_truth':'NOT_AVAILABLE; ambiguity sets are retained, no unique four-field validation',
 'truth_scoring':'Two-field alternatives from raw allele sets; never truncate a G-group representative as a unique genotype',
 'scope':'Unchanged backend on actual extracted reads; kernel timing excludes full BAM prefilter, whose baselines are separate'}
(OUT/'call_comparison.json').write_text(json.dumps(report,indent=2))
perfout=OUT/'perf_8';spec={'out':str(perfout.relative_to(ROOT)),'argv':command(8,perfout/'typing'),
 'cwd':str(ROOT/'.envs/hla/opt/hla-la/src'),'cpus':8,'memory':'256g','path_prefix':str(ROOT/'.envs/hla/bin'),'scope':report['scope']}
specfile=OUT/'perf_spec.json';specfile.write_text(json.dumps(spec,indent=2))
subprocess.run([sys.executable,str(ROOT/'benchmarks/perf_replay.py'),str(specfile)],check=True)
print('C actual-read kernel scaling, G-call comparison, two-field truth assessment and function profile completed.',flush=True)
