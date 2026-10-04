"""One route per user systemd unit; preserve old runs and write explicit final state."""
import datetime,json,os,pathlib,shutil,subprocess,sys,time,traceback
ROOT=pathlib.Path(__file__).resolve().parents[1];route=sys.argv[1]
OUT=ROOT/'profiles/persistent_20260928'/route;OUT.mkdir(parents=True,exist_ok=True)
PY='/opt/anaconda3/bin/python3';observers=[]
def state(value):
 tmp=OUT/'status.json.tmp';tmp.write_text(json.dumps(value,indent=2));tmp.replace(OUT/'status.json')
def run(argv,env=None):
 started=time.time()
 with (OUT/'steps.jsonl').open('a') as f:f.write(json.dumps({'event':'start','epoch':started,'argv':argv})+'\n')
 print('RUN',json.dumps(argv),flush=True)
 result=subprocess.run(argv,cwd=ROOT,env=env)
 with (OUT/'steps.jsonl').open('a') as f:f.write(json.dumps({'event':'exit','epoch':time.time(),'exit_code':result.returncode,'argv':argv})+'\n')
 if result.returncode:raise RuntimeError('Step failed with '+str(result.returncode)+': '+str(argv))
def script(name,*args):run([PY,'-u',str(ROOT/'benchmarks'/name),*map(str,args)])
started=time.time();state({'state':'RUNNING','route':route,'pid':os.getpid(),'started_epoch':started})
exit_code=0
try:
 if route=='A':
  data=ROOT/'data/A_pplacer'
  if not (data/'r232.receipt.json').exists():
   run(['/usr/bin/aria2c','--continue=true','--auto-file-renaming=false','--file-allocation=none',
        '--max-connection-per-server=8','--split=8','--summary-interval=300','--console-log-level=warn',
        '--dir='+str(data),'--out=gtdbtk_r232_data.tar.gz',
        'https://data.ace.uq.edu.au/public/gtdb/data/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz'])
  script('run_a.py');script('post_a.py',os.getpid())
 elif route=='B':
  base=ROOT/'benchmarks/B_dmel_full_20260928';profile=ROOT/'profiles/B_dmel_full_20260928'
  if base.exists() or profile.exists():raise RuntimeError('Fresh uninterrupted baseline namespace already exists')
  base.mkdir()
  for source in (ROOT/'benchmarks/B_dmel').glob('dmel.*'):shutil.copy2(source,base/source.name)
  env=dict(os.environ,PATH=str(ROOT/'.envs/repeatmodeler/bin')+':/usr/bin:/bin',BLAST_USAGE_REPORT='false')
  run([PY,str(ROOT/'benchmarks/run_profile.py'),'--out',str(profile),'--cwd',str(base),'--trace','--',
       str(ROOT/'.envs/repeatmodeler/bin/RepeatModeler'),'-database',str(base/'dmel'),'-threads','16','-srand','20260928'],env=env)
  script('summarize_repeatmodeler.py')
  script('analyze_process_trace.py',profile/'process.trace','--out',profile/'analysis')
  script('summarize_resource_samples.py',profile/'samples.jsonl','--out',profile/'resource_summary.json')
 elif route=='C':
  script('run_c.py');script('post_c.py',os.getpid())
 elif route=='D':
  script('fetch_d_cohort.py',200)
  observers.append(subprocess.Popen([PY,str(ROOT/'benchmarks/observe_d_imports.py'),str(os.getpid())]))
  script('run_d.py');script('compare_d_exports.py',os.getpid())
 else:raise ValueError('Unknown route')
 state({'state':'COMPLETED_SCHEDULED_COMPUTATION','route':route,'started_epoch':started,'finished_epoch':time.time(),
        'scope':'Scheduled measurements completed; route decision still requires evidence review.'})
except BaseException as e:
 exit_code=1;traceback.print_exc()
 state({'state':'FAILED','route':route,'started_epoch':started,'finished_epoch':time.time(),'error':repr(e)})
finally:
 for observer in observers:
  if observer.poll() is None:observer.terminate()
 for command in [[PY,str(ROOT/'benchmarks/summarize_profiles.py')],
                 [PY,str(ROOT/'benchmarks/update_data_index.py'),'--sync-central']]:
  result=subprocess.run(command,cwd=ROOT)
  if result.returncode:print('FINALIZATION_FAILED',command,result.returncode,flush=True)
raise SystemExit(exit_code)
