"""Finish authorized A/C measurements in new namespaces; record durable exit receipts."""
import json,os,pathlib,subprocess,sys,time,traceback
ROOT=pathlib.Path(__file__).resolve().parents[1];route=sys.argv[1];out=ROOT/'profiles/closure_20261002'/route
out.mkdir(parents=True,exist_ok=True);started=time.time();env=dict(os.environ)
def receipt(state,**extra):
 p=out/'status.json.tmp';p.write_text(json.dumps(dict(state=state,route=route,started_epoch=started,updated_epoch=time.time(),**extra),indent=2));p.replace(out/'status.json')
def run(argv):
 print('RUN',json.dumps(argv),flush=True);subprocess.run(argv,cwd=ROOT,env=env,check=True)
receipt('RUNNING',pid=os.getpid())
try:
 if route=='A':
  base=ROOT/'benchmarks/A_W2_retained_20261002';base.mkdir(exist_ok=False)
  env.update(PATH=str(ROOT/'.envs/gtdbtk/bin')+':/usr/bin:/bin',GTDBTK_DATA_PATH=str(ROOT/'data/A_pplacer/reference'),APPOPT_A_PROFILE='profiles/A/W2_retained_20261002',APPOPT_A_SCALING='profiles/A_scaling_20261002')
  run([sys.executable,str(ROOT/'benchmarks/run_profile.py'),'--out',str(ROOT/env['APPOPT_A_PROFILE']),'--cwd',str(base),'--trace','--','gtdbtk','classify_wf','--genome_dir',str(ROOT/'data/A_pplacer/W2'),'--out_dir',str(base/'classify'),'--extension','fna','--cpus','16','--pplacer_cpus','8','--debug','--keep_intermediates'])
  run([sys.executable,str(ROOT/'benchmarks/post_a.py'),str(os.getpid())])
 elif route=='C':
  env['APPOPT_C_SCALING']='profiles/C_scaling_20261002'
  run([sys.executable,str(ROOT/'benchmarks/post_c.py'),str(os.getpid())])
 else:raise ValueError(route)
 receipt('COMPLETED_SCHEDULED_COMPUTATION')
except BaseException as e:
 traceback.print_exc();receipt('FAILED',error=repr(e));raise
