"""Run only the repaired C diagnostic and finalize on exit; keep old receipts immutable."""
import csv,json,os,pathlib,subprocess,sys,time,traceback
ROOT=pathlib.Path(__file__).resolve().parents[1];OUT=ROOT/'profiles/closure_20261002/C_recovery';started=time.time()
def state(value,**extra):
 p=OUT/'status.json.tmp';p.write_text(json.dumps(dict(state=value,pid=os.getpid(),started_epoch=started,updated_epoch=time.time(),**extra),indent=2));p.replace(OUT/'status.json')
def calls(path):
 with path.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 return sorted((r['Locus'],r['Chromosome'],r['Allele']) for r in rows)
state('RUNNING')
try:
 subprocess.run([sys.executable,str(ROOT/'benchmarks/perf_replay.py'),str(ROOT/'profiles/C_scaling_20261002/perf_retry_spec.json')],cwd=ROOT,check=True)
 base=ROOT/'benchmarks/C_20260928/typing/WGS_NA12878_t8/hla';new=ROOT/'profiles/C_scaling_20261002/perf_8_retry_20261002/typing/hla'
 compare={name:{'baseline_rows':len(calls(base/name)),'perf_rows':len(calls(new/name)),'allele_calls_exact':calls(base/name)==calls(new/name),'file_bytes_exact':(base/name).read_bytes()==(new/name).read_bytes()} for name in ['R1_bestguess.txt','R1_bestguess_G.txt']}
 (OUT/'call_comparison.json').write_text(json.dumps(compare,indent=2));state('COMPLETED_SCHEDULED_COMPUTATION')
except BaseException as e:
 traceback.print_exc();state('FAILED',error=repr(e));raise
finally:
 subprocess.run([sys.executable,str(ROOT/'benchmarks/finalize_closure.py')],cwd=ROOT,check=True)
