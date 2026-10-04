import ctypes,datetime,hashlib,json,os,pathlib,select,subprocess,tarfile
import psutil
ROOT=pathlib.Path(__file__).resolve().parents[1];DATA=ROOT/'data/A_pplacer';PROFILE=ROOT/'profiles/A';PROFILE.mkdir(exist_ok=True)
archive=DATA/'gtdbtk_r232_data.tar.gz'; marker=archive.with_name(archive.name+'.aria2')
pids=[]
for p in psutil.process_iter(['pid','name','cmdline']):
 try:
  if p.info['name']=='aria2c' and '--out=gtdbtk_r232_data.tar.gz' in (p.info['cmdline'] or []):pids.append(p.pid)
 except psutil.Error:pass
if len(pids)>1:raise RuntimeError('Multiple downloads; refuse ambiguous dependency')
if pids:
 print('WAIT_REFERENCE_PROCESS_EXIT',pids[0],flush=True)
 libc=ctypes.CDLL(None,use_errno=True)
 fd=libc.syscall(434,pids[0],0)
 if fd<0:raise OSError(ctypes.get_errno(),'pidfd_open')
 select.select([fd],[],[]);os.close(fd)
if marker.exists() or not archive.exists() or archive.stat().st_size!=60806405195:raise RuntimeError('Reference download incomplete')
receipt=DATA/'r232.receipt.json'
if not receipt.exists():
 md5=hashlib.md5();sha=hashlib.sha256()
 with archive.open('rb') as f:
  for block in iter(lambda:f.read(16*1024*1024),b''):md5.update(block);sha.update(block)
 if md5.hexdigest()!='25a59e0352b1fd150c589f56559767d4':raise RuntimeError('Official reference MD5 mismatch')
 receipt.write_text(json.dumps({'name':'GTDBTk_R232','path':str(archive.relative_to(ROOT)),'url':'https://data.ace.uq.edu.au/public/gtdb/data/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz','bytes':archive.stat().st_size,'md5':md5.hexdigest(),'sha256':sha.hexdigest(),'downloaded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'GTDB-Tk R232 reference package'},indent=2))
ref=DATA/'reference';ref.mkdir(exist_ok=True)
if not (ref/'.extracted').exists():
 with tarfile.open(archive) as t:first=t.next().name
 if first.split('/')[0] not in ['release232','r232']:raise RuntimeError('Unexpected archive root: '+first)
 with (PROFILE/'extract.log').open('w') as log:subprocess.run(['tar','-xzf',str(archive),'-C',str(ref),'--strip-components=1'],stdout=log,stderr=subprocess.STDOUT,check=True)
 (ref/'.extracted').write_text('official archive MD5 verified\n')
env=dict(os.environ,PATH=str(ROOT/'.envs/gtdbtk/bin')+':/usr/bin:/bin',GTDBTK_DATA_PATH=str(ref))
for label in ['W1','W2']:
 out=ROOT/'benchmarks'/f'A_{label}_8';out.mkdir(exist_ok=True)
 profile=PROFILE/f'{label}_8'
 if (profile/'status.json').exists():
  if json.loads((profile/'status.json').read_text()).get('exit_code')==0:continue
  raise RuntimeError('Existing failed run; use fresh output namespace')
 if (profile/'invocation.json').exists():raise RuntimeError('Interrupted output exists; use fresh output namespace')
 cmd=['/opt/anaconda3/bin/python3',str(ROOT/'benchmarks/run_profile.py'),'--out',str(profile),'--cwd',str(out),'--trace','--','gtdbtk','classify_wf','--genome_dir',str(DATA/label),'--out_dir',str(out/'classify'),'--extension','fna','--cpus','16','--pplacer_cpus','8','--debug']
 subprocess.run(cmd,env=env,check=True)
 print(label,'CLASSIFY_COMPLETE',flush=True)
print('TWO_WORKLOAD_BASELINES_COMPLETE; scaling/internal profile/equivalence assessment remain',flush=True)
