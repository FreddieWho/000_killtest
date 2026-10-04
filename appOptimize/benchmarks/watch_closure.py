"""Refresh reports on terminal controller events only; no periodic polling."""
import ctypes,json,os,pathlib,select,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];libc=ctypes.CDLL(None,use_errno=True)
fd=libc.inotify_init1(os.O_CLOEXEC)
if fd<0:raise OSError(ctypes.get_errno(),'inotify_init1')
routes=['A','C'];pidfds={};last=None
for r in routes:
 p=ROOT/'profiles/closure_20261002'/r
 if libc.inotify_add_watch(fd,os.fsencode(p),0x8|0x80)<0:raise OSError(ctypes.get_errno(),'inotify_add_watch')
 d=json.loads((p/'status.json').read_text())
 if d.get('state')=='RUNNING':
  pidfd=libc.syscall(434,d['pid'],0)
  if pidfd>=0:pidfds[r]=pidfd
try:
 while True:
  states={r:json.loads((ROOT/'profiles/closure_20261002'/r/'status.json').read_text())['state'] for r in routes}
  terminal={r:s for r,s in states.items() if s in ['FAILED','COMPLETED_SCHEDULED_COMPUTATION']}
  if terminal!=last:
   subprocess.run([sys.executable,str(ROOT/'benchmarks/finalize_closure.py')],cwd=ROOT,check=True);last=terminal
  if len(terminal)==len(routes):break
  pending=[pidfds[r] for r in routes if r not in terminal and r in pidfds]
  ready=select.select([fd,*pending],[],[])[0]
  if fd in ready:os.read(fd,65536)
  for r,pidfd in pidfds.items():
   if pidfd in ready:
    d=json.loads((ROOT/'profiles/closure_20261002'/r/'status.json').read_text())
    if d['state']=='RUNNING':raise RuntimeError('Controller exited without terminal receipt: '+r)
 print('Terminal events collected; scientific/report gaps remain explicit.',flush=True)
finally:
 os.close(fd)
 for pidfd in pidfds.values():os.close(pidfd)
