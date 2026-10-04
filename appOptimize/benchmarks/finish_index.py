"""Refresh completed external assets when existing acquisition-dependent controllers exit."""
import ctypes,json,os,pathlib,select,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1];libc=ctypes.CDLL(None,use_errno=True)
waiting={};events=[]
for arg in sys.argv[1:]:
 pid=int(arg);fd=libc.syscall(434,pid,0)
 if fd>=0:waiting[fd]=pid
while waiting:
 for fd in select.select(list(waiting),[],[])[0]:
  pid=waiting.pop(fd);os.close(fd)
  result=subprocess.run([sys.executable,str(ROOT/'benchmarks/update_data_index.py'),'--sync-central'])
  events.append({'controller_pid':pid,'observed_exit_epoch':time.time(),'index_update_exit_code':result.returncode})
  (ROOT/'profiles/index_completion_events.json').write_text(json.dumps(events,indent=2))
  if result.returncode:raise RuntimeError('Index update failed')
