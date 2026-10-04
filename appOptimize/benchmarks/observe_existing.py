import ctypes,json,os,pathlib,select,sys,time
import psutil
pid=int(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
libc=ctypes.CDLL(None,use_errno=True);fd=libc.syscall(434,pid,0)
if fd<0:raise OSError(ctypes.get_errno(),'pidfd_open')
start=time.time()
with (out/'samples.jsonl').open('w') as f:
 while True:
  rec={'epoch':time.time(),'elapsed_since_observer_s':time.time()-start,'processes':[]}
  try:procs=[psutil.Process(pid),*psutil.Process(pid).children(recursive=True)]
  except psutil.Error:procs=[]
  for p in procs:
   try:rec['processes'].append({'pid':p.pid,'name':p.name(),'cpu':p.cpu_times()._asdict(),'rss':p.memory_info().rss,'io':p.io_counters()._asdict(),'fds':p.num_fds(),'threads':p.num_threads()})
   except psutil.Error:pass
  f.write(json.dumps(rec)+'\n');f.flush()
  ready,_,_=select.select([fd],[],[],30)
  if ready:break
os.close(fd)
(out/'status.json').write_text(json.dumps({'status':'OBSERVED_PROCESS_EXIT','observed_interval_s':time.time()-start,'pid':pid},indent=2))
print('Observed process exit. Samples cover only the recorded observation interval.',flush=True)
