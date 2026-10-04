"""Profile an authorized workload; 30s process samples, 300s temp metadata samples."""
import argparse,json,os,pathlib,subprocess,time
import psutil
p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--cwd',required=True);p.add_argument('--trace',action='store_true');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
out=pathlib.Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
cmd=a.command[1:] if a.command[:1]==['--'] else a.command
cmd=['/usr/bin/time','-v','-o',str(out/'time.txt'),*cmd]
if a.trace:cmd=['strace','--seccomp-bpf','-f','-ttt','-T','-s','4096','-e','trace=process','-o',str(out/'process.trace'),*cmd]
start=time.time();(out/'invocation.json').write_text(json.dumps({'command':cmd,'cwd':a.cwd,'started_epoch':start,'sample_seconds':30,'storage_sample_seconds':300},indent=2))
with (out/'stdout.log').open('w') as log,(out/'samples.jsonl').open('w') as samples:
 proc=subprocess.Popen(cmd,cwd=a.cwd,stdout=log,stderr=subprocess.STDOUT)
 last_storage=0;last_memory={}
 while proc.poll() is None:
  record={'epoch':time.time(),'elapsed_s':time.time()-start,'processes':[]}
  try:children=[psutil.Process(proc.pid),*psutil.Process(proc.pid).children(recursive=True)]
  except psutil.Error:children=[]
  for child in children:
   try:
    with child.oneshot():
     entry={'pid':child.pid,'name':child.name(),'cpu':child.cpu_times()._asdict(),'rss':child.memory_info().rss,'io':child.io_counters()._asdict(),'status':child.status(),'fds':child.num_fds(),'threads':child.num_threads()}
     if time.time()-last_memory.get(child.pid,0)>=300:
      try:entry['memory_full_info']=child.memory_full_info()._asdict();last_memory[child.pid]=time.time()
      except psutil.Error:pass
     record['processes'].append(entry)
   except psutil.Error:pass
  if time.time()-last_storage>=300:
   count=0;size=0
   for root,_,files in os.walk(a.cwd):
    for name in files:
     try:size+=os.stat(os.path.join(root,name)).st_size;count+=1
     except FileNotFoundError:pass
   record.update(temp_files=count,temp_bytes=size);last_storage=time.time()
  samples.write(json.dumps(record)+'\n');samples.flush()
  try:proc.wait(timeout=30)
  except subprocess.TimeoutExpired:pass
status={'exit_code':proc.returncode,'wall_s':time.time()-start,'finished_epoch':time.time(),'status':'COMPLETE' if proc.returncode==0 else 'FAILED'}
tmp=out/'status.json.tmp';tmp.write_text(json.dumps(status,indent=2));tmp.replace(out/'status.json');print(json.dumps(status),flush=True)
raise SystemExit(proc.returncode)
