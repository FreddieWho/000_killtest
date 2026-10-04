"""Summarize sampled observations without claiming exact peaks or accounting for short-lived children."""
import argparse,json,pathlib
p=argparse.ArgumentParser();p.add_argument('samples',type=pathlib.Path);p.add_argument('--out',required=True,type=pathlib.Path);a=p.parse_args()
first={};last={};times=[];max_rss=0;max_fds=0;max_threads=0;max_temp_bytes=None;max_temp_files=None
with a.samples.open() as f:
 for line in f:
  try:r=json.loads(line)
  except json.JSONDecodeError:continue # live writer may not have finished its final line
  times.append(r['epoch']);rss=fds=threads=0
  for proc in r.get('processes',[]):
   pid=proc['pid'];cpu=proc['cpu']['user']+proc['cpu']['system']
   first.setdefault(pid,cpu);last[pid]=cpu
   rss+=proc['rss'];fds+=proc.get('fds',0);threads+=proc.get('threads',0)
  max_rss=max(max_rss,rss);max_fds=max(max_fds,fds);max_threads=max(max_threads,threads)
  if 'temp_bytes' in r:max_temp_bytes=max(max_temp_bytes or 0,r['temp_bytes'])
  if 'temp_files' in r:max_temp_files=max(max_temp_files or 0,r['temp_files'])
result={'source':str(a.samples),'sample_count':len(times),'first_epoch':min(times) if times else None,
 'last_epoch':max(times) if times else None,'observed_cpu_delta_s':sum(max(0,last[k]-first[k]) for k in last),
 'sampled_max_sum_rss_bytes':max_rss,'sampled_max_sum_fds':max_fds or None,'sampled_max_sum_threads':max_threads or None,
 'sampled_max_directory_bytes':max_temp_bytes,'sampled_max_directory_files':max_temp_files,
 'caveats':['CPU delta misses pre-observation work and unsampled short children.','RSS sums double-count shared pages; not PSS or exact tree peak.','FD sums count handles, not unique open files.','Directory size includes inputs and retained outputs; not pure temporary bytes.','All peaks are sampled lower bounds for their stated measurement.']}
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2));print(json.dumps(result))
