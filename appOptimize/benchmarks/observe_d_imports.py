"""Observe only this project's importer processes until its controller exits."""
import ctypes,json,os,pathlib,re,select,sys,time
import psutil
ROOT=pathlib.Path(__file__).resolve().parents[1]
pid=int(sys.argv[1]);libc=ctypes.CDLL(None,use_errno=True);fd=libc.syscall(434,pid,0)
if fd<0:raise OSError(ctypes.get_errno(),'pidfd_open')
started=time.time();last_memory={}
while True:
 for proc in psutil.process_iter(['pid','name']):
  if proc.info['name'] not in ['java','vcf2genomicsdb']:continue
  try:
   argv=proc.cmdline();command=' '.join(argv)
   if str(ROOT/'benchmarks/D_20260928') not in command:continue
   if proc.info['name']=='java' and 'GenomicsDBImport' not in argv:continue
   match=re.search(r'/N(50|200)_(native|gatk_native_reader|gatk_java)(?:/|\s|$)',command)
   if not match:continue
   stage='N'+match[1]+'_'+('native_import' if match[2]=='native' else match[2])
   cpu=proc.cpu_times();record={'epoch':time.time(),'elapsed_since_observer_s':time.time()-started,
    'processes':[{'pid':proc.pid,'name':proc.name(),'create_time':proc.create_time(),
     'cpu':cpu._asdict(),'rss':proc.memory_info().rss,'io':proc.io_counters()._asdict(),
     'fds':proc.num_fds(),'threads':proc.num_threads()}]}
   if time.time()-last_memory.get(proc.pid,0)>=300:
    record['memory_full_info']=proc.memory_full_info()._asdict();last_memory[proc.pid]=time.time()
   out=ROOT/'profiles/D_20260928'/stage;out.mkdir(exist_ok=True)
   with (out/'resource_samples.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
  except psutil.Error:continue
 if select.select([fd],[],[],30)[0]:break
os.close(fd)
(ROOT/'profiles/D_resource_observer_status.json').write_text(json.dumps({'controller_pid':pid,'started_epoch':started,'finished_epoch':time.time(),'status':'CONTROLLER_EXIT_OBSERVED'}))
