import ctypes,json,os,pathlib,select,subprocess,time,xml.etree.ElementTree as ET
ROOT=pathlib.Path(__file__).resolve().parents[1];DATA=ROOT/'data/D_genomicsdb';BASE=ROOT/'benchmarks/D_20260928';BASE.mkdir(exist_ok=True)
PROF=ROOT/'profiles/D_20260928';PROF.mkdir(exist_ok=True)
ns={'s':'http://s3.amazonaws.com/doc/2006-03-01/'}
ids=[x.text.rstrip('/').split('/')[-1] for x in ET.parse(ROOT/'sources/D_s3_listing.xml').findall('s:CommonPrefixes/s:Prefix',ns)][:200]
def wait_inputs(n):
 libc=ctypes.CDLL(None,use_errno=True);fd=libc.inotify_init1(os.O_CLOEXEC)
 if fd<0:raise OSError(ctypes.get_errno(),'inotify_init1')
 if libc.inotify_add_watch(fd,os.fsencode(DATA),0x8|0x80)<0:raise OSError(ctypes.get_errno(),'inotify_add_watch')
 print('WAIT_INPUTS',n,'filesystem completion events; no periodic polling',flush=True)
 try:
  while not all((DATA/(x+'.receipt.json')).exists() for x in ids[:n]):
   select.select([fd],[],[]);os.read(fd,65536)
 finally:os.close(fd)
 paths=[str(DATA/(x+'.chr22.g.vcf.gz')) for x in ids[:n]]
 (BASE/f'samples_{n}.txt').write_text('\n'.join(paths)+'\n')
 (BASE/f'samplemap_{n}.tsv').write_text(''.join(x+'\t'+p+'\n' for x,p in zip(ids[:n],paths)))
 return paths
def run(name,image,command,trace=False):
 out=PROF/name;out.mkdir(exist_ok=True);receipt=out/'status.json'
 if receipt.exists():
  old=json.loads(receipt.read_text())
  if old['exit_code']==0:return
  raise RuntimeError('Prior failed run needs a fresh output namespace: '+name)
 loader='/prof/ld-linux-x86-64.so.2'
 payload=['--library-path','/prof','/prof/time','-v','-o',str(out/'time.txt'),*command]
 args=['/opt/anaconda3/bin/python3',str(ROOT/'benchmarks/docker_exec.py'),'run','--rm','--user',f'{os.getuid()}:{os.getgid()}','--cpus','8','--memory','64g','-e','OMP_NUM_THREADS=8','-v',f'{ROOT}:{ROOT}','-v',f'{ROOT}/profiles/profiler_tools:/prof:ro','-w',str(ROOT),'--entrypoint',loader,image,*payload]
 (out/'invocation.json').write_text(json.dumps({'command':args,'storage':'local rotational sda; cache uncontrolled','image':image},indent=2))
 start=time.time()
 with (out/'stdout.log').open('w') as log:r=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT)
 tmp=out/'status.json.tmp'
 tmp.write_text(json.dumps({'exit_code':r.returncode,'launcher_wall_s':time.time()-start,'status':'COMPLETE' if r.returncode==0 else 'FAILED'},indent=2));tmp.replace(receipt)
 print(name,'exit',r.returncode,'launcher wall',time.time()-start,flush=True)
 if r.returncode:raise RuntimeError('Run failed: '+name)
(BASE/'chr22.intervals').write_text('chr22\n')
for n in [50,200]:
 wait_inputs(n)
 workspace=BASE/f'N{n}_native'
 run(f'N{n}_native_init','ghcr.io/genomicsdb/genomicsdb:v1.5.4',['vcf2genomicsdb_init','-w',str(workspace),'-s',str(BASE/f'samples_{n}.txt'),'-i',str(BASE/'chr22.intervals')])
 run(f'N{n}_native_import','ghcr.io/genomicsdb/genomicsdb:v1.5.4',['vcf2genomicsdb',str(workspace/'loader.json')])
 for mode in ['java','native_reader']:
  workspace=BASE/f'N{n}_gatk_{mode}'
  cmd=['/gatk/gatk','--java-options','-Xmx12g -XX:ActiveProcessorCount=8','GenomicsDBImport','--genomicsdb-workspace-path',str(workspace),'--sample-name-map',str(BASE/f'samplemap_{n}.tsv'),'-L','chr22','--batch-size','50','--reader-threads','1','--tmp-dir',str(BASE/f'tmp_{n}_{mode}')]
  pathlib.Path(cmd[-1]).mkdir(exist_ok=True)
  if mode=='native_reader':cmd+=['--bypass-feature-reader','true']
  run(f'N{n}_gatk_{mode}','broadinstitute/gatk:4.6.2.0',cmd)
print('PRIMARY_BASELINES_COMPLETE',flush=True)
