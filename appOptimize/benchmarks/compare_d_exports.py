"""Same-reader, five-region import equivalence spot check; not full-cohort equivalence."""
import ctypes,hashlib,json,os,pathlib,select,subprocess,sys
import pysam
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/D_20260928';PROF=ROOT/'profiles/D_20260928'
OUT=ROOT/'profiles/D_export_comparison';OUT.mkdir(exist_ok=True)
libc=ctypes.CDLL(None,use_errno=True);controller=libc.syscall(434,int(sys.argv[1]),0)
if controller<0:raise OSError(ctypes.get_errno(),'pidfd_open')
intervals=[f'chr22:{start}-{start+1000}' for start in [20000000,25000000,30000000,40000000,50000000]]
reference=(ROOT/'data/D_genomicsdb/reference.fa').resolve()
def wait_receipt(path):
 while not path.exists():
  parent=path.parent
  while not parent.exists():parent=parent.parent
  fd=libc.inotify_init1(os.O_CLOEXEC)
  if fd<0:raise OSError(ctypes.get_errno(),'inotify')
  if libc.inotify_add_watch(fd,os.fsencode(parent),0x8|0x80|0x100)<0:raise OSError(ctypes.get_errno(),'watch')
  try:
   if path.exists():break
   ready=select.select([fd,controller],[],[])[0]
   if controller in ready and not path.exists():raise RuntimeError('Controller exited before '+str(path))
   if fd in ready:os.read(fd,65536)
  finally:os.close(fd)
 result=json.loads(path.read_text())
 if result['exit_code']!=0:raise RuntimeError('Import failed: '+str(path))
def canonical(path):
 with pysam.VariantFile(str(path)) as vcf:
  samples=sorted(vcf.header.samples);records=[]
  for rec in vcf:
   records.append({'chrom':rec.chrom,'pos':rec.pos,'stop':rec.stop,'id':rec.id,'ref':rec.ref,
    'alts':rec.alts,'qual':rec.qual,'filter':sorted(rec.filter),'info':dict(rec.info),
    'samples':{s:{'fields':dict(rec.samples[s]),'phased':rec.samples[s].phased} for s in samples}})
 records.sort(key=lambda r:(r['chrom'],r['pos'],r['ref'],r['alts'] or ()))
 payload=json.dumps({'samples':samples,'records':records},sort_keys=True,allow_nan=False)
 return {'samples':len(samples),'records':len(records),'canonical_sha256':hashlib.sha256(payload.encode()).hexdigest()}
for n in [50,200]:
 results={}
 for mode in ['native','gatk_java','gatk_native_reader']:
  profile=f'N{n}_'+('native_import' if mode=='native' else mode)
  wait_receipt(PROF/profile/'status.json')
  out=OUT/f'N{n}_{mode}';out.mkdir(exist_ok=True);vcf=out/'export.vcf'
  if not (out/'status.json').exists():
   command=[sys.executable,str(ROOT/'benchmarks/docker_exec.py'),'run','--rm','--network','none','--user',f'{os.getuid()}:{os.getgid()}',
    '--cpus','2','--memory','8g','-v',f'{ROOT}:{ROOT}:ro','-v',f'{out}:{out}',
    '-v',f'{reference}:/reference.fa:ro','-v',f'{reference}.fai:/reference.fa.fai:ro',
    '-v',f'{ROOT}/data/D_genomicsdb/reference.dict:/reference.dict:ro',
    '--entrypoint','/gatk/gatk','broadinstitute/gatk:4.6.2.0','--java-options','-Xmx4g -XX:ActiveProcessorCount=2',
    'SelectVariants','-R','/reference.fa','-V','gendb://'+str(BASE/f'N{n}_{mode}'),'-O',str(vcf)]
   for interval in intervals:command+=['-L',interval]
   (out/'invocation.json').write_text(json.dumps({'command':command,'scope':'Five regions only; no performance inference'},indent=2))
   with (out/'stdout.log').open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
   (out/'status.json').write_text(json.dumps({'exit_code':r.returncode}))
   if r.returncode:raise RuntimeError('Export failed: '+str(out))
  elif json.loads((out/'status.json').read_text())['exit_code']!=0:raise RuntimeError('Prior export failed; use fresh namespace')
  results[mode]=canonical(vcf)
  print(n,mode,results[mode],flush=True)
 report={'n':n,'intervals':intervals,'results':results,
  'canonical_match':len({r['canonical_sha256'] for r in results.values()})==1,
  'scope':'Same GATK reader; five 1001bp regions; not whole-chromosome or genotype-caller equivalence'}
 (OUT/f'N{n}_comparison.json').write_text(json.dumps(report,indent=2))
os.close(controller)
