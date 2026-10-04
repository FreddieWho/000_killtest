import ctypes,datetime,hashlib,json,os,pathlib,select,subprocess,tarfile
ROOT=pathlib.Path(__file__).resolve().parents[1];DATA=ROOT/'data/C_hla';BASE=ROOT/'benchmarks/C_20260928';BASE.mkdir(exist_ok=True)
ENV=ROOT/'.envs/hla';env=dict(os.environ,PATH=str(ENV/'bin')+':/usr/bin:/bin')
archive=DATA/'PRG_MHC_GRCh38_withIMGT.tar.gz'
if archive.with_name(archive.name+'.aria2').exists():raise RuntimeError('Graph download not complete')
md5=hashlib.md5();sha=hashlib.sha256()
with archive.open('rb') as f:
 for block in iter(lambda:f.read(16*1024*1024),b''):md5.update(block);sha.update(block)
if md5.hexdigest()!='525a8aa0c7f357bf29fe2c75ef1d477d':raise RuntimeError('Official graph MD5 mismatch')
(DATA/'graph.receipt.json').write_text(json.dumps({'name':'HLA_LA_graph','path':str(archive.relative_to(ROOT)),'url':'https://zenodo.org/records/19336310/files/PRG_MHC_GRCh38_withIMGT.tar.gz?download=1','bytes':archive.stat().st_size,'md5':md5.hexdigest(),'sha256':sha.hexdigest(),'downloaded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Official HLA-LA graph reference'},indent=2))
graphs=DATA/'graphs';graphs.mkdir(exist_ok=True);graph=graphs/'PRG_MHC_GRCh38_withIMGT'
if not (graphs/'.extracted').exists():
 with tarfile.open(archive) as t:t.extractall(graphs,filter='data')
 (graphs/'.extracted').write_text('official MD5 verified\n')
def run(name,command):
 out=ROOT/'profiles/C'/name
 if (out/'status.json').exists():
  if json.loads((out/'status.json').read_text()).get('exit_code')==0:return
  raise RuntimeError('Existing failed run receipt: '+name)
 subprocess.run(['/opt/anaconda3/bin/python3',str(ROOT/'benchmarks/run_profile.py'),'--out',str(out),'--cwd',str(BASE),'--trace','--',*command],env=env,check=True)
def wait_files(files):
 libc=ctypes.CDLL(None,use_errno=True);fd=libc.inotify_init1(os.O_CLOEXEC)
 if fd<0:raise OSError(ctypes.get_errno(),'inotify')
 libc.inotify_add_watch(fd,os.fsencode(DATA),0x8|0x80)
 try:
  while not all(p.exists() for p in files):select.select([fd],[],[]);os.read(fd,65536)
 finally:os.close(fd)
if not (graph/'serializedGRAPH').exists():
 run('prepare_graph',[str(ENV/'opt/hla-la/bin/HLA-LA'),'--action','prepareGraph','--PRG_graph_dir',str(graph)])
for kind in ['WGS','WES']:
 print('WAIT_INPUT_RECEIPTS',kind,flush=True)
 pattern=kind+'_NA12878.'
 inputs=[p for p in DATA.glob(pattern+'*.bam') if '.mapped.' in p.name or '.unmapped.' in p.name]
 if len(inputs)!=2:
  suffix='low_coverage' if kind=='WGS' else 'exome'
  inputs=[DATA/f'{kind}_NA12878.{mode}.ILLUMINA.bwa.CEU.{suffix}.20121211.bam' for mode in ['mapped','unmapped']]
 wait_files([p.with_name(p.name+'.receipt.json') for p in inputs])
 headers=[]
 for p in inputs:
  header=subprocess.check_output(['samtools','view','-H',str(p)],env=env,text=True)
  headers.append([line for line in header.splitlines() if line.startswith('@SQ')])
 if headers[0]!=headers[1]:raise RuntimeError('Mapped/unmapped reference dictionaries differ')
 merged=BASE/(kind+'.full.bam')
 if not merged.exists():run(kind+'_merge',['samtools','merge','-@','4','-c','-p','-o',str(merged),*[str(p) for p in inputs]])
 if not merged.with_name(merged.name+'.bai').exists():
  subprocess.run(['samtools','index','-@','4',str(merged)],env=env,check=True)
 working=BASE/'typing';working.mkdir(exist_ok=True)
 run(kind+'_typing_8',['HLA-LA.pl','--BAM',str(merged),'--graph',graph.name,'--customGraphDir',str(graphs),'--sampleID',kind+'_NA12878_t8','--workingDir',str(working),'--maxThreads','8'])
print('WGS_WES_BASELINES_COMPLETE; RNA HLA-HD remains blocked until software is supplied',flush=True)
