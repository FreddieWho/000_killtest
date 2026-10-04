import csv,datetime,hashlib,json,pathlib,subprocess,xml.etree.ElementTree as ET
ROOT=pathlib.Path(__file__).resolve().parents[1];OUT=ROOT/'data/C_hla';OUT.mkdir(exist_ok=True)
items=[]
for label in ['WGS','WES']:
 tree=ET.parse(ROOT/f'sources/C_{label}_phase3_listing.xml')
 for row in tree.findall('{*}Contents'):
  key=row.find('{*}Key').text
  if ('.mapped.' in key or '.unmapped.' in key) and key.endswith(('.bam','.bam.bai')):
   name=label+'_'+key.rsplit('/',1)[-1];items.append((name,'https://1000genomes.s3.amazonaws.com/'+key,int(row.find('{*}Size').text)))
with (ROOT/'sources/C_RNA_metadata.tsv').open() as f:
 row=next(r for r in csv.DictReader(f,delimiter='\t') if r['run_accession']=='ERR356372')
 for url,size in zip(row['fastq_ftp'].split(';'),row['fastq_bytes'].split(';')):items.append((url.rsplit('/',1)[-1],'https://'+url,int(size)))
manifest=OUT/'download.aria2.txt';manifest.write_text(''.join(url+'\n out='+name+'\n' for name,url,_ in items))
status=subprocess.run(['aria2c','--continue=true','--max-concurrent-downloads=3','--max-connection-per-server=4','--split=4','--file-allocation=none','--summary-interval=300','--console-log-level=warn','--dir='+str(OUT),'--input-file='+str(manifest)]).returncode
for name,url,size in items:
 p=OUT/name
 if not p.exists() or p.stat().st_size!=size or p.with_name(p.name+'.aria2').exists():print('NOT_READY',name,flush=True);continue
 rec={'name':name,'path':str(p.relative_to(ROOT)),'url':url,'bytes':size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'downloaded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Real NA12878 HLA '+('RNA-seq' if name.startswith('ERR') else name[:3])+' workload'}
 p.with_name(p.name+'.receipt.json').write_text(json.dumps(rec,indent=2));print('READY',name,flush=True)
raise SystemExit(status)
