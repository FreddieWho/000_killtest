import concurrent.futures as cf
import datetime, hashlib, json, pathlib, urllib.request, xml.etree.ElementTree as ET
import os, sys
os.environ["CURL_CA_BUNDLE"]="/etc/ssl/certs/ca-certificates.crt"
os.environ["SSL_CERT_FILE"]="/etc/ssl/certs/ca-certificates.crt"
import pysam
ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/'data/D_genomicsdb'; OUT.mkdir(parents=True,exist_ok=True)
ns={'s':'http://s3.amazonaws.com/doc/2006-03-01/'}
xml=ET.parse(ROOT/'sources/D_s3_listing.xml')
prefixes=[x.text for x in xml.findall('s:CommonPrefixes/s:Prefix',ns)][:int(sys.argv[1]) if len(sys.argv)>1 else 200]
def run(prefix):
 sample=prefix.rstrip('/').split('/')[-1]
 url='https://1000genomes-dragen.s3.us-west-2.amazonaws.com/'+prefix+sample+'.hard-filtered.gvcf.gz'
 receipt=OUT/(sample+'.receipt.json')
 if receipt.exists():return json.loads(receipt.read_text())
 idx=OUT/(sample+'.source.tbi'); dest=OUT/(sample+'.chr22.g.vcf.gz'); temp=dest.with_suffix('.partial')
 if not idx.exists(): urllib.request.urlretrieve(url+'.tbi',idx)
 n=0
 with pysam.TabixFile(url,index=str(idx)) as source,pysam.BGZFile(str(temp),'w') as output:
  for line in source.header:output.write((line+'\n').encode())
  for line in source.fetch('chr22'):
   output.write((line+'\n').encode()); n+=1
 if n==0:raise ValueError('Empty chromosome for '+sample)
 temp.rename(dest);pysam.tabix_index(str(dest),preset='vcf',force=True)
 rec={'name':sample,'path':str(dest.relative_to(ROOT)),'url':url,'region':'chr22','records':n,'bytes':dest.stat().st_size,'sha256':hashlib.file_digest(dest.open('rb'),'sha256').hexdigest(),'downloaded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Real 1000 Genomes DRAGEN gVCF cohort import benchmark; region extracted without changing records'}
 receipt.write_text(json.dumps(rec,indent=2));return rec
with cf.ThreadPoolExecutor(max_workers=4) as pool, (ROOT/'profiles/D_download_receipts.jsonl').open('a') as log:
 jobs={pool.submit(run,p):p for p in prefixes}
 for job in cf.as_completed(jobs):
  try:r=job.result();log.write(json.dumps(r)+'\n');log.flush();print(r['name'],r['records'],r['bytes'],flush=True)
  except Exception as e:log.write(json.dumps({'prefix':jobs[job],'error':repr(e)})+'\n');log.flush();print('ERROR',jobs[job],repr(e),flush=True)
rows=[json.loads(p.read_text()) for p in sorted(OUT.glob('*.receipt.json'))]
for n in [50,200]:
 if len(rows)>=n:
  (OUT/f'samplemap_{n}.tsv').write_text(''.join(r['name']+'\t'+str(ROOT/r['path'])+'\n' for r in rows[:n]))
print('COMPLETE',len(rows),flush=True)
expected={p.rstrip('/').split('/')[-1] for p in prefixes}
missing=sorted(expected-{r['name'] for r in rows})
if missing:
 print('INCOMPLETE_REQUIRED_SAMPLES',json.dumps(missing),flush=True)
 raise SystemExit(1)
