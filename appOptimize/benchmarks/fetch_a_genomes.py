import concurrent.futures as cf,datetime,gzip,hashlib,json,pathlib,shutil,urllib.request,urllib.parse
ROOT=pathlib.Path(__file__).resolve().parents[1]; OUT=ROOT/'data/A_pplacer';(OUT/'W1').mkdir(exist_ok=True);(OUT/'W2').mkdir(exist_ok=True)
raw=json.loads((ROOT/'sources/A_uncultured_MAGs.json').read_text())['reports']
rows=[r for r in raw if 1000000<=int(r.get('assembly_stats',{}).get('total_sequence_length',0))<=10000000][:100]
(ROOT/'sources/A_selected_MAGs.json').write_text(json.dumps(rows,indent=2))
def get(r,group):
 acc=r['accession'];name=r['assembly_info']['assembly_name'].replace(' ','_');digits=acc.split('_')[1].split('.')[0];base=acc+'_'+name
 url='https://ftp.ncbi.nlm.nih.gov/genomes/all/'+acc[:3]+'/'+digits[:3]+'/'+digits[3:6]+'/'+digits[6:9]+'/'+urllib.parse.quote(base)+'/'+urllib.parse.quote(base)+'_genomic.fna.gz'
 dest=OUT/group/(acc+'.fna.gz');receipt=dest.with_suffix('.receipt.json')
 if receipt.exists():return json.loads(receipt.read_text())
 urllib.request.urlretrieve(url,dest);n=0;length=0
 with gzip.open(dest,'rt') as source,(OUT/group/(acc+'.fna')).open('w') as out:
  for line in source:
   out.write(line)
   if line.startswith('>'):n+=1
   else:length+=len(line.strip())
 rec={'name':acc,'path':str(dest.relative_to(ROOT)),'url':url,'bytes':dest.stat().st_size,'sha256':hashlib.file_digest(dest.open('rb'),'sha256').hexdigest(),'sequences':n,'bases':length,'downloaded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'GTDB-Tk '+group+' real genome workload','assembly_release':r.get('assembly_info',{}).get('release_date'),'organism':r.get('organism',{}).get('organism_name')}
 receipt.write_text(json.dumps(rec,indent=2));return rec
print(get({'accession':'GCF_000005845.2','assembly_info':{'assembly_name':'ASM584v2'},'organism':{'organism_name':'Escherichia coli K-12 MG1655'}},'W1'),flush=True)
print(get(rows[0],'W2'),flush=True)
with cf.ThreadPoolExecutor(max_workers=4) as pool:
 jobs={pool.submit(get,r,'W2'):r['accession'] for r in rows[1:]}
 for j in cf.as_completed(jobs):
  try:r=j.result();print(r['name'],r['bases'],flush=True)
  except Exception as e:print('ERROR',jobs[j],repr(e),flush=True)
print('FINISHED',flush=True)
