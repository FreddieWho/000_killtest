"""Test DP metadata causality on five exported intervals without altering original arrays."""
import collections,hashlib,json,os,pathlib,shutil,subprocess,sys
import pysam
ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/'profiles/D_schema_diagnostic_20261002';SHADOW=ROOT/'benchmarks/D_schema_diagnostic_20261002'
OUT.mkdir(exist_ok=False);SHADOW.mkdir(exist_ok=False)
source=ROOT/'benchmarks/D_20260928/N200_native'
for p in source.iterdir():
 if p.is_dir():(SHADOW/p.name).symlink_to(p,target_is_directory=True)
 else:shutil.copy2(p,SHADOW/p.name)
p=SHADOW/'vidmap.json';doc=json.loads(p.read_text());before=None
for field in doc['fields']:
 if field['name']=='DP':before=dict(field);field.pop('VCF_field_combine_operation',None)
p.write_text(json.dumps(doc,indent=2))
cmd=json.loads((ROOT/'profiles/D_schema_diagnostic_20260928/invocation.json').read_text())['command']
cmd=[x.replace('profiles/D_schema_diagnostic_20260928','profiles/D_schema_diagnostic_20261002').replace('benchmarks/D_schema_diagnostic_20260928','benchmarks/D_schema_diagnostic_20261002') for x in cmd]
cmd=[sys.executable,str(ROOT/'benchmarks/docker_exec.py'),*cmd[1:]]
(OUT/'invocation.json').write_text(json.dumps({'command':cmd,'only_semantic_edit':{'DP_before':before,'removed':'VCF_field_combine_operation'}},indent=2))
with (OUT/'stdout.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
(OUT/'status.json').write_text(json.dumps({'exit_code':r.returncode}))
if r.returncode:raise SystemExit(r.returncode)
def records(path):
 with pysam.VariantFile(str(path)) as f:
  names=sorted(f.header.samples)
  rs=[{'chrom':r.chrom,'pos':r.pos,'stop':r.stop,'id':r.id,'ref':r.ref,'alts':r.alts,'qual':r.qual,'filter':sorted(r.filter),'info':dict(r.info),'samples':{s:{'fields':dict(r.samples[s]),'phased':r.samples[s].phased} for s in names}} for r in f]
 rs.sort(key=lambda r:(r['chrom'],r['pos'],r['ref'],r['alts'] or ()))
 return {'samples':names,'records':rs}
raw=records(ROOT/'profiles/D_export_comparison/N200_native/export.vcf');gatk=records(ROOT/'profiles/D_export_comparison/N200_gatk_java/export.vcf');fixed=records(OUT/'export.vcf')
counts=collections.Counter();details=[]
assert raw['samples']==gatk['samples'] and len(raw['records'])==len(gatk['records'])
for a,b in zip(raw['records'],gatk['records']):
 assert (a['chrom'],a['pos'],a['ref'],a['alts'])==(b['chrom'],b['pos'],b['ref'],b['alts'])
 changed=[k for k in a if a[k]!=b[k]]
 if changed:
  for k in changed:counts[k]+=1
  details.append({'pos':a['pos'],'fields':changed,'native_info':a['info'],'gatk_info':b['info']})
report={'samples':len(gatk['samples']),'records':len(gatk['records']),'original_differing_records':len(details),'changed_top_level_fields':dict(counts),'all_original_differences_only_INFO_DP':all(d['fields']==['info'] and {k:v for k,v in d['native_info'].items() if k!='DP'}=={k:v for k,v in d['gatk_info'].items() if k!='DP'} for d in details),'DP_metadata_normalized_exact_match':fixed==gatk,'canonical_sha256':hashlib.sha256(json.dumps(fixed,sort_keys=True,allow_nan=False).encode()).hexdigest(),'scope':'Five 1001bp chr22 intervals; original arrays unchanged; not whole-chromosome equivalence','differences':details}
(OUT/'comparison.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='differences'}))
