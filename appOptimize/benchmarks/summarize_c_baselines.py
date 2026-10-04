"""Baseline truth and process phases, preserving ambiguity and first-use index cost."""
import collections,csv,json,pathlib,re
ROOT=pathlib.Path(__file__).resolve().parents[1]
truth=next(r for r in csv.DictReader((ROOT/'data/C_hla/truth_20181129.tsv').open(),delimiter='\t') if r['Sample ID']=='NA12878')
report={}
for assay in ['WGS','WES']:
 p=ROOT/'benchmarks/C_20260928/typing'/(assay+'_NA12878_t8')/'hla/R1_bestguess.txt'
 rows=list(csv.DictReader(p.open(),delimiter='\t'));details=[]
 for locus in ['A','B','C','DQB1','DRB1']:
  expected=sorted([truth[f'HLA-{locus} 1'],truth[f'HLA-{locus} 2']]);calls=[r for r in rows if r['Locus'].removeprefix('HLA-')==locus]
  alternatives=[sorted(set(re.findall(r'(?<![\d:])(\d+:\d+)(?::\d+)*(?:[A-Z])?',r['Allele']))) for r in calls]
  confirmed=[a[0] for a in alternatives if len(a)==1];observed=sorted(confirmed) if len(confirmed)==2 else None
  compatible=max(sum(e[i] in alternatives[i] for i in [0,1]) for e in [expected,expected[::-1]]) if len(alternatives)==2 else 0
  details.append({'locus':locus,'truth':expected,'reported':[r['Allele'] for r in calls],'two_field_alternatives':alternatives,'confirmed_alleles':sum((collections.Counter(expected)&collections.Counter(confirmed)).values()),'compatible_alleles':compatible,'unambiguous_genotype_match':observed==expected,'coverage':[float(r['AverageCoverage']) for r in calls]})
 base=ROOT/'profiles/C'/(assay+'_typing_8');time=json.loads((base/'status.json').read_text())['wall_s'];intervals=json.loads((base/'analysis/exec_intervals.json').read_text())['completed'];bwas=[r for r in intervals if pathlib.Path(r['executable']).name=='bwa'];backend=next(r for r in intervals if pathlib.Path(r['executable']).name=='HLA-LA' and '--action' in r['argv']);index_s=sum(r['wall_s'] for r in bwas if 'index' in r['argv']);mem_s=sum(r['wall_s'] for r in bwas if 'mem' in r['argv'])
 report[assay]={'typing_wall_s':time,'backend_wall_s':backend['wall_s'],'first_use_bwa_index_wall_s':index_s,'bwa_mem_wall_s':mem_s,'outside_backend_s':time-backend['wall_s'],'backend_other_s':backend['wall_s']-index_s-mem_s,'five_locus_truth':details,'allele_denominator':10,'confirmed_alleles':sum(d['confirmed_alleles'] for d in details),'compatible_alleles':sum(d['compatible_alleles'] for d in details),'unambiguous_matching_genotypes':sum(d['unambiguous_genotype_match'] for d in details)}
report['limits']=['WGS historical low-coverage, not 30x clinical WGS','Compatible ambiguity is not confirmed allele accuracy','No unique four-field truth available','Backend other combines graph, inference, sorting and other work; not a pure DP measurement','Original WGS includes first-use shared BWA reference index; replay reuses index and cannot be compared as a cold-end-to-end speedup','BAM merge input preparation and graph preparation are separately measured']
(ROOT/'profiles/closure_20261002/C_baseline_evidence.json').write_text(json.dumps(report,indent=2));print(json.dumps({a:{k:v for k,v in d.items() if k!='five_locus_truth'} for a,d in report.items() if a!='limits'}))
