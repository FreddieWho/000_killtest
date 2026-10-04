"""Summarize immutable baselines; distinguish measured quantities from conditional ceilings."""
import collections,csv,hashlib,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1];OUT=ROOT/'profiles/closure_20261002';OUT.mkdir(exist_ok=True)
def timefile(path):
 d={}
 for line in path.read_text().splitlines():
  if ': ' in line:
   k,v=line.strip().split(': ',1);d[k]=v
 wall=next(v for k,v in d.items() if k.startswith('Elapsed (wall clock)'))
 return {'wall_s':sum(float(v)*60**i for i,v in enumerate(reversed(wall.split(':')))), 'user_s':float(d['User time (seconds)']),'system_s':float(d['System time (seconds)']),'peak_rss_kib':int(d['Maximum resident set size (kbytes)'])}
rows=[r for r in csv.DictReader((ROOT/'benchmarks/B_round_metrics.tsv').open(),delimiter='\t') if r['genome']=='dmel_full'];b=timefile(ROOT/'profiles/B_dmel_full_20260928/time.txt');comparison=sum(float(r['comparison_s']) for r in rows if r['comparison_s']!='NA')
phases={k:sum(float(r[k]) for r in rows if r[k]!='NA') for k in ['round_s','extraction_s','sampling_label_s','comparison_s','repeatscout_s','recon_s','instance_s','refinement_s']}
# Round-2+ extraction is contained in sampling_label_s; do not double count it.
components={'sampling_label_includes_TRF_and_TE_mask':phases['sampling_label_s'],'round1_sequence_extraction':float(rows[0]['extraction_s']),'similarity_comparison':comparison,'RepeatScout':phases['repeatscout_s'],'RECON_cluster':phases['recon_s'],'instance_gathering':phases['instance_s'],'refinement':phases['refinement_s']}
components['unattributed_including_final_classification_parse_IO_scheduling']=b['wall_s']-sum(components.values())
lib=ROOT/'benchmarks/B_dmel_full_20260928/dmel-families.fa';headers=[l[1:] for l in lib.read_text().splitlines() if l.startswith('>')]
b.update(rounds=rows,phase_labels_s=phases,additive_log_partition_s=components,comparison_fraction=comparison/b['wall_s'],conditional_ceiling_if_comparison_unchanged=b['wall_s']/comparison,mean_CPU_cores=(b['user_s']+b['system_s'])/b['wall_s'],final_families=len(headers),family_class_counts=dict(collections.Counter(h.split('#',1)[1].split()[0] if '#' in h else 'UNLABELLED' for h in headers)),library_sha256=hashlib.sha256(lib.read_bytes()).hexdigest(),resource_samples=json.loads((ROOT/'profiles/B_dmel_full_20260928/resource_summary.json').read_text()),limits=['Comparison includes its own orchestration and output; 1.17x is conditional, not an unconditional bound on every possible change.','Parse, IO, scheduling and other lack separable direct wall-time measures; residual is not labelled pure IO.','No architecture prototype, so pairwise family identity/remask equivalence is not claimed.','CPU sampling of two batches supports alignment-dominated interpretation, not all batches as a theorem.'])
(OUT/'B_evidence.json').write_text(json.dumps(b,indent=2))
a_rows=list(csv.DictReader((ROOT/'benchmarks/A_W2_8/classify/gtdbtk.bac120.summary.tsv').open(),delimiter='\t'));execs=json.loads((ROOT/'profiles/A/W2_8/analysis/exec_intervals.json').read_text())['completed'];calls=[r for r in execs if pathlib.Path(r['executable']).name=='pplacer' and '-o' in r['argv']];a=timefile(ROOT/'profiles/A/W2_8/time.txt');segments=sorted((r['start_epoch'],r['end_epoch']) for r in calls);merged=[]
for start,end in segments:
 if merged and start<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],end)
 else:merged.append([start,end])
placement_wall=sum(y-x for x,y in merged);a.update(input_genomes=len(a_rows),classification_methods=dict(collections.Counter(r['classification_method'] for r in a_rows)),placed_genomes=[r['user_genome'] for r in a_rows if r['classification_method']!='ani_screen'],pplacer_calls=len(calls),pplacer_union_wall_s=placement_wall,pplacer_fraction=placement_wall/a['wall_s'],conditional_speedup_with_3x_pplacer=1/(1-placement_wall/a['wall_s']+placement_wall/a['wall_s']/3),old_zero_count_invalid_reason='GTDB-Tk removed intermediate jplace; trace contains 9 placement calls and summary has 54 topology/RED genomes',calls=calls)
(OUT/'A_baseline_evidence.json').write_text(json.dumps(a,indent=2))
d={}
for n in [50,200]:
 entries={mode:timefile(ROOT/f'profiles/D_20260928/N{n}_{mode}/time.txt') for mode in ['native_init','native_import','gatk_java','gatk_native_reader']}
 entries['native_init_plus_import_s']=entries['native_init']['wall_s']+entries['native_import']['wall_s'];entries['Java_over_native_total']=entries['gatk_java']['wall_s']/entries['native_init_plus_import_s'];entries['Java_over_bypass']=entries['gatk_java']['wall_s']/entries['gatk_native_reader']['wall_s'];entries['Java_over_bypass_peak_RSS']=entries['gatk_java']['peak_rss_kib']/entries['gatk_native_reader']['peak_rss_kib'];d[str(n)]=entries
(OUT/'D_evidence.json').write_text(json.dumps(d,indent=2))
print(json.dumps({'B_comparison_fraction':b['comparison_fraction'],'B_ceiling':b['conditional_ceiling_if_comparison_unchanged'],'B_families':len(headers),'A_placed':len(a['placed_genomes']),'A_fraction':a['pplacer_fraction'],'D':d},indent=2))
