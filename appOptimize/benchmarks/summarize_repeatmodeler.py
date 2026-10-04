"""Preserve RepeatModeler log phase labels without relabelling masking as I/O."""
import csv, json, pathlib, re
ROOT = pathlib.Path(__file__).resolve().parents[1]
fields = ['genome','round','state','round_s','extraction_s','trf_s','te_mask_s',
          'sampling_label_s','comparison_s','repeatscout_s','recon_s','instance_s',
          'refinement_s','families','comparison_fraction','noncomparison_removal_ceiling',
          'measurement_scope']
patterns = {'extraction_s':'Sequence extraction', 'trf_s':'TRFMask time',
            'te_mask_s':'TE Masking time','sampling_label_s':'Sampling Time:',
            'comparison_s':'Comparison Time:','repeatscout_s':'- RepeatScout:',
            'recon_s':'RECON Elapsed:', 'instance_s':'Instance Gathering:',
            'refinement_s':'Family Refinement:', 'round_s':'Round Time:'}
rows=[]
for genome, path in [('dmel','profiles/B_dmel.log'),('rice','profiles/B_rice/stdout.log'),('dmel_full','profiles/B_dmel_full_20260928/stdout.log')]:
    if not (ROOT/path).exists():continue
    current=None
    for line in (ROOT/path).read_text(errors='replace').splitlines():
        m=re.search(r'RepeatModeler Round #\s*(\d+)',line)
        if m:
            current={'genome':genome,'round':int(m[1]),'state':'INCOMPLETE',
                     'measurement_scope':'baseline_log' if genome!='rice' else 'initial_trace_contamination_detach_boundary_required'}
            rows.append(current)
        if current is None:
            continue
        t=re.search(r'(\d+):(\d+):(\d+) \(hh:mm:ss\)',line)
        if not t:
            continue
        seconds=int(t[1])*3600+int(t[2])*60+int(t[3])
        for field, label in patterns.items():
            if label in line:
                current[field]=current.get(field,0)+seconds
        if 'Round Time:' in line:
            current['state']='ROUND_COMPLETE'
            f=re.search(r'(\d+) families discovered',line)
            if f:current['families']=int(f[1])
            if current.get('comparison_s') and current.get('round_s'):
                current['comparison_fraction']=current['comparison_s']/current['round_s']
                current['noncomparison_removal_ceiling']=current['round_s']/current['comparison_s']
rice_first=next((r for r in rows if r['genome']=='rice' and r['round']==1 and r.get('round_s')),None)
detach=ROOT/'profiles/B_rice/tracer_detach.json';invocation=ROOT/'profiles/B_rice/invocation.json'
if rice_first and detach.exists() and invocation.exists():
    started=json.loads(invocation.read_text())['started_epoch'];stopped=json.loads(detach.read_text())['epoch']
    earliest_round1_end=started+rice_first['round_s']
    if earliest_round1_end-stopped>2:
        for row in rows:
            if row['genome']=='rice' and row['round']>1:row['measurement_scope']='post_detach_stage_log_not_full_workflow'
        (ROOT/'profiles/B_rice/post_detach_round_boundary.json').write_text(json.dumps({
            'workflow_started_epoch':started,'round1_duration_s':rice_first['round_s'],
            'round1_end_lower_bound_epoch':earliest_round1_end,'tracer_detached_epoch':stopped,
            'margin_s':earliest_round1_end-stopped,
            'inference':'Round 1 starts no earlier than workflow start; its end lower bound is after detach, so rounds 2 onward are post-detach.'},indent=2))
out=ROOT/'benchmarks/B_round_metrics.tsv'
with out.open('w') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
    for row in rows:w.writerow({k:row.get(k,'NA') for k in fields})
print('Saved',len(rows),'rounds; incomplete rounds excluded from ceiling estimates.')
