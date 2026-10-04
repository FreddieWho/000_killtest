"""Extract completed exec intervals; overlapping process times are not additive wall time."""
import argparse, csv, json, pathlib, re

p = argparse.ArgumentParser()
p.add_argument('trace', type=pathlib.Path)
p.add_argument('--out', required=True, type=pathlib.Path)
a = p.parse_args()
active, pending, rows = {}, {}, []
line_re = re.compile(r'^(\d+)\s+(\d+\.\d+)\s+(.*)$')
exec_re = re.compile(r'^execve\(("(?:\\.|[^"\\])*"),\s*(\[.*\]),\s*0x')
last_epoch = None

def finish(pid, epoch, reason):
    if pid in active:
        row = active.pop(pid)
        row.update(end_epoch=epoch, wall_s=epoch-row['start_epoch'], end_reason=reason)
        rows.append(row)

with a.trace.open(errors='replace') as f:
    for line in f:
        m = line_re.match(line)
        if not m:
            continue
        pid, epoch, body = int(m[1]), float(m[2]), m[3]
        last_epoch = epoch
        ex = exec_re.match(body)
        if ex:
            try:
                executable, argv = json.loads(ex[1]), json.loads(ex[2])
            except json.JSONDecodeError:
                continue  # truncated argv is not reconstructed or executed
            candidate = dict(pid=pid, executable=executable, argv=argv, start_epoch=epoch)
            if '<unfinished ...>' in body:
                pending[pid] = candidate
            elif re.search(r'\)\s+= 0(?:\s|$)', body):
                finish(pid, epoch, 'next_exec')
                active[pid] = candidate
        elif '<... execve resumed>' in body:
            candidate = pending.pop(pid, None)
            if candidate and re.search(r'= 0(?:\s|$)', body):
                finish(pid, candidate['start_epoch'], 'next_exec')
                active[pid] = candidate
        elif body.startswith('+++ exited with') or body.startswith('+++ killed by'):
            finish(pid, epoch, body)

a.out.mkdir(parents=True, exist_ok=True)
with (a.out/'exec_intervals.json').open('w') as f:
    json.dump({'source':str(a.trace), 'completed':rows, 'open_intervals':list(active.values()),
               'last_trace_epoch':last_epoch, 'caveat':'Exec wall spans include child waits. Do not sum overlapping spans as workflow wall time.'}, f, indent=2)
with (a.out/'exec_intervals.tsv').open('w') as f:
    fields = ['pid','executable','start_epoch','end_epoch','wall_s','end_reason','argv']
    writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    for row in rows:
        writer.writerow(dict(row, argv=json.dumps(row['argv'])))
print(json.dumps({'completed_exec_intervals':len(rows),'open_intervals':len(active),'output':str(a.out)}))
