"""Profile an explicit argv in an isolated container, writing only a fresh output directory.

Usage: python perf_replay.py spec.json
Spec: {"out":"profiles/replay_name", "argv":[...], "path_prefix":".../bin", "cpus":8}
Inputs and the executable must already be available under the project root.
This is a diagnostic replay; its elapsed time is not an uninstrumented speed result.
"""
import json, os, pathlib, subprocess, sys, time

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = json.loads(pathlib.Path(sys.argv[1]).read_text())
out = (ROOT / spec['out']).resolve()
if ROOT not in out.parents or out == ROOT or out.exists():
    raise ValueError('Output must be a fresh directory below the project root')
argv = spec['argv']
if not argv or not all(isinstance(x, str) for x in argv):
    raise ValueError('argv must be a nonempty string array')
out.mkdir(parents=True)
image = spec.get('image','broadinstitute/gatk:4.6.2.0')
docker = [sys.executable,str(ROOT/'benchmarks/docker_exec.py')]
loader = '/prof/ld-linux-x86-64.so.2'
base = docker + ['run','--rm','--user','0:0','--group-add',str(os.getgid()),'--read-only','--network','none','--cap-drop','ALL',
        '--cap-add','SYS_ADMIN','--security-opt','seccomp='+str(ROOT/'profiles/profiler_tools/perf-seccomp.json'),'--pids-limit','512',
        '--cpus',str(spec.get('cpus',8)),'--memory',spec.get('memory','256g'),
        '--tmpfs','/tmp:rw,nosuid,nodev,size=2g',
        '-v',f'{ROOT}:{ROOT}:ro','-v',f'{out}:{out}:rw',
        '-v',f'{ROOT}/profiles/profiler_tools:/prof:ro','-w',spec.get('cwd',str(out)),
        '-e','PATH='+spec.get('path_prefix','')+':/usr/bin:/bin']
for key, value in spec.get('env',{}).items():
    base += ['-e', key+'='+str(value)]
base += ['--entrypoint',loader,image,'--library-path','/prof']
cmd = base + ['/prof/time','-v','-o',str(out/'time.txt'),loader,'--library-path','/prof',
              '/prof/perf','record','--no-bpf-event','-e','cpu-clock:u','-F','49',
              '-o',str(out/'perf.data'),'--',*argv]
(out/'invocation.json').write_text(json.dumps({'spec':spec,'command':cmd,
    'measurement_scope':'diagnostic_replay_49Hz_not_clean_baseline',
    'started_epoch':time.time()},indent=2))
subprocess.run(docker+['run','--rm','--user','0:0','--network','none','--cap-drop','ALL','--cap-add','CHOWN',
                '-v',f'{out}:/output','--entrypoint','chown',image,f'0:{os.getgid()}','/output'],check=True)
with (out/'stdout.log').open('w') as log:
    result = subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
if (out/'perf.data').exists():
    with (out/'perf_report.txt').open('w') as report:
        subprocess.run(base+['/prof/perf','report','--stdio','--no-children','--percent-limit','0.5',
                            '-i',str(out/'perf.data')],stdout=report,stderr=subprocess.STDOUT)
    with (out/'perf_header.txt').open('w') as header:
        subprocess.run(base+['/prof/perf','report','--header-only','-i',str(out/'perf.data')],
                       stdout=header,stderr=subprocess.STDOUT)
# Docker owns only this dedicated output tree; restore ownership for subsequent analysis.
subprocess.run(docker+['run','--rm','--network','none','--cap-drop','ALL','--cap-add','CHOWN',
                '-v',f'{out}:/output','--entrypoint','chown',image,'-R',f'{os.getuid()}:{os.getgid()}','/output'],check=True)
(out/'status.json').write_text(json.dumps({'exit_code':result.returncode,'finished_epoch':time.time(),
    'status':'COMPLETE' if result.returncode==0 else 'FAILED'},indent=2))
raise SystemExit(result.returncode)
