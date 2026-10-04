import sys,os,json,platform,importlib.metadata as m
from pathlib import Path
packages=['numpy','scipy','pandas','anndata','scanpy','liana','mudata','cellphonedb','fastccc','statsmodels','patsy','networkx','scikit-learn','h5py','matplotlib','tabulate','numba','docrep','geosketch','fbpca','numpy-groupies']
d={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'logical_cpus':os.cpu_count(),'allowed_cpus':len(os.sched_getaffinity(0)),'packages':{},'thread_variables':{v:os.environ.get(v) for v in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_CACHE_DIR','MPLCONFIGDIR','LD_LIBRARY_PATH']}}
for p in packages:
 try:d['packages'][p]=m.version(p)
 except m.PackageNotFoundError:d['packages'][p]='NOT_FOUND'
Path('results/environment.json').write_text(json.dumps(d,indent=2))
print(json.dumps(d['packages'],indent=2))
