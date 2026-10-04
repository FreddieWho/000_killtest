"""Same two resources, scores and paired fits in separate fresh processes."""
import argparse, json, resource, time
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from patsy import dmatrix
from cci_core import aggregate, communication, fit_linear
from streaming_input import read_lr_counts

p = argparse.ArgumentParser()
p.add_argument('mode', choices=['full', 'stream'])
p.add_argument('out')
p.add_argument('--chunk-size', type=int, default=4096)
args = p.parse_args()
out = Path(args.out); out.mkdir(parents=True, exist_ok=False)
root = Path('results/kang_run_03')
names = ['CellChat', 'LIANA_CellPhoneDB']
resources = [pd.read_csv(root / f'{n}_mapped_resource.tsv', sep='\t') for n in names]
t = time.perf_counter()
if args.mode == 'full':
    a = ad.read_h5ad('data/raw/kang_2018.h5ad')
    counts, raw_obs, genes, lib = a.X, a.obs, a.var_names, None
else:
    counts, raw_obs, genes, lib = read_lr_counts('data/raw/kang_2018.h5ad', resources, args.chunk_size)
obs = pd.DataFrame({'patient': raw_obs.replicate.astype(str), 'condition': raw_obs.label.astype(str), 'celltype': raw_obs.cell_type.astype(str)})
obs['sample'] = obs.patient + '__' + obs.condition
agg = aggregate(counts, obs, library_size=lib)
checks, fits = [], {}
for name, r in zip(names, resources):
    result = communication(agg, genes, r, obs)
    v = result['scores']
    samples = result['samples']
    good = np.isfinite(v).all(axis=0) & (np.nanstd(v, axis=0) > 1e-12)
    meta = pd.DataFrame({'patient': [s.split('__')[0] for s in samples], 'condition': [s.split('__')[1] for s in samples]})
    x = dmatrix('~ condition + patient', meta, return_type='dataframe')
    b, se, pv, df = fit_linear(v[:, good], x, x.columns.get_loc('condition[T.stim]'))
    fits[name] = dict(scores=v, support=result['support'], good=good, beta=b, se=se, p=pv)
elapsed = time.perf_counter() - t
peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
# Validation and disk output excluded from the measured computational interval.
for name, values in fits.items():
    with np.load(root / f'{name}_communication.npz') as old:
        np.testing.assert_allclose(values['scores'], old['scores'], rtol=1e-12, atol=1e-12)
        np.testing.assert_array_equal(values['support'], old['support'])
    np.savez_compressed(out / f'{name}.npz', **values)
    checks.append(name)
receipt = dict(mode=args.mode, seconds=elapsed, peak_RSS_MiB=peak, cells=len(obs), genes=len(genes), resources=checks, frozen_scores_and_support='PASS', chunk_size=args.chunk_size, status='DONE')
(out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt))
