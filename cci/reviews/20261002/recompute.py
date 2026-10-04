"""Independent review of frozen results; does not rerun or replace research runs."""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import r2_score

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
selected = json.loads((ROOT / 'results/selected_runs.json').read_text())
programs = pd.read_csv(ROOT / selected['kang'] / 'programs.tsv', sep='\t')
result = {'scope': 'Frozen-input audit, not independent biological validation', 'ic': {}}
graph_audit = {'checked': 0, 'mismatched_size_or_certificate': 0}
for name, frame in programs.groupby('database'):
    z = np.load(ROOT / selected['kang'] / f'{name}_communication.npz')
    res = pd.read_csv(ROOT / selected['kang'] / f'{name}_mapped_resource.tsv', sep='\t')
    features = pd.read_csv(ROOT / selected['kang'] / f'{name}_features.tsv', sep='\t')
    labels, samples, k = features.sender.unique().tolist(), z['samples'].tolist(), len(res)
    simple = ~(res.ligand_subunits.str.contains('|', regex=False) | res.receptor_subunits.str.contains('|', regex=False))
    for row in frame.itertuples():
        offset = (labels.index(row.sender) * len(labels) + labels.index(row.receiver)) * k
        active = z['support'][samples.index(row.sample), offset:offset + k] & simple.to_numpy()
        edges = [('L:' + r.ligand_subunits, 'R:' + r.receptor_subunits) for r in res[active].itertuples()]
        graph = nx.Graph(edges)
        # Independent reconstruction and a general maximum matching algorithm.
        optimum = len(nx.max_weight_matching(graph, maxcardinality=True))
        targets = set(str(row.targets).split(';'))
        valid = optimum == row.ic100 and all(a in targets or b in targets for a, b in edges)
        graph_audit['checked'] += 1
        graph_audit['mismatched_size_or_certificate'] += int(not valid)
result['simple_graph_audit'] = graph_audit
for name, g in programs.groupby('database'):
    g = g[g['count'] > 0].copy()
    x, y = np.log1p(g['count'].to_numpy()), g.ic100.to_numpy()
    pred = np.zeros(len(g))
    for donor in g.patient.unique():
        test = (g.patient == donor).to_numpy()
        model = IsotonicRegression(out_of_bounds='clip').fit(x[~test], y[~test])
        pred[test] = model.predict(x[test])
    pairs = g.sender.astype(str) + '>' + g.receiver.astype(str)
    centered_y = y - pd.Series(y, index=g.index).groupby(pairs).transform('mean').to_numpy()
    centered_pred = pred - pd.Series(pred, index=g.index).groupby(pairs).transform('mean').to_numpy()
    result['ic'][name] = {
        'n_programs': len(g), 'heldout_donor_R2': r2_score(y, pred),
        'per_donor_R2': {str(d): r2_score(y[g.patient == d], pred[g.patient == d]) for d in g.patient.unique()},
        'rounded_prediction_exact_accuracy': float(np.mean(np.rint(pred) == y)),
        'MAE': float(np.mean(abs(pred - y))),
        'pair_centered_R2_descriptive': r2_score(centered_y, centered_pred),
        'within_pair_fraction_of_variance': float(np.var(centered_y) / np.var(y)),
        'centering_note': 'Post-hoc pooled pair means: descriptive decomposition only, not a newly fitted or independently validated model.',
    }
result['null_table_means'] = {}
for label, path in [('pooled_full', selected['pooled_full']), ('prototype', selected['kang'])] + [(Path(p).name, p) for p in selected['native_analyses']]:
    frame = pd.read_csv(ROOT / path / 'null_replicates.tsv', sep='\t')
    result['null_table_means'][label] = frame.groupby('method')[['FPR', 'FDP']].mean().to_dict(orient='index')
z = np.load(ROOT / selected['kang'] / 'CellChat_communication.npz')
samples, scores = z['samples'].tolist(), z['scores']
donors = sorted({s.split('__')[0] for s in samples})
good = np.isfinite(scores).all(axis=0) & (np.nanstd(scores, axis=0) > 1e-12)
diff = scores[[samples.index(d + '__stim') for d in donors]][:, good] - scores[[samples.index(d + '__ctrl') for d in donors]][:, good]
signs = np.array(list(itertools.product([-1., 1.], repeat=len(donors))))
means = abs(signs @ diff / len(donors))
rounded = stats.rankdata(-np.round(means, 12), axis=0, method='max') / len(signs)
unrounded = stats.rankdata(-means, axis=0, method='max') / len(signs)
result['exact_p_rounding_changed_entries'] = int((rounded != unrounded).sum())
targets = np.random.default_rng(20260928).choice(diff.shape[1], diff.shape[1] // 10, replace=False)
result['prototype_injection'] = {'n_targets': len(targets), 'zero_SD_targets': int((np.std(diff[:, targets], axis=0, ddof=1) == 0).sum())}
(OUT / 'metrics.json').write_text(json.dumps(result, indent=2))
print(json.dumps({'status': 'DONE', 'output': str(OUT / 'metrics.json'), 'ic_R2': {k: v['heldout_donor_R2'] for k, v in result['ic'].items()}}))
