"""Bounded input memory; retain full-library denominators before LR selection."""
import anndata as ad
import numpy as np
from scipy import sparse


def read_lr_counts(path, resources, chunk_size=4096):
    if chunk_size < 1:
        raise ValueError('chunk_size must be positive')
    a = ad.read_h5ad(path, backed='r')
    try:
        if not a.var_names.is_unique:
            raise ValueError('Aggregate duplicate gene symbols before LR selection')
        needed = set()
        for r in resources:
            for col in ['ligand_subunits', 'receptor_subunits']:
                for value in r[col]:
                    needed.update(str(value).split('|'))
        keep = np.flatnonzero(a.var_names.isin(needed))
        parts, libraries = [], []
        for start in range(0, a.n_obs, chunk_size):
            block = sparse.csr_matrix(a.X[start:start + chunk_size])
            # Sum in float64, as in the full-matrix normalization path.
            libraries.append(np.asarray(block.sum(axis=1, dtype=np.float64)).ravel())
            parts.append(block[:, keep])
        return sparse.vstack(parts, format='csr'), a.obs.copy(), a.var_names[keep].copy(), np.concatenate(libraries)
    finally:
        a.file.close()
