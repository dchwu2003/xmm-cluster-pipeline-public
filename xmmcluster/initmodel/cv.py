"""Leave-one-cluster-out evaluation of the initial-kT model.

All fits of one cluster are held out together, and every prediction comes from a model that never saw
that cluster, so the score measures what happens for a new cluster. Baselines: the catalogue kT (one
value per cluster, the same for every region) and the median kT of the training clusters.
"""
from __future__ import annotations

import math

import numpy as np

from xmmcluster.initmodel.features import matrix, weights

N_TREES = 300
SEED = 0


def model(n_trees: int = N_TREES, seed: int = SEED):
    """Extremely randomised trees: robust, no tuning, fast."""
    from sklearn.ensemble import ExtraTreesRegressor
    return ExtraTreesRegressor(n_estimators=n_trees, min_samples_leaf=3, max_features=0.5, n_jobs=-1,
                               random_state=seed)


def cluster_folds(groups: np.ndarray, n_folds: int = 0):
    """(train, test) index arrays; each cluster is in exactly one test fold and never in its training set.
    ``n_folds`` 0 means one fold per cluster (leave-one-cluster-out)."""
    clusters = sorted(set(groups))
    n = n_folds or len(clusters)
    fold_of = {c: i % n for i, c in enumerate(clusters)}
    f = np.array([fold_of[g] for g in groups])
    for k in range(n):
        yield np.where(f != k)[0], np.where(f == k)[0]


def evaluate(rows: list[dict], n_folds: int = 0, n_trees: int = N_TREES, seed: int = SEED,
             with_abundance: bool = False) -> list[dict]:
    X = matrix(rows)
    y = np.log(np.array([r['kT'] for r in rows]))
    w = weights(rows)
    zl = np.array([r['Z'] for r in rows])
    groups = np.array([r['cluster'] for r in rows])
    pred, const, predZ = (np.full(len(rows), np.nan) for _ in range(3))
    for tr, te in cluster_folds(groups, n_folds):
        pred[te] = model(n_trees, seed).fit(X[tr], y[tr], sample_weight=w[tr]).predict(X[te])
        const[te] = np.median(y[tr])
        if with_abundance:
            ok = tr[np.isfinite(zl[tr])]
            predZ[te] = model(n_trees, seed).fit(X[ok], zl[ok], sample_weight=w[ok]).predict(X[te])
    out = []
    for i, r in enumerate(rows):
        out.append(dict(cluster=r['cluster'], kT=r['kT'], kT_pred=float(np.exp(pred[i])),
                        kT_const=float(np.exp(const[i])),
                        kT_catalog=r['kT_catalog'] if r['kT_catalog'] == r['kT_catalog'] else None,
                        Z=r['Z'] if r['Z'] == r['Z'] else None,
                        Z_pred=float(predZ[i]) if with_abundance else None))
    return out


def scores(ev: list[dict], col: str) -> dict:
    """Median |ln(kT_start / kT_fit)| and the fractions within 10 % / 20 % and beyond 50 %."""
    d = np.array([abs(math.log(e[col] / e['kT'])) for e in ev if e[col]])
    if not len(d):
        return {}
    return dict(n=int(len(d)), median_abs_dlnkT=round(float(np.median(d)), 3),
                within_10pct=round(float(np.mean(d < math.log(1.1))), 3),
                within_20pct=round(float(np.mean(d < math.log(1.2))), 3),
                beyond_50pct=round(float(np.mean(d > math.log(1.5))), 3))


def summary(ev: list[dict]) -> dict:
    """Prediction and baselines on all spectra, and prediction vs catalogue on the spectra that have one."""
    with_cat = [e for e in ev if e['kT_catalog']]
    return dict(n_spectra=len(ev), n_clusters=len({e['cluster'] for e in ev}),
                prediction=scores(ev, 'kT_pred'), median_training_kT=scores(ev, 'kT_const'),
                n_catalog=len(with_cat), n_catalog_clusters=len({e['cluster'] for e in with_cat}),
                prediction_on_catalog=scores(with_cat, 'kT_pred'), catalog=scores(with_cat, 'kT_catalog'))


def by_kT_bin(ev: list[dict], edges=(0, 1, 2, 4, 8, 12, 99)) -> list[dict]:
    out = []
    for lo, hi in zip(edges, edges[1:]):
        sub = [e for e in ev if lo <= e['kT'] < hi]
        if sub:
            out.append(dict(lo=lo, hi=hi, n=len(sub), n_clusters=len({e['cluster'] for e in sub}),
                            prediction=scores(sub, 'kT_pred').get('median_abs_dlnkT'),
                            catalog=scores(sub, 'kT_catalog').get('median_abs_dlnkT')))
    return out
