import math

import numpy as np

from xmmcluster.initmodel import cv, data, features


def check_folds(groups, n_folds):
    tested = np.zeros(len(groups), int)
    for tr, te in cv.cluster_folds(groups, n_folds):
        assert not set(groups[tr]) & set(groups[te]), 'a cluster is in both train and test'
        assert len(tr) + len(te) == len(groups)
        tested[te] += 1
    assert (tested == 1).all(), 'every row is tested exactly once'


def test_folds_never_share_a_cluster():
    groups = np.array(['a', 'a', 'b', 'c', 'c', 'c', 'd', 'e', 'e'])
    check_folds(groups, 0)
    check_folds(groups, 2)
    assert len(list(cv.cluster_folds(groups, 0))) == 5


def test_folds_on_bundled_data():
    rows = data.load()
    groups = np.array([r['cluster'] for r in rows])
    check_folds(groups, 0)
    check_folds(groups, 10)


def synthetic(n_clusters=12, per=6, seed=1):
    rng = np.random.default_rng(seed)
    bands = ['0p5-0p8', '0p8-1p2', '1p2-2', '2-3p5', '3p5-6', '6-10']
    rows = []
    for c in range(n_clusters):
        kT = float(np.exp(rng.uniform(0, 2.3)))
        for _ in range(per):
            r = dict(cluster=f'c{c}', kT=kT * float(np.exp(rng.normal(0, 0.05))), e_kT=0.05 * kT, Z=0.3,
                     z=0.1, nh=0.02, band_lo=0.5, has_mos=1, has_pn=1, net_counts=1e4, area_arcmin2=20.0,
                     r_in=0.0, r_out=300.0, kT_catalog=kT * 1.1)
            for i, b in enumerate(bands):        # harder spectra for hotter clusters
                r[f'lf_{b}'] = -4 - 0.3 * i / math.log(1 + kT) + rng.normal(0, 0.01)
                r[f'lsky_{b}'] = r[f'lf_{b}'] - 1
                r[f'lnet_{b}'] = r[f'lf_{b}']
            rows.append(r)
    return rows


def test_evaluate_is_deterministic_and_learns():
    rows = synthetic()
    a = cv.evaluate(rows, n_trees=30)
    b = cv.evaluate(rows, n_trees=30)
    assert np.allclose([x['kT_pred'] for x in a], [x['kT_pred'] for x in b], rtol=1e-12, atol=0)
    s = cv.summary(a)
    assert s['n_clusters'] == 12 and s['n_catalog'] == len(rows)
    assert s['prediction']['median_abs_dlnkT'] < s['median_training_kT']['median_abs_dlnkT']


def test_feature_matrix_shape():
    rows = synthetic(2, 2)
    X = features.matrix(rows)
    assert X.shape == (4, len(features.feature_names(rows)))
    assert np.isfinite(X).all()
