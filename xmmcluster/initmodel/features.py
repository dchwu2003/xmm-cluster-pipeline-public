"""Feature matrix of the initial model (no fit needed for any of these)."""
from __future__ import annotations

import math

import numpy as np

REF_BAND = '1p2-2'


def bands(rows: list[dict]) -> list[str]:
    return [k[3:] for k in rows[0] if k.startswith('lf_')]


def feature_names(rows: list[dict]) -> list[str]:
    b = bands(rows)
    names = ['z', 'nh', 'band_lo', 'has_mos', 'has_pn', 'log_counts', 'log_area', 'r_in_frac']
    names += [f'hr_{x}' for x in b] + [f'hrnet_{x}' for x in b] + [f'skyfrac_{x}' for x in b]
    return names + ['lsb']


def matrix(rows: list[dict]) -> np.ndarray:
    """Band ratios against 1.2-2 keV of the source (sky included) and of source - sky, sky / source per
    band, surface brightness, counts, area, z, N_H. A missing sky or band becomes the sentinel -9."""
    b = bands(rows)
    X = []
    for r in rows:
        v = [r['z'], r['nh'], r['band_lo'], r['has_mos'], r['has_pn'],
             math.log10(max(r['net_counts'], 1.0)), math.log10(max(r['area_arcmin2'], 1e-3)),
             r['r_in'] / r['r_out'] if r['r_out'] and r['r_out'] == r['r_out'] else math.nan]
        v += [r[f'lf_{x}'] - r[f'lf_{REF_BAND}'] for x in b]
        v += [r[f'lnet_{x}'] - r[f'lnet_{REF_BAND}'] for x in b]
        v += [r[f'lsky_{x}'] - r[f'lf_{x}'] for x in b]
        v += [r[f'lf_{REF_BAND}']]
        X.append(v)
    return np.nan_to_num(np.asarray(X, float), nan=-9.0)


def weights(rows: list[dict], floor: float = 0.05) -> np.ndarray:
    """1 / (sigma_lnkT^2 + floor^2): precise fits count more, but no single fit dominates."""
    rel = np.array([r['e_kT'] / r['kT'] for r in rows])
    return 1.0 / (rel ** 2 + floor ** 2)
