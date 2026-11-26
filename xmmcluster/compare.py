"""Per-paper comparison of reproduced values with the published ones.

Input: ``data/summary/reproduction_points.csv``, one row per compared value (paper, cluster, region,
quantity, ours, published, ratio = ours / published, and the verdict of the automatic data-quality
checks). Two kinds of values are kept in the table but left out of every statistic:

* clusters whose data the automatic checks rejected (``autofit_verdict == 'rejected'``: no trusted
  fit, e.g. a flared observation). The verdict comes from data-quality rules (``xmmcluster.rules``)
  and never sees the published value;
* temperatures at the upper limit of the fit (within 2 % of 30 keV): not measurements.

Both rules are applied whatever the ratio is, so a bad ratio cannot be removed by hand.
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

from xmmcluster.rules import KT_LIMITS

ROOT = Path(__file__).resolve().parents[1]
POINTS = ROOT / 'data' / 'summary' / 'reproduction_points.csv'
KT_MAX = KT_LIMITS[1]

# paper key -> (citation, telescope, what is compared)
PAPERS = {
    'zhang2008': ('Zhang et al. 2008 (LoCuSS)', 'XMM', 'kT in 0.2-0.5 r500'),
    'planck2011': ('Planck Collaboration 2011 XI', 'XMM', 'T_X in 0.15-0.75 R500'),
    'planckval': ('Planck Collaboration 2011 IX, 2012 I, 2013 IV', 'XMM', 'T_X in 0.15-0.75 R500'),
    'lovisari2015': ('Lovisari et al. 2015', 'XMM', 'group kT'),
    'lovisari2020': ('Lovisari et al. 2020', 'XMM', 'kT in 0-1 and 0.15-1 R500'),
    'rossetti2024': ('Rossetti et al. 2024 (CHEX-MATE)', 'XMM', 'kT profiles (annuli)'),
    'johnson2009': ('Johnson et al. 2009 (2dXGS)', 'XMM', 'group kT in 0.1-0.3 r500'),
    'pratt2009': ('Pratt et al. 2009 (REXCESS)', 'XMM', 'kT in 0.15-1 R500'),
    'osullivan2017': ("O'Sullivan et al. 2017 (CLoGS)", 'XMM', 'group system kT'),
    'xgap2025': ('Eckert et al. 2025 (X-GAP)', 'XMM', 'group kT in 0-300 kpc'),
    'schellenberger2015': ('Schellenberger et al. 2015 (HIFLUGCS)', 'XMM', 'EPIC kT in 0.7-7 keV, r < 3.5 arcmin'),
    'bulbul2019': ('Bulbul et al. 2019 (SPT-SZ)', 'XMM', 'MOS kT in < R500 and 0.15-1 R500'),
}


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def left_out(p: dict) -> str:
    """Why a value is left out of the statistics ('' if it is kept). Never looks at the published value."""
    if p.get('autofit_verdict') == 'rejected':
        return 'autofit rejected'
    if p['quantity'] == 'kT' and _f(p['ours']) is not None and _f(p['ours']) >= 0.98 * KT_MAX:
        return 'kT at fit limit'
    return ''


def load(path: Path = POINTS) -> list[dict]:
    out = []
    for r in csv.DictReader(open(path)):
        r = dict(r)
        for k in ('ours', 'e_ours', 'published', 'e_published', 'ratio'):
            r[k] = _f(r[k])
        r['left_out'] = left_out(r)
        out.append(r)
    return out


def quantile(xs: list[float], q: float) -> float:
    xs = sorted(xs)
    i = q * (len(xs) - 1)
    lo = int(i)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)


def within(ratios: list[float], tol: float) -> float:
    return sum(abs(x - 1) <= tol + 1e-12 for x in ratios) / len(ratios)


def per_paper(pts: list[dict]) -> list[dict]:
    out = []
    for paper in sorted({p['paper'] for p in pts}, key=lambda k: list(PAPERS).index(k) if k in PAPERS else 99):
        cite, tel, what = PAPERS.get(paper, (paper, 'XMM', ''))
        allp = [p for p in pts if p['paper'] == paper]
        ps = [p for p in allp if not p['left_out']]
        if not ps:
            continue
        r = [p['ratio'] for p in ps]
        out.append(dict(paper=paper, citation=cite, telescope=tel, compared=what, quantity=ps[0]['quantity'],
                        n_clusters=len({p['cluster'] for p in ps}), n_values=len(ps),
                        n_left_out=len(allp) - len(ps),
                        median_ratio=round(statistics.median(r), 3), p16=round(quantile(r, 0.16), 3),
                        p84=round(quantile(r, 0.84), 3),
                        within_10pct=round(within(r, 0.10), 2), within_20pct=round(within(r, 0.20), 2)))
    return out


def headline(pts: list[dict], quantity: str = 'kT') -> dict:
    """All kept values of one quantity together."""
    ps = [p for p in pts if p['quantity'] == quantity and not p['left_out']]
    r = [p['ratio'] for p in ps]
    return dict(n_papers=len({p['paper'] for p in ps}), n_values=len(r),
                n_pairs=len({(p['paper'], p['cluster']) for p in ps}), n_clusters=len({p['cluster'] for p in ps}),
                n_left_out=sum(1 for p in pts if p['quantity'] == quantity and p['left_out']),
                median_ratio=round(statistics.median(r), 3), within_10pct=round(within(r, 0.10), 3),
                within_20pct=round(within(r, 0.20), 3))


def write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def figure(pts: list[dict], papers: list[dict], out: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    order = sorted(papers, key=lambda p: p['median_ratio'])
    fig, ax = plt.subplots(figsize=(7.5, 0.32 * len(order) + 1.3))
    for i, p in enumerate(order):
        r = [x['ratio'] for x in pts if x['paper'] == p['paper'] and not x['left_out']]
        ax.scatter(r, [i] * len(r), s=8, alpha=0.25 if len(r) > 30 else 0.6, color='#2a6fb0', linewidths=0)
        ax.plot([p['p16'], p['p84']], [i, i], color='#2a6fb0', lw=2)
        ax.plot(p['median_ratio'], i, 'o', color='white', mec='#2a6fb0', ms=6, mew=1.6)
    ax.axvline(1, color='0.3', lw=0.8)
    ax.axvspan(0.9, 1.1, color='0.92', zorder=0)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{p['citation']}  (n={p['n_values']})" for p in order])
    ax.set_xlim(0.5, 1.6)
    ax.set_xlabel('ours / published kT\n(circle: median, bar: 16-84 %, shaded: ±10 %)')
    ax.set_title('Reproduced temperatures per paper', loc='left')
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
