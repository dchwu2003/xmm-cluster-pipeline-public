"""Initial-kT prediction from band fluxes, every cluster held out in turn.

    python -m demo.initmodel [--folds 0] [--trees 300] [--figure docs/figures/initmodel.png]

Trains extra trees on data/initmodel/spectra.csv with leave-one-cluster-out cross-validation
(``--folds 0``: one fold per cluster; ``--folds 10``: ten folds of whole clusters, faster), prints the
median error of the predicted starting kT against the fitted kT, next to the catalogue-kT baseline,
and writes the predicted-vs-fitted figure. Seeds are fixed, so the numbers repeat.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from xmmcluster.initmodel import cv, data

ROOT = Path(__file__).resolve().parents[1]


def figure(ev: list[dict], out: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullLocator
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    ax.scatter([e['kT'] for e in ev], [e['kT_pred'] for e in ev], s=5, alpha=0.35, color='#2a6fb0', linewidths=0)
    lo, hi = 0.3, 25
    ax.plot([lo, hi], [lo, hi], color='0.3', lw=0.8)
    for f in (0.8, 1.2):
        ax.plot([lo, hi], [lo * f, hi * f], color='0.6', lw=0.6, ls='--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10]))
        axis.set_major_formatter(FormatStrFormatter('%g'))
        axis.set_minor_locator(NullLocator())
    ax.set_xlabel('fitted kT (keV)')
    ax.set_ylabel('predicted starting kT (keV)')
    ax.set_title(f"Initial kT, cluster held out\n({len(ev)} spectra of {len({e['cluster'] for e in ev})} clusters; "
                 'dashed: ±20 %)', loc='left')
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def pct(s: dict) -> str:
    return (f"median error {s['median_abs_dlnkT']:.1%}, {s['within_10pct']:.0%} within 10 %, "
            f"{s['within_20pct']:.0%} within 20 %, {s['beyond_50pct']:.0%} beyond 50 %")


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--folds', type=int, default=0, help='0: leave-one-cluster-out (default)')
    ap.add_argument('--trees', type=int, default=cv.N_TREES)
    ap.add_argument('--figure', type=Path, default=ROOT / 'docs' / 'figures' / 'initmodel.png')
    a = ap.parse_args(argv)
    t0 = time.time()
    rows = data.load()
    ev = cv.evaluate(rows, n_folds=a.folds, n_trees=a.trees)
    s = cv.summary(ev)
    folds = a.folds or s['n_clusters']
    print(f"{s['n_spectra']} spectra of {s['n_clusters']} clusters, {folds} folds of whole clusters, {a.trees} trees")
    print(f"  prediction:          {pct(s['prediction'])}")
    print(f"  median training kT:  {pct(s['median_training_kT'])}")
    print(f"On the {s['n_catalog']} spectra ({s['n_catalog_clusters']} clusters) with a catalogue kT:")
    print(f"  prediction:          {pct(s['prediction_on_catalog'])}")
    print(f"  catalogue kT:        {pct(s['catalog'])}")
    print('Median error by fitted kT (prediction / catalogue):')
    for b in cv.by_kT_bin(ev):
        c = f"{b['catalog']:.1%}" if b['catalog'] is not None else '-'
        print(f"  {b['lo']:>2}-{b['hi']:<2} keV  {b['n']:5d} spectra {b['n_clusters']:4d} clusters  "
              f"{b['prediction']:.1%} / {c}")
    figure(ev, a.figure)
    print(f'figure: {a.figure.relative_to(ROOT) if a.figure.is_relative_to(ROOT) else a.figure}'
          f'   ({time.time() - t0:.0f} s)')
    s['folds'], s['trees'], s['bins'] = folds, a.trees, cv.by_kT_bin(ev)
    return s


if __name__ == '__main__':
    main()
