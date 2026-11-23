"""Reproduced temperatures against the published ones, per paper.

    python -m demo.reproductions [--figure docs/figures/reproductions.png]

Reads data/summary/reproduction_points.csv, prints one line per paper and the total, writes
data/summary/reproduction_papers.csv and the per-paper figure.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from xmmcluster import compare

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--figure', type=Path, default=ROOT / 'docs' / 'figures' / 'reproductions.png')
    a = ap.parse_args(argv)
    pts = compare.load()
    papers = compare.per_paper(pts)
    compare.write(ROOT / 'data' / 'summary' / 'reproduction_papers.csv', papers)
    print(f"{'paper':44s} {'clusters':>8s} {'values':>6s} {'median':>7s} {'16-84 %':>12s} {'<=10 %':>7s} {'<=20 %':>7s}")
    for p in papers:
        print(f"{p['citation']:44s} {p['n_clusters']:8d} {p['n_values']:6d} {p['median_ratio']:7.3f} "
              f"{p['p16']:5.2f}-{p['p84']:<5.2f}  {p['within_10pct']:6.0%} {p['within_20pct']:7.0%}")
    h = compare.headline(pts)
    print(f"\nAll temperatures: {h['n_papers']} papers, {h['n_values']} values of {h['n_pairs']} paper-cluster pairs "
          f"({h['n_clusters']} clusters); median ours/published {h['median_ratio']:.3f}, "
          f"{h['within_10pct']:.0%} within 10 %, {h['within_20pct']:.0%} within 20 %. "
          f"Left out by the data-quality rules: {h['n_left_out']} values.")
    compare.figure(pts, papers, a.figure)
    print(f'figure: {a.figure.relative_to(ROOT) if a.figure.is_relative_to(ROOT) else a.figure}')
    return h


if __name__ == '__main__':
    main()
