"""Write every result number of README.md from the bundled CSVs (nothing in those blocks is typed by hand).

    python -m demo.readme_numbers           # rewrites the blocks between <!-- numbers:NAME --> markers
    python -m demo.readme_numbers --check   # exit 1 if README.md is out of date (used by the tests)

Runs both demos (the initial model with its full leave-one-cluster-out setting) and refreshes the
figures in docs/figures.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from demo import initmodel, reproductions
from xmmcluster import compare

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / 'README.md'
SNAPSHOT = ROOT / 'data' / 'SNAPSHOT'


def pc(x: float) -> str:
    return f'{x:.0%}'


def blocks() -> dict[str, str]:
    h = reproductions.main([])
    papers = compare.per_paper(compare.load())
    s = initmodel.main([])
    p, pcat, cat, const = s['prediction'], s['prediction_on_catalog'], s['catalog'], s['median_training_kT']
    snap = SNAPSHOT.read_text().strip()
    out = {}
    out['snapshot'] = f'Snapshot of the data: {snap}.'
    out['headline'] = (
        f"- **{h['n_papers']} papers reproduced**: {h['n_values']:,} temperatures of {h['n_pairs']} paper-cluster pairs "
        f"({h['n_clusters']} clusters). Median ratio ours / published **{h['median_ratio']:.2f}**; "
        f"**{pc(h['within_10pct'])}** within 10 % and **{pc(h['within_20pct'])}** within 20 %. "
        f"{h['n_left_out']} further values are left out by the data-quality rules (rejected data, or kT at the fit limit), "
        f"whatever their ratio.\n"
        f"- **Initial-model prototype**: from band fluxes measured before any fit, extra trees predict the starting "
        f"kT of a cluster they never saw ({s['n_spectra']:,} spectra of {s['n_clusters']} clusters, leave-one-cluster-out) "
        f"with a median error of **{p['median_abs_dlnkT']:.1%}**; {pc(p['within_20pct'])} are within 20 %. "
        f"On the {s['n_catalog']:,} spectra ({s['n_catalog_clusters']} clusters) that have a catalogue temperature, "
        f"the prediction's median error is {pcat['median_abs_dlnkT']:.1%} against {cat['median_abs_dlnkT']:.1%} "
        f"for the catalogue kT; the median kT of the training clusters is off by {const['median_abs_dlnkT']:.1%}.")
    rows = ['| paper | compared | clusters | values | median ratio | 16-84 % | within 10 % | within 20 % |',
            '|---|---|---|---|---|---|---|---|']
    for q in papers:
        rows.append(f"| {q['citation']} | {q['compared']} | {q['n_clusters']} | {q['n_values']} | {q['median_ratio']:.2f} "
                    f"| {q['p16']:.2f}-{q['p84']:.2f} | {pc(q['within_10pct'])} | {pc(q['within_20pct'])} |")
    out['papers'] = '\n'.join(rows)
    rows = ['| fitted kT | spectra | clusters | prediction | catalogue kT |', '|---|---|---|---|---|']
    for b in s['bins']:
        hi = f"{b['hi']} keV" if b['hi'] < 99 else 'keV and above'
        c = f"{b['catalog']:.1%}" if b['catalog'] is not None else '-'
        rows.append(f"| {b['lo']}{'-' if b['hi'] < 99 else ' '}{hi} | {b['n']:,} | {b['n_clusters']} | {b['prediction']:.1%} | {c} |")
    out['bins'] = '\n'.join(rows)
    return out


def render(text: str, new: dict[str, str]) -> str:
    for name, body in new.items():
        # the block body never crosses another marker, so an empty block cannot swallow the next one
        pat = re.compile(rf'(<!-- numbers:{name} -->\n)(?:(?!<!-- numbers:|<!-- /numbers -->).)*(<!-- /numbers -->)', re.S)
        if len(pat.findall(text)) != 1:
            raise SystemExit(f'README.md must have exactly one block {name}')
        text = pat.sub(lambda m: m.group(1) + body + '\n' + m.group(2), text)
    return text


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args(argv)
    old = README.read_text()
    new = render(old, blocks())
    if a.check:
        if new != old:
            print('README.md numbers are out of date: run python -m demo.readme_numbers', file=sys.stderr)
            return 1
        print('README.md numbers are up to date')
        return 0
    README.write_text(new)
    print('README.md updated')
    return 0


if __name__ == '__main__':
    sys.exit(main())
