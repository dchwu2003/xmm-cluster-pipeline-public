"""Does the starting temperature matter? Summary of replayed fits.

A replay refits a logged 1T fit unchanged except for its starting ICM kT (the logged value, the
model's prediction, the catalogue kT, 1 keV, 10 keV). The replays run on the analysis server; this
module only summarises their output rows (record_id, start, stat, kT, seconds). The replay data are
not bundled: they were run on clusters outside this repository's data scope.
"""
from __future__ import annotations

import statistics


def summarize(rows: list[dict]) -> list[dict]:
    """Per start value: how often the fit ends at the best minimum found from any start (statistic
    within 1 of the lowest), how often kT ends > 5 % from the best-minimum kT, and the wall time."""
    rows = [r for r in rows if r.get('stat') not in ('', None, 'None')]
    by: dict[str, list[dict]] = {}
    for r in rows:
        by.setdefault(r['record_id'], []).append(r)
    out = []
    for st in sorted({r['start'] for r in rows}):
        hit = miss_kT = n = 0
        secs = []
        for rs in by.values():
            mine = [r for r in rs if r['start'] == st]
            if not mine or len(rs) < 2:
                continue
            best = min(rs, key=lambda r: float(r['stat']))
            m = mine[0]
            n += 1
            hit += float(m['stat']) - float(best['stat']) < 1.0
            miss_kT += abs(float(m['kT']) / float(best['kT']) - 1) > 0.05
            secs.append(float(m['seconds']))
        out.append(dict(start=st, n=n, at_best_minimum=round(hit / n, 3) if n else None,
                        kT_off_5pct=round(miss_kT / n, 3) if n else None,
                        median_seconds=round(statistics.median(secs), 1) if secs else None))
    return out
