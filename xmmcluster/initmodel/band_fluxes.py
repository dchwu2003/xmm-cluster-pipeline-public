"""Band photon fluxes of the source and sky spectra of one logged fit (reference implementation).

This is how the ``lf_`` / ``lsky_`` / ``lnet_`` columns of ``data/initmodel/spectra.csv`` were
measured. It reads the XMM spectra, backgrounds and responses of a fit record, so it runs where the
spectra are (needs ``astropy``); the demos do not use it.
For each fit record the XMM spectra are split into source regions and the sky annulus
(``spectra[i]['region'] == 'outer'``). For each group of spectra and each band the
QPB-subtracted counts are converted to a photon flux per arcmin^2,

    F_band = sum(net counts) / sum(exposure * <ARF>_band * area_arcmin2) / band width,

with <ARF>_band the mean effective area over the band. MOS and pn are combined in this way,
so the features do not depend on which cameras survived the screening. The features are
what an initial-model generator can measure before any fit:

    lf_<band>     log10 source flux per arcmin^2
    lsky_<band>   log10 sky-annulus flux per arcmin^2 (the same cameras)
    lnet_<band>   log10 (source - sky) flux per arcmin^2 (a crude double subtraction)
    net_counts, area_arcmin2, has_mos, has_pn, qpb_fraction

Spectra are cached by path, so shared spectra are read once.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

BANDS = ((0.5, 0.8), (0.8, 1.2), (1.2, 2.0), (2.0, 3.5), (3.5, 6.0), (6.0, 10.0))
BAND_NAMES = [f'{lo:g}-{hi:g}'.replace('.', 'p') for lo, hi in BANDS]


def band_key(i: int) -> str:
    return BAND_NAMES[i]


def data_groups(xcm_text: str) -> list[dict]:
    """XMM data groups of an .xcm in order (RASS left out): data, background, arf paths."""
    groups, bkg, arf = [], {}, {}
    for line in xcm_text.splitlines():
        w = line.split()
        if len(w) >= 3 and w[0] == 'data':
            groups.append(dict(index=int(w[1].split(':')[1]), data=w[2]))
        elif len(w) >= 3 and w[0] == 'backgrnd':
            bkg[int(w[1])] = w[2]
        elif len(w) >= 3 and w[0] == 'arf' and w[1].startswith('1:'):
            arf[int(w[1].split(':')[1])] = w[2]
    out = []
    for g in groups:
        if g['data'].endswith('rass.pi'):
            continue
        g['bkg'] = bkg.get(g['index'])
        g['arf'] = arf.get(g['index'])
        out.append(g)
    return out


_cache: dict = {}


def spectrum_bands(path: Path, bkg: Path | None, arf: Path | None, rmf: Path | None) -> dict | None:
    """Net counts and exposure * <ARF> per band for one spectrum; None if a file is missing."""
    key = (str(path), str(bkg), str(arf))
    if key in _cache:
        return _cache[key]
    from astropy.io import fits
    import numpy as np
    try:
        with fits.open(path) as h:
            sp = h['SPECTRUM']
            counts = np.asarray(sp.data['COUNTS'], float)
            exp = float(sp.header['EXPOSURE'])
            bs = float(sp.header.get('BACKSCAL', 1.0))
            resp = sp.header.get('RESPFILE')
        if rmf is None and resp and resp.upper() != 'NONE':
            rmf = path.parent / resp
        if rmf is None:
            rmf = path.with_name(path.name.split('-')[0] + '.rmf')
        with fits.open(rmf) as h:
            eb = h['EBOUNDS'].data
            emid = 0.5 * (np.asarray(eb['E_MIN'], float) + np.asarray(eb['E_MAX'], float))
        net = counts.copy()
        if bkg is not None and bkg.exists():
            with fits.open(bkg) as h:
                b = h[1]
                bc = np.asarray(b.data['COUNTS'] if 'COUNTS' in b.columns.names
                                else b.data['RATE'] * float(b.header['EXPOSURE']), float)
                scale = exp / float(b.header['EXPOSURE']) * bs / float(b.header.get('BACKSCAL', 1.0))
            net = counts - bc * scale
        with fits.open(arf) as h:
            a = h['SPECRESP'].data
            ae = 0.5 * (np.asarray(a['ENERG_LO'], float) + np.asarray(a['ENERG_HI'], float))
            area = np.asarray(a['SPECRESP'], float)
    except (OSError, KeyError, TypeError, ValueError) as e:
        print(f'skip {path}: {e}', file=sys.stderr)
        _cache[key] = None
        return None
    out = dict(exposure=exp, total=float(counts.sum()))
    for i, (lo, hi) in enumerate(BANDS):
        m = (emid >= lo) & (emid < hi)
        ma = (ae >= lo) & (ae < hi)
        out[f'net_{i}'] = float(net[m].sum())
        out[f'expa_{i}'] = exp * float(area[ma].mean())
    _cache[key] = out
    return out


def record_features(rec: dict) -> dict | None:
    base = Path(rec['remote_dir'])
    groups = data_groups(rec.get('xcm_text', ''))
    spectra = rec.get('spectra') or []
    if not groups or len(groups) != len(spectra):
        return None
    acc = {r: {f'{k}_{i}': 0.0 for i in range(len(BANDS)) for k in ('net', 'expa')} for r in ('src', 'sky')}
    total = dict(src=0.0, sky=0.0)
    area = dict(src=0.0, sky=0.0)
    inst = set()
    qpb = []
    for g, s in zip(groups, spectra):
        reg = 'sky' if s.get('region') == 'outer' else 'src'
        p = (base / g['data']).resolve()
        b = (base / g['bkg']).resolve() if g.get('bkg') else None
        a = (base / g['arf']).resolve() if g.get('arf') else p.with_name(p.name.split('-')[0] + '.arf')
        sb = spectrum_bands(p, b, a, None)
        if sb is None:
            return None
        ar = float(s.get('area_arcmin2') or 0) or 1.0
        for i in range(len(BANDS)):
            acc[reg][f'net_{i}'] += sb[f'net_{i}']
            acc[reg][f'expa_{i}'] += sb[f'expa_{i}'] * ar
        if reg == 'src':
            area['src'] += ar
            total['src'] += sb['total']
            inst.add('pn' if s.get('instrument') == 'pn' else 'mos')
            if s.get('qpb_fraction') is not None:
                qpb.append(float(s['qpb_fraction']))
        else:
            total['sky'] += sb['total']
    if total['src'] <= 0:
        return None
    row = dict(record_id=rec['record_id'], net_counts=round(sum(acc['src'][f'net_{i}'] for i in range(len(BANDS))), 1),
               area_arcmin2=round(area['src'] / max(1, len([s for s in spectra if s.get('region') != 'outer'])), 3),
               has_mos=int('mos' in inst), has_pn=int('pn' in inst),
               qpb_fraction=round(sum(qpb) / len(qpb), 4) if qpb else '')
    for i, (lo, hi) in enumerate(BANDS):
        w = hi - lo
        f = {}
        for reg in ('src', 'sky'):
            ea = acc[reg][f'expa_{i}']
            f[reg] = acc[reg][f'net_{i}'] / ea / w if ea > 0 else float('nan')
        k = band_key(i)
        row[f'lf_{k}'] = _log(f['src'])
        row[f'lsky_{k}'] = _log(f['sky'])
        row[f'lnet_{k}'] = _log(f['src'] - f['sky']) if not math.isnan(f['sky']) else row[f'lf_{k}']
    return row


def _log(x: float):
    return round(math.log10(x), 5) if x == x and x > 0 else ''
