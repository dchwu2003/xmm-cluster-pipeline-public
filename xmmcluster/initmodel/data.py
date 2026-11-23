"""The bundled training table: one row per fitted spectrum.

``data/initmodel/spectra.csv`` columns:

* ``cluster``, ``paper``, ``region``: which fit (refits of the same spectrum, e.g. with and without a
  RASS sky spectrum, are merged into one row with their geometric-mean kT; ``n_fits`` counts them);
* labels ``kT`` (keV), ``e_kT``, ``Z`` (solar, empty when it was fixed): the fitted 1T values;
* inputs measured before any fit: ``z``, ``nh`` (1e22 cm^-2), ``band_lo`` (lowest band edge, keV),
  ``has_mos``, ``has_pn``, ``net_counts``, ``area_arcmin2``, ``r_in`` / ``r_out`` (arcsec), and per
  band (0.5-0.8, 0.8-1.2, 1.2-2, 2-3.5, 3.5-6, 6-10 keV) the log10 photon flux per arcmin^2 of the
  source (``lf_``), of the sky annulus (``lsky_``) and of their difference (``lnet_``);
* baseline ``kT_catalog``: kT from the MCXC M500 (Piffaretti et al. 2011) through the M-T relation of
  Arnaud et al. 2005, for clusters in MCXC (``catalog`` names the entry); empty otherwise.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPECTRA = ROOT / 'data' / 'initmodel' / 'spectra.csv'
TEXT = ('cluster', 'paper', 'region', 'catalog')


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def load(path: Path = SPECTRA) -> list[dict]:
    rows = []
    for r in csv.DictReader(open(path)):
        rows.append({k: (v if k in TEXT else _num(v)) for k, v in r.items()})
    return rows
