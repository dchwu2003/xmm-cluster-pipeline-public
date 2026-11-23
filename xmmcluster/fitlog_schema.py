"""Schema of the fit log: one JSON record per spectral fit, kept as training data.

Every fit the pipeline runs is stored with what was fitted, the complete initial model, every final
parameter, the model choices and the published value it is compared with. The published value
(``target``) is stored for the comparison only; no fit, rule or model input reads it.

``validate(record)`` returns a list of problems ('' -> valid). The fit log itself is not part of this
repository (it holds server paths); ``data/initmodel/spectra.csv`` is derived from it.
"""
from __future__ import annotations

import re

SCHEMA = 2          # 2 added model_choice; schema-1 records are still valid
SCHEMAS = (1, 2)

# field -> (type or tuple of types, required, meaning)
FIELDS = {
    'schema': ((int, str), True, 'schema version'),
    'record_id': (str, True, '16 hex digits: hash of paper, cluster, label, run, command and inputs'),
    'paper': (str, True, "the paper reproduced ('autofit' for the automatic fit)"),
    'cluster': (str, True, 'cluster or group name'),
    'label': (str, True, 'fit variant, e.g. r015-r500_1T_xmm3'),
    'obsids': (list, True, 'XMM ObsIDs (10 digits)'),
    'region': (dict, True, 'shape, centre, radii (arcsec and r500), r500 source'),
    'redshift': ((int, float), True, 'held fixed unless the rules free it'),
    'nh_1e22': ((int, float, type(None)), True, 'Galactic N_H (1e22 cm^-2), held fixed'),
    'spectra': (list, True, 'per spectrum: instrument, camera, ObsID, counts, area_arcmin2, exposure_s, qpb_fraction'),
    'xcm_text': (str, True, 'the complete initial XSPEC model file'),
    'initial_params': (list, True, 'parsed initial parameters with their physical role'),
    'model': (dict, True, 'model expressions (main and named models such as the soft protons) and generator'),
    'settings': (dict, False, 'XSPEC settings: statistic, method, abundance table, cross sections'),
    'ignores': (list, False, 'ignored channel ranges per data group'),
    'extra_commands': (list, False, 'XSPEC commands run after loading the model'),
    'model_choice': (dict, True, 'ICM model, number of temperatures, band, background, sky, N_H source, ...'),
    'converged': (bool, True, 'the fit converged'),
    'result': (dict, True, 'statistic, dof, bins and every final parameter'),
    'summary': (dict, True, 'kT1, kT1_sigma, abundance, ... of the ICM (kT, e_kT for single-region drivers)'),
    'target': (dict, False, 'published values compared with; never used to steer a fit'),
    'run_id': (str, True, 'run identifier (the logged command stream of that run)'),
    'command_id': (str, False, 'logged command that ran the fit'),
    'seconds': ((int, float), False, 'wall time of the fit'),
    'time_utc': (str, True, 'ISO 8601 time of the fit'),
    'notes': (str, False, 'free text'),
}
ROLES = {'icm', 'lhb', 'galactic_halo', 'cxb', 'absorber', 'scale', 'softp'}   # plus instr_<line>
RESULT_KEYS = ('dof', 'params')
_HEX16 = re.compile(r'^[0-9a-f]{16}$')
_OBSID = re.compile(r'^\d{10}$')


def validate(rec: dict) -> list[str]:
    """Problems of one record; an empty list means the record follows the schema."""
    out = []
    for name, (typ, required, _) in FIELDS.items():
        if name not in rec:
            if required:
                out.append(f'missing {name}')
            continue
        if not isinstance(rec[name], typ) or (typ is int and isinstance(rec[name], bool)):
            out.append(f'{name}: expected {typ}, got {type(rec[name]).__name__}')
    for name in rec:
        if name not in FIELDS and name != 'remote_dir':
            out.append(f'unknown field {name}')
    if out:
        return out
    if int(rec['schema']) not in SCHEMAS:
        out.append(f'schema {rec["schema"]} not in {SCHEMAS}')
    if not _HEX16.match(rec['record_id']):
        out.append('record_id: not 16 hex digits')
    out += [f'obsid {o!r}: not 10 digits' for o in rec['obsids'] if not _OBSID.match(str(o))]
    for k in RESULT_KEYS:
        if k not in rec['result']:
            out.append(f'result: missing {k}')
    if rec['converged'] and not {'kT', 'kT1'} & set(rec['summary']):
        out.append('summary: no ICM temperature (kT or kT1)')
    for p in rec['initial_params']:
        if p.get('role') and p['role'] not in ROLES and not p['role'].startswith('instr_'):
            out.append(f"initial_params: unknown role {p['role']!r}")
    for s in rec['spectra']:
        if not {'instrument', 'counts'} <= set(s):
            out.append('spectra: entry without instrument or counts')
    return out


def training_label(rec: dict) -> dict | None:
    """The (kT, sigma, Z) label of a record for the initial model, or None if the fit is not usable:
    converged, 0.2 < kT < 25 keV, kT / sigma >= 3; Z only when it was free."""
    s = rec.get('summary') or {}
    try:
        kT = float(s.get('kT1', s.get('kT')))
        e = float(s.get('kT1_sigma', s.get('e_kT')))
    except (TypeError, ValueError):
        return None
    if not rec.get('converged') or not 0.2 < kT < 25 or not e > 0 or kT / e < 3:
        return None
    Z = s.get('abundance', s.get('Z'))
    return dict(kT=kT, e_kT=e, Z=None if s.get('abundance_fixed') else Z)
