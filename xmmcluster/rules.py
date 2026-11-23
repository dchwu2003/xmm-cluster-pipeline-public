"""Decision rules of the automatic fit, each taken from a reproduction.

Pure functions without I/O, so every rule is unit tested offline. The numbers are on the pipeline's
scale (QPB = filter-wheel-closed spectra scaled by corner rates; sky and soft protons fitted jointly
with an outer annulus). None of these rules sees a published value.
"""
from __future__ import annotations

import math

# --- scaling relations (initial region before the first fit) --------------------------------
OM = 0.3            # flat LCDM, H0 = 70 (as in the reproductions)
# Arnaud, Pointecouteau & Pratt 2005 (A&A 441, 893), whole sample, h70:
#   h(z) R500 = 1104 kpc (kT / 5 keV)^0.57,   h(z) M500 = 3.84e14 Msun (kT / 5 keV)^1.71
R500_5KEV_KPC, R500_SLOPE = 1104.0, 0.57
M500_5KEV, M500_SLOPE = 3.84e14, 1.71
R500_ITER_TOL = 0.10    # refit when the r500 implied by the fitted kT moves by more than this
R500_ITER_MAX = 2

# --- spectral bands (Pratt 2002, Maughan 2008, Leccardi 2008a) ------------------------------
HOT_KT = 3.0                          # above: hot-cluster bands
BAND_HOT = {'mos': 0.7, 'pn': 1.0}    # MOS kT stable above 0.6-0.8 keV; pn below 1 keV biases kT low
BAND_COOL = {'mos': 0.5, 'pn': 0.5}   # groups: Fe-L needs the soft band (Lovisari 2015 0.5 keV)
BAND_HI = {'mos': 10.0, 'pn': 10.0}
BAND_HI_SP = {'mos': 12.0, 'pn': 14.0}  # soft-proton excess: let the band it dominates pin the index
RASS_BAND = (0.1, 2.0)

# --- screening (Leccardi 2008a, Martino 2014, Zhang 2008) -----------------------------------
MIN_EXPOSURE_S = 3000.0     # A1413 pn: 152 s after filtering
RSB_DROP = 2.0              # R_SB (6-12 keV, FWC-normalised) above this: drop the exposure (paper scale 1.4)
RSB_FLAG = 1.66             # ... above this: flag (paper scale 1.15)
HARD_RATIO_HIBAND = 1.15    # observed / scaled QPB at 10-12 (MOS), 10-14 keV (pn): fit to 12/14 keV
NH_MISMATCH = 2.0           # an input N_H off the 21 cm value by more than x2 is a typo (seen for A1763, A1914)

# --- regions (Martino 2014, Lovisari 2015, Zhang 2008) ---------------------------------------
POSITION_MISMATCH_ARCMIN = 1.0  # catalogue vs literature X-ray centroid: flag
CENTROID_RADIUS_ARCMIN = 4.0  # 0.5-2 keV MOS centroid within 4' (Leccardi 2008a stand-in for large-radius isophotes)
CENTROID_MAX_ARCSEC = 60.0    # larger shifts are not trusted: keep the catalogue position and flag
                              # (A267: flared MOS pulled the 4' centroid 226" off the cluster)
CORE_R500 = 0.15            # core excision radius; the whole-region fit is kept as a bracket
FOV_ARCSEC = 600.0          # main region capped at 10 arcmin (an outer annulus has to fit inside the FOV)
BKG_ARCSEC = (720.0, 870.0)  # outer annulus 12-14.5 arcmin (own ICM; constrains sky and soft protons)
BKG_GAP_ARCSEC = 60.0
REDSHIFT_FREE_COUNTS = 50000  # main-region counts above which z is free (Leccardi 2008b: A1763 7% off);
                              # hot clusters only (Fe-K): a 1 keV group ran to the -7% limit (NGC 4325)
REDSHIFT_FREE_FRAC = 0.07

# --- quality flags after the fit (Zhang 2008, Maughan 2008, Leccardi 2008a/b, Martino 2014) -----
REDCHI_MULTI_T = 1.4        # mergers A781, A2259, A2409: a 1T fit of a multi-T region
REDCHI_BAD = 2.0            # flared A781 80 ks: biased kT with a tiny error
PN_MOS_RANGE = (0.8, 1.2)
SP_INDEX_LO = 0.1           # lower limit of the soft-proton index in fitgen.lovisari
KT_LIMITS = (0.1, 30.0)
Z_PEGGED = 0.005
Z_HIGH = 1.5                # solar: no cluster region is this metal rich in a 1T fit
MAX_REL_ERR_ITER = 0.2      # r500 is only moved by fits with kT known to better than 20%

# --- AtomDB (AtomDB version test) -------------------------------------------------------------
ATOMDB_SENSITIVE_KT = 2.0
# fitted kT with AtomDB 1.3.1 / with 3.0.9 (same spectra): ~1 keV 0.92, 1.7 keV 0.96, >= 2.3 keV 1
ATOMDB131_SCALE = ((1.0, 0.92), (1.7, 0.96), (2.3, 1.0))


def ez(z: float, om: float = OM) -> float:
    return math.sqrt(om * (1 + z) ** 3 + 1 - om)


def r500_kpc_from_kT(kT: float, z: float) -> float:
    return R500_5KEV_KPC * (kT / 5.0) ** R500_SLOPE / ez(z)


def kT_from_m500(m500_msun: float, z: float) -> float:
    return 5.0 * (m500_msun * ez(z) / M500_5KEV) ** (1 / M500_SLOPE)


def bands(kT: float, hiband: bool = False) -> dict:
    lo = BAND_HOT if kT >= HOT_KT else BAND_COOL
    hi = BAND_HI_SP if hiband else BAND_HI
    b = {cam: (lo[k], hi[k]) for cam, k in (('mos1', 'mos'), ('mos2', 'mos'), ('pn', 'pn'))}
    b['rass'] = RASS_BAND
    return b


def regions(r500_arcsec: float) -> tuple[list[dict], tuple[float, float], list[str]]:
    """Main regions (core-excised primary + whole bracket), the outer annulus and flags."""
    flags = []
    outer = r500_arcsec
    if outer > FOV_ARCSEC:
        outer = FOV_ARCSEC
        flags.append('region_truncated_at_10arcmin')
        flags.append('outer_annulus_contains_icm')  # nearby group: the 12-14.5' annulus is not source-free
    main = [dict(tag='r015-r500', inner=CORE_R500 * r500_arcsec, outer=outer, primary=True),
            dict(tag='r0-r500', inner=0.0, outer=outer, primary=False)]
    bkg = (max(BKG_ARCSEC[0], outer + BKG_GAP_ARCSEC), BKG_ARCSEC[1])
    if bkg[1] - bkg[0] < BKG_GAP_ARCSEC:
        flags.append('outer_annulus_thin')
    return main, bkg, flags


def screen_spectrum(summary: dict, rsb: float | None = None, hard_ratio: float | None = None) -> tuple[bool, list[str]]:
    """Keep one camera/exposure? ``summary`` is the native_region summary entry."""
    if 'error' in summary or not summary.get('counts'):
        return False, ['no_products']
    if (summary.get('exposure_s') or 0) < MIN_EXPOSURE_S:
        return False, ['short_exposure']
    if summary.get('qpb_counts_in_region') is None:
        return False, ['no_corner_data']
    if rsb is not None and rsb > RSB_DROP:
        return False, ['soft_protons_rsb_drop']
    flags = []
    if rsb is not None and rsb > RSB_FLAG:
        flags.append('soft_protons_rsb_flag')
    if hard_ratio is not None and hard_ratio > HARD_RATIO_HIBAND:
        flags.append('hard_excess')
    return True, flags


def use_hiband(hard_ratios: list[float]) -> bool:
    """Fit to 12/14 keV when the median hard-band ratio of the kept spectra is high."""
    r = sorted(v for v in hard_ratios if v is not None)
    return bool(r) and r[len(r) // 2] > HARD_RATIO_HIBAND


def nh_choice(nh_dl: float | None, nh_input: float | None) -> tuple[float | None, list[str]]:
    """N_H fixed to the 21 cm value (free N_H runs ~1e20 low and biases kT up)."""
    if nh_dl is None:
        return nh_input, ['nh_from_input']
    if nh_input and not (1 / NH_MISMATCH <= nh_input / nh_dl <= NH_MISMATCH):
        return nh_dl, ['nh_input_mismatch']
    return nh_dl, []


def redshift_range(z: float, counts: int, kT: float) -> tuple[float, float] | None:
    if counts >= REDSHIFT_FREE_COUNTS and kT >= HOT_KT:
        return (z * (1 - REDSHIFT_FREE_FRAC), z * (1 + REDSHIFT_FREE_FRAC))
    return None


def atomdb131_scale(kT: float) -> float:
    """Factor from an AtomDB 3.0.9 kT to the AtomDB 1.3.x kT of the same spectrum."""
    pts = ATOMDB131_SCALE
    if kT <= pts[0][0]:
        return pts[0][1]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if kT <= x1:
            return y0 + (y1 - y0) * (kT - x0) / (x1 - x0)
    return 1.0


def use_centroid(shift_arcsec: float) -> bool:
    return shift_arcsec <= CENTROID_MAX_ARCSEC


def trusted(row: dict) -> bool:
    """Good enough to move r500 (A267: a pn-only fit off the cluster gave 1.3 keV, Z 3.3)."""
    kT, e = row.get('kT'), row.get('e_kT')
    bad = {'bad_fit', 'kT_at_limit', 'Z_pegged_0', 'Z_high', 'redshift_at_limit'}
    return bool(kT and e and e / kT < MAX_REL_ERR_ITER) and not bad & set(qa_flags(row))


def qa_flags(row: dict) -> list[str]:
    """Flags of one fit row (kT, Z, chi2, dof, sp_index_min)."""
    f = []
    kT, Z = row.get('kT'), row.get('Z')
    if kT is None:
        return ['fit_failed']
    red = row['chi2'] / row['dof'] if row.get('dof') else None
    if red is not None and red > REDCHI_BAD:
        f.append('bad_fit')
    elif red is not None and red > REDCHI_MULTI_T:
        f.append('multi_T')
    if kT <= KT_LIMITS[0] * 1.02 or kT >= KT_LIMITS[1] * 0.98:
        f.append('kT_at_limit')
    if Z is not None and Z <= Z_PEGGED:
        f.append('Z_pegged_0')
    if Z is not None and Z > Z_HIGH:
        f.append('Z_high')
    zr, zf = row.get('redshift_range'), row.get('redshift_fit')
    if zr and zf is not None and (zf <= zr[0] * 1.0001 or zf >= zr[1] * 0.9999):
        f.append('redshift_at_limit')
    sp = row.get('sp_index_min')
    if sp is not None and sp <= SP_INDEX_LO * 1.001:
        f.append('sp_index_at_limit')
    if kT < ATOMDB_SENSITIVE_KT:
        f.append('atomdb_sensitive')
    return f


def camera_consistency(rows: list[dict]) -> list[str]:
    """pn vs MOS kT of the same region."""
    by = {r['fit']: r.get('kT') for r in rows}
    if by.get('pn') and by.get('mos'):
        q = by['pn'] / by['mos']
        if not PN_MOS_RANGE[0] <= q <= PN_MOS_RANGE[1]:
            return [f'pn_mos_{q:.2f}']
    return []


# --- verdict per cluster -----------------------------------------------------------------------
CHECK_PLAN = ('all_mos_dropped', 'centroid_rejected', 'sky_unconstrained', 'r500_not_converged',
              'r500_not_iterated_untrusted_fit')
CHECK_FIT = ('multi_T', 'bad_fit', 'pn_mos_', 'redshift_at_limit')


def verdict(plan_flags: list[str], best: dict | None) -> str:
    """'rejected' (no trusted primary fit), 'check' (a flag a person should look at) or 'ok'."""
    if not best or not trusted(best):
        return 'rejected'
    fit = (best.get('flags') or '').split(';')
    if any(f in CHECK_PLAN for f in plan_flags) or any(f.startswith(CHECK_FIT) for f in fit if f):
        return 'check'
    return 'ok'
