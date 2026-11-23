import math

from xmmcluster import rules
from xmmcluster.initmodel import replay


def test_scaling_relations_are_consistent():
    # M500 from kT and back (Arnaud et al. 2005)
    m = rules.M500_5KEV * (7.0 / 5.0) ** rules.M500_SLOPE / rules.ez(0.2)
    assert math.isclose(rules.kT_from_m500(m, 0.2), 7.0, rel_tol=1e-9)
    assert math.isclose(rules.r500_kpc_from_kT(5.0, 0.0), rules.R500_5KEV_KPC)


def test_bands_depend_on_temperature():
    assert rules.bands(1.0)['pn'][0] == 0.5      # groups need the Fe-L band
    assert rules.bands(6.0)['pn'][0] == 1.0      # pn below 1 keV biases hot-cluster kT


def test_quality_flags_and_verdict():
    good = dict(kT=5.0, e_kT=0.2, Z=0.3, chi2=1100, dof=1000)
    assert rules.qa_flags(good) == []
    assert 'kT_at_limit' in rules.qa_flags(dict(good, kT=29.8))
    assert 'bad_fit' in rules.qa_flags(dict(good, chi2=2500))
    assert rules.verdict([], good) == 'ok'
    assert rules.verdict(['centroid_rejected'], good) == 'check'
    assert rules.verdict([], dict(good, e_kT=2.0)) == 'rejected'   # 40 % error: not trusted
    assert rules.verdict([], None) == 'rejected'


def test_replay_summary():
    rows = [dict(record_id='a', start='logged', stat='100', kT='5.0', seconds='10'),
            dict(record_id='a', start='1keV', stat='130', kT='1.2', seconds='30'),
            dict(record_id='b', start='logged', stat='50', kT='2.0', seconds='12'),
            dict(record_id='b', start='1keV', stat='50.2', kT='2.02', seconds='20')]
    s = {r['start']: r for r in replay.summarize(rows)}
    assert s['logged']['at_best_minimum'] == 1.0
    assert s['1keV']['at_best_minimum'] == 0.5 and s['1keV']['kT_off_5pct'] == 0.5
