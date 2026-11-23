import copy

from xmmcluster import fitlog_schema as fs

RECORD = dict(
    schema=2, record_id='0123456789abcdef', paper='lovisari2020', cluster='A2029', label='r015-r500_1T_xmm3',
    obsids=['0111270201'], region={'shape': 'annulus', 'inner_r500': 0.15, 'outer_r500': 1.0},
    redshift=0.0767, nh_1e22=0.0325,
    spectra=[{'instrument': 'mos1', 'camera': 'mos1S001', 'obsid': '0111270201', 'counts': 91048,
              'area_arcmin2': 27.4, 'exposure_s': 34550.6, 'qpb_fraction': 0.011}],
    model={'expressions': {'main': 'constant*apec'}}, xcm_text='statistic chi\nmodel constant*apec\n', initial_params=[{'index': 2, 'par': 'kT', 'value': 8.0, 'role': 'icm'}],
    extra_commands=['fit 500'], model_choice={'icm_component': 'apec', 'n_icm_temperatures': 1},
    converged=True, result={'dof': 1200, 'params': []}, summary={'kT1': 8.3, 'kT1_sigma': 0.2, 'abundance': 0.35},
    target={'kT': 8.2}, run_id='lovisari2020_20261121T103059', time_utc='2026-11-21T10:30:59-05:00')


def test_valid_record():
    assert fs.validate(RECORD) == []


def test_problems_are_reported():
    r = copy.deepcopy(RECORD)
    del r['summary']
    assert 'missing summary' in fs.validate(r)
    r = copy.deepcopy(RECORD)
    r['obsids'] = ['12345']
    assert any('not 10 digits' in p for p in fs.validate(r))
    r = copy.deepcopy(RECORD)
    r['initial_params'][0]['role'] = 'mystery'
    assert any('unknown role' in p for p in fs.validate(r))
    r = copy.deepcopy(RECORD)
    r['converged'] = 'yes'
    assert any(p.startswith('converged') for p in fs.validate(r))
    r = copy.deepcopy(RECORD)
    r['extra'] = 1
    assert 'unknown field extra' in fs.validate(r)


def test_training_label():
    assert fs.training_label(RECORD) == {'kT': 8.3, 'e_kT': 0.2, 'Z': 0.35}
    r = copy.deepcopy(RECORD)
    r['summary']['kT1_sigma'] = 4.0            # kT / sigma < 3: too uncertain to learn from
    assert fs.training_label(r) is None
    r = copy.deepcopy(RECORD)
    r['converged'] = False
    assert fs.training_label(r) is None
    r = copy.deepcopy(RECORD)
    r['summary']['abundance_fixed'] = True
    assert fs.training_label(r)['Z'] is None
