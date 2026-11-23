import math

from xmmcluster import compare


def pt(paper='p', cluster='c', ours=5.0, published=5.0, verdict='ok', quantity='kT'):
    p = dict(paper=paper, cluster=cluster, region='r', quantity=quantity, ours=ours, e_ours=None,
             published=published, e_published=None, ratio=ours / published, autofit_verdict=verdict)
    p['left_out'] = compare.left_out(p)
    return p


def test_left_out_never_depends_on_the_published_value():
    for published in (1.0, 5.0, 100.0):
        assert compare.left_out(pt(published=published)) == ''
        assert compare.left_out(pt(published=published, verdict='rejected')) == 'autofit rejected'
        assert compare.left_out(pt(ours=29.5, published=published)) == 'kT at fit limit'
    assert compare.left_out(pt(ours=29.0)) == ''
    assert compare.left_out(pt(ours=29.5, quantity='Z')) == ''


def test_per_paper_statistics():
    pts = [pt(cluster=f'c{i}', ours=r, published=1.0) for i, r in enumerate((0.8, 0.95, 1.0, 1.05, 1.3))]
    pts.append(pt(cluster='bad', ours=2.0, published=1.0, verdict='rejected'))
    (p,) = compare.per_paper(pts)
    assert p['n_values'] == 5 and p['n_clusters'] == 5 and p['n_left_out'] == 1
    assert p['median_ratio'] == 1.0
    assert p['within_10pct'] == 0.6          # 0.95, 1.0, 1.05
    assert p['within_20pct'] == 0.8          # and 0.8
    h = compare.headline(pts)
    assert h['n_values'] == 5 and h['n_pairs'] == 5 and h['n_left_out'] == 1


def test_quantile_interpolates():
    assert compare.quantile([1, 2, 3, 4, 5], 0.5) == 3
    assert math.isclose(compare.quantile([0, 10], 0.16), 1.6)


def test_bundled_points_are_consistent():
    pts = compare.load()
    assert len(pts) > 500
    for p in pts:
        assert p['quantity'] == 'kT'
        assert math.isclose(p['ratio'], p['ours'] / p['published'], rel_tol=2e-3), p
        assert p['paper'] in compare.PAPERS
