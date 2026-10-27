"""Lock the American Gut diet->gut->BMI mediation harness (example/american_gut_mediation): on the
synthetic demo (planted partial mediation diet -> higher diversity -> lower BMI), the harness must
recover a negative, significant indirect effect through Shannon and a proportion mediated strictly
between 0 and 1."""
import os
import sys
import pytest

pytest.importorskip("statsmodels")

EX = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "example", "american_gut_mediation"))


def _load():
    # Several examples ship modules named analyze.py / make_demo.py; force this one's to load even
    # if another example's are cached from an earlier test in the same pytest session.
    sys.path.insert(0, EX)
    for m in ("analyze", "make_demo"):
        sys.modules.pop(m, None)
    import analyze
    import make_demo
    return analyze, make_demo


def test_demo_recovers_partial_mediation():
    analyze, make_demo = _load()
    r = analyze.run(make_demo.make(), n_boot=200, label="test")
    assert r["outcome_type"] == "continuous"
    assert r["path_a_x_to_m"] > 0          # better diet -> higher diversity
    assert r["path_b_m_to_y"] < 0          # higher diversity -> lower BMI
    assert r["indirect_effect_nie"] < 0    # net indirect effect is negative
    assert r["nie_ci"][1] < 0              # NIE 95% CI excludes 0 (significant mediation)
    assert 0.0 < r["proportion_mediated"] < 1.0   # partial mediation


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
