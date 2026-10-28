"""Lock the within-person diet -> gut-function harness (example/diet_function_timeseries): on the
synthetic demo, the planted personalized fiber -> SCFA-capacity coupling must be 'consistent' (pooled
within-person r>0, significant, majority same-sign) while the null glycogen axis must NOT be."""
import os
import sys
import pytest

pytest.importorskip("statsmodels")
pytest.importorskip("scipy")

EX = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "example", "diet_function_timeseries"))


def _by_axis():
    # Several examples ship analyze.py / make_demo.py; force this one's even if another is cached.
    sys.path.insert(0, EX)
    for m in ("analyze", "make_demo"):
        sys.modules.pop(m, None)
    import analyze
    res = analyze.run_demo(verbose=False)
    return {r["axis"]: r for r in res}


def test_scfa_couples_within_person():
    by = _by_axis()
    assert "scfa" in by
    scfa = by["scfa"]
    assert scfa["k"] >= 10                       # pooled across many people
    assert scfa["pooled_r"] > 0                  # planted positive fiber -> SCFA coupling
    assert scfa["p_value"] < 0.05                # significant pooled effect
    assert scfa["frac_same_sign"] >= 0.6         # most people share the direction
    assert scfa["consistent"] is True


def test_null_axis_not_consistent():
    by = _by_axis()
    assert "glycogen_starch" in by
    assert by["glycogen_starch"]["consistent"] is False   # no planted effect -> must not be flagged


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
