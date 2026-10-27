"""Lock the cross-cohort function<->cardiometabolic reproducibility harness
(example/gut_function_cardiometabolic): on the synthetic demo, the planted SCFA-depletion-in-T2D
signal must REPLICATE (consistent direction across all 3 cohorts + significant pooled OR<1) while the
null glycogen axis must NOT be flagged as replicating."""
import os
import sys
import pytest

pytest.importorskip("statsmodels")

EX = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "example", "gut_function_cardiometabolic"))


def _by_axis():
    # Several examples ship modules named analyze.py / make_demo.py; force this one's to load even
    # if another example's are cached from an earlier test in the same pytest session.
    sys.path.insert(0, EX)
    for m in ("analyze", "make_demo"):
        sys.modules.pop(m, None)
    import analyze
    res = analyze.run_demo(verbose=False)
    return {r["axis"]: r for r in res}


def test_scfa_signal_replicates():
    by = _by_axis()
    assert "scfa" in by, "SCFA axis should be modelable in the demo cohorts"
    scfa = by["scfa"]
    assert scfa["k"] == 3                       # pooled across all three cohorts
    assert scfa["pooled_or"] < 1.0              # lower SCFA capacity in T2D (planted direction)
    assert scfa["p_value"] < 0.05               # significant pooled estimate
    assert scfa["replicates"] is True           # same direction in every cohort AND pooled-significant
    assert all(o < 1.0 for _, o, _, _ in scfa["per_cohort"])


def test_null_axis_not_flagged():
    by = _by_axis()
    assert "glycogen_starch" in by
    assert by["glycogen_starch"]["replicates"] is False   # no planted effect -> must not replicate


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
