"""Verification tests for pipeflow_bench.

Each test pins a distinct property:
- laminar friction is exactly 64/Re and used below transition;
- the Colebrook root satisfies the implicit equation to the requested
  tolerance (residual check, not a copy of the solver loop);
- smooth-pipe limit: Colebrook -> Blasius within a few percent at
  moderate Re (documented approximation band);
- rough-pipe limit: f becomes independent of Re at high Re (fully rough
  regime), matching the analytic fully-rough value
  1/sqrt(f) = -2 log10(eps/3.7D) within tolerance;
- Haaland agrees with Colebrook within its documented ~2-3% accuracy;
- monotonicity: f increases with roughness and decreases with Re (smooth);
- regime flags incl. the transition band;
- dimensional consistency of the pressure loss against a hand computation;
- input validation (nonpositive/NaN everywhere).
"""

import numpy as np
import pytest

from pipeflow_bench import (
    RE_LAMINAR_MAX,
    RE_TURBULENT_MIN,
    blasius_friction,
    colebrook_friction,
    colebrook_residual,
    haaland_friction,
    laminar_friction,
    pressure_loss,
    reynolds,
    solve_pipe,
)


def test_laminar_friction_exact_and_regime():
    assert laminar_friction(1600.0) == pytest.approx(64.0 / 1600.0)
    # v=0.02 m/s in a 50 mm pipe of water: Re ~ 998 (deep laminar).
    result = solve_pipe(rho=998.0, mu=1.0e-3, v=0.02, D=0.05, L=20.0,
                        eps=1.5e-5)
    assert result["regime"] == "laminar"
    assert result["friction_factor"] == pytest.approx(
        64.0 / result["Re"], rel=1e-14
    )
    assert result["converged"]


def test_colebrook_root_satisfies_equation():
    for Re in (4000.0, 1.0e4, 1.0e5, 1.0e6):
        for rel_rough in (0.0, 1e-5, 1e-3, 0.01):
            f, history, converged = colebrook_friction(Re, rel_rough)
            assert converged
            assert abs(colebrook_residual(f, Re, rel_rough)) <= 1e-12
            # The final residual is smaller here; individual bisection
            # residual magnitudes need not decrease monotonically.
            assert abs(history[-1]) <= abs(history[0]) + 1e-15


def test_smooth_limit_approaches_blasius():
    # At moderate Re on a smooth pipe, Colebrook and Blasius agree within
    # the documented approximation band of Blasius (a few percent).
    for Re in (1.0e4, 5.0e4):
        f_c, _, _ = colebrook_friction(Re, 0.0)
        f_b = blasius_friction(Re)
        assert abs(f_c - f_b) / f_b < 0.05


def test_fully_rough_limit_matches_analytic():
    rel_rough = 0.01
    f_analytic = (1.0 / (-2.0 * np.log10(rel_rough / 3.7))) ** 2
    # At Re = 1e8/1e9 the viscous term 2.51/(Re sqrt(f)) still shifts the
    # root by ~1.6e-5 relative per decade - tolerances reflect that
    # physical offset, which shrinks as Re grows (measured 1.6e-6 at 1e9).
    f8, _, _ = colebrook_friction(1.0e8, rel_rough)
    assert f8 == pytest.approx(f_analytic, rel=5e-4)
    f9, _, _ = colebrook_friction(1.0e9, rel_rough)
    assert f9 == pytest.approx(f_analytic, rel=1e-5)
    assert abs(f9 - f_analytic) < abs(f8 - f_analytic)  # converging to it


def test_haaland_accuracy_band():
    for Re in (1.0e4, 1.0e5, 1.0e6):
        for rel_rough in (1e-5, 1e-3):
            f_c, _, _ = colebrook_friction(Re, rel_rough)
            f_h = haaland_friction(Re, rel_rough)
            assert abs(f_h - f_c) / f_c < 0.03


def test_friction_monotonicity():
    f1, _, _ = colebrook_friction(1.0e5, 1e-4)
    f2, _, _ = colebrook_friction(1.0e5, 1e-3)
    assert f2 > f1  # rougher -> higher f
    f3, _, _ = colebrook_friction(1.0e6, 1e-4)
    assert f3 < f1  # smooth-ish: higher Re -> lower f


def test_transition_band_flagged():
    Re_band = 3000.0
    result = solve_pipe(rho=998.0, mu=998.0 * 0.05 * 0.05 / Re_band, v=0.05,
                        D=0.05, L=1.0, eps=1e-5)
    assert RE_LAMINAR_MAX <= result["Re"] < RE_TURBULENT_MIN
    assert result["regime"] == "transition-band"


def test_pressure_loss_hand_computation():
    # f=0.02, L=100 m, D=0.1 m, rho=1000, v=2 m/s:
    # dp = 0.02 * 1000 * (1000*4/2) = 0.02*1000*2000 = 40000 Pa.
    assert pressure_loss(0.02, 100.0, 0.1, 1000.0, 2.0) == pytest.approx(4.0e4)
    dp = solve_pipe(rho=1000.0, mu=1e-3, v=2.0, D=0.1, L=100.0, eps=4.5e-5)
    expected = dp["friction_factor"] * (100.0 / 0.1) * (1000.0 * 4.0 / 2.0)
    assert dp["dp_Pa"] == pytest.approx(expected, rel=1e-14)


def test_input_validation():
    with pytest.raises(ValueError, match="mu"):
        reynolds(998.0, 1.0, 0.05, -1e-3)
    with pytest.raises(ValueError, match="Re must be finite positive"):
        laminar_friction(float("nan"))
    with pytest.raises(ValueError, match="Colebrook requires"):
        colebrook_friction(2000.0, 1e-4)
    with pytest.raises(ValueError, match="roughness"):
        colebrook_friction(1.0e5, 0.5)
    with pytest.raises(ValueError, match="L"):
        pressure_loss(0.02, -1.0, 0.1, 1000.0, 2.0)


def test_requested_residual_is_not_replaced_by_width_tolerance():
    # Previously converged at |g|=1.031374e-9 > requested 1e-9.
    f, history, converged = colebrook_friction(1e7, 0, tol=1e-9)
    assert converged
    assert abs(colebrook_residual(f, 1e7, 0)) <= 1e-9
    assert history[-1] == colebrook_residual(f, 1e7, 0)


@pytest.mark.parametrize("velocity", [0.01, 0.03, 0.10])
@pytest.mark.parametrize("roughness", [float("nan"), float("inf"), -1e-5, 0.02])
def test_invalid_roughness_rejected_before_every_regime(velocity, roughness):
    with pytest.raises(ValueError, match="roughness"):
        solve_pipe(1000, .001, velocity, .1, 10, roughness)


@pytest.mark.parametrize("tol", [0, -1, float("nan"), float("inf"), True])
def test_invalid_tolerance_rejected(tol):
    with pytest.raises(ValueError, match="tol"):
        colebrook_friction(1e5, 0, tol=tol)
    with pytest.raises(ValueError, match="tol"):
        solve_pipe(1000, .001, .01, .1, 10, 0, tol=tol)


@pytest.mark.parametrize("budget", [0, -1, 1.5, True])
def test_invalid_iteration_budget_rejected(budget):
    with pytest.raises(ValueError, match="max_iter"):
        colebrook_friction(1e5, 0, max_iter=budget)


def test_budget_exhaustion_and_stagnation_do_not_claim_convergence():
    f, history, converged = colebrook_friction(1e5, 1e-4, max_iter=1)
    assert not converged and len(history) == 1
    assert history[-1] == colebrook_residual(f, 1e5, 1e-4)
    f, history, converged = colebrook_friction(1e7, 0, tol=1e-30)
    assert (abs(history[-1]) <= 1e-30) == converged
    assert len(history) < 200  # exact success or finite precision stagnation


@pytest.mark.parametrize("bad", [0, -1, float("nan"), float("inf")])
def test_scalar_and_array_laminar_validation_match(bad):
    for value in (bad, [1000, bad]):
        with pytest.raises(ValueError, match="Re"):
            laminar_friction(value)
    with pytest.raises(ValueError, match="f"):
        pressure_loss(bad, 10, .1, 1000, 1)


def test_nonrepresentable_derived_results_raise():
    with pytest.raises(FloatingPointError, match="Re"):
        reynolds(1e308, 1e308, .1, 1)
    with pytest.raises(FloatingPointError, match="pressure"):
        pressure_loss(.02, 100, .1, 1e308, 1e308)


def test_recorded_trace_contains_actual_bracketed_evaluations():
    from pipeflow_bench import colebrook_trace
    trace = colebrook_trace(4000, 1e-4)
    f, history, converged = colebrook_friction(4000, 1e-4)
    assert trace["converged"] == converged
    assert trace["friction_factor"] == f
    assert trace["residual_history"] == history
    widths = []
    for row in trace["iterations"]:
        assert row["lo"] <= row["f"] <= row["hi"]
        assert row["residual"] == colebrook_residual(row["f"], 4000, 1e-4)
        assert colebrook_residual(row["lo"], 4000, 1e-4)*colebrook_residual(row["hi"],4000,1e-4) <= 0
        widths.append(row["hi"]-row["lo"])
    assert np.all(np.diff(widths) < 0)


def test_reference_approximations_reject_outside_stated_domains():
    with pytest.raises(ValueError, match="Blasius"):
        blasius_friction([1e4, 1e7])
    with pytest.raises(ValueError, match="Haaland"):
        haaland_friction(3000, 0)


def test_churchill_laminar_limit_is_64_over_Re():
    from pipeflow_bench import churchill_friction

    for re_ in (10.0, 100.0, 1000.0):
        assert churchill_friction(re_, 0.0) == pytest.approx(64.0 / re_, rel=2e-3)


def test_churchill_tracks_colebrook_in_turbulent_regime():
    from pipeflow_bench import churchill_friction

    worst = 0.0
    for re_ in (4000.0, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8):
        for eps_d in (0.0, 1.0e-5, 1.0e-4, 1.0e-3, 5.0e-3, 0.05):
            f_ch = churchill_friction(re_, eps_d)
            f_cb, _, converged = colebrook_friction(re_, eps_d)
            assert converged
            worst = max(worst, abs(f_ch - f_cb) / f_cb)
    assert worst < 0.04, f"Churchill vs Colebrook worst deviation {worst:.3%}"


def test_churchill_monotone_decreasing_in_Re_smooth():
    from pipeflow_bench import churchill_friction

    res = np.geomspace(4000.0, 1.0e7, 50)
    f = churchill_friction(res, 0.0)
    assert np.all(np.diff(f) < 0)


def test_churchill_input_validation():
    from pipeflow_bench import churchill_friction

    with pytest.raises(ValueError):
        churchill_friction(0.0, 0.0)
    with pytest.raises(ValueError):
        churchill_friction(1.0e5, -1.0)
    with pytest.raises(ValueError):
        churchill_friction(1.0e5, 0.2)
    with pytest.raises(ValueError):
        churchill_friction(np.nan, 0.0)
