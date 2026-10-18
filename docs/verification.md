# Verification evidence

Reproduced by `python -m pytest tests -q` (39 test cases) and
`python examples/run_verification.py` (figures + GIF + CSV/JSON under
`docs/results/`). Verified on Windows 11, Python 3.13.5, NumPy 2.1.3; CI
on Ubuntu.

## What each test establishes

| Test | Property established |
| --- | --- |
| `test_laminar_friction_exact_and_regime` | f = 64/Re exactly and used below Re = 2300. |
| `test_colebrook_root_satisfies_equation` | The returned root satisfies the implicit equation to ≤ 1e-12 residual across a Re × roughness grid (residual check, not a mirror of the solver loop). |
| `test_smooth_limit_approaches_blasius` | Smooth-pipe Colebrook matches Blasius within 5% in its validity range. |
| `test_fully_rough_limit_matches_analytic` | Approaches the analytic fully-rough value with the viscous term's physical offset quantified (1.6e-5 relative at Re = 1e8, shrinking to 1.6e-6 at 1e9). |
| `test_haaland_accuracy_band` | Haaland within 3% of Colebrook across the sweep. |
| `test_friction_monotonicity` | f increases with roughness; decreases with Re (smooth-ish). |
| `test_transition_band_flagged` | 2300 ≤ Re < 4000 results carry the explicit band flag. |
| `test_pressure_loss_hand_computation` | Darcy-Weisbach matches a hand computation exactly (4.0e4 Pa case). |
| `test_input_validation` | Nonpositive/NaN physical inputs and out-of-domain Colebrook calls raise with clear messages. |

## Recorded run (`docs/results/`)

- `friction_chart.png/.svg`: regime chart with laminar line, Colebrook
  curves for four roughnesses, Blasius reference; companion
  `explicit_comparison.csv` with the measured Haaland maximum relative
  difference **1.34%**.
- `residual_history.png/.csv`: bisection residual convergence to ≤ 1e-12.
- `pipeline_sweep.png/.csv`: Δp per 100 m vs velocity for water in a
  50 mm pipe (ε = 45 µm) with regime split; laminar-to-turbulent transition
  included explicitly, with a separate low-speed panel.
- `root_convergence.gif`: 43 frames of the recorded bisection convergence
  on the implicit equation (residual history animated; converged root
  marked).
- `run_summary.json`: metrics + runtime context. Operating point
  v = 1.5 m/s: Re = 74,850, f = 0.0225314, Δp = 50.59 kPa per 100 m,
  converged.

## Limits of this evidence

- Standard correlation accuracy: Colebrook itself is an empirical
  correlation; application-dependent physical error is not assessed — our tests
  verify the solver against the equation and analytic limits, not the
  correlation against nature.
- Transition band (2300–4000) has no exact physics; results there are
  flagged approximations.
- Single-phase incompressible fully developed flow only.


## New regressions and saved solver states

The additional 30 collected cases cover the reproduced width-stop failure
at Re=1e7, NaN/Inf/negative/out-of-domain roughness in all three regimes,
invalid tolerances and budgets, consistent scalar/array laminar validation,
invalid friction inputs, nonrepresentable derived values, comparison-domain
guards and failed-budget/stagnation flags. A trace check independently
evaluates every recorded midpoint and both bracket-end signs.

`bisection_trace.csv` stores all **43 actual solver probes**, including
pre-probe lo/hi, f and signed residual. The animated root is
**0.0400084312335583**, final |g| **2.1e-13**. The GIF is 1,100×600 px,
white with no grids, and shows the full residual history and shrinking
bracket. It represents iteration count, not physical time.

`explicit_comparison.csv` now includes each residual, probe count and
convergence flag; its 30 converged roots meet 1e-12. `stopping_regression.csv`
records the Re=1e7 smooth-pipe solve at tolerances 1e-9, 1e-12 and 1e-30
(the last residual happens to evaluate to exactly zero in this run), plus
a Re=4000, roughness=1e-5 solve that correctly returns False at tolerance
1e-30 on this runtime. The pressure sweep has **100 operating points**, with
laminar, transition-band and turbulent labels preserved.

The low residuals establish numerical solution of the chosen empirical
equation. No experiment, pressure-pipe test or universal uncertainty
interval is claimed by this evidence.
