# Verification evidence

Reproduced by `python -m pytest tests -q` (9 tests) and
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
  visible.
- `root_convergence.gif`: 41 frames of the recorded bisection convergence
  on the implicit equation (residual history animated; converged root
  marked).
- `run_summary.json`: metrics + runtime context. Operating point
  v = 1.5 m/s: Re = 74,850, f = 0.0221, Δp = 50.59 kPa per 100 m,
  converged.

## Limits of this evidence

- Standard correlation accuracy: Colebrook itself is an empirical
  correlation (~±5% against experiment in its fit range) — our tests
  verify the solver against the equation and analytic limits, not the
  correlation against nature.
- Transition band (2300–4000) has no exact physics; results there are
  flagged approximations.
- Single-phase incompressible fully developed flow only.
