# 🚰 Pipeflow Bench

**Friction, pressure loss, and implicit-solver diagnostics: an original pipe-flow tool — exact laminar law, Colebrook-White solved by bracketed bisection with residual histories, and an explicit transition-band flag — with the implicit root verified against its own equation to 1e-12, the analytic fully-rough limit, and hand-computed pressure losses.**

[![tests](https://github.com/khalil-official/pipeflow-bench/actions/workflows/ci.yml/badge.svg)](https://github.com/khalil-official/pipeflow-bench/actions/workflows/ci.yml)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![python](https://img.shields.io/badge/python-3.12%20%7C%203.13-blue.svg)](requirements.txt)

<div align="center">

**▶ Watch the implicit solver converge:**

<img src="docs/results/root_convergence.gif" alt="Animated bisection convergence on the implicit Colebrook equation" width="560">

<sub>41 frames of the recorded bisection residual history on the implicit
Colebrook–White equation (Re = 4·10³, ε/D = 10⁻⁴); the dashed line marks
the converged root. Produced by
<code>examples/run_verification.py</code>.</sub>

</div>

## 🎯 Purpose and scope

A transparent pipe-flow pressure-loss tool built around an honestly solved
implicit friction relation: the Colebrook equation is solved by bracketed
bisection with visible residual histories and explicit regime/transition-band
flags, and the explicit Haaland approximation is carried as a comparison
reference — never as the solver. Scope: single-phase incompressible fully
developed flow in circular pipes. No minor losses, compressibility, or
non-circular ducts.

## 🧮 Model, method, and diagnostics

- Laminar (Re < 2300): exact f = 64/Re.
- Turbulent (Re ≥ 4000): Colebrook–White by bisection on
  g(f) = 1/√f + 2 log₁₀(ε/3.7D + 2.51/(Re√f)); bracket validated; residual
  ≤ 1e-12; full residual history returned.
- Transition band (2300–4000): flagged, turbulent correlation applied with
  the band recorded — a stated modeling choice, never silent.
- Δp = f (L/D) ρv²/2; ṁ = ρvπD²/4; finite-positive validation throughout.

Full details: [docs/model.md](docs/model.md); API:
[docs/api.md](docs/api.md).

## ⚡ Quick start

```bash
python -m pip install -r requirements-dev.txt
python -m pytest tests -q                       # 9 verification tests
python examples/run_verification.py             # figures + GIF + raw data
```

Verified: Python 3.13.5, NumPy 2.1.3, Windows 11; CI on Ubuntu. Tests
< 1 s; results script ≈ 5 s.

## 🧱 Architecture and API

```
pipeflow_bench/
└── colebrook.py        # regime logic, implicit solver, references, pipeline
tests/test_colebrook.py # 9 verification tests
examples/run_verification.py
docs/ model.md  api.md  verification.md  results/
```

```python
from pipeflow_bench import solve_pipe

result = solve_pipe(rho=998.0, mu=1.0e-3, v=1.5, D=0.05, L=100.0, eps=4.5e-5)
print(result["Re"], result["regime"], result["friction_factor"], result["dp_Pa"])
# 74850.0 turbulent 0.0221... 50594.17...
```

## 📊 Results and metrics

| Case | Metric | Value |
| --- | --- | --- |
| Colebrook root vs its own equation (Re 4e3–1e7, ε/D 0–1e-2) | max residual | **≤ 1e-12** |
| Fully-rough limit vs analytic | relative offset at Re=1e8 / 1e9 | 1.6e-5 / 1.6e-6 (viscous-term physics quantified) |
| Haaland vs Colebrook (recorded sweep) | max relative difference | **1.34%** |
| Blasius vs smooth Colebrook | relative difference | < 5% |
| Darcy-Weisbach vs hand computation | error | exact (4.0e4 Pa case) |
| Operating point v=1.5 m/s, D=50 mm | Re / f / Δp per 100 m | 74,850 / 0.0221 / 50.59 kPa |

Raw data: [docs/results/](docs/results) — `explicit_comparison.csv`,
`residual_history.csv`, `pipeline_sweep.csv`, `run_summary.json`; cover
[friction_chart.png](docs/results/friction_chart.png),
[pipeline_sweep.png](docs/results/pipeline_sweep.png),
[residual_history.png](docs/results/residual_history.png), animation
`root_convergence.gif`.

## ✅ Verification

Nine tests: equation-satisfaction at the returned root, analytic
fully-rough limit, Blasius/Haaland accuracy bands, monotonicity, transition
band flagging, exact hand pressure loss, input validation. Evidence and
limits: [docs/verification.md](docs/verification.md). Not claimed:
correlation-vs-nature accuracy (Colebrook itself is empirical ±5%);
anything outside fully developed incompressible pipe flow.

## ⚠️ Limitations

- Empirical correlations: solver correctness ≠ physical accuracy of the
  correlation.
- Transition band results are flagged approximations.
- No minor losses / fittings / networks.

## 📚 References

- S. E. Haaland, "Simple and explicit formulas for the friction factor in
  turbulent pipe flow" (1983).
- C. F. Colebrook, "Turbulent flow in pipes" (1939).
- Moody chart literature (L. F. Moody, 1944).

Dependencies: NumPy (BSD-3), Matplotlib (PSF-based), Pillow (HPND), pytest
(MIT). This repository's own code is original and MIT-licensed; see
[LICENSE](LICENSE).
