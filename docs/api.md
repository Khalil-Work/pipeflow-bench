# API reference

Module `pipeflow_bench.colebrook` (re-exported at package level).

- `reynolds(rho, v, D, mu) -> Re` — validated.
- `laminar_friction(Re) -> 64/Re` — scalar (validated) or array.
- `colebrook_residual(f, Re, rel_rough) -> g(f)` — the implicit equation's
  residual at a trial f.
- `colebrook_friction(Re, rel_rough, tol=1e-12, max_iter=200)
  -> (f, residual_history, converged)` — bracketed bisection; requires
  Re ≥ 4000 and 0 ≤ ε/D < 0.1 (ValueError otherwise; FloatingPointError if
  not bracketed).
- `haaland_friction(Re, rel_rough)`, `blasius_friction(Re)` — explicit
  comparison references (not used by the solver).
- `pressure_loss(f, L, D, rho, v) -> Δp [Pa]` — Darcy-Weisbach.
- `solve_pipe(rho, mu, v, D, L, eps, tol=1e-12) -> dict` — full pipeline:
  `Re`, `regime` ('laminar'/'turbulent'/'transition-band'),
  `friction_factor`, `dp_Pa`, `mass_flow_kg_s`, `residual_history`,
  `converged`, `relative_roughness`.

## Example (executed)

```python
from pipeflow_bench import solve_pipe

result = solve_pipe(rho=998.0, mu=1.0e-3, v=1.5, D=0.05, L=100.0, eps=4.5e-5)
print(result["Re"], result["regime"], result["friction_factor"],
      result["dp_Pa"], result["converged"])
# 74850.0 turbulent 0.0221... 50594.17... True
```
