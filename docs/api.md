# API reference

Public imports are re-exported by `pipeflow_bench`. Physical pipeline
inputs are finite positive real scalars in SI; roughness eps may be zero.
Booleans are rejected by scalar pipeline/solver validation. Invalid inputs
raise `ValueError`; nonrepresentable derived results raise
`FloatingPointError` instead of returning infinities or zeros.

| Function | Return and contract |
| --- | --- |
| `reynolds(rho, v, D, mu)` | Re=ρvD/μ; validates inputs and finite positive result |
| `laminar_friction(Re)` | 64/Re; scalar or array with every element positive/finite |
| `colebrook_residual(f, Re, rel_rough)` | Signed g(f); broadcasts arrays; f/Re positive and roughness in [0,0.1); rejects nonfinite residuals |
| `colebrook_friction(Re, rel_rough, tol=1e-12, max_iter=200)` | `(f, residual_history, converged)`; scalar Re>=4000 and 0<=roughness<0.1 |
| `colebrook_trace(Re, rel_rough, tol=1e-12, max_iter=200)` | Dictionary with `friction_factor`, `residual_history`, `converged`, and actual `iterations` records |
| `haaland_friction(Re, rel_rough)` | Explicit turbulent comparison, Re>=4000, roughness [0,0.1); broadcasts arrays |
| `blasius_friction(Re)` | Smooth-pipe comparison, restricted to 4000<=Re<=1e5; scalar or array |
| `pressure_loss(f, L, D, rho, v)` | Darcy-Weisbach Δp [Pa]; f and every physical scalar positive/finite |
| `solve_pipe(rho, mu, v, D, L, eps, tol=1e-12)` | Regime-specific friction, pressure loss and mass flow with diagnostics |

## Stopping and diagnostics

`tol` is a finite positive absolute tolerance on the dimensionless
Colebrook residual; `max_iter` is an integer >=1 (not a boolean).
`converged=True` guarantees abs(history[-1])<=tol. Width reduction alone
never counts as convergence. If the budget runs out or no distinct
floating-point midpoint remains, return the last evaluated f and `False`.
Every trace record stores `iteration`, `lo`, `hi`, `f`, `residual`, with
the bracket recorded **before** the probe. The final residual corresponds
to the returned f, including unsuccessful solves.

`solve_pipe` returns `Re`, `regime`, `friction_factor`, `dp_Pa`,
`mass_flow_kg_s`, `residual_history`, `converged`, `relative_roughness`.
It validates eps/D<0.1 and tolerance in all regimes. For Re<2300,
the direct laminar law has the diagnostic history `[0.0]`; this is not
a Colebrook residual. For 2300<=Re<4000 the result is explicitly flagged
`transition-band` and uses the turbulent equation as a modeling choice.
If a root solve returns False, pressure loss is only a provisional value;
callers must inspect the flag before accepting it.

## Worked operating point

```python
from pipeflow_bench import solve_pipe, colebrook_trace

result = solve_pipe(rho=998, mu=.001, v=1.5, D=.05, L=100, eps=45e-6)
assert result["converged"]
print(result["Re"], result["friction_factor"], result["dp_Pa"])
# 74850.0 0.022531362454543326 50594.17439167704

trace = colebrook_trace(4000, 1e-4)
print(trace["iterations"][-1])  # actual final bracket, midpoint and residual
```

Arrays are supported for equation/reference helpers as stated, not for
the scalar root solver or complete pipe operating point. These friction
factors use the **Darcy** convention (four times a Fanning factor).
