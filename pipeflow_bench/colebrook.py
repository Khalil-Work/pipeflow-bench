"""Pipeflow Bench: pipe-flow pressure loss with laminar/turbulent handling.

Friction factor by regime:
- Laminar (Re < 2300):  f = 64 / Re                       (exact)
- Turbulent (Re >= 4000): Colebrook-White implicit equation, solved by
  bisection on
      1/sqrt(f) = -2 log10( eps/(3.7 D) + 2.51/(Re sqrt(f)) )
  with a documented transition band 2300 <= Re < 4000 flagged (results in
  the band are reported with a band flag, never silently interpolated).

Pressure loss (Darcy-Weisbach):  dp = f (L/D) (rho v^2 / 2)   [Pa]
Mass flow from velocity:         mdot = rho v A              [kg/s]

Reference explicit approximations (for comparison, never as the solver):
- Haaland (1983):  1/sqrt(f) = -1.8 log10( (eps/3.7D)^1.11 + 6.9/Re )
- Blasius (smooth, Re < 1e5): f = 0.316 Re^-0.25

Dimension checks: Re = rho v D / mu [-]; eps/D relative roughness [-].
All symbols SI. The Colebrook residual is monitored: the returned root
satisfies |g(f)| <= tol with g the equation above; the residual history is
returned for diagnostics.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "reynolds",
    "laminar_friction",
    "colebrook_residual",
    "colebrook_friction",
    "haaland_friction",
    "blasius_friction",
    "pressure_loss",
    "solve_pipe",
    "RE_LAMINAR_MAX",
    "RE_TURBULENT_MIN",
]

RE_LAMINAR_MAX = 2300.0
RE_TURBULENT_MIN = 4000.0


def reynolds(rho, v, D, mu):
    """Reynolds number rho*v*D/mu [-] with validation."""
    for name, value in (("rho", rho), ("v", v), ("D", D), ("mu", mu)):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite positive, got {value}")
    return rho * v * D / mu


def laminar_friction(Re):
    if np.ndim(Re) == 0:
        if not np.isfinite(Re) or Re <= 0.0:
            raise ValueError(f"Re must be finite positive, got {Re}")
        return 64.0 / Re
    return 64.0 / np.asarray(Re, dtype=float)


def colebrook_residual(f, Re, rel_rough):
    """g(f) = 1/sqrt(f) + 2 log10( eps/(3.7D) + 2.51/(Re sqrt(f)) )."""
    return 1.0 / np.sqrt(f) + 2.0 * np.log10(
        rel_rough / 3.7 + 2.51 / (Re * np.sqrt(f))
    )


def _colebrook_bisect(Re, rel_rough, tol, max_iter):
    """Bracketed bisection on the Colebrook equation (no Re-domain guard)."""
    lo, hi = 1e-6, 0.2
    g_lo, g_hi = colebrook_residual(lo, Re, rel_rough), colebrook_residual(hi, Re, rel_rough)
    if g_lo * g_hi > 0.0:
        raise FloatingPointError(
            f"Colebrook root not bracketed for Re={Re}, eps/D={rel_rough}"
        )
    history = []
    f = 0.5 * (lo + hi)
    converged = False
    for _ in range(max_iter):
        f = 0.5 * (lo + hi)
        g = colebrook_residual(f, Re, rel_rough)
        history.append(g)
        if abs(g) <= tol or (hi - lo) <= tol * f:
            converged = True
            break
        if g * g_lo > 0.0:
            lo, g_lo = f, g
        else:
            hi = f
    return f, history, converged


def colebrook_friction(Re, rel_rough, tol=1e-12, max_iter=200):
    """Solve Colebrook-White for f by bracketed bisection.

    Valid for Re >= RE_TURBULENT_MIN and 0 <= rel_rough < 0.1 (documented
    domain). Returns (f, residual_history, converged_flag).
    """
    if not np.isfinite(Re) or Re < RE_TURBULENT_MIN:
        raise ValueError(
            f"Colebrook requires Re >= {RE_TURBULENT_MIN}, got {Re} "
            "(use laminar_friction below transition or flag the band)"
        )
    if not (np.isfinite(rel_rough) and 0.0 <= rel_rough < 0.1):
        raise ValueError(f"relative roughness must be in [0, 0.1), got {rel_rough}")
    return _colebrook_bisect(Re, rel_rough, tol, max_iter)


def haaland_friction(Re, rel_rough):
    """Explicit Haaland approximation (comparison reference)."""
    return (1.0 / (-1.8 * np.log10((rel_rough / 3.7) ** 1.11 + 6.9 / Re))) ** 2


def blasius_friction(Re):
    """Smooth-pipe Blasius approximation (Re < 1e5 comparison reference)."""
    return 0.316 / Re**0.25


def pressure_loss(f, L, D, rho, v):
    """Darcy-Weisbach pressure loss [Pa]."""
    for name, value in (("L", L), ("D", D), ("rho", rho), ("v", v)):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite positive, got {value}")
    return f * (L / D) * (rho * v * v / 2.0)


def solve_pipe(rho, mu, v, D, L, eps, tol=1e-12):
    """Full pipeline: regime detection -> friction -> pressure loss.

    Returns a dict with Re, regime ('laminar'/'turbulent'/'transition-band'),
    friction factor (Colebrook in the band, flagged), dp [Pa], mdot [kg/s],
    residual history, and convergence flag. The transition band uses the
    turbulent correlation with an explicit band warning recorded - it is a
    modeling choice, stated, not hidden.
    """
    Re = reynolds(rho, v, D, mu)
    rel_rough = eps / D
    if Re < RE_LAMINAR_MAX:
        f = laminar_friction(Re)
        regime = "laminar"
        history, converged = [0.0], True
    elif Re >= RE_TURBULENT_MIN:
        f, history, converged = colebrook_friction(Re, rel_rough, tol=tol)
        regime = "turbulent"
    else:
        # Transition band: the turbulent correlation is used WITH an
        # explicit band flag - a stated modeling choice, not hidden
        # interpolation. The private bisection has no Re-domain guard.
        f, history, converged = _colebrook_bisect(Re, rel_rough, tol, 200)
        regime = "transition-band"
    dp = pressure_loss(f, L, D, rho, v)
    return {
        "Re": Re,
        "regime": regime,
        "friction_factor": f,
        "dp_Pa": dp,
        "mass_flow_kg_s": rho * v * np.pi * D * D / 4.0,
        "residual_history": history,
        "converged": converged,
        "relative_roughness": rel_rough,
    }
