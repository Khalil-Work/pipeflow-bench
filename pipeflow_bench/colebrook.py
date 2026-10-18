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
All symbols SI. A true convergence flag guarantees |g(f)| <= tol with g
the equation above. A false flag returns the final provisional point;
the residual history is returned for diagnostics.
"""

from __future__ import annotations

import math
from numbers import Real
import numpy as np

__all__ = [
    "reynolds",
    "laminar_friction",
    "colebrook_residual",
    "colebrook_friction",
    "colebrook_trace",
    "haaland_friction",
    "blasius_friction",
    "pressure_loss",
    "solve_pipe",
    "RE_LAMINAR_MAX",
    "RE_TURBULENT_MIN",
]

RE_LAMINAR_MAX = 2300.0
RE_TURBULENT_MIN = 4000.0


def _positive(name, value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive scalar")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite positive, got {value}")
    return value


def _array_positive(name, value):
    array = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(array) | (array <= 0)):
        raise ValueError(f"{name} must be finite positive")
    return array


def _finite_result(name, value):
    if np.any(~np.isfinite(value)) or np.any(np.asarray(value) <= 0):
        raise FloatingPointError(f"derived {name} must be finite positive; rescale inputs")
    return float(value) if np.ndim(value) == 0 else value


def _roughness(value):
    rough = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(rough) | (rough < 0) | (rough >= 0.1)):
        raise ValueError("relative roughness must be finite and in [0, 0.1)")
    return rough


def _controls(tol, max_iter):
    tol = _positive("tol", tol)
    if isinstance(max_iter, (bool, np.bool_)) or not isinstance(max_iter, (int, np.integer)) or max_iter < 1:
        raise ValueError("max_iter must be an integer >= 1")
    return tol, int(max_iter)


def reynolds(rho, v, D, mu):
    """Reynolds number rho*v*D/mu [-]; finite positive scalar inputs/output."""
    rho, v, D, mu = (_positive(n, x) for n, x in (("rho",rho),("v",v),("D",D),("mu",mu)))
    return _finite_result("Re", rho * v * D / mu)


def laminar_friction(Re):
    """64/Re; validate every element consistently for scalar or array input."""
    re = _array_positive("Re", Re)
    with np.errstate(over="ignore", divide="ignore", under="ignore"):
        return _finite_result("laminar friction", 64.0/re)


def colebrook_residual(f, Re, rel_rough):
    """g(f), broadcast over positive f/Re and roughness in [0,0.1).

    Re>0 is allowed here to evaluate the flagged transition-band equation.
    The public turbulent solver separately enforces Re>=4000.
    """
    f, re = _array_positive("f", f), _array_positive("Re", Re)
    rough = _roughness(rel_rough)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        root = np.sqrt(f)
        value = 1/root + 2*np.log10(rough/3.7 + (2.51/re)/root)
    if np.any(~np.isfinite(value)):
        raise FloatingPointError("Colebrook residual is not finite")
    return float(value) if np.ndim(value) == 0 else value


def _colebrook_bisect(Re, rel_rough, tol, max_iter, trace=None):
    """Bisection; convergence means abs(g)<=tol, never width alone."""
    lo, hi = 1e-6, 0.2
    g_lo = colebrook_residual(lo, Re, rel_rough)
    g_hi = colebrook_residual(hi, Re, rel_rough)
    if np.signbit(g_lo) == np.signbit(g_hi):
        raise FloatingPointError(f"Colebrook root not bracketed for Re={Re}, eps/D={rel_rough}")
    history = []
    for iteration in range(1, max_iter+1):
        f = lo+(hi-lo)/2
        g = colebrook_residual(f, Re, rel_rough)
        history.append(g)
        if trace is not None:
            trace.append({"iteration":iteration, "lo":lo, "hi":hi, "f":f, "residual":g})
        if abs(g) <= tol:
            return f, history, True
        if f == lo or f == hi:
            break  # floating-point stagnation; do not claim tolerance was met
        if np.signbit(g) == np.signbit(g_lo):
            lo, g_lo = f, g
        else:
            hi = f
    return f, history, False


def _turbulent_inputs(Re, rel_rough, tol, max_iter):
    Re = _positive("Re", Re)
    if Re < RE_TURBULENT_MIN:
        raise ValueError(f"Colebrook requires Re >= {RE_TURBULENT_MIN}, got {Re}")
    rough = _roughness(rel_rough)
    if rough.ndim != 0:
        raise ValueError("relative roughness must be a scalar for the solver")
    tol, max_iter = _controls(tol, max_iter)
    return Re, float(rough), tol, max_iter


def colebrook_friction(Re, rel_rough, tol=1e-12, max_iter=200):
    """Return (f, residual_history, converged) for Re>=4000, 0<=eps/D<0.1.

    A true flag guarantees abs(residual_history[-1])<=tol. Exhausted budgets
    or floating-point stagnation return the last evaluated point with False.
    """
    return _colebrook_bisect(*_turbulent_inputs(Re, rel_rough, tol, max_iter))


def colebrook_trace(Re, rel_rough, tol=1e-12, max_iter=200):
    """Same solve, plus actual pre-probe brackets and evaluated midpoint data."""
    trace = []
    f, history, converged = _colebrook_bisect(*_turbulent_inputs(Re, rel_rough, tol, max_iter), trace=trace)
    return {"friction_factor":f, "residual_history":history,
            "converged":converged, "iterations":trace}


def haaland_friction(Re, rel_rough):
    """Explicit turbulent comparison; Re>=4000, 0<=eps/D<0.1."""
    re = _array_positive("Re", Re)
    if np.any(re < RE_TURBULENT_MIN):
        raise ValueError("Haaland comparison requires Re >= 4000")
    rough = _roughness(rel_rough)
    value = 1 / (-1.8*np.log10((rough/3.7)**1.11 + 6.9/re))**2
    return _finite_result("Haaland friction", value)


def blasius_friction(Re):
    """Smooth-pipe comparison, restricted to 4000<=Re<=1e5."""
    re = _array_positive("Re", Re)
    if np.any((re < 4000) | (re > 1e5)):
        raise ValueError("Blasius comparison requires 4000 <= Re <= 1e5")
    return _finite_result("Blasius friction", 0.316/re**0.25)


def pressure_loss(f, L, D, rho, v):
    """Darcy-Weisbach pressure loss [Pa]; positive finite f and SI inputs."""
    f, L, D, rho, v = (_positive(n,x) for n,x in (("f",f),("L",L),("D",D),("rho",rho),("v",v)))
    return _finite_result("pressure loss", f*(L/D)*(rho*v*v/2))


def solve_pipe(rho, mu, v, D, L, eps, tol=1e-12):
    """Regime -> friction -> pressure loss; invalid roughness raises in ALL regimes.

    Returns Re, regime, friction, pressure loss, mass flow, residual history,
    convergence flag and relative roughness. The transition band explicitly
    uses the turbulent correlation as a flagged modeling choice. Laminar
    flow is direct; its diagnostic residual [0] is not a Colebrook evaluation.
    """
    rho, mu, v, D, L = (_positive(n,x) for n,x in (("rho",rho),("mu",mu),("v",v),("D",D),("L",L)))
    if isinstance(eps, (bool, np.bool_)) or not isinstance(eps, Real) or not math.isfinite(float(eps)) or eps < 0:
        raise ValueError("eps roughness must be a finite non-negative scalar")
    tol, _ = _controls(tol, 200)
    rel_rough = float(_roughness(float(eps)/D))
    Re = reynolds(rho, v, D, mu)
    if Re < RE_LAMINAR_MAX:
        f = laminar_friction(Re)
        regime, history, converged = "laminar", [0.0], True
    else:
        f, history, converged = _colebrook_bisect(Re, rel_rough, tol, 200)
        regime = "turbulent" if Re >= RE_TURBULENT_MIN else "transition-band"
    dp = pressure_loss(f, L, D, rho, v)
    mass_flow = _finite_result("mass flow", rho*v*np.pi*D*D/4)
    return {"Re":Re, "regime":regime, "friction_factor":f, "dp_Pa":dp,
            "mass_flow_kg_s":mass_flow, "residual_history":history,
            "converged":converged, "relative_roughness":rel_rough}
