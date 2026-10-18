# Model and method

## Problem

Darcy-Weisbach pressure loss in a circular pipe with friction factor by
flow regime.

## Relations implemented

**Reynolds number** Re = ρ v D / μ [-] (fluid density ρ [kg/m³], mean
velocity v [m/s], diameter D [m], dynamic viscosity μ [Pa·s]).

**Laminar regime (Re < 2300):** exact Hagen–Poiseuille friction factor

    f = 64 / Re.

**Turbulent regime (Re ≥ 4000):** the Colebrook–White implicit equation

    1/√f = −2 log₁₀( ε/(3.7 D) + 2.51/(Re √f) ),

solved by bracketed bisection on g(f) = 1/√f + 2 log₁₀(ε/(3.7D) +
2.51/(Re√f)) over f ∈ [1e-6, 0.2] (bracket validity checked; root residual
≤ the requested absolute tolerance is the only success criterion).
Default tol=1e-12; the recorded Re=4000, ε/D=1e-4 case takes 43 probes.
Residual histories and optional actual bracket/midpoint traces are returned.
Budget exhaustion or floating-point stagnation returns a false flag.

**Transition band (2300 ≤ Re < 4000):** flagged `"transition-band"`; the
turbulent correlation is applied WITH the band recorded — a stated modeling
choice, never silent interpolation.

**Pressure loss (Darcy-Weisbach):**

    Δp = f (L/D) (ρ v²/2)  [Pa],  ṁ = ρ v π D²/4 [kg/s].

## Explicit references (comparisons, never the solver)

- Haaland (1983): 1/√f = −1.8 log₁₀( (ε/3.7D)^1.11 + 6.9/Re ) — measured
  maximum relative difference vs Colebrook across the recorded sweep:
  **1.34%** (documented accuracy band of the approximation).
- Blasius (smooth pipes, Re < 1e5): f = 0.316 Re^−0.25 — agrees with
  Colebrook within 5% in its validity range.

## Validation and domain checks

Finite-positive validation of every physical input (ρ, v, D, μ, L);
relative roughness domain 0 ≤ ε/D < 0.1; Colebrook's public contract
requires Re ≥ 4000 (the band path uses the internal bisection directly and
flags the result). Laminar `laminar_friction` also accepts arrays for
charting; invalid elements raise consistently for scalars and arrays.
Roughness and tolerance are checked before regime selection, including
laminar flow. Derived Re, pressure loss and mass flow must remain finite
and positive; overflow/underflow raises rather than silently propagating.

## Limits

Single-phase incompressible fully developed flow; no minor losses, no
compressibility, no non-circular ducts. The transition band has no exact
correlation — the band flag exists precisely so users know when physics is
uncertain.


## Symbols, units and algorithm

| Symbol | Meaning | SI unit |
| --- | --- | --- |
| ρ | Density | kg/m³ |
| μ | Dynamic viscosity | Pa·s |
| v | Mean axial speed | m/s |
| D | Circular internal diameter | m |
| L | Pipe length | m |
| ε | Equivalent wall roughness | m |
| f | Darcy friction factor | dimensionless |
| Re, ε/D, g | Reynolds number, relative roughness, residual | dimensionless |
| Δp, ṁ | Frictional pressure loss, mass rate | Pa, kg/s |

Bisection first checks opposite residual signs on [1e-6,0.2]. It evaluates
one midpoint, records its pre-probe bracket and residual, and stops only
if abs(g)<=tol. Otherwise the sign identifies the half containing a root.
If no distinct floating midpoint exists, the solver reports nonconvergence.
Shrinking brackets do not imply monotonically decreasing residual
magnitudes: a probe can lie farther from the root than the previous probe.
The same loop drives the animation; no illustrative or fabricated probe
coordinates are used. A residual tolerance verifies the equation, not an
absolute friction-factor error bound or correlation accuracy.

The friction-chart comparison uses true Reynolds-number coordinates from
the saved 30-case CSV. The Blasius line ends at Re=1e5. The 100-point
pressure sweep includes low velocities in all three regimes and a separate
low-speed view; transition-band points are never labeled turbulent.

The model assumes steady, single-phase, incompressible, fully developed
flow in a straight circular pipe. It excludes entry effects, fittings,
elevation, compressibility, multiphase flow, heat transfer and networks.
