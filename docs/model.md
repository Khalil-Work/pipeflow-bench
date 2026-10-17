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
≤ 1e-12 or relative width ≤ tol). Iteration count ~60 for full double
precision; residual histories returned for diagnostics.

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
charting; scalar misuse still raises.

## Limits

Single-phase incompressible fully developed flow; no minor losses, no
compressibility, no non-circular ducts. The transition band has no exact
correlation — the band flag exists precisely so users know when physics is
uncertain.
