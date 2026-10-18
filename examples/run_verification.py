"""Produce verification figures, animation, and raw data for Pipeflow Bench.

Cases:
1. Friction-factor chart: Moody-type curves f(Re) for several eps/D with the
   laminar line, plus the Haaland/Blasius comparison references.
2. Implicit-vs-explicit: |Colebrook - Haaland| relative difference across a
   parameter sweep (the honest accuracy of the explicit approximation).
3. Solver diagnostics: bisection residual histories for representative
   operating points.
4. Pressure-loss pipeline sweep: dp vs velocity for a documented pipe.

Outputs under docs/results/: cover figure, sweep CSVs, a actual-iteration GIF
(convergence progression: the animating bracket/residual is a real
parameter-progression of solver states), run_summary.json.
"""

from __future__ import annotations

import json
import pathlib
import platform
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as mpl_animation
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from pipeflow_bench import (
    blasius_friction,
    colebrook_friction,
    colebrook_trace,
    colebrook_residual,
    haaland_friction,
    laminar_friction,
    solve_pipe,
)

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "docs" / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

BLUE, RED, GREEN, GRAY = "#1f4e79", "#b23a48", "#0b6e4f", "#555555"
plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
    }
)


def make_cover():
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 4.4))
    Re_lam = np.linspace(600.0, 2300.0, 50)
    Re_turb = np.logspace(np.log10(4000.0), 8.0, 200)
    ax[0].loglog(Re_lam, laminar_friction(Re_lam), "-", color=GRAY, lw=1.6,
                 label="laminar 64/Re")
    colors = ["#1f4e79", "#0b6e4f", "#8a6d1f", "#b23a48"]
    for rel_rough, color in zip((0.0, 1e-4, 1e-3, 1e-2), colors):
        f = [colebrook_friction(Re, rel_rough)[0] for Re in Re_turb]
        ax[0].loglog(Re_turb, f, "-", color=color, lw=1.6,
                     label=f"Colebrook ε/D = {rel_rough:g}")
    re_blasius = Re_turb[Re_turb <= 1e5]
    ax[0].loglog(re_blasius, blasius_friction(re_blasius), ":", color=GRAY, lw=1.2,
                 label="Blasius (smooth ref.)")
    ax[0].set_xlabel("Reynolds number Re [-]")
    ax[0].set_ylabel("Darcy friction factor f [-]")
    ax[0].set_title("Friction factor by regime")
    ax[0].legend(frameon=False, fontsize=7)
    ax[0].set_ylim(0.005, 0.12)

    # Implicit vs explicit comparison.
    rows = ["Re,rel_rough,f_colebrook,f_haaland,rel_diff,final_residual,iterations,converged"]
    max_diff = 0.0
    for Re in (4.0e3, 1.0e4, 5.0e4, 1.0e5, 1.0e6, 1.0e7):
        for rel_rough in (0.0, 1e-5, 1e-4, 1e-3, 1e-2):
            f_c, history, converged = colebrook_friction(Re, rel_rough)
            if not converged:
                raise RuntimeError("comparison root failed requested residual tolerance")
            f_h = haaland_friction(Re, rel_rough)
            diff = abs(f_h - f_c) / f_c
            max_diff = max(max_diff, diff)
            rows.append(f"{Re:.6e},{rel_rough:.6e},{f_c:.12e},{f_h:.12e},{diff:.9e},"
                        f"{history[-1]:.12e},{len(history)},{converged}")
    (RESULTS / "explicit_comparison.csv").write_text("\n".join(rows), encoding="utf-8")
    ax[1].semilogx([float(r.split(",")[0]) for r in rows[1:]], [float(r.split(",")[4]) for r in rows[1:]],
                   "o", ms=3, color=RED)
    ax[1].set_xlabel("Reynolds number Re [-]")
    ax[1].set_ylabel("|f_Haaland − f_Colebrook| / f_Colebrook [-]")
    ax[1].set_title(f"Implicit solver vs explicit approximation\n(max rel. diff {max_diff:.1%})")
    fig.suptitle("Pipeflow Bench — Colebrook implicit solver and references", fontsize=13)
    fig.tight_layout()
    fig.savefig(RESULTS / "friction_chart.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "friction_chart.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return max_diff


def make_residual_history_figure():
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    rows = ["iteration,residual"]
    histories = {}
    for Re, rel_rough, color in ((4.0e3, 1e-4, BLUE), (1.0e5, 1e-3, GREEN),
                                 (1.0e7, 1e-2, RED)):
        _, history, _ = colebrook_friction(Re, rel_rough)
        histories[(Re, rel_rough)] = history
        ax.semilogy(np.arange(1, len(history) + 1), np.abs(history), "o-",
                    ms=3, color=color, lw=1.2,
                    label=f"Re = {Re:.0e}, ε/D = {rel_rough:g}")
    ax.set_xlabel("bisection iteration [-]")
    ax.set_ylabel("|Colebrook residual g(f)| [-]")
    ax.set_title("Implicit-solver residual convergence")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "residual_history.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "residual_history.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)
    key = (4.0e3, 1e-4)
    for i, g in enumerate(histories[key]):
        rows.append(f"{i + 1},{abs(g):.9e}")
    (RESULTS / "residual_history.csv").write_text("\n".join(rows), encoding="utf-8")
    return histories


def make_pipeline_sweep():
    velocities = np.r_[np.linspace(0.005, 0.08, 40), np.linspace(0.09, 3.0, 60)]
    rows = ["v_m_s,Re,regime,friction_factor,dp_Pa,dp_kPa_per_100m"]
    for v in velocities:
        r = solve_pipe(rho=998.0, mu=1.0e-3, v=v, D=0.05, L=100.0, eps=4.5e-5)
        rows.append(
            f"{v:.6f},{r['Re']:.9e},{r['regime']},{r['friction_factor']:.9e},"
            f"{r['dp_Pa']:.9e},{r['dp_Pa'] / 1000.0:.9e}"
        )
    (RESULTS / "pipeline_sweep.csv").write_text("\n".join(rows), encoding="utf-8")

    fig, axes = plt.subplots(1,2,figsize=(11,4.6))
    data = np.genfromtxt(RESULTS / "pipeline_sweep.csv", delimiter=",", names=True,
                         dtype=None, encoding="utf-8")
    for i, ax in enumerate(axes):
        for regime, color, label in (("laminar", BLUE, "laminar: 64/Re"),
                                     ("transition-band", "#B47B28", "transition: flagged Colebrook"),
                                     ("turbulent", RED, "turbulent: Colebrook")):
            selected = data["regime"] == regime
            if i == 1:
                selected &= data["v_m_s"] <= 0.08
            ax.plot(data["v_m_s"][selected], data["dp_kPa_per_100m"][selected],
                    "o", ms=3, color=color, label=label)
        ax.set_xlabel("Mean velocity [m/s]")
        ax.set_ylabel("Pressure loss [kPa per 100 m]")
        ax.set_title("Full operating range" if i == 0 else "Low-speed regime detail")
        ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Water · D = 50 mm · roughness = 45 µm · length = 100 m",fontsize=12)
    fig.tight_layout()
    fig.savefig(RESULTS / "pipeline_sweep.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "pipeline_sweep.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_gif():
    trace = colebrook_trace(4000., 1e-4)
    if not trace["converged"]:
        raise RuntimeError("animation solve failed requested residual tolerance")
    records = trace["iterations"]
    header = "iteration,lo,hi,f,residual\n"
    (RESULTS / "bisection_trace.csv").write_text(header + "\n".join(
        ",".join(str(row[key]) for key in ("iteration","lo","hi","f","residual"))
        for row in records), encoding="utf-8")
    fig, axes = plt.subplots(1,2,figsize=(11,6),dpi=100)
    fig.subplots_adjust(left=0.08,right=0.97,bottom=0.18,top=0.80,wspace=0.30)
    fs = np.linspace(1e-4,0.2,1000)
    gs = colebrook_residual(fs,4000.,1e-4)
    residuals = np.abs(trace["residual_history"])
    footer = fig.text(.5,.065,"",ha="center",fontsize=10)
    writer = mpl_animation.PillowWriter(fps=8)
    with writer.saving(fig, RESULTS / "root_convergence.gif", dpi=100):
        for i,row in enumerate(records):
            left, right = axes
            left.clear()
            right.clear()
            left.plot(fs,gs,color=GRAY,lw=1.5,label="Colebrook equation")
            left.axhline(0,color="#999999",lw=0.7)
            left.axvline(trace["friction_factor"],color=GREEN,ls="--",lw=1,label="Final root")
            left.plot(row["f"],row["residual"],"o",color=RED,ms=6,label="Evaluated midpoint")
            left.plot([row["lo"],row["hi"]],[-4.2,-4.2],"|-",color="#B47B28",lw=2,ms=10,label="Current bracket")
            left.set(xlim=(0,0.2),ylim=(-5,10),xlabel="Darcy friction factor f",ylabel="Signed residual g(f)")
            left.legend(loc="upper right",fontsize=8,frameon=False)
            right.semilogy(np.arange(1,i+2),np.maximum(residuals[:i+1],1e-15),"o-",ms=3,color=BLUE)
            right.axhline(1e-12,color=RED,ls="--",lw=1,label="Requested tolerance 1e−12")
            right.set(xlim=(0,len(records)+1),ylim=(1e-15,1e3),xlabel="Actual bisection iteration",ylabel="Absolute residual |g(f)|")
            right.legend(loc="upper right",fontsize=8,frameon=False)
            fig.suptitle(f"Colebrook bisection · iteration {i+1}/{len(records)}\nRe = 4,000 · relative roughness = 0.0001",fontsize=12)
            footer.set_text(f"f = {row['f']:.12f}    |g| = {abs(row['residual']):.3e}    bracket width = {row['hi']-row['lo']:.3e}")
            writer.grab_frame()
    plt.close(fig)
    return len(records)


def main():
    started = time.perf_counter()
    max_diff = make_cover()
    histories = make_residual_history_figure()
    make_pipeline_sweep()
    frames = make_gif()

    # Representative operating point diagnostics for the summary.
    op = solve_pipe(rho=998.0, mu=1.0e-3, v=1.5, D=0.05, L=100.0, eps=4.5e-5)
    comparison = np.genfromtxt(RESULTS / "explicit_comparison.csv", delimiter=",", names=True)
    tolerance_rows = []
    for re, rough, tol in ((1e7, 0., 1e-9), (1e7, 0., 1e-12),
                           (1e7, 0., 1e-30), (4000., 1e-5, 1e-30)):
        f, hist, converged = colebrook_friction(re, rough, tol=tol)
        tolerance_rows.append(f"{re},{rough},{tol},{f},{hist[-1]},{len(hist)},{converged}")
    (RESULTS / "stopping_regression.csv").write_text(
        "Re,relative_roughness,tolerance,f,final_residual,iterations,converged\n"+"\n".join(tolerance_rows),encoding="utf-8")
    summary = {
        "wall_time_s":time.perf_counter()-started,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "system": platform.platform(),
        },
        "results": {
            "max_haaland_rel_diff": max_diff,
            "comparison_cases":len(comparison),
            "max_comparison_residual":float(np.max(np.abs(comparison["final_residual"]))),
            "operating_point": {
                "v_m_s": 1.5, "D_m": 0.05, "eps_m": 4.5e-5,
                "Re": op["Re"], "regime": op["regime"],
                "f": op["friction_factor"],
                "dp_Pa_per_100m": op["dp_Pa"],
                "converged": op["converged"],
                "final_residual": op["residual_history"][-1],
            },
        },
        "animation": {"file": "root_convergence.gif", "frames": frames,
                      "content": "bisection convergence on the implicit "
                                 "Colebrook equation (recorded residual "
                                 "history, brackets and evaluated midpoint data)"},
    }
    (RESULTS / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
