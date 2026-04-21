"""
bgk_three_operator_stability_scan.py

Stability / positivity scan for the BGK-inspired three-operator experiment.

Purpose
-------
This script refines the Phase 2 study by scanning a list of dt values for a
fixed stiff regime, typically eps = 1e-3, and recording:

- stable / unstable
- positive / not positive
- minimum value of f over the run
- final minimum of f
- mass drift
- L1 error vs fine reference
- max L1 growth ratio

Model
-----
We split

    f_t + v f_x + F f_v = ( M[f] - f ) / eps

into three operators:
    A: f_t = - v f_x
    B: f_t = ( M[f] - f ) / eps
    C: f_t = - F f_v

using pythOS fractional_step(...) with built-in 3-operator alpha tables.

Recommended use
---------------
Start with:
    eps = 1e-3
    methods = Strang-3, OS32_7op_minLEM-3, PP3_4A-3, Yoshida-3

Examples
--------
python experiments/bgk_three_operator_stability_scan.py
python experiments/bgk_three_operator_stability_scan.py --eps 1e-3 --T 0.2
python experiments/bgk_three_operator_stability_scan.py --dts 1e-2,8e-3,6e-3,5e-3,4e-3,3e-3,2.5e-3,2e-3,1.5e-3,1e-3
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Sequence

import numpy as np

# ---------------------------------------------------------------------
# Robust repo-root imports
# ---------------------------------------------------------------------
THIS = Path(__file__).resolve()
REPO = THIS.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from fractional_step import fractional_step  # noqa: E402


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def parse_csv_floats(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_csv_strings(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def periodic_shift_1d(values: np.ndarray, shift_cells: float) -> np.ndarray:
    """
    Periodic linear interpolation for a 1D array.
    Positive shift_cells means evaluate values(x - shift).
    """
    n = values.shape[0]
    idx = np.arange(n, dtype=float)
    src = (idx - shift_cells) % n
    i0 = np.floor(src).astype(int)
    i1 = (i0 + 1) % n
    theta = src - i0
    return (1.0 - theta) * values[i0] + theta * values[i1]


def shift_along_x_periodic(f: np.ndarray, shift_cells_per_v: np.ndarray) -> np.ndarray:
    Nx, Nv = f.shape
    out = np.empty_like(f)
    for j in range(Nv):
        out[:, j] = periodic_shift_1d(f[:, j], shift_cells_per_v[j])
    return out


def shift_along_v_periodic(f: np.ndarray, shift_cells: float) -> np.ndarray:
    Nx, Nv = f.shape
    out = np.empty_like(f)
    for i in range(Nx):
        out[i, :] = periodic_shift_1d(f[i, :], shift_cells)
    return out


# ---------------------------------------------------------------------
# Moments / Maxwellian
# ---------------------------------------------------------------------
def compute_moments(
    f: np.ndarray,
    v: np.ndarray,
    dv: float,
    rho_floor: float = 1e-14,
    T_floor: float = 1e-12,
):
    rho = np.sum(f, axis=1) * dv
    rho_safe = np.maximum(rho, rho_floor)

    mom1 = np.sum(f * v[None, :], axis=1) * dv
    u = mom1 / rho_safe

    centered = v[None, :] - u[:, None]
    energy = np.sum(f * centered * centered, axis=1) * dv
    T = np.maximum(energy / rho_safe, T_floor)

    return rho, u, T


def maxwellian_from_moments(
    rho: np.ndarray,
    u: np.ndarray,
    T: np.ndarray,
    v: np.ndarray,
    rho_floor: float = 1e-14,
    T_floor: float = 1e-12,
) -> np.ndarray:
    rho_safe = np.maximum(rho, rho_floor)
    T_safe = np.maximum(T, T_floor)
    pref = rho_safe[:, None] / np.sqrt(2.0 * np.pi * T_safe[:, None])
    z = (v[None, :] - u[:, None]) ** 2 / (2.0 * T_safe[:, None])
    return pref * np.exp(-z)


# ---------------------------------------------------------------------
# Initial condition
# ---------------------------------------------------------------------
def initial_condition(x: np.ndarray, v: np.ndarray, profile: str = "mixed") -> np.ndarray:
    X = x[:, None]
    V = v[None, :]

    if profile == "mixed":
        rho = 1.0 + 0.2 * np.sin(2.0 * np.pi * X)
        u = 0.5 * np.cos(2.0 * np.pi * X)
        T = 0.8 + 0.1 * np.sin(4.0 * np.pi * X)

        M1 = rho / np.sqrt(2.0 * np.pi * T) * np.exp(-((V - u) ** 2) / (2.0 * T))
        M2 = 0.35 * rho / np.sqrt(2.0 * np.pi * T) * np.exp(-((V + 0.7 * u) ** 2) / (2.0 * (1.25 * T)))
        f0 = 0.75 * M1 + 0.25 * M2
        return np.maximum(f0, 1e-14)

    if profile == "simple":
        rho = 1.0 + 0.1 * np.sin(2.0 * np.pi * X)
        u = 0.25 * np.cos(2.0 * np.pi * X)
        T = 1.0 + 0.05 * np.sin(2.0 * np.pi * X)
        f0 = rho / np.sqrt(2.0 * np.pi * T) * np.exp(-((V - u) ** 2) / (2.0 * T))
        return np.maximum(f0, 1e-14)

    raise ValueError("profile must be 'mixed' or 'simple'")


# ---------------------------------------------------------------------
# Problem definition
# ---------------------------------------------------------------------
@dataclass
class BGKPhase2Problem:
    x: np.ndarray
    v: np.ndarray
    dx: float
    dv: float
    force: float
    eps: float


def build_problem(
    Nx: int,
    Nv: int,
    xmax: float,
    vmax: float,
    force: float,
    eps: float,
) -> BGKPhase2Problem:
    x = np.linspace(0.0, xmax, Nx, endpoint=False)
    v = np.linspace(-vmax, vmax, Nv, endpoint=False)
    dx = xmax / Nx
    dv = 2.0 * vmax / Nv
    return BGKPhase2Problem(x=x, v=v, dx=dx, dv=dv, force=force, eps=eps)


# ---------------------------------------------------------------------
# Analytic subflows
# ---------------------------------------------------------------------
def make_transport_flow(problem: BGKPhase2Problem) -> Callable[[float, float, np.ndarray], np.ndarray]:
    v = problem.v
    dx = problem.dx

    def flow(_t: float, h: float, f: np.ndarray) -> np.ndarray:
        shift_cells = (v * h) / dx
        return shift_along_x_periodic(f, shift_cells)

    return flow


def make_force_flow(problem: BGKPhase2Problem) -> Callable[[float, float, np.ndarray], np.ndarray]:
    dv = problem.dv
    F = problem.force

    def flow(_t: float, h: float, f: np.ndarray) -> np.ndarray:
        shift_cells = (F * h) / dv
        return shift_along_v_periodic(f, shift_cells)

    return flow


def make_bgk_relaxation_flow(problem: BGKPhase2Problem) -> Callable[[float, float, np.ndarray], np.ndarray]:
    v = problem.v
    dv = problem.dv
    eps = problem.eps

    def flow(_t: float, h: float, f: np.ndarray) -> np.ndarray:
        rho, u, T = compute_moments(f, v, dv)
        M = maxwellian_from_moments(rho, u, T, v)
        E = np.exp(-h / eps)
        return E * f + (1.0 - E) * M

    return flow


# ---------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------
def total_mass(f: np.ndarray, dx: float, dv: float) -> float:
    return float(np.sum(f) * dx * dv)


def l1_norm(f: np.ndarray, dx: float, dv: float) -> float:
    return float(np.sum(np.abs(f)) * dx * dv)


@dataclass
class RunMetrics:
    stable: bool
    positive: bool
    err_l1: float
    min_f_over_run: float
    final_min_f: float
    final_mass: float
    mass_rel_drift: float
    max_l1_ratio: float
    exploded_at_step: int
    dt_used: float
    steps: int


def run_method(
    method_name: str,
    dt: float,
    T: float,
    f0: np.ndarray,
    problem: BGKPhase2Problem,
    reference: np.ndarray,
    positivity_tol: float = -1e-12,
    explosion_factor: float = 1e8,
) -> RunMetrics:
    flows = [
        make_transport_flow(problem),
        make_bgk_relaxation_flow(problem),
        make_force_flow(problem),
    ]

    steps = max(1, int(round(T / dt)))
    dt = T / steps

    f = np.array(f0, dtype=float, copy=True)
    m0 = total_mass(f0, problem.dx, problem.dv)
    l10 = l1_norm(f0, problem.dx, problem.dv)

    min_f_over_run = float(np.min(f))
    max_l1_ratio = 1.0
    stable = True
    exploded_at_step = -1
    t = 0.0

    for n in range(steps):
        f_new = fractional_step(
            functions=flows,
            delta_t=dt,
            initial_y=f,
            initial_t=t,
            final_t=t + dt,
            alpha=method_name,
            methods={(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"},
        )
        f_new = np.array(f_new, dtype=float, copy=False)

        if not np.all(np.isfinite(f_new)):
            stable = False
            exploded_at_step = n + 1
            f = f_new
            break

        min_f_over_run = min(min_f_over_run, float(np.min(f_new)))

        l1r = l1_norm(f_new, problem.dx, problem.dv) / max(l10, 1e-30)
        max_l1_ratio = max(max_l1_ratio, l1r)

        if l1r > explosion_factor:
            stable = False
            exploded_at_step = n + 1
            f = f_new
            break

        f = f_new
        t += dt

    if stable and np.all(np.isfinite(f)):
        err_l1 = l1_norm(f - reference, problem.dx, problem.dv)
        final_min_f = float(np.min(f))
        final_mass = total_mass(f, problem.dx, problem.dv)
    else:
        err_l1 = float("inf")
        final_min_f = float("nan")
        final_mass = float("nan")

    mass_rel_drift = abs(final_mass - m0) / max(abs(m0), 1e-30) if np.isfinite(final_mass) else float("inf")
    positive = bool(np.isfinite(final_min_f) and final_min_f >= positivity_tol)

    return RunMetrics(
        stable=stable,
        positive=positive,
        err_l1=err_l1,
        min_f_over_run=min_f_over_run,
        final_min_f=final_min_f,
        final_mass=final_mass,
        mass_rel_drift=mass_rel_drift,
        max_l1_ratio=max_l1_ratio,
        exploded_at_step=exploded_at_step,
        dt_used=dt,
        steps=steps,
    )


def compute_reference(
    ref_method: str,
    ref_dt: float,
    T: float,
    f0: np.ndarray,
    problem: BGKPhase2Problem,
) -> np.ndarray:
    flows = [
        make_transport_flow(problem),
        make_bgk_relaxation_flow(problem),
        make_force_flow(problem),
    ]
    f_ref = fractional_step(
        functions=flows,
        delta_t=ref_dt,
        initial_y=f0,
        initial_t=0.0,
        final_t=T,
        alpha=ref_method,
        methods={(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"},
    )
    return np.array(f_ref, dtype=float)


# ---------------------------------------------------------------------
# Scan runner
# ---------------------------------------------------------------------
def run_scan(
    eps: float,
    T: float,
    dts: Sequence[float],
    methods: Sequence[str],
    problem_base: Dict[str, float],
    profile: str,
    ref_method: str,
    ref_refine: int,
    outdir: Path,
    positivity_tol: float,
    explosion_factor: float,
) -> None:
    problem = build_problem(
        Nx=int(problem_base["Nx"]),
        Nv=int(problem_base["Nv"]),
        xmax=float(problem_base["xmax"]),
        vmax=float(problem_base["vmax"]),
        force=float(problem_base["force"]),
        eps=eps,
    )

    f0 = initial_condition(problem.x, problem.v, profile=profile)

    dt_min = min(dts)
    ref_dt = dt_min / ref_refine
    reference = compute_reference(ref_method, ref_dt, T, f0, problem)

    outdir.mkdir(parents=True, exist_ok=True)
    detail_path = outdir / f"bgk_three_operator_scan_eps{str(eps).replace('.', 'p').replace('-', 'm')}_detail.csv"
    summary_path = outdir / f"bgk_three_operator_scan_eps{str(eps).replace('.', 'p').replace('-', 'm')}_summary.csv"

    detail_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []

    print("\n" + "=" * 92)
    print(f"BGK THREE-OPERATOR STABILITY SCAN | eps = {eps:.3e} | T = {T:g}")
    print("=" * 92)
    print(f"Reference method = {ref_method} | ref_dt = {ref_dt:.4e}")

    for method_name in methods:
        print(f"\nMethod: {method_name}")

        largest_stable_dt = float("nan")
        largest_positive_dt = float("nan")
        stable_count = 0
        positive_count = 0

        for dt in dts:
            metrics = run_method(
                method_name=method_name,
                dt=dt,
                T=T,
                f0=f0,
                problem=problem,
                reference=reference,
                positivity_tol=positivity_tol,
                explosion_factor=explosion_factor,
            )

            print(
                f"  dt_req={dt:10.4e} | "
                f"dt_used={metrics.dt_used:10.4e} | "
                f"stable={metrics.stable!s:5s} | "
                f"positive={metrics.positive!s:5s} | "
                f"err={metrics.err_l1:12.5e} | "
                f"min_f_run={metrics.min_f_over_run:12.5e} | "
                f"mass_drift={metrics.mass_rel_drift:10.4e}"
            )

            detail_rows.append(
                {
                    "eps": eps,
                    "T": T,
                    "method": method_name,
                    "dt_requested": dt,
                    "dt_used": metrics.dt_used,
                    "steps": metrics.steps,
                    "ref_method": ref_method,
                    "ref_dt": ref_dt,
                    "Nx": int(problem_base["Nx"]),
                    "Nv": int(problem_base["Nv"]),
                    "xmax": float(problem_base["xmax"]),
                    "vmax": float(problem_base["vmax"]),
                    "force": float(problem_base["force"]),
                    "profile": profile,
                    "stable": int(metrics.stable),
                    "positive": int(metrics.positive),
                    "err_l1": metrics.err_l1,
                    "min_f_over_run": metrics.min_f_over_run,
                    "final_min_f": metrics.final_min_f,
                    "final_mass": metrics.final_mass,
                    "mass_rel_drift": metrics.mass_rel_drift,
                    "max_l1_ratio": metrics.max_l1_ratio,
                    "exploded_at_step": metrics.exploded_at_step,
                }
            )

            if metrics.stable:
                stable_count += 1
                if np.isnan(largest_stable_dt) or metrics.dt_used > largest_stable_dt:
                    largest_stable_dt = metrics.dt_used

            if metrics.stable and metrics.positive:
                positive_count += 1
                if np.isnan(largest_positive_dt) or metrics.dt_used > largest_positive_dt:
                    largest_positive_dt = metrics.dt_used

        summary_rows.append(
            {
                "eps": eps,
                "method": method_name,
                "stable_count": stable_count,
                "positive_count": positive_count,
                "num_dts_tested": len(dts),
                "largest_stable_dt": largest_stable_dt,
                "largest_positive_dt": largest_positive_dt,
            }
        )

    with detail_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "eps",
                "T",
                "method",
                "dt_requested",
                "dt_used",
                "steps",
                "ref_method",
                "ref_dt",
                "Nx",
                "Nv",
                "xmax",
                "vmax",
                "force",
                "profile",
                "stable",
                "positive",
                "err_l1",
                "min_f_over_run",
                "final_min_f",
                "final_mass",
                "mass_rel_drift",
                "max_l1_ratio",
                "exploded_at_step",
            ],
        )
        writer.writeheader()
        writer.writerows(detail_rows)

    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "eps",
                "method",
                "stable_count",
                "positive_count",
                "num_dts_tested",
                "largest_stable_dt",
                "largest_positive_dt",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print("\nSummary:")
    for row in summary_rows:
        print(
            f"  {str(row['method']):20s} | "
            f"largest_stable_dt={row['largest_stable_dt']} | "
            f"largest_positive_dt={row['largest_positive_dt']}"
        )

    print(f"\nSaved detail CSV:  {detail_path}")
    print(f"Saved summary CSV: {summary_path}")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="BGK-inspired three-operator stability and positivity scan.")

    ap.add_argument("--eps", type=float, default=1e-3)
    ap.add_argument("--dts", type=str, default="1e-2,8e-3,6e-3,5e-3,4e-3,3e-3,2.5e-3,2e-3,1.5e-3,1e-3")
    ap.add_argument(
        "--methods",
        type=str,
        default="Strang-3,OS32_7op_minLEM-3,PP3_4A-3,Yoshida-3",
    )
    ap.add_argument("--reference-method", type=str, default="Yoshida-3")
    ap.add_argument("--reference-refine", type=int, default=8)
    ap.add_argument("--T", type=float, default=0.2)

    ap.add_argument("--Nx", type=int, default=40)
    ap.add_argument("--Nv", type=int, default=80)
    ap.add_argument("--xmax", type=float, default=1.0)
    ap.add_argument("--vmax", type=float, default=8.0)
    ap.add_argument("--force", type=float, default=0.5)

    ap.add_argument("--profile", type=str, choices=["mixed", "simple"], default="mixed")
    ap.add_argument("--positivity-tol", type=float, default=-1e-12)
    ap.add_argument("--explosion-factor", type=float, default=1e8)
    ap.add_argument("--outdir", type=str, default="experiments/bgk_three_operator_scan_outputs")

    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    dts = parse_csv_floats(args.dts)
    methods = parse_csv_strings(args.methods)

    problem_base = {
        "Nx": args.Nx,
        "Nv": args.Nv,
        "xmax": args.xmax,
        "vmax": args.vmax,
        "force": args.force,
    }

    outdir = REPO / args.outdir

    print("Running BGK three-operator stability scan with:")
    print(f"  eps              = {args.eps}")
    print(f"  dts              = {dts}")
    print(f"  methods          = {methods}")
    print(f"  reference_method = {args.reference_method}")
    print(f"  reference_refine = {args.reference_refine}")
    print(f"  T                = {args.T}")
    print(f"  Nx, Nv           = {args.Nx}, {args.Nv}")
    print(f"  xmax, vmax       = {args.xmax}, {args.vmax}")
    print(f"  force            = {args.force}")
    print(f"  profile          = {args.profile}")
    print(f"  positivity_tol   = {args.positivity_tol}")
    print(f"  explosion_factor = {args.explosion_factor}")
    print(f"  outdir           = {outdir}")

    run_scan(
        eps=args.eps,
        T=args.T,
        dts=dts,
        methods=methods,
        problem_base=problem_base,
        profile=args.profile,
        ref_method=args.reference_method,
        ref_refine=args.reference_refine,
        outdir=outdir,
        positivity_tol=args.positivity_tol,
        explosion_factor=args.explosion_factor,
    )


if __name__ == "__main__":
    main()