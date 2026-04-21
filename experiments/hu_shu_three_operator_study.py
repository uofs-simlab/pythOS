#!/usr/bin/env python3
"""
hu_shu_three_operator_study.py

Screening study for built-in pythOS three-operator splitting methods on the
Hu–Shu-style BGK problem.

Purpose
-------
Compare shortlisted three-operator methods on the same BGK setting used by the
current Hu–Shu driver, but without overwriting the existing script.

We use a three-operator split:
    A : transport in x
    B : exact homogeneous BGK relaxation
    C : identity / null operator

Why include C?
--------------
The goal here is to test the built-in three-operator coefficient tables
(Strang-3, OS32_7op_minLEM-3, PP3_4A-3, Yoshida-3) on the Hu–Shu/BGK problem
with minimal additional modeling assumptions. Using a null third operator is a
clean screening step before attempting a more elaborate genuine three-way kinetic
decomposition.

Outputs
-------
For each eps and dt, the script records:
- L2 error in f vs a fine reference
- L1 errors in rho, u, T
- min(f), max(f)
- positivity flag
- mass drift
- runtime

Examples
--------
python experiments/hu_shu_three_operator_study.py --case table
python experiments/hu_shu_three_operator_study.py --case mixed --dt 2.5e-4
python experiments/hu_shu_three_operator_study.py --case table --eps-list 1e0,1e-1,1e-2,1e-3 --dt 2.5e-4
"""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
from typing import Callable, Dict, List

import numpy as np
import sys

# ---------------------------------------------------------------------
# Repository-local imports
# ---------------------------------------------------------------------
THIS = Path(__file__).resolve()
REPO = THIS.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from fractional_step import fractional_step

# Reuse existing Hu–Shu BGK building blocks from the current driver.
from experiments.hu_shu_expk2_bgk import (
    Grid,
    make_grid,
    moments,
    maxwellian,
    collision_frequency,
    initial_condition,
    constant_eps,
    mixed_regime_eps,
    get_transport_rhs,
    l1_x,
    l2_xv,
)


ArrayLike = np.ndarray


# ---------------------------------------------------------------------
# Exact / analytic subflows
# ---------------------------------------------------------------------
def bgk_homogeneous_exact(
    g: ArrayLike,
    dt_collision: float,
    grid: Grid,
    eps_x: np.ndarray,
    eta_model: str,
) -> ArrayLike:
    """
    Exact homogeneous BGK solve:
      d_t f = (eta/eps) (M[f] - f)

    Since rho,u,T are conserved under the homogeneous BGK flow, M[g] remains
    fixed over the substep, giving:
      f(h) = exp(-eta h / eps) g + (1 - exp(-eta h / eps)) M[g]
    """
    rho, u, T = moments(g, grid)
    M = maxwellian(rho, u, T, grid)
    eta = collision_frequency(rho, T, model=eta_model)
    alpha = np.exp(-eta * dt_collision / np.maximum(eps_x, 1e-14))
    return alpha[:, None] * g + (1.0 - alpha)[:, None] * M


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


def transport_exact_shift(
    f: ArrayLike,
    dt: float,
    grid: Grid,
) -> ArrayLike:
    """
    Analytic transport subflow for:
      f_t + v f_x = 0
    i.e.
      f(x, v, t + dt) = f(x - v dt, v, t)
    using periodic interpolation in x.
    """
    out = np.empty_like(f)
    shift_cells = (grid.v * dt) / grid.dx
    for j in range(grid.Nv):
        out[:, j] = periodic_shift_1d(f[:, j], shift_cells[j])
    return out


def identity_flow(
    f: ArrayLike,
    dt: float,
    grid: Grid,
) -> ArrayLike:
    _ = dt, grid
    return np.array(f, copy=True)


def make_A_flow(grid: Grid) -> Callable[[float, float, ArrayLike], ArrayLike]:
    def flow(_t: float, h: float, f: ArrayLike) -> ArrayLike:
        return transport_exact_shift(f, h, grid)
    return flow


def make_B_flow(grid: Grid, eps_x: np.ndarray, eta_model: str) -> Callable[[float, float, ArrayLike], ArrayLike]:
    def flow(_t: float, h: float, f: ArrayLike) -> ArrayLike:
        return bgk_homogeneous_exact(f, h, grid, eps_x, eta_model)
    return flow


def make_C_flow(grid: Grid) -> Callable[[float, float, ArrayLike], ArrayLike]:
    def flow(_t: float, h: float, f: ArrayLike) -> ArrayLike:
        return identity_flow(f, h, grid)
    return flow


# ---------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------
def total_mass(f: ArrayLike, grid: Grid) -> float:
    return float(np.sum(f) * grid.dx * grid.dv)


def compare_solution(f: ArrayLike, fref: ArrayLike, grid: Grid) -> Dict[str, float]:
    rho, u, T = moments(f, grid)
    rho_ref, u_ref, T_ref = moments(fref, grid)

    return {
        "f_L2": l2_xv(f, fref, grid.dx, grid.dv),
        "rho_L1": l1_x(rho, rho_ref, grid.dx),
        "u_L1": l1_x(u, u_ref, grid.dx),
        "T_L1": l1_x(T, T_ref, grid.dx),
        "min_f": float(np.min(f)),
        "max_f": float(np.max(f)),
    }


def positivity_flag(min_f: float, tol: float = -1e-12) -> bool:
    return bool(min_f >= tol)


# ---------------------------------------------------------------------
# Time evolution
# ---------------------------------------------------------------------
def evolve_three_operator(
    f0: ArrayLike,
    Tfinal: float,
    dt: float,
    grid: Grid,
    eps_x: np.ndarray,
    alpha_name: str,
    eta_model: str,
) -> ArrayLike:
    f = np.array(f0, copy=True)
    t = 0.0

    flows = [
        make_A_flow(grid),
        make_B_flow(grid, eps_x, eta_model),
        make_C_flow(grid),
    ]

    while t < Tfinal - 1e-15:
        dt_step = min(dt, Tfinal - t)
        f = fractional_step(
            functions=flows,
            delta_t=dt_step,
            initial_y=f,
            initial_t=t,
            final_t=t + dt_step,
            alpha=alpha_name,
            methods={(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"},
        )
        f = np.array(f, dtype=float, copy=False)
        t += dt_step

    return f


# ---------------------------------------------------------------------
# Experiment cases
# ---------------------------------------------------------------------
def run_table(args: argparse.Namespace) -> None:
    grid = make_grid(args.Nx, args.Nv, xmax=args.xmax, vmax=args.vmax)
    f0 = initial_condition(grid, args.init)

    eps_list = [float(s.strip()) for s in args.eps_list.split(",")]
    methods = [s.strip() for s in args.methods.split(",")]

    dt = args.dt
    dt_ref = dt / args.refine_factor

    outdir = REPO / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "hu_shu_three_operator_table.csv"

    rows: List[Dict[str, object]] = []

    print("=" * 84)
    print("HU–SHU THREE-OPERATOR EPSILON TABLE STUDY")
    print(f"methods          = {methods}")
    print(f"dt               = {dt:.3e}")
    print(f"reference method = {args.reference_method}")
    print(f"reference dt     = {dt_ref:.3e}")
    print(f"Tfinal           = {args.T:.3f}")
    print("=" * 84)

    for eps_val in eps_list:
        eps_x = constant_eps(grid, eps_val)

        fref = evolve_three_operator(
            f0=f0,
            Tfinal=args.T,
            dt=dt_ref,
            grid=grid,
            eps_x=eps_x,
            alpha_name=args.reference_method,
            eta_model=args.eta_model,
        )

        ref_mass = total_mass(fref, grid)

        for method_name in methods:
            t0 = time.perf_counter()
            f = evolve_three_operator(
                f0=f0,
                Tfinal=args.T,
                dt=dt,
                grid=grid,
                eps_x=eps_x,
                alpha_name=method_name,
                eta_model=args.eta_model,
            )
            elapsed = time.perf_counter() - t0

            errs = compare_solution(f, fref, grid)
            mass = total_mass(f, grid)
            mass_drift = abs(mass - ref_mass) / max(abs(ref_mass), 1e-30)
            positive = positivity_flag(errs["min_f"], tol=args.positivity_tol)

            print(
                f"eps={eps_val:10.3e} | "
                f"method={method_name:18s} | "
                f"f_L2={errs['f_L2']:10.4e} | "
                f"rho_L1={errs['rho_L1']:10.4e} | "
                f"min(f)={errs['min_f']:10.4e} | "
                f"positive={positive!s:5s} | "
                f"runtime={elapsed:8.3f}s"
            )

            rows.append(
                {
                    "case": "table",
                    "eps": eps_val,
                    "method": method_name,
                    "dt": dt,
                    "reference_method": args.reference_method,
                    "reference_dt": dt_ref,
                    "T": args.T,
                    "Nx": args.Nx,
                    "Nv": args.Nv,
                    "eta_model": args.eta_model,
                    "init": args.init,
                    "f_L2": errs["f_L2"],
                    "rho_L1": errs["rho_L1"],
                    "u_L1": errs["u_L1"],
                    "T_L1": errs["T_L1"],
                    "min_f": errs["min_f"],
                    "max_f": errs["max_f"],
                    "positive": int(positive),
                    "mass_rel_drift_vs_ref": mass_drift,
                    "runtime_sec": elapsed,
                }
            )

    with csv_path.open("w", newline="") as fcsv:
        writer = csv.DictWriter(
            fcsv,
            fieldnames=[
                "case",
                "eps",
                "method",
                "dt",
                "reference_method",
                "reference_dt",
                "T",
                "Nx",
                "Nv",
                "eta_model",
                "init",
                "f_L2",
                "rho_L1",
                "u_L1",
                "T_L1",
                "min_f",
                "max_f",
                "positive",
                "mass_rel_drift_vs_ref",
                "runtime_sec",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved: {csv_path}")


def run_mixed(args: argparse.Namespace) -> None:
    grid = make_grid(args.Nx, args.Nv, xmax=args.xmax, vmax=args.vmax)
    f0 = initial_condition(grid, "mixed")
    eps_x = mixed_regime_eps(grid, eps0=args.eps0)
    methods = [s.strip() for s in args.methods.split(",")]

    dt = args.dt
    dt_ref = dt / args.refine_factor

    outdir = REPO / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "hu_shu_three_operator_mixed.csv"

    rows: List[Dict[str, object]] = []

    print("=" * 84)
    print("HU–SHU THREE-OPERATOR MIXED-REGIME STUDY")
    print(f"methods          = {methods}")
    print(f"dt               = {dt:.3e}")
    print(f"reference method = {args.reference_method}")
    print(f"reference dt     = {dt_ref:.3e}")
    print(f"Tfinal           = {args.T:.3f}")
    print(f"eps0             = {args.eps0:.3e}")
    print("=" * 84)

    fref = evolve_three_operator(
        f0=f0,
        Tfinal=args.T,
        dt=dt_ref,
        grid=grid,
        eps_x=eps_x,
        alpha_name=args.reference_method,
        eta_model=args.eta_model,
    )

    ref_mass = total_mass(fref, grid)

    for method_name in methods:
        t0 = time.perf_counter()
        f = evolve_three_operator(
            f0=f0,
            Tfinal=args.T,
            dt=dt,
            grid=grid,
            eps_x=eps_x,
            alpha_name=method_name,
            eta_model=args.eta_model,
        )
        elapsed = time.perf_counter() - t0

        errs = compare_solution(f, fref, grid)
        mass = total_mass(f, grid)
        mass_drift = abs(mass - ref_mass) / max(abs(ref_mass), 1e-30)
        positive = positivity_flag(errs["min_f"], tol=args.positivity_tol)

        rho, u, T = moments(f, grid)
        rho_ref, u_ref, T_ref = moments(fref, grid)

        print(
            f"method={method_name:18s} | "
            f"f_L2={errs['f_L2']:10.4e} | "
            f"rho_L1={errs['rho_L1']:10.4e} | "
            f"min(f)={errs['min_f']:10.4e} | "
            f"positive={positive!s:5s} | "
            f"runtime={elapsed:8.3f}s"
        )

        rows.append(
            {
                "case": "mixed",
                "method": method_name,
                "dt": dt,
                "reference_method": args.reference_method,
                "reference_dt": dt_ref,
                "T": args.T,
                "Nx": args.Nx,
                "Nv": args.Nv,
                "eta_model": args.eta_model,
                "eps0": args.eps0,
                "f_L2": errs["f_L2"],
                "rho_L1": errs["rho_L1"],
                "u_L1": errs["u_L1"],
                "T_L1": errs["T_L1"],
                "min_f": errs["min_f"],
                "max_f": errs["max_f"],
                "positive": int(positive),
                "mass_rel_drift_vs_ref": mass_drift,
                "rho_max_abs_err": float(np.max(np.abs(rho - rho_ref))),
                "u_max_abs_err": float(np.max(np.abs(u - u_ref))),
                "T_max_abs_err": float(np.max(np.abs(T - T_ref))),
                "runtime_sec": elapsed,
            }
        )

    with csv_path.open("w", newline="") as fcsv:
        writer = csv.DictWriter(
            fcsv,
            fieldnames=[
                "case",
                "method",
                "dt",
                "reference_method",
                "reference_dt",
                "T",
                "Nx",
                "Nv",
                "eta_model",
                "eps0",
                "f_L2",
                "rho_L1",
                "u_L1",
                "T_L1",
                "min_f",
                "max_f",
                "positive",
                "mass_rel_drift_vs_ref",
                "rho_max_abs_err",
                "u_max_abs_err",
                "T_max_abs_err",
                "runtime_sec",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved: {csv_path}")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Three-operator screening study on the Hu–Shu BGK problem.")

    p.add_argument("--case", choices=["table", "mixed"], default="table")
    p.add_argument(
        "--methods",
        type=str,
        default="Strang-3,OS32_7op_minLEM-3,PP3_4A-3,Yoshida-3",
        help="Comma-separated pythOS alpha names.",
    )
    p.add_argument(
        "--reference-method",
        type=str,
        default="Yoshida-3",
        help="Reference built-in alpha name.",
    )

    p.add_argument("--eta-model", choices=["constant1", "rho", "rho_sqrtT"], default="constant1")

    p.add_argument("--Nx", type=int, default=80)
    p.add_argument("--Nv", type=int, default=120)
    p.add_argument("--xmax", type=float, default=2.0)
    p.add_argument("--vmax", type=float, default=15.0)
    p.add_argument("--T", type=float, default=0.1)
    p.add_argument("--dt", type=float, default=2.5e-4)
    p.add_argument("--refine-factor", type=int, default=4)

    p.add_argument("--eps-list", type=str, default="1e0,1e-1,1e-2,1e-3")
    p.add_argument("--eps0", type=float, default=1e-5)
    p.add_argument("--init", choices=["accuracy", "paper", "mixed", "simple"], default="accuracy")

    p.add_argument("--positivity-tol", type=float, default=-1e-12)
    p.add_argument("--outdir", type=str, default="experiments/hu_shu_three_operator_outputs")

    return p


def main() -> None:
    args = build_parser().parse_args()

    print("Running Hu–Shu three-operator study with:")
    print(f"  case             = {args.case}")
    print(f"  methods          = {args.methods}")
    print(f"  reference_method = {args.reference_method}")
    print(f"  eta_model        = {args.eta_model}")
    print(f"  Nx, Nv           = {args.Nx}, {args.Nv}")
    print(f"  xmax, vmax       = {args.xmax}, {args.vmax}")
    print(f"  T                = {args.T}")
    print(f"  dt               = {args.dt}")
    print(f"  refine_factor    = {args.refine_factor}")
    print()

    if args.case == "table":
        run_table(args)
    elif args.case == "mixed":
        run_mixed(args)
    else:
        raise ValueError(f"Unknown case: {args.case}")


if __name__ == "__main__":
    main()