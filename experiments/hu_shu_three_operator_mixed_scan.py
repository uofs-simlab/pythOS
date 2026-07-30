#!/usr/bin/env python3
"""
hu_shu_three_operator_mixed_scan.py

Focused timestep scan for the mixed-regime Hu–Shu three-operator screening study.

Purpose
-------
Starting from the mixed-regime failure observed at dt = 2.5e-4, this script
scans smaller timesteps to determine when each shortlisted built-in three-
operator method becomes usable again.

We keep the same three-operator screening decomposition used in
hu_shu_three_operator_study.py:
    A : transport in x
    B : exact homogeneous BGK relaxation
    C : identity / null operator

Outputs
-------
For each method and dt:
- L2 error in f vs reference
- L1 errors in rho, u, T
- min(f), max(f)
- positivity flag
- finite/nonfinite flag
- mass drift vs reference
- runtime

Examples
--------
python experiments/hu_shu_three_operator_mixed_scan.py
python experiments/hu_shu_three_operator_mixed_scan.py --dts 2.5e-4,1.25e-4,6.25e-5,3.125e-5
python experiments/hu_shu_three_operator_mixed_scan.py --reference-method Strang-3 --refine-factor 8
"""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
from typing import Dict, List

import numpy as np
import sys

# ---------------------------------------------------------------------
# Robust repo-root imports
# ---------------------------------------------------------------------
THIS = Path(__file__).resolve()
REPO = THIS.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from experiments.hu_shu_three_operator_study import (
    evolve_three_operator,
    compare_solution,
    total_mass,
    positivity_flag,
)
from experiments.hu_shu_expk2_bgk import (
    make_grid,
    initial_condition,
    mixed_regime_eps,
)

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def parse_csv_floats(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_csv_strings(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


# ---------------------------------------------------------------------
# Main mixed scan
# ---------------------------------------------------------------------
def run_scan(args: argparse.Namespace) -> None:
    grid = make_grid(args.Nx, args.Nv, xmax=args.xmax, vmax=args.vmax)
    f0 = initial_condition(grid, "mixed")
    eps_x = mixed_regime_eps(grid, eps0=args.eps0)

    methods = parse_csv_strings(args.methods)
    dts = parse_csv_floats(args.dts)

    outdir = REPO / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "hu_shu_three_operator_mixed_scan.csv"

    rows: List[Dict[str, object]] = []

    dt_min = min(dts)
    dt_ref = dt_min / args.refine_factor

    print("=" * 92)
    print("HU–SHU THREE-OPERATOR MIXED-REGIME TIMESTEP SCAN")
    print(f"methods          = {methods}")
    print(f"dts              = {dts}")
    print(f"reference method = {args.reference_method}")
    print(f"reference dt     = {dt_ref:.6e}")
    print(f"Tfinal           = {args.T:.6f}")
    print(f"eps0             = {args.eps0:.3e}")
    print("=" * 92)

    # One fine reference for the mixed regime
    print("\nComputing reference solution...")
    t_ref0 = time.perf_counter()
    fref = evolve_three_operator(
        f0=f0,
        Tfinal=args.T,
        dt=dt_ref,
        grid=grid,
        eps_x=eps_x,
        alpha_name=args.reference_method,
        eta_model=args.eta_model,
    )
    ref_elapsed = time.perf_counter() - t_ref0
    ref_mass = total_mass(fref, grid)
    print(f"Reference done in {ref_elapsed:.3f}s")

    # Scan all methods / dts
    for method_name in methods:
        print(f"\nMethod: {method_name}")

        for dt in dts:
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

            finite = bool(np.all(np.isfinite(f)))

            if finite:
                errs = compare_solution(f, fref, grid)
                mass = total_mass(f, grid)
                mass_drift = abs(mass - ref_mass) / max(abs(ref_mass), 1e-30)
                positive = positivity_flag(errs["min_f"], tol=args.positivity_tol)
            else:
                errs = {
                    "f_L2": float("inf"),
                    "rho_L1": float("inf"),
                    "u_L1": float("inf"),
                    "T_L1": float("inf"),
                    "min_f": float("nan"),
                    "max_f": float("nan"),
                }
                mass = float("nan")
                mass_drift = float("inf")
                positive = False

            usable = bool(
                finite
                and positive
                and np.isfinite(errs["f_L2"])
                and errs["f_L2"] <= args.usable_f_l2_threshold
            )

            print(
                f"  dt={dt:10.4e} | "
                f"finite={finite!s:5s} | "
                f"positive={positive!s:5s} | "
                f"usable={usable!s:5s} | "
                f"f_L2={errs['f_L2']:12.5e} | "
                f"rho_L1={errs['rho_L1']:12.5e} | "
                f"min(f)={errs['min_f']:12.5e} | "
                f"runtime={elapsed:8.3f}s"
            )

            rows.append(
                {
                    "case": "mixed_scan",
                    "method": method_name,
                    "dt": dt,
                    "reference_method": args.reference_method,
                    "reference_dt": dt_ref,
                    "T": args.T,
                    "Nx": args.Nx,
                    "Nv": args.Nv,
                    "eta_model": args.eta_model,
                    "eps0": args.eps0,
                    "finite": int(finite),
                    "positive": int(positive),
                    "usable": int(usable),
                    "f_L2": errs["f_L2"],
                    "rho_L1": errs["rho_L1"],
                    "u_L1": errs["u_L1"],
                    "T_L1": errs["T_L1"],
                    "min_f": errs["min_f"],
                    "max_f": errs["max_f"],
                    "final_mass": mass,
                    "mass_rel_drift_vs_ref": mass_drift,
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
                "finite",
                "positive",
                "usable",
                "f_L2",
                "rho_L1",
                "u_L1",
                "T_L1",
                "min_f",
                "max_f",
                "final_mass",
                "mass_rel_drift_vs_ref",
                "runtime_sec",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved: {csv_path}")

    # Console summary
    print("\nSummary by method:")
    for method_name in methods:
        subset = [r for r in rows if r["method"] == method_name]
        largest_finite = None
        largest_positive = None
        largest_usable = None

        for r in subset:
            dt = float(r["dt"])
            if int(r["finite"]) == 1:
                largest_finite = dt if largest_finite is None else max(largest_finite, dt)
            if int(r["positive"]) == 1:
                largest_positive = dt if largest_positive is None else max(largest_positive, dt)
            if int(r["usable"]) == 1:
                largest_usable = dt if largest_usable is None else max(largest_usable, dt)

        print(
            f"  {method_name:20s} | "
            f"largest_finite_dt={largest_finite} | "
            f"largest_positive_dt={largest_positive} | "
            f"largest_usable_dt={largest_usable}"
        )


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Mixed-regime timestep scan for Hu–Shu three-operator methods."
    )

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
    p.add_argument(
        "--dts",
        type=str,
        default="2.5e-4,1.25e-4,6.25e-5,3.125e-5",
        help="Comma-separated dt values to scan.",
    )

    p.add_argument("--eta-model", choices=["constant1", "rho", "rho_sqrtT"], default="constant1")

    p.add_argument("--Nx", type=int, default=80)
    p.add_argument("--Nv", type=int, default=120)
    p.add_argument("--xmax", type=float, default=2.0)
    p.add_argument("--vmax", type=float, default=15.0)
    p.add_argument("--T", type=float, default=0.1)
    p.add_argument("--refine-factor", type=int, default=4)

    p.add_argument("--eps0", type=float, default=1e-5)
    p.add_argument("--positivity-tol", type=float, default=-1e-12)
    p.add_argument(
        "--usable-f-l2-threshold",
        type=float,
        default=1e-1,
        help="Simple threshold used for a practical 'usable' flag.",
    )
    p.add_argument("--outdir", type=str, default="experiments/hu_shu_three_operator_outputs")

    return p


def main() -> None:
    args = build_parser().parse_args()

    print("Running Hu–Shu mixed-regime three-operator scan with:")
    print(f"  methods              = {args.methods}")
    print(f"  reference_method     = {args.reference_method}")
    print(f"  dts                  = {args.dts}")
    print(f"  eta_model            = {args.eta_model}")
    print(f"  Nx, Nv               = {args.Nx}, {args.Nv}")
    print(f"  xmax, vmax           = {args.xmax}, {args.vmax}")
    print(f"  T                    = {args.T}")
    print(f"  refine_factor        = {args.refine_factor}")
    print(f"  eps0                 = {args.eps0}")
    print(f"  usable_f_l2_threshold= {args.usable_f_l2_threshold}")
    print()

    run_scan(args)


if __name__ == "__main__":
    main()