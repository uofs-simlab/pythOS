"""
abc_higher_order_stability_scan.py

Phase 1C stability-window scan for higher-order 3-operator splitting methods.

Purpose
-------
For each method and each stiffness regime, scan over a list of dt values and
record whether the method remains stable under a chosen criterion.

This is designed to answer:
- Which methods have the largest stable timestep window?
- Do higher-order methods lose robustness earlier than Strang-type methods?
- How does that depend on the problem mode (decay vs growth) and stiffness?

Model
-----
We solve y' = (A + B + C) y on R^2 with analytic subflows for A, B, C, using
pythOS fractional_step(...) and built-in alpha tables.

Modes
-----
1. decay:
   stiff dissipative B
2. growth:
   one expanding and one contracting stiff direction

Outputs
-------
For each (mode, Bmult, method, dt), writes:
- stable flag
- final error (if stable)
- max ||y_n|| / ||y_0||
- final ||y_T|| / ||y_0||
- max step ratio ||y_{n+1}|| / ||y_n||
- exploded_at_step
- largest stable dt per method/Bmult summary

Examples
--------
python experiments/abc_higher_order_stability_scan.py
python experiments/abc_higher_order_stability_scan.py --mode decay --Bmults 20,50,100
python experiments/abc_higher_order_stability_scan.py --mode growth --Bmults 20,50
python experiments/abc_higher_order_stability_scan.py --dtmax 2e-2 --levels 10
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------
# Repository-local imports
# ---------------------------------------------------------------------
THIS = Path(__file__).resolve()
REPO = THIS.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from fractional_step import fractional_step  # noqa: E402

try:
    from scipy.linalg import expm  # type: ignore
except Exception:
    expm = None


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------
def parse_csv_floats(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_csv_strings(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def norm2(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def safe_ratio(a: float, b: float) -> float:
    if b == 0.0:
        return float("inf") if a != 0.0 else 1.0
    return a / b


def expm_fallback(M: np.ndarray) -> np.ndarray:
    w, V = np.linalg.eig(M)
    Vinv = np.linalg.inv(V)
    E = V @ np.diag(np.exp(w)) @ Vinv
    if np.max(np.abs(E.imag)) < 1e-12:
        return E.real
    return E


def mat_expm(M: np.ndarray) -> np.ndarray:
    if expm is not None:
        return expm(M)
    return expm_fallback(M)


def geometric_dt_list(dtmin: float, dtmax: float, levels: int) -> List[float]:
    if levels < 2:
        return [dtmax]
    vals = np.geomspace(dtmax, dtmin, num=levels)
    return [float(v) for v in vals]


# ---------------------------------------------------------------------
# Problem definition
# ---------------------------------------------------------------------
def make_abc_matrices(Bmult: float, mode: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    A = np.array(
        [
            [0.0, 2.0],
            [-2.0, 0.0],
        ],
        dtype=float,
    )

    if mode == "decay":
        B = np.array(
            [
                [-1.0, 0.0],
                [0.0, -0.05],
            ],
            dtype=float,
        ) * Bmult
    elif mode == "growth":
        B = np.array(
            [
                [1.0, 0.0],
                [0.0, -0.10],
            ],
            dtype=float,
        ) * Bmult
    else:
        raise ValueError("mode must be 'decay' or 'growth'")

    C = np.array(
        [
            [0.0, 1.5],
            [0.4, 0.0],
        ],
        dtype=float,
    )

    return A, B, C


def exact_solution(y0: np.ndarray, T: float, A: np.ndarray, B: np.ndarray, C: np.ndarray) -> np.ndarray:
    return mat_expm(T * (A + B + C)) @ y0


def make_linear_flow(L: np.ndarray) -> Callable[[float, float, np.ndarray], np.ndarray]:
    def flow(_t: float, h: float, y: np.ndarray) -> np.ndarray:
        return mat_expm(h * L) @ y
    return flow


# ---------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------
@dataclass
class RunMetrics:
    err: float
    stable: bool
    max_norm_ratio: float
    final_norm_ratio: float
    max_step_ratio: float
    final_norm: float
    exploded_at_step: int
    steps: int
    dt_used: float


def advance_one_step(
    y: np.ndarray,
    t: float,
    dt: float,
    flows: List[Callable[[float, float, np.ndarray], np.ndarray]],
    alpha_name: str,
) -> np.ndarray:
    y1 = fractional_step(
        functions=flows,
        delta_t=dt,
        initial_y=np.array(y, dtype=float),
        initial_t=t,
        final_t=t + dt,
        alpha=alpha_name,
        methods={(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"},
    )
    return np.array(y1, dtype=np.complex128 if np.iscomplexobj(y1) else float)


def simulate_method(
    alpha_name: str,
    dt: float,
    T: float,
    y0: np.ndarray,
    flows: List[Callable[[float, float, np.ndarray], np.ndarray]],
    y_exact: np.ndarray,
    explosion_factor: float,
) -> RunMetrics:
    steps = max(1, int(round(T / dt)))
    dt = T / steps

    y = np.array(y0, dtype=float)
    t = 0.0

    n0 = norm2(y0)
    max_norm_ratio = 1.0
    max_step_ratio = 1.0
    stable = True
    exploded_at_step = -1

    for n in range(steps):
        y_new = advance_one_step(y, t, dt, flows, alpha_name)

        if not np.all(np.isfinite(y_new)):
            stable = False
            exploded_at_step = n + 1
            y = y_new
            break

        ny = norm2(np.asarray(y, dtype=np.complex128))
        ny_new = norm2(np.asarray(y_new, dtype=np.complex128))

        if ny > 0.0:
            max_step_ratio = max(max_step_ratio, ny_new / ny)
        else:
            max_step_ratio = float("inf")

        if n0 > 0.0:
            max_norm_ratio = max(max_norm_ratio, ny_new / n0)

        if ny_new > explosion_factor * max(1.0, n0):
            stable = False
            exploded_at_step = n + 1
            y = y_new
            break

        y = np.array(y_new, copy=True)
        t += dt

    final_norm = norm2(np.asarray(y, dtype=np.complex128))
    final_norm_ratio = safe_ratio(final_norm, n0)

    if stable and np.all(np.isfinite(y)):
        err = norm2(np.asarray(y, dtype=np.complex128) - np.asarray(y_exact, dtype=np.complex128))
    else:
        err = float("inf")

    return RunMetrics(
        err=err,
        stable=stable,
        max_norm_ratio=max_norm_ratio,
        final_norm_ratio=final_norm_ratio,
        max_step_ratio=max_step_ratio,
        final_norm=final_norm,
        exploded_at_step=exploded_at_step,
        steps=steps,
        dt_used=dt,
    )


# ---------------------------------------------------------------------
# Scan runner
# ---------------------------------------------------------------------
def run_scan(
    mode: str,
    Bmults: Sequence[float],
    dts: Sequence[float],
    methods: Sequence[str],
    T: float,
    y0: np.ndarray,
    outdir: Path,
    explosion_factor: float,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    detail_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []

    for Bmult in Bmults:
        A, B, C = make_abc_matrices(Bmult, mode)
        flows = [make_linear_flow(A), make_linear_flow(B), make_linear_flow(C)]
        y_exact = exact_solution(y0, T, A, B, C)

        print("\n" + "=" * 80)
        print(f"STABILITY SCAN | mode = {mode} | Bmult = {Bmult:g} | T = {T:.6g}")
        print("=" * 80)

        for method_name in methods:
            print(f"\nMethod: {method_name}")

            largest_stable_dt = float("nan")
            smallest_unstable_dt = float("nan")
            stable_count = 0

            for dt in dts:
                metrics = simulate_method(
                    alpha_name=method_name,
                    dt=dt,
                    T=T,
                    y0=y0,
                    flows=flows,
                    y_exact=y_exact,
                    explosion_factor=explosion_factor,
                )

                print(
                    f"  dt_req={dt:10.4e} | "
                    f"dt_used={metrics.dt_used:10.4e} | "
                    f"stable={metrics.stable!s:5s} | "
                    f"err={metrics.err:12.5e} | "
                    f"max||y||/||y0||={metrics.max_norm_ratio:10.4e} | "
                    f"max step ratio={metrics.max_step_ratio:10.4e}"
                )

                detail_rows.append(
                    {
                        "mode": mode,
                        "Bmult": Bmult,
                        "method": method_name,
                        "dt_requested": dt,
                        "dt_used": metrics.dt_used,
                        "steps": metrics.steps,
                        "T": T,
                        "y0_0": float(y0[0]),
                        "y0_1": float(y0[1]),
                        "stable": int(metrics.stable),
                        "err": metrics.err,
                        "max_norm_ratio": metrics.max_norm_ratio,
                        "final_norm_ratio": metrics.final_norm_ratio,
                        "max_step_ratio": metrics.max_step_ratio,
                        "final_norm": metrics.final_norm,
                        "exploded_at_step": metrics.exploded_at_step,
                    }
                )

                if metrics.stable:
                    stable_count += 1
                    if np.isnan(largest_stable_dt) or metrics.dt_used > largest_stable_dt:
                        largest_stable_dt = metrics.dt_used
                else:
                    if np.isnan(smallest_unstable_dt) or metrics.dt_used < smallest_unstable_dt:
                        smallest_unstable_dt = metrics.dt_used

            summary_rows.append(
                {
                    "mode": mode,
                    "Bmult": Bmult,
                    "method": method_name,
                    "stable_count": stable_count,
                    "num_dts_tested": len(dts),
                    "largest_stable_dt": largest_stable_dt,
                    "smallest_unstable_dt": smallest_unstable_dt,
                }
            )

        print("\nLargest stable dt by method:")
        for row in [r for r in summary_rows if r["mode"] == mode and r["Bmult"] == Bmult]:
            print(
                f"  {str(row['method']):20s} -> "
                f"largest_stable_dt = {row['largest_stable_dt']!s}"
            )

    detail_path = outdir / f"abc_stability_scan_{mode}_detail.csv"
    summary_path = outdir / f"abc_stability_scan_{mode}_summary.csv"

    with detail_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "mode",
                "Bmult",
                "method",
                "dt_requested",
                "dt_used",
                "steps",
                "T",
                "y0_0",
                "y0_1",
                "stable",
                "err",
                "max_norm_ratio",
                "final_norm_ratio",
                "max_step_ratio",
                "final_norm",
                "exploded_at_step",
            ],
        )
        writer.writeheader()
        writer.writerows(detail_rows)

    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "mode",
                "Bmult",
                "method",
                "stable_count",
                "num_dts_tested",
                "largest_stable_dt",
                "smallest_unstable_dt",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\nSaved detail CSV:  {detail_path}")
    print(f"Saved summary CSV: {summary_path}")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Stability-window scan for higher-order ABC splitting methods."
    )

    parser.add_argument(
        "--mode",
        type=str,
        choices=["decay", "growth"],
        default="growth",
        help="Which regime to scan.",
    )
    parser.add_argument(
        "--Bmults",
        type=str,
        default="20,50,100",
        help="Comma-separated stiffness multipliers.",
    )
    parser.add_argument(
        "--dts",
        type=str,
        default="",
        help="Optional explicit comma-separated dt list. If omitted, use geometric grid from dtmax..dtmin.",
    )
    parser.add_argument(
        "--dtmax",
        type=float,
        default=2e-2,
        help="Largest requested dt when --dts is not provided.",
    )
    parser.add_argument(
        "--dtmin",
        type=float,
        default=3.125e-4,
        help="Smallest requested dt when --dts is not provided.",
    )
    parser.add_argument(
        "--levels",
        type=int,
        default=7,
        help="Number of geometric dt levels when --dts is not provided.",
    )
    parser.add_argument(
        "--methods",
        type=str,
        default="Strang-3,OS32_7op_minLEM-3,PP3_4A-3,Yoshida-3",
        help="Comma-separated pythOS alpha names to compare.",
    )
    parser.add_argument(
        "--T",
        type=float,
        default=0.5,
        help="Final time.",
    )
    parser.add_argument(
        "--y0",
        type=str,
        default="1.0,-0.5",
        help="Initial condition as 'y0,y1'.",
    )
    parser.add_argument(
        "--explosion_factor",
        type=float,
        default=1e10,
        help="Declare unstable if ||y|| exceeds explosion_factor * max(1, ||y0||).",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="experiments/abc_higher_order_stability_outputs",
        help="Directory for CSV outputs.",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.dts.strip():
        dts = parse_csv_floats(args.dts)
    else:
        dts = geometric_dt_list(args.dtmin, args.dtmax, args.levels)

    Bmults = parse_csv_floats(args.Bmults)
    methods = parse_csv_strings(args.methods)
    y0_vals = parse_csv_floats(args.y0)

    if len(y0_vals) != 2:
        raise ValueError("--y0 must have exactly two components, e.g. --y0 1.0,-0.5")

    y0 = np.array(y0_vals, dtype=float)
    outdir = REPO / args.outdir

    print("Running higher-order stability scan with:")
    print(f"  mode             = {args.mode}")
    print(f"  Bmults           = {Bmults}")
    print(f"  dts              = {dts}")
    print(f"  methods          = {methods}")
    print(f"  T                = {args.T}")
    print(f"  y0               = {y0}")
    print(f"  explosion_factor = {args.explosion_factor}")
    print(f"  outdir           = {outdir}")

    run_scan(
        mode=args.mode,
        Bmults=Bmults,
        dts=dts,
        methods=methods,
        T=args.T,
        y0=y0,
        outdir=outdir,
        explosion_factor=args.explosion_factor,
    )


if __name__ == "__main__":
    main()