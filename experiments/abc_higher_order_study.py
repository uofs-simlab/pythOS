"""
abc_higher_order_study.py

Phase 1 higher-order 3-operator splitting study on a small linear ABC test.

What it does
------------
- Defines a linear 2x2 system y' = (A + B + C) y with tunable stiffness via Bmult.
- Uses pythOS fractional_step(...) with ANALYTIC subflows for A, B, C.
- Compares several built-in 3-operator splitting methods already available in
  fractional_step.py, including higher-order candidates.
- For each (method, dt), computes:
    * final-time error vs exact expm solution
    * observed order from a dt sweep
    * simple stability / growth diagnostics:
        - max ||y_n|| / ||y_0||
        - final ||y_T|| / ||y_0||
        - max step ratio ||y_{n+1}|| / ||y_n||
- Saves one CSV per Bmult regime.

Examples
--------
python experiments/abc_higher_order_study.py
python experiments/abc_higher_order_study.py --Bmults 1,20,50,100
python experiments/abc_higher_order_study.py --dts 1e-2,5e-3,2.5e-3,1.25e-3
python experiments/abc_higher_order_study.py --methods Godunov-3,Strang-3,Yoshida-3,PP3_4A-3
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
    """
    Small-matrix expm fallback using eigendecomposition.
    Fine for this 2x2 experiment.
    """
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


def estimate_order(errs: Sequence[float], dts: Sequence[float]) -> float:
    """
    Least-squares fit to log(err) = p log(dt) + c, ignoring non-finite values.
    """
    xs: List[float] = []
    ys: List[float] = []
    for e, dt in zip(errs, dts):
        if np.isfinite(e) and e > 0.0 and dt > 0.0:
            xs.append(math.log(dt))
            ys.append(math.log(e))

    if len(xs) < 2:
        return float("nan")

    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)
    A = np.vstack([x, np.ones_like(x)]).T
    p, _ = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(p)


# ---------------------------------------------------------------------
# Problem definition
# ---------------------------------------------------------------------
def make_abc_matrices(Bmult: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    A, B, C are chosen so that:
    - operators do not commute
    - B controls stiffness / damping scale
    """
    A = np.array(
        [
            [0.0, 1.0],
            [-1.0, 0.0],
        ],
        dtype=float,
    )  # rotation-like

    B = np.array(
        [
            [-1.0, 0.0],
            [0.0, -0.2],
        ],
        dtype=float,
    ) * Bmult  # damping / stiffness

    C = np.array(
        [
            [0.0, 0.6],
            [0.0, 0.0],
        ],
        dtype=float,
    )  # shear-like

    return A, B, C


def exact_solution(y0: np.ndarray, T: float, A: np.ndarray, B: np.ndarray, C: np.ndarray) -> np.ndarray:
    return mat_expm(T * (A + B + C)) @ y0


# ---------------------------------------------------------------------
# Analytic subflows for fractional_step(..., methods={(i,): "ANALYTIC"})
# ---------------------------------------------------------------------
def make_linear_flow(L: np.ndarray) -> Callable[[float, float, np.ndarray], np.ndarray]:
    """
    fractional_step's analytic callback expects signature:
        flow(current_t, delta_t, y)
    """
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


def advance_one_step(
    y: np.ndarray,
    t: float,
    dt: float,
    flows: List[Callable[[float, float, np.ndarray], np.ndarray]],
    alpha_name: str,
) -> np.ndarray:
    """
    Advance exactly one global splitting step of size dt using fractional_step.
    """
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
    explosion_factor: float = 1e12,
) -> RunMetrics:
    """
    Repeatedly calls fractional_step one global step at a time so that we can
    track growth diagnostics across the trajectory.
    """
    steps = max(1, int(round(T / dt)))
    dt = T / steps  # force exact final time hit

    y = np.array(y0, dtype=float)
    t = 0.0

    n0 = norm2(y0)
    max_norm_ratio = 1.0
    max_step_ratio = 1.0
    stable = True

    for _ in range(steps):
        y_new = advance_one_step(y, t, dt, flows, alpha_name)

        if not np.all(np.isfinite(y_new)):
            stable = False
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
    )


# ---------------------------------------------------------------------
# Main study
# ---------------------------------------------------------------------
def run_study(
    Bmults: Sequence[float],
    dts: Sequence[float],
    methods: Sequence[str],
    T: float,
    y0: np.ndarray,
    outdir: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    for Bmult in Bmults:
        A, B, C = make_abc_matrices(Bmult)
        flows = [make_linear_flow(A), make_linear_flow(B), make_linear_flow(C)]
        y_exact = exact_solution(y0, T, A, B, C)

        rows: List[Dict[str, object]] = []
        err_by_method: Dict[str, List[float]] = {m: [] for m in methods}

        print("\n" + "=" * 72)
        print(f"HIGHER-ORDER ABC STUDY | Bmult = {Bmult:g} | T = {T:.6g}")
        print("=" * 72)

        for method_name in methods:
            print(f"\nMethod: {method_name}")

            for dt in dts:
                metrics = simulate_method(
                    alpha_name=method_name,
                    dt=dt,
                    T=T,
                    y0=y0,
                    flows=flows,
                    y_exact=y_exact,
                )
                err_by_method[method_name].append(metrics.err)

                print(
                    f"  dt={dt:10.4e} | "
                    f"err={metrics.err:12.5e} | "
                    f"stable={metrics.stable!s:5s} | "
                    f"max||y||/||y0||={metrics.max_norm_ratio:10.4e} | "
                    f"max step ratio={metrics.max_step_ratio:10.4e}"
                )

                rows.append(
                    {
                        "Bmult": Bmult,
                        "method": method_name,
                        "dt": dt,
                        "T": T,
                        "y0_0": float(y0[0]),
                        "y0_1": float(y0[1]),
                        "err": metrics.err,
                        "stable": int(metrics.stable),
                        "max_norm_ratio": metrics.max_norm_ratio,
                        "final_norm_ratio": metrics.final_norm_ratio,
                        "max_step_ratio": metrics.max_step_ratio,
                        "final_norm": metrics.final_norm,
                    }
                )

        p_by_method: Dict[str, float] = {
            m: estimate_order(err_by_method[m], dts) for m in methods
        }

        print("\nObserved orders:")
        for m in methods:
            print(f"  {m:20s} -> p ≈ {p_by_method[m]:.4f}")

        for row in rows:
            row["observed_order_fit"] = p_by_method[str(row["method"])]

        csv_path = outdir / f"abc_higher_order_Bmult{str(Bmult).replace('.', 'p')}.csv"
        with csv_path.open("w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "Bmult",
                    "method",
                    "dt",
                    "T",
                    "y0_0",
                    "y0_1",
                    "err",
                    "stable",
                    "max_norm_ratio",
                    "final_norm_ratio",
                    "max_step_ratio",
                    "final_norm",
                    "observed_order_fit",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)

        print(f"\nSaved: {csv_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Higher-order ABC splitting study using pythOS fractional_step.")

    parser.add_argument(
        "--Bmults",
        type=str,
        default="1,20,50,100",
        help="Comma-separated stiffness multipliers for B.",
    )
    parser.add_argument(
        "--dts",
        type=str,
        default="1e-2,5e-3,2.5e-3,1.25e-3",
        help="Comma-separated time steps for the order study.",
    )
    parser.add_argument(
        "--methods",
        type=str,
        default="Godunov-3,Strang-3,StrangACB-3,StrangBAC-3,StrangBCA-3,StrangCAB-3,StrangCBA-3,OS32_7op_minLEM-3,Yoshida-3,PP3_4A-3",
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
        "--outdir",
        type=str,
        default="experiments/abc_higher_order_outputs",
        help="Directory where CSV files will be saved.",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    Bmults = parse_csv_floats(args.Bmults)
    dts = parse_csv_floats(args.dts)
    methods = parse_csv_strings(args.methods)
    y0_vals = parse_csv_floats(args.y0)

    if len(y0_vals) != 2:
        raise ValueError("--y0 must have exactly two components, e.g. --y0 1.0,-0.5")

    y0 = np.array(y0_vals, dtype=float)
    outdir = REPO / args.outdir

    print("Running higher-order ABC study with:")
    print(f"  Bmults  = {Bmults}")
    print(f"  dts     = {dts}")
    print(f"  methods = {methods}")
    print(f"  T       = {args.T}")
    print(f"  y0      = {y0}")
    print(f"  outdir  = {outdir}")

    run_study(
        Bmults=Bmults,
        dts=dts,
        methods=methods,
        T=args.T,
        y0=y0,
        outdir=outdir,
    )


if __name__ == "__main__":
    main()