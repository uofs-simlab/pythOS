"""
Adaptive and fixed-step higher-order splitting study for a three-operator ABC test problem.

The script supports adaptive integration through step doubling, fixed-step
comparisons for the same methods, and CSV summaries suitable for plots or
tables.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple, Any

import numpy as np

try:
    from experiments.utils.adaptive_splitting import (
        AdaptiveOptions,
        adaptive_integrate,
        l2_norm,
    )
except ModuleNotFoundError:
    import sys

    THIS_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.dirname(THIS_DIR)
    if REPO_ROOT not in sys.path:
        sys.path.insert(0, REPO_ROOT)

    from experiments.utils.adaptive_splitting import (
        AdaptiveOptions,
        adaptive_integrate,
        l2_norm,
    )


# =============================================================================
# Problem definition
# =============================================================================

@dataclass
class ABCProblem:
    Bmult: float

    def full_rhs(self, y: np.ndarray, t: float) -> np.ndarray:
        y1, y2, y3 = y

        Ay = np.array([
            y2,
            -y1,
            0.0,
        ])

        lam = self.Bmult
        By = np.array([
            0.0,
            -lam * y2,
            -0.5 * lam * y3,
        ])

        Cy = np.array([
            -0.2 * y1 * y3,
            0.1 * y1 * y1,
            0.3 * y1 - 0.1 * y3,
        ])

        return Ay + By + Cy

    def phi_A(self, y: np.ndarray, dt: float) -> np.ndarray:
        y1, y2, y3 = y
        c = math.cos(dt)
        s = math.sin(dt)
        return np.array([
            c * y1 + s * y2,
            -s * y1 + c * y2,
            y3,
        ], dtype=float)

    def phi_B(self, y: np.ndarray, dt: float) -> np.ndarray:
        y1, y2, y3 = y
        lam = self.Bmult
        return np.array([
            y1,
            math.exp(-lam * dt) * y2,
            math.exp(-0.5 * lam * dt) * y3,
        ], dtype=float)

    def phi_C(self, y: np.ndarray, dt: float) -> np.ndarray:
        def rhs_c(z: np.ndarray) -> np.ndarray:
            z1, _, z3 = z
            return np.array([
                -0.2 * z1 * z3,
                0.1 * z1 * z1,
                0.3 * z1 - 0.1 * z3,
            ], dtype=float)

        z = y.copy()
        k1 = rhs_c(z)
        k2 = rhs_c(z + 0.5 * dt * k1)
        k3 = rhs_c(z + 0.5 * dt * k2)
        k4 = rhs_c(z + dt * k3)
        return z + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def make_abc_problem(Bmult: float) -> Tuple[ABCProblem, np.ndarray]:
    problem = ABCProblem(Bmult=Bmult)
    y0 = np.array([1.0, -0.25, 0.5], dtype=float)
    return problem, y0


# =============================================================================
# Splitting methods
# =============================================================================

def strang_abc_step(problem: ABCProblem, y: np.ndarray, t: float, dt: float) -> np.ndarray:
    z = problem.phi_A(y, 0.5 * dt)
    z = problem.phi_B(z, 0.5 * dt)
    z = problem.phi_C(z, dt)
    z = problem.phi_B(z, 0.5 * dt)
    z = problem.phi_A(z, 0.5 * dt)
    return z


def yoshida4_from_strang(
    base2: Callable[[np.ndarray, float, float], np.ndarray],
) -> Callable[[np.ndarray, float, float], np.ndarray]:
    cbrt2 = 2.0 ** (1.0 / 3.0)
    w1 = 1.0 / (2.0 - cbrt2)
    w0 = -cbrt2 / (2.0 - cbrt2)

    def step(y: np.ndarray, t: float, dt: float) -> np.ndarray:
        z = base2(y, t, w1 * dt)
        z = base2(z, t + w1 * dt, w0 * dt)
        z = base2(z, t + (w1 + w0) * dt, w1 * dt)
        return z

    return step


def make_base_stepper(problem: ABCProblem, method: str) -> Tuple[Callable[[np.ndarray, float, float], np.ndarray], int]:
    method = method.lower()

    if method == "strang2":
        return (lambda y, t, dt: strang_abc_step(problem, y, t, dt), 2)

    if method == "yoshida4":
        base2 = lambda y, t, dt: strang_abc_step(problem, y, t, dt)
        return yoshida4_from_strang(base2), 4

    raise ValueError(f"Unknown method '{method}'. Supported: strang2, yoshida4")


# =============================================================================
# Reference solver
# =============================================================================

def rk4_step(rhs: Callable[[np.ndarray, float], np.ndarray], y: np.ndarray, t: float, dt: float) -> np.ndarray:
    k1 = rhs(y, t)
    k2 = rhs(y + 0.5 * dt * k1, t + 0.5 * dt)
    k3 = rhs(y + 0.5 * dt * k2, t + 0.5 * dt)
    k4 = rhs(y + dt * k3, t + dt)
    return y + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


def integrate_reference(
    rhs: Callable[[np.ndarray, float], np.ndarray],
    y0: np.ndarray,
    T: float,
    dt_ref: float,
) -> np.ndarray:
    t = 0.0
    y = y0.copy()
    while t < T - 1.0e-15:
        dt = min(dt_ref, T - t)
        y = rk4_step(rhs, y, t, dt)
        t += dt
    return y


# =============================================================================
# Fixed-step solver
# =============================================================================

def integrate_fixed_step(
    stepper: Callable[[np.ndarray, float, float], np.ndarray],
    y0: np.ndarray,
    T: float,
    dt: float,
) -> Tuple[np.ndarray, int, float]:
    t = 0.0
    y = y0.copy()
    n_steps = 0
    tic = time.perf_counter()

    while t < T - 1.0e-15:
        h = min(dt, T - t)
        y = stepper(y, t, h)
        t += h
        n_steps += 1

    wall = time.perf_counter() - tic
    return y, n_steps, wall


# =============================================================================
# IO helpers
# =============================================================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def save_adaptive_history_csv(path: str, result) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "attempt_index",
            "t",
            "dt_try",
            "dt_suggested",
            "accepted",
            "err_est",
            "cpu_seconds",
        ])
        for i, rec in enumerate(result.records):
            writer.writerow([
                i,
                rec.t,
                rec.dt_try,
                rec.dt_suggested,
                int(rec.accepted),
                rec.err_est,
                rec.cpu_seconds,
            ])


def append_summary_csv(path: str, row: Dict[str, Any]) -> None:
    file_exists = os.path.exists(path)
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def parse_float_list(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


# =============================================================================
# Study runners
# =============================================================================

def run_one_adaptive_case(
    Bmult: float,
    method: str,
    T: float,
    dt0: float,
    atol: float,
    rtol: float,
    dt_ref: float,
    outdir: str,
) -> Dict[str, Any]:
    problem, y0 = make_abc_problem(Bmult=Bmult)
    stepper, order = make_base_stepper(problem, method=method)

    opts = AdaptiveOptions(
        order=order,
        atol=atol,
        rtol=rtol,
        dt_min=1.0e-12,
        dt_max=max(dt0, T),
        safety=0.9,
        growth_max=2.0,
        shrink_min=0.2,
        max_reject=20,
        use_weighted_rms=True,
    )

    result = adaptive_integrate(
        stepper=stepper,
        y0=y0,
        t0=0.0,
        tfinal=T,
        dt0=dt0,
        opts=opts,
    )

    y_ref = integrate_reference(problem.full_rhs, y0, T, dt_ref=dt_ref)
    err_global = l2_norm(result.final_state - y_ref)

    history_name = (
        f"abc_history_mode-adaptive_method-{method}_B{Bmult:g}_"
        f"atol{atol:.0e}_rtol{rtol:.0e}.csv"
    )
    save_adaptive_history_csv(os.path.join(outdir, history_name), result)

    return {
        "mode": "adaptive",
        "method": method,
        "order": order,
        "Bmult": Bmult,
        "T": T,
        "dt0": dt0,
        "fixed_dt": "",
        "atol": atol,
        "rtol": rtol,
        "dt_ref": dt_ref,
        "success": int(result.success),
        "message": result.message,
        "n_steps": result.n_accept,
        "n_accept": result.n_accept,
        "n_reject": result.n_reject,
        "dt_min_used": result.dt_min_used,
        "dt_max_used": result.dt_max_used,
        "dt_avg_used": result.dt_avg_used,
        "wall_seconds": result.wall_seconds,
        "global_l2_error": err_global,
    }


def run_one_fixed_case(
    Bmult: float,
    method: str,
    T: float,
    fixed_dt: float,
    dt_ref: float,
) -> Dict[str, Any]:
    problem, y0 = make_abc_problem(Bmult=Bmult)
    stepper, order = make_base_stepper(problem, method=method)

    y_final, n_steps, wall = integrate_fixed_step(
        stepper=stepper,
        y0=y0,
        T=T,
        dt=fixed_dt,
    )

    y_ref = integrate_reference(problem.full_rhs, y0, T, dt_ref=dt_ref)
    err_global = l2_norm(y_final - y_ref)

    return {
        "mode": "fixed",
        "method": method,
        "order": order,
        "Bmult": Bmult,
        "T": T,
        "dt0": "",
        "fixed_dt": fixed_dt,
        "atol": "",
        "rtol": "",
        "dt_ref": dt_ref,
        "success": 1,
        "message": "Fixed-step integration completed successfully.",
        "n_steps": n_steps,
        "n_accept": n_steps,
        "n_reject": 0,
        "dt_min_used": fixed_dt,
        "dt_max_used": fixed_dt,
        "dt_avg_used": fixed_dt,
        "wall_seconds": wall,
        "global_l2_error": err_global,
    }


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="ABC higher-order splitting study: adaptive vs fixed.")
    parser.add_argument("--mode", type=str, default="adaptive",
                        choices=["adaptive", "fixed", "both"],
                        help="Run adaptive, fixed, or both.")
    parser.add_argument("--method", type=str, default="yoshida4",
                        choices=["strang2", "yoshida4"],
                        help="Base splitting method.")
    parser.add_argument("--Bmults", type=str, default="100,200,500",
                        help="Comma-separated stiffness multipliers.")
    parser.add_argument("--T", type=float, default=0.5,
                        help="Final time.")
    parser.add_argument("--dt0", type=float, default=1.0e-2,
                        help="Initial adaptive step size guess.")
    parser.add_argument("--atols", type=str, default="1e-6,1e-8",
                        help="Comma-separated absolute tolerances for adaptive mode.")
    parser.add_argument("--rtols", type=str, default="1e-4,1e-6",
                        help="Comma-separated relative tolerances for adaptive mode.")
    parser.add_argument("--fixed-dts", type=str, default="1e-2,5e-3,2e-3,1e-3,5e-4",
                        help="Comma-separated fixed time steps for fixed mode.")
    parser.add_argument("--dt-ref", type=float, default=1.0e-6,
                        help="Reference RK4 step size.")
    parser.add_argument("--outdir", type=str, default="experiments/outputs/abc_adaptive_higher_order",
                        help="Output directory.")
    parser.add_argument("--overwrite-summary", action="store_true",
                        help="Overwrite summary CSV if it already exists.")
    args = parser.parse_args()

    ensure_dir(args.outdir)

    Bmults = parse_float_list(args.Bmults)
    atols = parse_float_list(args.atols)
    rtols = parse_float_list(args.rtols)
    fixed_dts = parse_float_list(args.fixed_dts)

    summary_path = os.path.join(args.outdir, "abc_higher_order_comparison_summary.csv")
    if args.overwrite_summary and os.path.exists(summary_path):
        os.remove(summary_path)

    print("=" * 72)
    print("ABC HIGHER-ORDER SPLITTING COMPARISON STUDY")
    print(f"mode      = {args.mode}")
    print(f"method    = {args.method}")
    print(f"Bmults    = {Bmults}")
    print(f"T         = {args.T}")
    print(f"dt0       = {args.dt0}")
    print(f"atols     = {atols}")
    print(f"rtols     = {rtols}")
    print(f"fixed_dts = {fixed_dts}")
    print(f"dt_ref    = {args.dt_ref}")
    print(f"outdir    = {args.outdir}")
    print("=" * 72)

    if args.mode in ("adaptive", "both"):
        for Bmult in Bmults:
            for atol in atols:
                for rtol in rtols:
                    print(
                        f"\nAdaptive: method={args.method}, "
                        f"Bmult={Bmult:g}, atol={atol:.0e}, rtol={rtol:.0e}"
                    )

                    row = run_one_adaptive_case(
                        Bmult=Bmult,
                        method=args.method,
                        T=args.T,
                        dt0=args.dt0,
                        atol=atol,
                        rtol=rtol,
                        dt_ref=args.dt_ref,
                        outdir=args.outdir,
                    )
                    append_summary_csv(summary_path, row)

                    print(
                        f"  success={row['success']} "
                        f"accept={row['n_accept']} reject={row['n_reject']} "
                        f"dt[min,avg,max]=({row['dt_min_used']:.3e}, "
                        f"{row['dt_avg_used']:.3e}, {row['dt_max_used']:.3e}) "
                        f"error={row['global_l2_error']:.6e} "
                        f"wall={row['wall_seconds']:.3f}s"
                    )

    if args.mode in ("fixed", "both"):
        for Bmult in Bmults:
            for fixed_dt in fixed_dts:
                print(
                    f"\nFixed: method={args.method}, "
                    f"Bmult={Bmult:g}, dt={fixed_dt:.3e}"
                )

                row = run_one_fixed_case(
                    Bmult=Bmult,
                    method=args.method,
                    T=args.T,
                    fixed_dt=fixed_dt,
                    dt_ref=args.dt_ref,
                )
                append_summary_csv(summary_path, row)

                print(
                    f"  steps={row['n_steps']} "
                    f"error={row['global_l2_error']:.6e} "
                    f"wall={row['wall_seconds']:.3f}s"
                )

    print("\nDone.")
    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()