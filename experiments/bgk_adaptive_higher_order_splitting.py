"""
Adaptive higher-order splitting study for a BGK-style kinetic test problem.

This script reuses the generic adaptive step-doubling driver, records BGK-specific
diagnostics such as min(f), and compares adaptive and fixed-step runs.

The file is structured so it can be connected to an existing three-operator BGK
implementation through make_bgk_problem(...), make_bgk_stepper(...), and
make_reference_solver(...).
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

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


Array = np.ndarray


# =============================================================================
# Problem wrapper
# =============================================================================

@dataclass
class BGKProblem:
    """
    Generic wrapper for a BGK-type problem.

    Required capabilities:
    - initial state y0
    - optional reference integration
    - one-step higher-order splitting method through make_bgk_stepper(...)
    - diagnostics on state, especially min(f)

    State convention:
    - by default we assume y is an ndarray containing f over phase space
    - if your code uses a custom object, adapt diagnostics_fn and norm_fn
    """
    case: str
    Nx: int
    Nv: int
    xmax: float
    vmax: float
    eps0: float
    transport: str
    eta_model: str
    init: str

    def initial_state(self) -> Array:
        """
        Placeholder initial state used for repository structure testing.
        """
        x = np.linspace(-self.xmax, self.xmax, self.Nx, endpoint=False)
        v = np.linspace(-self.vmax, self.vmax, self.Nv, endpoint=False)
        X, V = np.meshgrid(x, v, indexing="ij")

        rho = 1.0 + 0.2 * np.sin(np.pi * X / max(self.xmax, 1.0))
        u = 0.3 * np.cos(np.pi * X / max(self.xmax, 1.0))
        T = 1.0 + 0.1 * np.sin(2.0 * np.pi * X / max(self.xmax, 1.0))

        f = rho / np.sqrt(2.0 * np.pi * T) * np.exp(-0.5 * ((V - u) ** 2) / T)
        return f.astype(float)

    def diagnostics(self, y: Array) -> Dict[str, Any]:
        """
        BGK-specific diagnostics stored per attempted adaptive step.
        """
        arr = np.asarray(y, dtype=float)
        return {
            "min_f": float(np.min(arr)),
            "max_f": float(np.max(arr)),
            "l2_state": float(l2_norm(arr)),
        }

    def error_norm(self, diff: Array) -> float:
        """
        Error norm used when not using weighted RMS.
        """
        return float(l2_norm(diff))


# =============================================================================
# Hook points to existing BGK implementation
# =============================================================================

def make_bgk_problem(
    case: str,
    Nx: int,
    Nv: int,
    xmax: float,
    vmax: float,
    eps0: float,
    transport: str,
    eta_model: str,
    init: str,
) -> BGKProblem:
    """
    Construct the BGK problem.
    """
    return BGKProblem(
        case=case,
        Nx=Nx,
        Nv=Nv,
        xmax=xmax,
        vmax=vmax,
        eps0=eps0,
        transport=transport,
        eta_model=eta_model,
        init=init,
    )


def make_bgk_stepper(problem: BGKProblem, method: str) -> Tuple[Callable[[Array, float, float], Array], int]:
    """
    Return (stepper, nominal_order).

    This is the main hook for an existing higher-order BGK splitting implementation.

    Expected signature:
        stepper(y, t, dt) -> y_next

    Supported placeholder methods:
    - strang2
    - yoshida4

    Note:
    The placeholder below is only a structural stand-in for the project-specific
    BGK operator splitting routine.
    """
    method = method.lower()

    def demo_collision_relax(y: Array, dt: float) -> Array:
        # mild relaxation toward local copy; structural placeholder only
        lam = max(problem.eps0, 1.0e-12)
        return np.exp(-dt / lam) * y + (1.0 - np.exp(-dt / lam)) * np.maximum(y, 0.0)

    def demo_transport(y: Array, dt: float) -> Array:
        # simple smoothing/roll placeholder; replace with actual transport solve
        z = np.asarray(y, dtype=float)
        rolled = np.roll(z, shift=1, axis=0)
        return z - 0.1 * dt * (z - rolled)

    def demo_aux(y: Array, dt: float) -> Array:
        # placeholder third operator
        z = np.asarray(y, dtype=float)
        return np.maximum(z * (1.0 - 0.02 * dt), 0.0)

    def strang_step(y: Array, t: float, dt: float) -> Array:
        z = demo_transport(y, 0.5 * dt)
        z = demo_collision_relax(z, 0.5 * dt)
        z = demo_aux(z, dt)
        z = demo_collision_relax(z, 0.5 * dt)
        z = demo_transport(z, 0.5 * dt)
        return z

    if method == "strang2":
        return strang_step, 2

    if method == "yoshida4":
        cbrt2 = 2.0 ** (1.0 / 3.0)
        w1 = 1.0 / (2.0 - cbrt2)
        w0 = -cbrt2 / (2.0 - cbrt2)

        def yoshida4_step(y: Array, t: float, dt: float) -> Array:
            z = strang_step(y, t, w1 * dt)
            z = strang_step(z, t + w1 * dt, w0 * dt)
            z = strang_step(z, t + (w1 + w0) * dt, w1 * dt)
            return z

        return yoshida4_step, 4

    raise ValueError(f"Unknown method '{method}'. Supported: strang2, yoshida4")


def make_reference_solver(problem: BGKProblem, reference_method: str) -> Callable[[Array, float, float], Array]:
    """
    Return a one-step reference integrator step(y, t, dt) -> y_next.
    """
    reference_method = reference_method.lower()

    def ref_step(y: Array, t: float, dt: float) -> Array:
        z = np.asarray(y, dtype=float)
        z = np.roll(z, shift=1, axis=0) * (0.02 * dt) + z * (1.0 - 0.02 * dt)
        z = np.maximum(z * np.exp(-0.2 * dt), 0.0)
        return z

    if reference_method in ("rk4", "exprk2", "strang2", "yoshida4"):
        return ref_step

    raise ValueError(f"Unknown reference method '{reference_method}'.")


# =============================================================================
# Integration helpers
# =============================================================================

def integrate_fixed_step(
    stepper: Callable[[Array, float, float], Array],
    y0: Array,
    T: float,
    dt: float,
) -> Tuple[Array, int, float]:
    t = 0.0
    y = np.array(y0, copy=True)
    n_steps = 0
    tic = time.perf_counter()

    while t < T - 1.0e-15:
        h = min(dt, T - t)
        y = stepper(y, t, h)
        t += h
        n_steps += 1

    wall = time.perf_counter() - tic
    return y, n_steps, wall


def integrate_reference(
    ref_stepper: Callable[[Array, float, float], Array],
    y0: Array,
    T: float,
    dt_ref: float,
) -> Array:
    t = 0.0
    y = np.array(y0, copy=True)

    while t < T - 1.0e-15:
        h = min(dt_ref, T - t)
        y = ref_stepper(y, t, h)
        t += h

    return y


# =============================================================================
# CSV helpers
# =============================================================================

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def append_summary_csv(path: str, row: Dict[str, Any]) -> None:
    exists = os.path.exists(path)
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


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
            "min_f",
            "max_f",
            "l2_state",
        ])
        for i, rec in enumerate(result.records):
            extra = rec.extra or {}
            writer.writerow([
                i,
                rec.t,
                rec.dt_try,
                rec.dt_suggested,
                int(rec.accepted),
                rec.err_est,
                rec.cpu_seconds,
                extra.get("min_f", ""),
                extra.get("max_f", ""),
                extra.get("l2_state", ""),
            ])


def parse_float_list(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


# =============================================================================
# Case runners
# =============================================================================

def run_one_adaptive_case(
    case: str,
    method: str,
    reference_method: str,
    T: float,
    dt0: float,
    atol: float,
    rtol: float,
    dt_ref: float,
    Nx: int,
    Nv: int,
    xmax: float,
    vmax: float,
    eps0: float,
    transport: str,
    eta_model: str,
    init: str,
    outdir: str,
) -> Dict[str, Any]:
    problem = make_bgk_problem(
        case=case,
        Nx=Nx,
        Nv=Nv,
        xmax=xmax,
        vmax=vmax,
        eps0=eps0,
        transport=transport,
        eta_model=eta_model,
        init=init,
    )
    y0 = problem.initial_state()

    stepper, order = make_bgk_stepper(problem, method=method)
    ref_stepper = make_reference_solver(problem, reference_method=reference_method)

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
        diagnostics_fn=problem.diagnostics,
    )

    y_ref = integrate_reference(ref_stepper, y0, T=T, dt_ref=dt_ref)
    err_global = l2_norm(result.final_state - y_ref)

    final_diag = problem.diagnostics(result.final_state)

    history_name = (
        f"bgk_history_mode-adaptive_method-{method}_case-{case}_"
        f"eps0-{eps0:.0e}_atol-{atol:.0e}_rtol-{rtol:.0e}.csv"
    )
    save_adaptive_history_csv(os.path.join(outdir, history_name), result)

    return {
        "mode": "adaptive",
        "case": case,
        "method": method,
        "reference_method": reference_method,
        "order": order,
        "T": T,
        "dt0": dt0,
        "fixed_dt": "",
        "atol": atol,
        "rtol": rtol,
        "dt_ref": dt_ref,
        "Nx": Nx,
        "Nv": Nv,
        "xmax": xmax,
        "vmax": vmax,
        "eps0": eps0,
        "transport": transport,
        "eta_model": eta_model,
        "init": init,
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
        "final_min_f": final_diag["min_f"],
        "final_max_f": final_diag["max_f"],
        "final_l2_state": final_diag["l2_state"],
    }


def run_one_fixed_case(
    case: str,
    method: str,
    reference_method: str,
    T: float,
    fixed_dt: float,
    dt_ref: float,
    Nx: int,
    Nv: int,
    xmax: float,
    vmax: float,
    eps0: float,
    transport: str,
    eta_model: str,
    init: str,
) -> Dict[str, Any]:
    problem = make_bgk_problem(
        case=case,
        Nx=Nx,
        Nv=Nv,
        xmax=xmax,
        vmax=vmax,
        eps0=eps0,
        transport=transport,
        eta_model=eta_model,
        init=init,
    )
    y0 = problem.initial_state()

    stepper, order = make_bgk_stepper(problem, method=method)
    ref_stepper = make_reference_solver(problem, reference_method=reference_method)

    y_final, n_steps, wall = integrate_fixed_step(stepper, y0, T=T, dt=fixed_dt)
    y_ref = integrate_reference(ref_stepper, y0, T=T, dt_ref=dt_ref)
    err_global = l2_norm(y_final - y_ref)
    final_diag = problem.diagnostics(y_final)

    return {
        "mode": "fixed",
        "case": case,
        "method": method,
        "reference_method": reference_method,
        "order": order,
        "T": T,
        "dt0": "",
        "fixed_dt": fixed_dt,
        "atol": "",
        "rtol": "",
        "dt_ref": dt_ref,
        "Nx": Nx,
        "Nv": Nv,
        "xmax": xmax,
        "vmax": vmax,
        "eps0": eps0,
        "transport": transport,
        "eta_model": eta_model,
        "init": init,
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
        "final_min_f": final_diag["min_f"],
        "final_max_f": final_diag["max_f"],
        "final_l2_state": final_diag["l2_state"],
    }


# =============================================================================
# Command-line interface
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Adaptive BGK higher-order splitting study.")
    parser.add_argument("--mode", type=str, default="adaptive",
                        choices=["adaptive", "fixed", "both"])
    parser.add_argument("--case", type=str, default="mixed",
                        help="BGK case label, e.g. mixed, smooth, stiff.")
    parser.add_argument("--method", type=str, default="yoshida4",
                        choices=["strang2", "yoshida4"])
    parser.add_argument("--reference-method", type=str, default="rk4",
                        help="Reference integrator label.")
    parser.add_argument("--T", type=float, default=0.5)
    parser.add_argument("--dt0", type=float, default=1.0e-3)
    parser.add_argument("--atols", type=str, default="1e-6,1e-8")
    parser.add_argument("--rtols", type=str, default="1e-4,1e-6")
    parser.add_argument("--fixed-dts", type=str, default="2e-3,1e-3,5e-4")
    parser.add_argument("--dt-ref", type=float, default=1.0e-5)

    parser.add_argument("--Nx", type=int, default=40)
    parser.add_argument("--Nv", type=int, default=100)
    parser.add_argument("--xmax", type=float, default=2.0)
    parser.add_argument("--vmax", type=float, default=15.0)
    parser.add_argument("--eps0-list", type=str, default="1e-2,1e-3,1e-4",
                        help="Comma-separated stiffness / regime parameters.")
    parser.add_argument("--transport", type=str, default="muscl2")
    parser.add_argument("--eta-model", type=str, default="constant1")
    parser.add_argument("--init", type=str, default="mixed")

    parser.add_argument("--outdir", type=str,
                        default="experiments/outputs/bgk_adaptive_higher_order")
    parser.add_argument("--overwrite-summary", action="store_true")
    args = parser.parse_args()

    ensure_dir(args.outdir)

    atols = parse_float_list(args.atols)
    rtols = parse_float_list(args.rtols)
    fixed_dts = parse_float_list(args.fixed_dts)
    eps0_list = parse_float_list(args.eps0_list)

    summary_path = os.path.join(args.outdir, "bgk_higher_order_comparison_summary.csv")
    if args.overwrite_summary and os.path.exists(summary_path):
        os.remove(summary_path)

    print("=" * 72)
    print("BGK HIGHER-ORDER SPLITTING COMPARISON STUDY")
    print(f"mode             = {args.mode}")
    print(f"case             = {args.case}")
    print(f"method           = {args.method}")
    print(f"reference_method = {args.reference_method}")
    print(f"T                = {args.T}")
    print(f"dt0              = {args.dt0}")
    print(f"atols            = {atols}")
    print(f"rtols            = {rtols}")
    print(f"fixed_dts        = {fixed_dts}")
    print(f"dt_ref           = {args.dt_ref}")
    print(f"Nx, Nv           = {args.Nx}, {args.Nv}")
    print(f"xmax, vmax       = {args.xmax}, {args.vmax}")
    print(f"eps0_list        = {eps0_list}")
    print(f"transport        = {args.transport}")
    print(f"eta_model        = {args.eta_model}")
    print(f"init             = {args.init}")
    print(f"outdir           = {args.outdir}")
    print("=" * 72)

    if args.mode in ("adaptive", "both"):
        for eps0 in eps0_list:
            for atol in atols:
                for rtol in rtols:
                    print(
                        f"\nAdaptive: case={args.case}, method={args.method}, "
                        f"eps0={eps0:.0e}, atol={atol:.0e}, rtol={rtol:.0e}"
                    )

                    row = run_one_adaptive_case(
                        case=args.case,
                        method=args.method,
                        reference_method=args.reference_method,
                        T=args.T,
                        dt0=args.dt0,
                        atol=atol,
                        rtol=rtol,
                        dt_ref=args.dt_ref,
                        Nx=args.Nx,
                        Nv=args.Nv,
                        xmax=args.xmax,
                        vmax=args.vmax,
                        eps0=eps0,
                        transport=args.transport,
                        eta_model=args.eta_model,
                        init=args.init,
                        outdir=args.outdir,
                    )
                    append_summary_csv(summary_path, row)

                    print(
                        f"  success={row['success']} "
                        f"accept={row['n_accept']} reject={row['n_reject']} "
                        f"dt[min,avg,max]=({row['dt_min_used']:.3e}, "
                        f"{row['dt_avg_used']:.3e}, {row['dt_max_used']:.3e}) "
                        f"error={row['global_l2_error']:.6e} "
                        f"min(f)={row['final_min_f']:.6e} "
                        f"wall={row['wall_seconds']:.3f}s"
                    )

    if args.mode in ("fixed", "both"):
        for eps0 in eps0_list:
            for fixed_dt in fixed_dts:
                print(
                    f"\nFixed: case={args.case}, method={args.method}, "
                    f"eps0={eps0:.0e}, dt={fixed_dt:.3e}"
                )

                row = run_one_fixed_case(
                    case=args.case,
                    method=args.method,
                    reference_method=args.reference_method,
                    T=args.T,
                    fixed_dt=fixed_dt,
                    dt_ref=args.dt_ref,
                    Nx=args.Nx,
                    Nv=args.Nv,
                    xmax=args.xmax,
                    vmax=args.vmax,
                    eps0=eps0,
                    transport=args.transport,
                    eta_model=args.eta_model,
                    init=args.init,
                )
                append_summary_csv(summary_path, row)

                print(
                    f"  steps={row['n_steps']} "
                    f"error={row['global_l2_error']:.6e} "
                    f"min(f)={row['final_min_f']:.6e} "
                    f"wall={row['wall_seconds']:.3f}s"
                )

    print("\nDone.")
    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()