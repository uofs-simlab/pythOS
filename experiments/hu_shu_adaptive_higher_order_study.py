"""
Adaptive higher-order study for the Hu–Shu BGK setting.

This script is designed to live alongside the validated
experiments/hu_shu_expk2_bgk.py baseline rather than replace it. It uses the
exported helpers from that file and runtime signature inspection so that calls
remain compatible with the available function signatures.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import os
import time
import traceback
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

try:
    import experiments.hu_shu_expk2_bgk as hs
except Exception as exc:
    raise ImportError(
        "Could not import experiments.hu_shu_expk2_bgk. "
        "Make sure you are running from the repo root or with the repo on PYTHONPATH."
    ) from exc

try:
    from fractional_step import fractional_step
except Exception:
    fractional_step = None

Array = np.ndarray


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def parse_float_list(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


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
            "rho_min",
            "u_min",
            "T_min",
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
                extra.get("rho_min", ""),
                extra.get("u_min", ""),
                extra.get("T_min", ""),
            ])


def call_with_known_kwargs(func: Callable, pool: Dict[str, Any]) -> Any:
    sig = inspect.signature(func)
    kwargs = {}
    accepts_var_kw = False
    for name, param in sig.parameters.items():
        if param.kind == inspect.Parameter.VAR_KEYWORD:
            accepts_var_kw = True
        elif param.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY):
            if name in pool:
                kwargs[name] = pool[name]
    if accepts_var_kw:
        return func(**pool)
    return func(**kwargs)


def extract_attr(obj: Any, candidates: List[str], default: Any = None) -> Any:
    for name in candidates:
        if hasattr(obj, name):
            return getattr(obj, name)
    return default


@dataclass
class HuShuProblem:
    case: str
    Nx: int
    Nv: int
    xmax: float
    vmax: float
    eps0: float
    transport: str
    eta_model: str
    init: str

    grid: Any = None
    x: Optional[Array] = None
    v: Optional[Array] = None
    f0: Optional[Array] = None
    eps_x: Optional[Array] = None

    def __post_init__(self) -> None:
        if fractional_step is None:
            raise ImportError("Could not import pythOS fractional_step. Make sure the pythOS environment is set correctly.")

        if not hasattr(hs, "make_grid"):
            raise AttributeError("hu_shu_expk2_bgk.py does not expose make_grid.")

        self.grid = call_with_known_kwargs(
            hs.make_grid,
            {
                "Nx": self.Nx,
                "Nv": self.Nv,
                "xmax": self.xmax,
                "vmax": self.vmax,
                "x_min": -self.xmax,
                "x_max": self.xmax,
                "v_min": -self.vmax,
                "v_max": self.vmax,
            },
        )

        self.x = extract_attr(self.grid, ["x", "xc", "x_grid", "x_nodes"])
        self.v = extract_attr(self.grid, ["v", "vc", "v_grid", "v_nodes"])

        if self.x is None or self.v is None:
            self.x = np.linspace(-self.xmax, self.xmax, self.Nx, endpoint=False)
            self.v = np.linspace(-self.vmax, self.vmax, self.Nv, endpoint=False)

        if not hasattr(hs, "initial_condition"):
            raise AttributeError("hu_shu_expk2_bgk.py does not expose initial_condition.")

        self.f0 = call_with_known_kwargs(
            hs.initial_condition,
            {
                "grid": self.grid,
                "x": self.x,
                "v": self.v,
                "Nx": self.Nx,
                "Nv": self.Nv,
                "kind": self.init,
                "init": self.init,
                "init_name": self.init,
                "case": self.case,
                "eps0": self.eps0,
                "eta_model": self.eta_model,
            },
        )

        if self.case == "mixed":
            if hasattr(hs, "mixed_regime_eps"):
                self.eps_x = call_with_known_kwargs(
                    hs.mixed_regime_eps,
                    {
                        "grid": self.grid,
                        "x": self.x,
                        "eps0": self.eps0,
                        "eps": self.eps0,
                    },
                )
            else:
                self.eps_x = self.eps0 * np.ones_like(self.x)
        else:
            self.eps_x = self.eps0 * np.ones_like(self.x)

    def diagnostics(self, y: Array) -> Dict[str, Any]:
        out: Dict[str, Any] = {"min_f": float(np.min(y))}
        try:
            rho, u, T = call_with_known_kwargs(
                hs.moments,
                {
                    "f": y,
                    "v": self.v,
                    "grid": self.grid,
                },
            )
            out["rho_min"] = float(np.min(rho))
            out["u_min"] = float(np.min(u))
            out["T_min"] = float(np.min(T))
        except Exception:
            pass
        return out

    def macro_errors(self, y: Array, y_ref: Array) -> Tuple[float, float, float]:
        try:
            rho, u, T = call_with_known_kwargs(hs.moments, {"f": y, "v": self.v, "grid": self.grid})
            rho_ref, u_ref, T_ref = call_with_known_kwargs(hs.moments, {"f": y_ref, "v": self.v, "grid": self.grid})
            return (
                float(np.mean(np.abs(rho - rho_ref))),
                float(np.mean(np.abs(u - u_ref))),
                float(np.mean(np.abs(T - T_ref))),
            )
        except Exception:
            return (float("nan"), float("nan"), float("nan"))

    def full_rhs(self, y: Array, t: float) -> Array:
        if not hasattr(hs, "full_rhs"):
            raise AttributeError("hu_shu_expk2_bgk.py does not expose full_rhs.")
        return call_with_known_kwargs(
            hs.full_rhs,
            {
                "f": y,
                "y": y,
                "t": t,
                "grid": self.grid,
                "x": self.x,
                "v": self.v,
                "eps0": self.eps0,
                "eps": self.eps0,
                "eps_x": self.eps_x,
                "case": self.case,
                "transport": self.transport,
                "transport_name": self.transport,
                "eta_model": self.eta_model,
            },
        )

    def exprk2_step(self, y: Array, t: float, dt: float) -> Array:
        if not hasattr(hs, "hu_shu_exprk2_step"):
            raise AttributeError("hu_shu_expk2_bgk.py does not expose hu_shu_exprk2_step.")
        return call_with_known_kwargs(
            hs.hu_shu_exprk2_step,
            {
                "f": y,
                "y": y,
                "t": t,
                "dt": dt,
                "grid": self.grid,
                "x": self.x,
                "v": self.v,
                "eps0": self.eps0,
                "eps": self.eps0,
                "eps_x": self.eps_x,
                "case": self.case,
                "transport": self.transport,
                "transport_name": self.transport,
                "eta_model": self.eta_model,
            },
        )

    def transport_rhs(self, y: Array) -> Array:
        if not hasattr(hs, "get_transport_rhs"):
            raise AttributeError("hu_shu_expk2_bgk.py does not expose get_transport_rhs.")
        rhs_fun = call_with_known_kwargs(
            hs.get_transport_rhs,
            {
                "name": self.transport,
                "transport": self.transport,
                "transport_name": self.transport,
                "grid": self.grid,
            },
        )
        return call_with_known_kwargs(rhs_fun, {"f": y, "grid": self.grid})

    def transport_flow(self) -> Callable[[float, float, Array], Array]:
        def _flow(t: float, dt: float, y: Array) -> Array:
            rhs = self.transport_rhs(y)
            return y + dt * rhs
        return _flow

    def collision_flow(self) -> Callable[[float, float, Array], Array]:
        if not hasattr(hs, "bgk_homogeneous_exact"):
            raise AttributeError("hu_shu_expk2_bgk.py does not expose bgk_homogeneous_exact.")
        if not hasattr(hs, "moments"):
            raise AttributeError("hu_shu_expk2_bgk.py does not expose moments.")
        if not hasattr(hs, "maxwellian"):
            raise AttributeError("hu_shu_expk2_bgk.py does not expose maxwellian.")

        def _flow(t: float, dt: float, y: Array) -> Array:
            # Build local equilibrium g from the current state y
            rho, u, T = call_with_known_kwargs(
                hs.moments,
                {
                    "f": y,
                    "v": self.v,
                    "grid": self.grid,
                },
            )

            g = call_with_known_kwargs(
                hs.maxwellian,
                {
                    "rho": rho,
                    "u": u,
                    "T": T,
                    "v": self.v,
                    "grid": self.grid,
                },
            )

            return call_with_known_kwargs(
                hs.bgk_homogeneous_exact,
                {
                    "f": y,
                    "g": g,
                    "dt_collision": dt,
                    "dt": dt,
                    "grid": self.grid,
                    "x": self.x,
                    "v": self.v,
                    "eps0": self.eps0,
                    "eps": self.eps0,
                    "eps_x": self.eps_x,
                    "case": self.case,
                    "eta_model": self.eta_model,
                },
            )

        return _flow

    def null_third_flow(self) -> Callable[[float, float, Array], Array]:
        def _flow(t: float, dt: float, y: Array) -> Array:
            return np.array(y, copy=True)
        return _flow


def rk4_step(rhs: Callable[[Array, float], Array], y: Array, t: float, dt: float) -> Array:
    k1 = rhs(y, t)
    k2 = rhs(y + 0.5 * dt * k1, t + 0.5 * dt)
    k3 = rhs(y + 0.5 * dt * k2, t + 0.5 * dt)
    k4 = rhs(y + dt * k3, t + dt)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate_reference(rhs: Callable[[Array, float], Array], y0: Array, T: float, dt_ref: float) -> Array:
    y = np.array(y0, copy=True)
    t = 0.0
    while t < T - 1.0e-15:
        h = min(dt_ref, T - t)
        y = rk4_step(rhs, y, t, h)
        t += h
    return y


def make_three_operator_stepper(problem: HuShuProblem, method_name: str) -> Tuple[Callable[[Array, float, float], Array], int]:
    methods = {(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"}
    flows = [problem.transport_flow(), problem.collision_flow(), problem.null_third_flow()]
    order_guess = {
        "Strang-3": 2,
        "OS32_7op_minLEM-3": 2,
        "PP3_4A-3": 3,
        "Yoshida-3": 4,
    }.get(method_name, 2)

    def _step(y: Array, t: float, dt: float) -> Array:
        return fractional_step(
            functions=flows,
            delta_t=dt,
            initial_y=y,
            initial_t=t,
            final_t=t + dt,
            alpha=method_name,
            methods=methods,
        )

    return _step, order_guess


def integrate_fixed_step(stepper: Callable[[Array, float, float], Array], y0: Array, T: float, dt: float) -> Tuple[Array, int, float]:
    y = np.array(y0, copy=True)
    t = 0.0
    n_steps = 0
    tic = time.perf_counter()
    while t < T - 1.0e-15:
        h = min(dt, T - t)
        y = stepper(y, t, h)
        t += h
        n_steps += 1
    wall = time.perf_counter() - tic
    return y, n_steps, wall


def run_one_adaptive_case(case: str, method: str, T: float, dt0: float, atol: float, rtol: float,
                          dt_ref: float, Nx: int, Nv: int, xmax: float, vmax: float, eps0: float,
                          transport: str, eta_model: str, init: str, outdir: str) -> Dict[str, Any]:
    problem = HuShuProblem(case, Nx, Nv, xmax, vmax, eps0, transport, eta_model, init)
    y0 = np.array(problem.f0, copy=True)
    stepper, order = make_three_operator_stepper(problem, method)

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

    y_ref = integrate_reference(problem.full_rhs, y0, T, dt_ref)
    f_l2 = l2_norm(result.final_state - y_ref)
    rho_l1, u_l1, T_l1 = problem.macro_errors(result.final_state, y_ref)
    diag = problem.diagnostics(result.final_state)

    history_path = os.path.join(
        outdir,
        f"hu_shu_history_mode-adaptive_method-{method}_case-{case}_eps0-{eps0:.0e}_atol-{atol:.0e}_rtol-{rtol:.0e}.csv",
    )
    save_adaptive_history_csv(history_path, result)

    return {
        "mode": "adaptive",
        "case": case,
        "method": method,
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
        "f_l2": f_l2,
        "rho_l1": rho_l1,
        "u_l1": u_l1,
        "T_l1": T_l1,
        "min_f": diag.get("min_f", float("nan")),
        "rho_min": diag.get("rho_min", float("nan")),
        "u_min": diag.get("u_min", float("nan")),
        "T_min": diag.get("T_min", float("nan")),
    }


def run_one_fixed_case(case: str, method: str, T: float, fixed_dt: float, dt_ref: float,
                       Nx: int, Nv: int, xmax: float, vmax: float, eps0: float,
                       transport: str, eta_model: str, init: str) -> Dict[str, Any]:
    problem = HuShuProblem(case, Nx, Nv, xmax, vmax, eps0, transport, eta_model, init)
    y0 = np.array(problem.f0, copy=True)
    stepper, _ = make_three_operator_stepper(problem, method)

    y_final, n_steps, wall = integrate_fixed_step(stepper, y0, T, fixed_dt)
    y_ref = integrate_reference(problem.full_rhs, y0, T, dt_ref)
    f_l2 = l2_norm(y_final - y_ref)
    rho_l1, u_l1, T_l1 = problem.macro_errors(y_final, y_ref)
    diag = problem.diagnostics(y_final)

    return {
        "mode": "fixed",
        "case": case,
        "method": method,
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
        "f_l2": f_l2,
        "rho_l1": rho_l1,
        "u_l1": u_l1,
        "T_l1": T_l1,
        "min_f": diag.get("min_f", float("nan")),
        "rho_min": diag.get("rho_min", float("nan")),
        "u_min": diag.get("u_min", float("nan")),
        "T_min": diag.get("T_min", float("nan")),
    }


def run_exprk2_baseline(case: str, T: float, fixed_dt: float, dt_ref: float,
                        Nx: int, Nv: int, xmax: float, vmax: float, eps0: float,
                        transport: str, eta_model: str, init: str) -> Dict[str, Any]:
    problem = HuShuProblem(case, Nx, Nv, xmax, vmax, eps0, transport, eta_model, init)
    y0 = np.array(problem.f0, copy=True)

    y_final, n_steps, wall = integrate_fixed_step(problem.exprk2_step, y0, T, fixed_dt)
    y_ref = integrate_reference(problem.full_rhs, y0, T, dt_ref)
    f_l2 = l2_norm(y_final - y_ref)
    rho_l1, u_l1, T_l1 = problem.macro_errors(y_final, y_ref)
    diag = problem.diagnostics(y_final)

    return {
        "mode": "baseline",
        "case": case,
        "method": "ExpRK2",
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
        "message": "ExpRK2 baseline completed successfully.",
        "n_steps": n_steps,
        "n_accept": n_steps,
        "n_reject": 0,
        "dt_min_used": fixed_dt,
        "dt_max_used": fixed_dt,
        "dt_avg_used": fixed_dt,
        "wall_seconds": wall,
        "f_l2": f_l2,
        "rho_l1": rho_l1,
        "u_l1": u_l1,
        "T_l1": T_l1,
        "min_f": diag.get("min_f", float("nan")),
        "rho_min": diag.get("rho_min", float("nan")),
        "u_min": diag.get("u_min", float("nan")),
        "T_min": diag.get("T_min", float("nan")),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 3: Adaptive higher-order Hu–Shu / BGK study.")
    parser.add_argument("--mode", type=str, default="both", choices=["adaptive", "fixed", "baseline", "both"])
    parser.add_argument("--case", type=str, default="mixed", choices=["mixed", "uniform"])
    parser.add_argument("--methods", type=str, default="Strang-3,PP3_4A-3,Yoshida-3")
    parser.add_argument("--T", type=float, default=0.5)
    parser.add_argument("--dt0", type=float, default=1.0e-4)
    parser.add_argument("--baseline-dt", type=float, default=3.0e-5)
    parser.add_argument("--atols", type=str, default="1e-6,1e-8")
    parser.add_argument("--rtols", type=str, default="1e-4,1e-6")
    parser.add_argument("--fixed-dts", type=str, default="2.5e-4,1.25e-4,6.25e-5,3.125e-5")
    parser.add_argument("--dt-ref", type=float, default=1.0e-5)
    parser.add_argument("--Nx", type=int, default=40)
    parser.add_argument("--Nv", type=int, default=150)
    parser.add_argument("--xmax", type=float, default=2.0)
    parser.add_argument("--vmax", type=float, default=15.0)
    parser.add_argument("--eps0", type=float, default=1.0e-5)
    parser.add_argument("--transport", type=str, default="muscl2")
    parser.add_argument("--eta-model", type=str, default="constant1")
    parser.add_argument("--init", type=str, default="mixed")
    parser.add_argument("--outdir", type=str, default="experiments/outputs/hu_shu_adaptive_higher_order")
    parser.add_argument("--overwrite-summary", action="store_true")
    args = parser.parse_args()

    ensure_dir(args.outdir)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    atols = parse_float_list(args.atols)
    rtols = parse_float_list(args.rtols)
    fixed_dts = parse_float_list(args.fixed_dts)

    summary_path = os.path.join(args.outdir, "hu_shu_adaptive_higher_order_summary.csv")
    if args.overwrite_summary and os.path.exists(summary_path):
        os.remove(summary_path)

    print("=" * 72)
    print("PHASE 3: HU–SHU ADAPTIVE HIGHER-ORDER STUDY")
    print(f"mode        = {args.mode}")
    print(f"case        = {args.case}")
    print(f"methods     = {methods}")
    print(f"T           = {args.T}")
    print(f"dt0         = {args.dt0}")
    print(f"baseline_dt = {args.baseline_dt}")
    print(f"atols       = {atols}")
    print(f"rtols       = {rtols}")
    print(f"fixed_dts   = {fixed_dts}")
    print(f"dt_ref      = {args.dt_ref}")
    print(f"Nx, Nv      = {args.Nx}, {args.Nv}")
    print(f"xmax, vmax  = {args.xmax}, {args.vmax}")
    print(f"eps0        = {args.eps0}")
    print(f"transport   = {args.transport}")
    print(f"eta_model   = {args.eta_model}")
    print(f"init        = {args.init}")
    print(f"outdir      = {args.outdir}")
    print("=" * 72)

    if args.mode in ("baseline", "both"):
        print("\nBaseline: ExpRK2")
        row = run_exprk2_baseline(
            args.case, args.T, args.baseline_dt, args.dt_ref,
            args.Nx, args.Nv, args.xmax, args.vmax, args.eps0,
            args.transport, args.eta_model, args.init
        )
        append_summary_csv(summary_path, row)
        print(
            f"  steps={row['n_steps']} "
            f"f_l2={row['f_l2']:.6e} "
            f"rho_l1={row['rho_l1']:.6e} "
            f"u_l1={row['u_l1']:.6e} "
            f"T_l1={row['T_l1']:.6e} "
            f"min(f)={row['min_f']:.6e} "
            f"wall={row['wall_seconds']:.3f}s"
        )

    if args.mode in ("adaptive", "both"):
        for method in methods:
            for atol in atols:
                for rtol in rtols:
                    print(f"\nAdaptive: {method}, atol={atol:.0e}, rtol={rtol:.0e}")
                    try:
                        row = run_one_adaptive_case(
                            args.case, method, args.T, args.dt0, atol, rtol, args.dt_ref,
                            args.Nx, args.Nv, args.xmax, args.vmax, args.eps0,
                            args.transport, args.eta_model, args.init, args.outdir
                        )
                    except Exception as exc:
                        tb = traceback.format_exc()
                        row = {
                            "mode": "adaptive",
                            "case": args.case,
                            "method": method,
                            "T": args.T,
                            "dt0": args.dt0,
                            "fixed_dt": "",
                            "atol": atol,
                            "rtol": rtol,
                            "dt_ref": args.dt_ref,
                            "Nx": args.Nx,
                            "Nv": args.Nv,
                            "xmax": args.xmax,
                            "vmax": args.vmax,
                            "eps0": args.eps0,
                            "transport": args.transport,
                            "eta_model": args.eta_model,
                            "init": args.init,
                            "success": 0,
                            "message": tb,
                            "n_steps": 0,
                            "n_accept": 0,
                            "n_reject": 0,
                            "dt_min_used": float("nan"),
                            "dt_max_used": float("nan"),
                            "dt_avg_used": float("nan"),
                            "wall_seconds": 0.0,
                            "f_l2": float("inf"),
                            "rho_l1": float("inf"),
                            "u_l1": float("inf"),
                            "T_l1": float("inf"),
                            "min_f": float("nan"),
                            "rho_min": float("nan"),
                            "u_min": float("nan"),
                            "T_min": float("nan"),
                        }
                    append_summary_csv(summary_path, row)
                    if row["success"] == 0:
                        print(f"  FAILED: {row['message']}")
                    else:
                        print(
                            f"  success={row['success']} "
                            f"accept={row['n_accept']} reject={row['n_reject']} "
                            f"dt_avg={row['dt_avg_used']:.3e} "
                            f"f_l2={row['f_l2']:.6e} "
                            f"rho_l1={row['rho_l1']:.6e} "
                            f"u_l1={row['u_l1']:.6e} "
                            f"T_l1={row['T_l1']:.6e} "
                            f"min(f)={row['min_f']:.6e}"
                        )

    if args.mode in ("fixed", "both"):
        for method in methods:
            for fixed_dt in fixed_dts:
                print(f"\nFixed: {method}, dt={fixed_dt:.3e}")
                try:
                    row = run_one_fixed_case(
                        args.case, method, args.T, fixed_dt, args.dt_ref,
                        args.Nx, args.Nv, args.xmax, args.vmax, args.eps0,
                        args.transport, args.eta_model, args.init
                    )
                except Exception as exc:
                    row = {
                        "mode": "fixed", "case": args.case, "method": method, "T": args.T,
                        "dt0": "", "fixed_dt": fixed_dt, "atol": "", "rtol": "", "dt_ref": args.dt_ref,
                        "Nx": args.Nx, "Nv": args.Nv, "xmax": args.xmax, "vmax": args.vmax,
                        "eps0": args.eps0, "transport": args.transport, "eta_model": args.eta_model, "init": args.init,
                        "success": 0, "message": repr(exc), "n_steps": 0, "n_accept": 0, "n_reject": 0,
                        "dt_min_used": fixed_dt, "dt_max_used": fixed_dt, "dt_avg_used": fixed_dt,
                        "wall_seconds": 0.0, "f_l2": float("inf"), "rho_l1": float("inf"),
                        "u_l1": float("inf"), "T_l1": float("inf"), "min_f": float("nan"),
                        "rho_min": float("nan"), "u_min": float("nan"), "T_min": float("nan"),
                    }
                append_summary_csv(summary_path, row)
                if row["success"] == 0:
                    print(f"  FAILED: {row['message']}")
                else:
                    print(
                        f"  success={row['success']} "
                        f"steps={row['n_steps']} "
                        f"f_l2={row['f_l2']:.6e} "
                        f"rho_l1={row['rho_l1']:.6e} "
                        f"u_l1={row['u_l1']:.6e} "
                        f"T_l1={row['T_l1']:.6e} "
                        f"min(f)={row['min_f']:.6e}"
                    )

    print("\nDone.")
    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()
