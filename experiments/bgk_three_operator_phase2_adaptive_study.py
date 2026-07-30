from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

THIS = Path(__file__).resolve()
REPO = THIS.parents[1] if len(THIS.parents) > 1 else THIS.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from fractional_step import fractional_step  # noqa: E402


def parse_csv_floats(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_csv_strings(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def trapz_weights_uniform(n: int, h: float) -> np.ndarray:
    return np.ones(n, dtype=float) * h


def periodic_shift_1d(values: np.ndarray, shift_cells: float) -> np.ndarray:
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


def compute_moments(
    f: np.ndarray,
    v: np.ndarray,
    dv: float,
    rho_floor: float = 1e-14,
    T_floor: float = 1e-12,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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


@dataclass
class BGKPhase2Problem:
    x: np.ndarray
    v: np.ndarray
    dx: float
    dv: float
    force: float
    eps: float


def build_problem(Nx: int, Nv: int, xmax: float, vmax: float, force: float, eps: float) -> BGKPhase2Problem:
    x = np.linspace(0.0, xmax, Nx, endpoint=False)
    v = np.linspace(-vmax, vmax, Nv, endpoint=False)
    dx = xmax / Nx
    dv = 2.0 * vmax / Nv
    return BGKPhase2Problem(x=x, v=v, dx=dx, dv=dv, force=force, eps=eps)


def make_transport_flow(problem: BGKPhase2Problem):
    v = problem.v
    dx = problem.dx

    def flow(_t: float, h: float, f: np.ndarray) -> np.ndarray:
        shift_cells = (v * h) / dx
        return shift_along_x_periodic(f, shift_cells)

    return flow


def make_force_flow(problem: BGKPhase2Problem):
    dv = problem.dv
    F = problem.force

    def flow(_t: float, h: float, f: np.ndarray) -> np.ndarray:
        shift_cells = (F * h) / dv
        return shift_along_v_periodic(f, shift_cells)

    return flow


def make_bgk_relaxation_flow(problem: BGKPhase2Problem):
    v = problem.v
    dv = problem.dv
    eps = problem.eps

    def flow(_t: float, h: float, f: np.ndarray) -> np.ndarray:
        rho, u, T = compute_moments(f, v, dv)
        M = maxwellian_from_moments(rho, u, T, v)
        E = np.exp(-h / eps)
        return E * f + (1.0 - E) * M

    return flow


def total_mass(f: np.ndarray, dx: float, dv: float) -> float:
    return float(np.sum(f) * dx * dv)


def l1_norm(f: np.ndarray, dx: float, dv: float) -> float:
    return float(np.sum(np.abs(f)) * dx * dv)


def make_flows(problem: BGKPhase2Problem):
    return [
        make_transport_flow(problem),
        make_bgk_relaxation_flow(problem),
        make_force_flow(problem),
    ]


def bgk_split_step(method_name: str, problem: BGKPhase2Problem, f: np.ndarray, t: float, dt: float) -> np.ndarray:
    flows = make_flows(problem)
    f_new = fractional_step(
        functions=flows,
        delta_t=dt,
        initial_y=f,
        initial_t=t,
        final_t=t + dt,
        alpha=method_name,
        methods={(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"},
    )
    return np.array(f_new, dtype=float, copy=False)


def compute_reference(ref_method: str, ref_dt: float, T: float, f0: np.ndarray, problem: BGKPhase2Problem) -> np.ndarray:
    flows = make_flows(problem)
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


@dataclass
class AdaptiveOptions:
    order: int
    atol: float
    rtol: float
    dt_min: float = 1e-12
    dt_max: float = 1.0
    safety: float = 0.9
    growth_max: float = 2.0
    shrink_min: float = 0.2
    max_reject: int = 20
    positivity_tol: float = -1e-12
    explosion_factor: float = 1e8


@dataclass
class AttemptRecord:
    t: float
    dt_try: float
    dt_suggested: float
    accepted: bool
    err_est: float
    cpu_seconds: float
    min_f: float
    max_f: float
    l1_ratio: float
    stable: bool


@dataclass
class AdaptiveRunMetrics:
    err: float
    min_f_over_run: float
    final_min_f: float
    final_mass: float
    mass_rel_drift: float
    max_l1_ratio: float
    positive: bool
    stable: bool
    n_accept: int
    n_reject: int
    dt_min_used: float
    dt_avg_used: float
    dt_max_used: float
    wall_seconds: float
    message: str
    history: List[AttemptRecord]


def weighted_rms_error(y_high: np.ndarray, y_low: np.ndarray, atol: float, rtol: float) -> float:
    scale = atol + rtol * np.maximum(np.abs(y_high), np.abs(y_low))
    diff = (y_high - y_low) / scale
    return float(np.sqrt(np.mean(diff * diff)))


def propose_new_dt(err: float, dt: float, opts: AdaptiveOptions) -> float:
    if not np.isfinite(err):
        factor = opts.shrink_min
    elif err <= 0.0:
        factor = opts.growth_max
    else:
        factor = opts.safety * (err ** (-1.0 / (opts.order + 1.0)))
        factor = min(opts.growth_max, max(opts.shrink_min, factor))
    return min(opts.dt_max, max(opts.dt_min, dt * factor))


def bgk_diagnostics(f: np.ndarray, problem: BGKPhase2Problem, l10: float) -> Tuple[float, float, float, bool]:
    min_f = float(np.min(f))
    max_f = float(np.max(f))
    l1r = l1_norm(f, problem.dx, problem.dv) / max(l10, 1e-30)
    stable = bool(np.all(np.isfinite(f)) and l1r <= 1e8)
    return min_f, max_f, l1r, stable


def run_adaptive_method(
    method_name: str,
    order: int,
    dt0: float,
    T: float,
    f0: np.ndarray,
    problem: BGKPhase2Problem,
    reference: np.ndarray,
    opts: AdaptiveOptions,
) -> AdaptiveRunMetrics:
    f = np.array(f0, dtype=float, copy=True)
    t = 0.0
    dt = min(max(dt0, opts.dt_min), opts.dt_max)
    m0 = total_mass(f0, problem.dx, problem.dv)
    l10 = l1_norm(f0, problem.dx, problem.dv)
    min_f_over_run = float(np.min(f))
    max_l1_ratio = 1.0
    accepted_dts: List[float] = []
    history: List[AttemptRecord] = []
    n_accept = 0
    n_reject = 0
    start = time.perf_counter()

    while t < T - 1e-15:
        dt = min(dt, T - t)
        local_rejects = 0

        while True:
            tic = time.perf_counter()
            full = bgk_split_step(method_name, problem, f, t, dt)
            half = bgk_split_step(method_name, problem, f, t, 0.5 * dt)
            two_half = bgk_split_step(method_name, problem, half, t + 0.5 * dt, 0.5 * dt)
            cpu = time.perf_counter() - tic

            min_f, max_f, l1r, stable = bgk_diagnostics(two_half, problem, l10)
            err = weighted_rms_error(two_half, full, opts.atol, opts.rtol)
            if not np.isfinite(err):
                stable = False
            if min_f < opts.positivity_tol:
                stable = False
                err = float("inf")

            dt_new = propose_new_dt(err, dt, opts)
            accepted = bool(stable and err <= 1.0)

            history.append(AttemptRecord(
                t=t,
                dt_try=dt,
                dt_suggested=dt_new,
                accepted=accepted,
                err_est=err,
                cpu_seconds=cpu,
                min_f=min_f,
                max_f=max_f,
                l1_ratio=l1r,
                stable=stable,
            ))

            if accepted:
                f = two_half
                t += dt
                n_accept += 1
                accepted_dts.append(dt)
                min_f_over_run = min(min_f_over_run, min_f)
                max_l1_ratio = max(max_l1_ratio, l1r)
                dt = dt_new
                break

            n_reject += 1
            local_rejects += 1
            dt = dt_new
            if dt <= opts.dt_min + 1e-30:
                wall = time.perf_counter() - start
                final_min_f = float(np.min(f))
                final_mass = total_mass(f, problem.dx, problem.dv)
                mass_rel_drift = abs(final_mass - m0) / max(abs(m0), 1e-30)
                err_final = l1_norm(f - reference, problem.dx, problem.dv)
                return AdaptiveRunMetrics(
                    err=err_final,
                    min_f_over_run=min_f_over_run,
                    final_min_f=final_min_f,
                    final_mass=final_mass,
                    mass_rel_drift=mass_rel_drift,
                    max_l1_ratio=max_l1_ratio,
                    positive=bool(final_min_f >= opts.positivity_tol),
                    stable=False,
                    n_accept=n_accept,
                    n_reject=n_reject,
                    dt_min_used=min(accepted_dts) if accepted_dts else 0.0,
                    dt_avg_used=float(np.mean(accepted_dts)) if accepted_dts else 0.0,
                    dt_max_used=max(accepted_dts) if accepted_dts else 0.0,
                    wall_seconds=wall,
                    message="Step size reached dt_min during rejection.",
                    history=history,
                )
            if local_rejects >= opts.max_reject:
                wall = time.perf_counter() - start
                final_min_f = float(np.min(f))
                final_mass = total_mass(f, problem.dx, problem.dv)
                mass_rel_drift = abs(final_mass - m0) / max(abs(m0), 1e-30)
                err_final = l1_norm(f - reference, problem.dx, problem.dv)
                return AdaptiveRunMetrics(
                    err=err_final,
                    min_f_over_run=min_f_over_run,
                    final_min_f=final_min_f,
                    final_mass=final_mass,
                    mass_rel_drift=mass_rel_drift,
                    max_l1_ratio=max_l1_ratio,
                    positive=bool(final_min_f >= opts.positivity_tol),
                    stable=False,
                    n_accept=n_accept,
                    n_reject=n_reject,
                    dt_min_used=min(accepted_dts) if accepted_dts else 0.0,
                    dt_avg_used=float(np.mean(accepted_dts)) if accepted_dts else 0.0,
                    dt_max_used=max(accepted_dts) if accepted_dts else 0.0,
                    wall_seconds=wall,
                    message="Exceeded maximum rejects for one advance.",
                    history=history,
                )

    wall = time.perf_counter() - start
    final_min_f = float(np.min(f))
    final_mass = total_mass(f, problem.dx, problem.dv)
    mass_rel_drift = abs(final_mass - m0) / max(abs(m0), 1e-30)
    err_final = l1_norm(f - reference, problem.dx, problem.dv)
    return AdaptiveRunMetrics(
        err=err_final,
        min_f_over_run=min_f_over_run,
        final_min_f=final_min_f,
        final_mass=final_mass,
        mass_rel_drift=mass_rel_drift,
        max_l1_ratio=max_l1_ratio,
        positive=bool(final_min_f >= opts.positivity_tol),
        stable=True,
        n_accept=n_accept,
        n_reject=n_reject,
        dt_min_used=min(accepted_dts) if accepted_dts else 0.0,
        dt_avg_used=float(np.mean(accepted_dts)) if accepted_dts else 0.0,
        dt_max_used=max(accepted_dts) if accepted_dts else 0.0,
        wall_seconds=wall,
        message="Adaptive integration completed successfully.",
        history=history,
    )


def run_fixed_method(method_name: str, dt: float, T: float, f0: np.ndarray, problem: BGKPhase2Problem, reference: np.ndarray,
                     positivity_tol: float = -1e-12, explosion_factor: float = 1e8):
    steps = max(1, int(round(T / dt)))
    dt = T / steps
    f = np.array(f0, dtype=float, copy=True)
    t = 0.0
    m0 = total_mass(f0, problem.dx, problem.dv)
    l10 = l1_norm(f0, problem.dx, problem.dv)
    min_f_over_run = float(np.min(f))
    max_l1_ratio = 1.0
    stable = True
    tic = time.perf_counter()
    for _ in range(steps):
        f_new = bgk_split_step(method_name, problem, f, t, dt)
        if not np.all(np.isfinite(f_new)):
            stable = False
            f = f_new
            break
        min_f_over_run = min(min_f_over_run, float(np.min(f_new)))
        l1r = l1_norm(f_new, problem.dx, problem.dv) / max(l10, 1e-30)
        max_l1_ratio = max(max_l1_ratio, l1r)
        if l1r > explosion_factor:
            stable = False
            f = f_new
            break
        f = f_new
        t += dt
    wall = time.perf_counter() - tic
    final_min_f = float(np.min(f)) if np.all(np.isfinite(f)) else float("nan")
    final_mass = total_mass(f, problem.dx, problem.dv) if np.all(np.isfinite(f)) else float("nan")
    mass_rel_drift = abs(final_mass - m0) / max(abs(m0), 1e-30) if np.isfinite(final_mass) else float("inf")
    positive = bool(final_min_f >= positivity_tol) if np.isfinite(final_min_f) else False
    err = l1_norm(f - reference, problem.dx, problem.dv) if stable and np.all(np.isfinite(f)) else float("inf")
    return {
        "err": err,
        "min_f_over_run": min_f_over_run,
        "final_min_f": final_min_f,
        "final_mass": final_mass,
        "mass_rel_drift": mass_rel_drift,
        "max_l1_ratio": max_l1_ratio,
        "positive": positive,
        "stable": stable,
        "n_steps": steps,
        "wall_seconds": wall,
        "dt": dt,
    }


def method_order(method_name: str) -> int:
    name = method_name.lower()
    if "strang" in name:
        return 2
    return 3


def save_history_csv(path: Path, history: List[AttemptRecord]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "attempt_index", "t", "dt_try", "dt_suggested", "accepted", "err_est", "cpu_seconds",
            "min_f", "max_f", "l1_ratio", "stable"
        ])
        for i, rec in enumerate(history):
            writer.writerow([
                i, rec.t, rec.dt_try, rec.dt_suggested, int(rec.accepted), rec.err_est, rec.cpu_seconds,
                rec.min_f, rec.max_f, rec.l1_ratio, int(rec.stable)
            ])


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Adaptive BGK-inspired three-operator higher-order study.")
    ap.add_argument("--mode", type=str, default="adaptive", choices=["adaptive", "fixed", "both"])
    ap.add_argument("--eps-list", type=str, default="1e-1,1e-2,1e-3")
    ap.add_argument("--fixed-dts", type=str, default="1e-2,5e-3,2.5e-3")
    ap.add_argument("--methods", type=str, default="Strang-3,OS32_7op_minLEM-3,PP3_4A-3,Yoshida-3")
    ap.add_argument("--reference-method", type=str, default="Yoshida-3")
    ap.add_argument("--reference-refine", type=int, default=8)
    ap.add_argument("--T", type=float, default=0.2)
    ap.add_argument("--dt0", type=float, default=1e-3)
    ap.add_argument("--atols", type=str, default="1e-6,1e-8")
    ap.add_argument("--rtols", type=str, default="1e-4,1e-6")
    ap.add_argument("--Nx", type=int, default=40)
    ap.add_argument("--Nv", type=int, default=80)
    ap.add_argument("--xmax", type=float, default=1.0)
    ap.add_argument("--vmax", type=float, default=8.0)
    ap.add_argument("--force", type=float, default=0.5)
    ap.add_argument("--profile", type=str, choices=["mixed", "simple"], default="mixed")
    ap.add_argument("--outdir", type=str, default="experiments/bgk_three_operator_adaptive_outputs")
    ap.add_argument("--overwrite-summary", action="store_true")
    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    eps_list = parse_csv_floats(args.eps_list)
    fixed_dts = parse_csv_floats(args.fixed_dts)
    methods = parse_csv_strings(args.methods)
    atols = parse_csv_floats(args.atols)
    rtols = parse_csv_floats(args.rtols)

    outdir = REPO / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "bgk_three_operator_adaptive_comparison_summary.csv"
    if args.overwrite_summary and summary_path.exists():
        summary_path.unlink()

    print("Running adaptive BGK-inspired three-operator study with:")
    print(f"  mode             = {args.mode}")
    print(f"  eps_list         = {eps_list}")
    print(f"  fixed_dts        = {fixed_dts}")
    print(f"  methods          = {methods}")
    print(f"  reference_method = {args.reference_method}")
    print(f"  reference_refine = {args.reference_refine}")
    print(f"  T                = {args.T}")
    print(f"  dt0              = {args.dt0}")
    print(f"  atols            = {atols}")
    print(f"  rtols            = {rtols}")
    print(f"  Nx, Nv           = {args.Nx}, {args.Nv}")
    print(f"  xmax, vmax       = {args.xmax}, {args.vmax}")
    print(f"  force            = {args.force}")
    print(f"  profile          = {args.profile}")
    print(f"  outdir           = {outdir}")

    wrote_header = False
    for eps in eps_list:
        problem = build_problem(args.Nx, args.Nv, args.xmax, args.vmax, args.force, eps)
        f0 = initial_condition(problem.x, problem.v, profile=args.profile)
        dt_min = min(fixed_dts + [args.dt0])
        ref_dt = dt_min / args.reference_refine
        reference = compute_reference(args.reference_method, ref_dt, args.T, f0, problem)

        print("\n" + "=" * 88)
        print(f"REGIME eps = {eps:.3e} | T = {args.T:g} | ref_dt = {ref_dt:.4e}")
        print("=" * 88)

        rows: List[dict] = []
        for method_name in methods:
            if args.mode in ("adaptive", "both"):
                order = method_order(method_name)
                for atol in atols:
                    for rtol in rtols:
                        opts = AdaptiveOptions(order=order, atol=atol, rtol=rtol, dt_max=max(args.dt0, args.T))
                        metrics = run_adaptive_method(method_name, order, args.dt0, args.T, f0, problem, reference, opts)
                        print(
                            f"Adaptive | {method_name:18s} | atol={atol:.0e} rtol={rtol:.0e} | "
                            f"accept={metrics.n_accept:4d} reject={metrics.n_reject:4d} | "
                            f"dt[min,avg,max]=({metrics.dt_min_used:.3e}, {metrics.dt_avg_used:.3e}, {metrics.dt_max_used:.3e}) | "
                            f"err={metrics.err:.5e} | stable={metrics.stable!s:5s} | positive={metrics.positive!s:5s} | "
                            f"min_f={metrics.min_f_over_run:.5e}"
                        )
                        hist_name = f"bgk_phase2_history_eps{eps:.0e}_{method_name}_atol{atol:.0e}_rtol{rtol:.0e}.csv".replace("/", "-")
                        save_history_csv(outdir / hist_name, metrics.history)
                        rows.append({
                            "mode": "adaptive",
                            "eps": eps,
                            "T": args.T,
                            "method": method_name,
                            "fixed_dt": "",
                            "dt0": args.dt0,
                            "atol": atol,
                            "rtol": rtol,
                            "ref_method": args.reference_method,
                            "ref_dt": ref_dt,
                            "Nx": args.Nx,
                            "Nv": args.Nv,
                            "xmax": args.xmax,
                            "vmax": args.vmax,
                            "force": args.force,
                            "profile": args.profile,
                            "err_l1": metrics.err,
                            "stable": int(metrics.stable),
                            "positive": int(metrics.positive),
                            "min_f_over_run": metrics.min_f_over_run,
                            "final_min_f": metrics.final_min_f,
                            "final_mass": metrics.final_mass,
                            "mass_rel_drift": metrics.mass_rel_drift,
                            "max_l1_ratio": metrics.max_l1_ratio,
                            "n_steps": metrics.n_accept,
                            "n_accept": metrics.n_accept,
                            "n_reject": metrics.n_reject,
                            "dt_min_used": metrics.dt_min_used,
                            "dt_avg_used": metrics.dt_avg_used,
                            "dt_max_used": metrics.dt_max_used,
                            "wall_seconds": metrics.wall_seconds,
                            "message": metrics.message,
                        })

            if args.mode in ("fixed", "both"):
                for dt in fixed_dts:
                    metrics = run_fixed_method(method_name, dt, args.T, f0, problem, reference)
                    print(
                        f"Fixed    | {method_name:18s} | dt={metrics['dt']:.3e} | "
                        f"steps={metrics['n_steps']:4d} | err={metrics['err']:.5e} | stable={metrics['stable']!s:5s} | "
                        f"positive={metrics['positive']!s:5s} | min_f={metrics['min_f_over_run']:.5e}"
                    )
                    rows.append({
                        "mode": "fixed",
                        "eps": eps,
                        "T": args.T,
                        "method": method_name,
                        "fixed_dt": metrics["dt"],
                        "dt0": "",
                        "atol": "",
                        "rtol": "",
                        "ref_method": args.reference_method,
                        "ref_dt": ref_dt,
                        "Nx": args.Nx,
                        "Nv": args.Nv,
                        "xmax": args.xmax,
                        "vmax": args.vmax,
                        "force": args.force,
                        "profile": args.profile,
                        "err_l1": metrics["err"],
                        "stable": int(metrics["stable"]),
                        "positive": int(metrics["positive"]),
                        "min_f_over_run": metrics["min_f_over_run"],
                        "final_min_f": metrics["final_min_f"],
                        "final_mass": metrics["final_mass"],
                        "mass_rel_drift": metrics["mass_rel_drift"],
                        "max_l1_ratio": metrics["max_l1_ratio"],
                        "n_steps": metrics["n_steps"],
                        "n_accept": metrics["n_steps"],
                        "n_reject": 0,
                        "dt_min_used": metrics["dt"],
                        "dt_avg_used": metrics["dt"],
                        "dt_max_used": metrics["dt"],
                        "wall_seconds": metrics["wall_seconds"],
                        "message": "Fixed-step integration completed successfully.",
                    })

        with summary_path.open("a", newline="") as f:
            fieldnames = list(rows[0].keys()) if rows else []
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            if not wrote_header and fieldnames:
                writer.writeheader()
                wrote_header = True
            for row in rows:
                writer.writerow(row)

    print(f"\nSaved summary: {summary_path}")


if __name__ == "__main__":
    main()
