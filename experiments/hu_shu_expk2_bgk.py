#!/usr/bin/env python3
"""
hu_shu_expk2_bgk.py

Reconstructed BGK experiment driver using the actual Hu–Shu ExpRK2 scheme:

  f^(0)   = exp((a0 dt / eps) Q) f^n
  f^(1)   = exp((a1 dt / eps) Q) ( f^(0) + b1 dt T(f^(0)) )
  f^(2)   = f^(1) + b2 dt T(f^(1))
  f^(n+1) = exp((a2 dt / eps) Q) [ w f^(2) + (1-w) exp(((1-a2) dt / eps) Q) f^n ]

with the Hu–Shu coefficients:
  w = 1/2, b1 = b2 = 1, a0 = a1 = a2 = 1/3

For BGK:
  Q(f) = eta ( M[f] - f )

and the homogeneous solve is exact:
  exp(s Q) g = exp(-eta s) g + (1 - exp(-eta s)) M[g].

Paper references:
- scheme form: Eq. (2.6)
- recommended coefficients: Eqs. (2.23) and (2.26)
- BGK homogeneous solve: Eq. (4.5)
- accuracy test initial data: Eqs. (7.1)–(7.3)
- mixed-regime epsilon(x): Eq. (7.6)

Practical note:
- This script reconstructs the *time integrator* from the paper.
- The paper's spatial transport discretization used WENO5 + positivity limiter.
  Here we keep simpler upwind / MUSCL2 options for experimentation.
"""

import argparse
import math
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple, Union

import numpy as np


ArrayLike = np.ndarray


@dataclass
class Grid:
    Nx: int
    Nv: int
    xmax: float
    vmax: float
    x: np.ndarray
    v: np.ndarray
    dx: float
    dv: float


# -----------------------------------------------------------------------------
# Grid / moments / Maxwellian
# -----------------------------------------------------------------------------

def make_grid(Nx: int, Nv: int, xmax: float = 2.0, vmax: float = 15.0) -> Grid:
    dx = xmax / Nx
    dv = 2.0 * vmax / Nv
    x = (np.arange(Nx) + 0.5) * dx
    v = -vmax + (np.arange(Nv) + 0.5) * dv
    return Grid(Nx=Nx, Nv=Nv, xmax=xmax, vmax=vmax, x=x, v=v, dx=dx, dv=dv)


def moments(f: ArrayLike, grid: Grid, min_temperature: float = 1e-14) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    1D velocity moments:
      rho = ∫ f dv
      u   = (1/rho) ∫ v f dv
      T   = (1/rho) ∫ (v-u)^2 f dv
    """
    v = grid.v[None, :]
    rho = np.sum(f, axis=1) * grid.dv
    rho_safe = np.maximum(rho, 1e-14)

    mom = np.sum(f * v, axis=1) * grid.dv
    u = mom / rho_safe

    centered2 = (v - u[:, None]) ** 2
    T = np.sum(f * centered2, axis=1) * grid.dv / rho_safe
    T = np.maximum(T, min_temperature)
    return rho, u, T


def maxwellian(rho: np.ndarray, u: np.ndarray, T: np.ndarray, grid: Grid) -> ArrayLike:
    """
    1D Maxwellian:
      M = rho / sqrt(2 pi T) * exp(-(v-u)^2 / (2T))
    """
    v = grid.v[None, :]
    rho = rho[:, None]
    u = u[:, None]
    T = np.maximum(T[:, None], 1e-14)

    coef = rho / np.sqrt(2.0 * np.pi * T)
    expo = np.exp(-0.5 * (v - u) ** 2 / T)
    return coef * expo


def collision_frequency(rho: np.ndarray, T: np.ndarray, model: str = "constant1") -> np.ndarray:
    """
    BGK collision frequency eta(rho, T).

    The paper allows eta to be any positive function of rho and T.
    For a clean baseline, default to eta = 1.

    Options:
      - constant1 : eta = 1
      - rho       : eta = rho
      - rho_sqrtT : eta = rho * sqrt(T)
    """
    if model == "constant1":
        return np.ones_like(rho)
    if model == "rho":
        return np.maximum(rho, 1e-14)
    if model == "rho_sqrtT":
        return np.maximum(rho * np.sqrt(np.maximum(T, 1e-14)), 1e-14)
    raise ValueError(f"Unknown eta model: {model}")


# -----------------------------------------------------------------------------
# Paper-aligned initial data / epsilon profiles
# -----------------------------------------------------------------------------

def initial_condition_accuracy(grid: Grid) -> ArrayLike:
    """
    Paper accuracy-test initial data, Eqs. (7.1)-(7.2):
      f(0,x,v) = 0.5 M_{rho,u,T} + 0.3 M_{rho,-0.5u,T}
      rho = 1 + 0.2 sin(pi x), u = 1, T = 1 / (1 + 0.2 sin(pi x))
    """
    x = grid.x
    rho = 1.0 + 0.2 * np.sin(np.pi * x)
    u = np.ones_like(x)
    T = 1.0 / (1.0 + 0.2 * np.sin(np.pi * x))

    M1 = maxwellian(rho, u, T, grid)
    M2 = maxwellian(rho, -0.5 * u, T, grid)
    f0 = 0.5 * M1 + 0.3 * M2
    return np.maximum(f0, 1e-14)


def initial_condition_simple(grid: Grid) -> ArrayLike:
    """
    A fallback smooth non-equilibrium initial condition.
    """
    x = grid.x
    rho = 1.0 + 0.1 * np.sin(2.0 * np.pi * x / grid.xmax)
    u = 0.5 * np.sin(2.0 * np.pi * x / grid.xmax)
    T = 1.0 + 0.1 * np.cos(2.0 * np.pi * x / grid.xmax)
    M = maxwellian(rho, u, T, grid)

    v = grid.v[None, :]
    perturb = 0.05 * np.sin(2.0 * np.pi * x / grid.xmax)[:, None] * (v ** 2 - 1.0) * np.exp(-0.1 * v ** 2)
    f0 = M * (1.0 + perturb)
    return np.maximum(f0, 1e-14)


def initial_condition(grid: Grid, kind: str) -> ArrayLike:
    if kind in ("accuracy", "paper", "mixed"):
        return initial_condition_accuracy(grid)
    if kind == "simple":
        return initial_condition_simple(grid)
    raise ValueError(f"Unknown initial condition kind: {kind}")


def constant_eps(grid: Grid, eps_value: float) -> np.ndarray:
    return np.full(grid.Nx, float(eps_value))


def mixed_regime_eps(grid: Grid, eps0: float = 1e-5) -> np.ndarray:
    """
    Paper mixed-regime epsilon, Eq. (7.6):
      eps(x) = eps0 + tanh(1 - 11(x-1)) + tanh(1 + 11(x-1))
    """
    x = grid.x
    return eps0 + np.tanh(1.0 - 11.0 * (x - 1.0)) + np.tanh(1.0 + 11.0 * (x - 1.0))


# -----------------------------------------------------------------------------
# Transport operators
# -----------------------------------------------------------------------------

def minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    out = np.zeros_like(a)
    mask = (a * b) > 0.0
    out[mask] = np.sign(a[mask]) * np.minimum(np.abs(a[mask]), np.abs(b[mask]))
    return out


def mc_limiter(dl: np.ndarray, dr: np.ndarray) -> np.ndarray:
    return minmod(0.5 * (dl + dr), minmod(2.0 * dl, 2.0 * dr))


def transport_rhs_upwind(f: ArrayLike, grid: Grid) -> ArrayLike:
    """
    First-order upwind discretization of T(f) = -v d_x f.
    Periodic in x.
    """
    rhs = np.zeros_like(f)
    for j, vj in enumerate(grid.v):
        q = f[:, j]
        if vj >= 0.0:
            dfdx = (q - np.roll(q, 1)) / grid.dx
        else:
            dfdx = (np.roll(q, -1) - q) / grid.dx
        rhs[:, j] = -vj * dfdx
    return rhs


def transport_rhs_muscl2(f: ArrayLike, grid: Grid) -> ArrayLike:
    """
    Second-order MUSCL finite-volume discretization of T(f) = -v d_x f.
    Periodic in x.
    """
    rhs = np.zeros_like(f)

    for j, vj in enumerate(grid.v):
        q = f[:, j]

        qm1 = np.roll(q, 1)
        qp1 = np.roll(q, -1)

        dl = q - qm1
        dr = qp1 - q
        slope = mc_limiter(dl, dr)

        qL = q + 0.5 * slope
        qR = q - 0.5 * slope

        # interface i+1/2
        qL_ip = qL
        qR_ip = np.roll(qR, -1)

        if vj >= 0.0:
            flux_ip = vj * qL_ip
        else:
            flux_ip = vj * qR_ip

        flux_im = np.roll(flux_ip, 1)
        rhs[:, j] = -(flux_ip - flux_im) / grid.dx

    return rhs


def get_transport_rhs(name: str) -> Callable[[ArrayLike, Grid], ArrayLike]:
    if name == "upwind":
        return transport_rhs_upwind
    if name == "muscl2":
        return transport_rhs_muscl2
    raise ValueError(f"Unknown transport method: {name}")


# -----------------------------------------------------------------------------
# BGK pieces
# -----------------------------------------------------------------------------

def bgk_homogeneous_exact(g: ArrayLike, dt_collision: float, grid: Grid, eps_x: np.ndarray, eta_model: str) -> ArrayLike:
    """
    Exact homogeneous BGK solve for:
      d_t f = (1/eps) Q(f),  Q(f) = eta (M[f] - f)

    Since rho,u,T are conserved by the homogeneous BGK flow, M[g] and eta[g]
    stay constant over this substep, giving:
      exp((dt/eps)Q) g = e^{-eta dt / eps} g + (1 - e^{-eta dt / eps}) M[g]
    """
    rho, u, T = moments(g, grid)
    M = maxwellian(rho, u, T, grid)
    eta = collision_frequency(rho, T, model=eta_model)
    alpha = np.exp(-eta * dt_collision / np.maximum(eps_x, 1e-14))
    return alpha[:, None] * g + (1.0 - alpha)[:, None] * M


def collision_rhs(f: ArrayLike, grid: Grid, eps_x: np.ndarray, eta_model: str) -> ArrayLike:
    rho, u, T = moments(f, grid)
    eta = collision_frequency(rho, T, model=eta_model)
    M = maxwellian(rho, u, T, grid)
    return (eta[:, None] / np.maximum(eps_x, 1e-14)[:, None]) * (M - f)


def full_rhs(f: ArrayLike, grid: Grid, eps_x: np.ndarray, transport_name: str, eta_model: str) -> ArrayLike:
    T_rhs = get_transport_rhs(transport_name)
    return T_rhs(f, grid) + collision_rhs(f, grid, eps_x, eta_model)


# -----------------------------------------------------------------------------
# Hu–Shu ExpRK2 from paper Eq. (2.6) with Eqs. (2.23), (2.26)
# -----------------------------------------------------------------------------

def hu_shu_exprk2_step(
    f: ArrayLike,
    dt: float,
    grid: Grid,
    eps_x: np.ndarray,
    transport_name: str,
    eta_model: str = "constant1",
) -> ArrayLike:
    """
    Reconstructed Hu–Shu ExpRK2 BGK step.

    Paper coefficients:
      w = 1/2, b1 = b2 = 1, a0 = a1 = a2 = 1/3
    """
    a0 = 1.0 / 3.0
    a1 = 1.0 / 3.0
    a2 = 1.0 / 3.0
    b1 = 1.0
    b2 = 1.0
    w = 0.5

    T_rhs = get_transport_rhs(transport_name)

    # f^(0) = exp((a0 dt / eps)Q) f^n
    f0 = bgk_homogeneous_exact(f, a0 * dt, grid, eps_x, eta_model)

    # f^(1) = exp((a1 dt / eps)Q) ( f^(0) + b1 dt T(f^(0)) )
    g1 = f0 + b1 * dt * T_rhs(f0, grid)
    f1 = bgk_homogeneous_exact(g1, a1 * dt, grid, eps_x, eta_model)

    # f^(2) = f^(1) + b2 dt T(f^(1))
    f2 = f1 + b2 * dt * T_rhs(f1, grid)

    # fn+1 = exp((a2 dt / eps)Q) [ w f^(2) + (1-w) exp(((1-a2)dt/eps)Q) f^n ]
    tail = bgk_homogeneous_exact(f, (1.0 - a2) * dt, grid, eps_x, eta_model)
    combo = w * f2 + (1.0 - w) * tail
    fn1 = bgk_homogeneous_exact(combo, a2 * dt, grid, eps_x, eta_model)

    return fn1


# -----------------------------------------------------------------------------
# Explicit reference methods
# -----------------------------------------------------------------------------

def ssp_rk2_step(
    f: ArrayLike,
    dt: float,
    grid: Grid,
    eps_x: np.ndarray,
    transport_name: str,
    eta_model: str,
) -> ArrayLike:
    rhs1 = full_rhs(f, grid, eps_x, transport_name, eta_model)
    f1 = f + dt * rhs1
    rhs2 = full_rhs(f1, grid, eps_x, transport_name, eta_model)
    return 0.5 * f + 0.5 * (f1 + dt * rhs2)


def rk4_step(
    f: ArrayLike,
    dt: float,
    grid: Grid,
    eps_x: np.ndarray,
    transport_name: str,
    eta_model: str,
) -> ArrayLike:
    k1 = full_rhs(f, grid, eps_x, transport_name, eta_model)
    k2 = full_rhs(f + 0.5 * dt * k1, grid, eps_x, transport_name, eta_model)
    k3 = full_rhs(f + 0.5 * dt * k2, grid, eps_x, transport_name, eta_model)
    k4 = full_rhs(f + dt * k3, grid, eps_x, transport_name, eta_model)
    return f + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def advance_one_step(
    f: ArrayLike,
    dt: float,
    grid: Grid,
    eps_x: np.ndarray,
    method: str,
    transport_name: str,
    eta_model: str,
) -> ArrayLike:
    if method == "exprk2":
        return hu_shu_exprk2_step(f, dt, grid, eps_x, transport_name, eta_model)
    if method == "ssp_rk2":
        return ssp_rk2_step(f, dt, grid, eps_x, transport_name, eta_model)
    if method == "rk4":
        return rk4_step(f, dt, grid, eps_x, transport_name, eta_model)
    raise ValueError(f"Unknown method: {method}")


def evolve(
    f0: ArrayLike,
    Tfinal: float,
    dt: float,
    grid: Grid,
    eps_x: np.ndarray,
    method: str,
    transport_name: str,
    eta_model: str,
) -> ArrayLike:
    f = f0.copy()
    t = 0.0
    while t < Tfinal - 1e-15:
        dt_step = min(dt, Tfinal - t)
        f = advance_one_step(f, dt_step, grid, eps_x, method, transport_name, eta_model)
        t += dt_step
    return f


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

def l1_x(a: np.ndarray, b: np.ndarray, dx: float) -> float:
    return float(np.sum(np.abs(a - b)) * dx)


def l2_xv(a: np.ndarray, b: np.ndarray, dx: float, dv: float) -> float:
    return float(np.sqrt(np.sum((a - b) ** 2) * dx * dv))


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


def observed_orders(errors: List[float]) -> List[float]:
    out = [np.nan]
    for i in range(1, len(errors)):
        if errors[i] <= 0.0 or errors[i - 1] <= 0.0:
            out.append(np.nan)
        else:
            out.append(math.log(errors[i - 1] / errors[i], 2.0))
    return out


def estimate_transport_cfl_dt(grid: Grid) -> float:
    vmax = np.max(np.abs(grid.v))
    return grid.dx / max(vmax, 1e-14)


def print_run_summary(method: str, transport: str, dt: float, elapsed: float, errs: Dict[str, float]) -> None:
    print(f"method={method}, transport={transport}, dt={dt:.6e}")
    print(f"runtime   = {elapsed:.3f} s")
    print(f"L2 f      = {errs['f_L2']:.6e}")
    print(f"L1 rho    = {errs['rho_L1']:.6e}")
    print(f"L1 u      = {errs['u_L1']:.6e}")
    print(f"L1 T      = {errs['T_L1']:.6e}")
    print(f"min(f)    = {errs['min_f']:.6e}")
    print(f"max(f)    = {errs['max_f']:.6e}")
    if errs["min_f"] < -1e-12:
        print("WARNING: negative values detected.")
    print("-" * 72)


# -----------------------------------------------------------------------------
# Experiment cases
# -----------------------------------------------------------------------------

def run_timeorder(args: argparse.Namespace) -> None:
    grid = make_grid(args.Nx, args.Nv, xmax=args.xmax, vmax=args.vmax)
    f0 = initial_condition(grid, args.init)
    eps_x = constant_eps(grid, args.eps)

    dt_list = [float(s.strip()) for s in args.dt_list.split(",")]
    dt_ref = min(dt_list) / args.refine_factor

    cfl_dt = estimate_transport_cfl_dt(grid)
    print("=" * 72)
    print("TIME-ORDER STUDY")
    print(f"reference method = {args.reference_method}")
    print(f"reference dt     = {dt_ref:.6e}")
    print(f"eps              = {args.eps:.3e}")
    print(f"Tfinal           = {args.T:.3f}")
    print(f"heuristic dx/|v|max CFL dt ≈ {cfl_dt:.6e}")
    print("=" * 72)

    t0 = time.perf_counter()
    fref = evolve(f0, args.T, dt_ref, grid, eps_x, args.reference_method, args.transport, args.eta_model)
    tref = time.perf_counter() - t0
    print(f"Reference done in {tref:.3f} s")
    print("-" * 72)

    f_errors = []
    rho_errors = []
    u_errors = []
    T_errors = []

    for dt in dt_list:
        if dt > cfl_dt:
            print(f"WARNING: dt={dt:.3e} exceeds heuristic transport CFL {cfl_dt:.3e}")

        t1 = time.perf_counter()
        f = evolve(f0, args.T, dt, grid, eps_x, args.method, args.transport, args.eta_model)
        elapsed = time.perf_counter() - t1

        errs = compare_solution(f, fref, grid)
        print_run_summary(args.method, args.transport, dt, elapsed, errs)

        f_errors.append(errs["f_L2"])
        rho_errors.append(errs["rho_L1"])
        u_errors.append(errs["u_L1"])
        T_errors.append(errs["T_L1"])

    ord_f = observed_orders(f_errors)
    ord_rho = observed_orders(rho_errors)
    ord_u = observed_orders(u_errors)
    ord_T = observed_orders(T_errors)

    print("Observed orders")
    print(f"{'dt':>12} {'f_L2':>14} {'ord':>8} {'rho_L1':>14} {'ord':>8} {'u_L1':>14} {'ord':>8} {'T_L1':>14} {'ord':>8}")
    for i, dt in enumerate(dt_list):
        print(
            f"{dt:12.4e} "
            f"{f_errors[i]:14.6e} {ord_f[i]:8.3f} "
            f"{rho_errors[i]:14.6e} {ord_rho[i]:8.3f} "
            f"{u_errors[i]:14.6e} {ord_u[i]:8.3f} "
            f"{T_errors[i]:14.6e} {ord_T[i]:8.3f}"
        )


def run_table(args: argparse.Namespace) -> None:
    grid = make_grid(args.Nx, args.Nv, xmax=args.xmax, vmax=args.vmax)
    f0 = initial_condition(grid, args.init)
    eps_list = [float(s.strip()) for s in args.eps_list.split(",")]

    dt = args.dt
    dt_ref = dt / args.refine_factor

    print("=" * 72)
    print("EPSILON TABLE STUDY")
    print(f"method           = {args.method}")
    print(f"reference method = {args.reference_method}")
    print(f"dt               = {dt:.3e}")
    print(f"Tfinal           = {args.T:.3f}")
    print("=" * 72)
    print(f"{'eps':>12} {'f_L2':>14} {'rho_L1':>14} {'u_L1':>14} {'T_L1':>14} {'min(f)':>14} {'runtime(s)':>12}")

    for eps_val in eps_list:
        eps_x = constant_eps(grid, eps_val)
        fref = evolve(f0, args.T, dt_ref, grid, eps_x, args.reference_method, args.transport, args.eta_model)

        t0 = time.perf_counter()
        f = evolve(f0, args.T, dt, grid, eps_x, args.method, args.transport, args.eta_model)
        elapsed = time.perf_counter() - t0

        errs = compare_solution(f, fref, grid)
        print(
            f"{eps_val:12.4e} "
            f"{errs['f_L2']:14.6e} "
            f"{errs['rho_L1']:14.6e} "
            f"{errs['u_L1']:14.6e} "
            f"{errs['T_L1']:14.6e} "
            f"{errs['min_f']:14.6e} "
            f"{elapsed:12.3f}"
        )


def run_mixed(args: argparse.Namespace) -> None:
    grid = make_grid(args.Nx, args.Nv, xmax=args.xmax, vmax=args.vmax)
    f0 = initial_condition(grid, "mixed")
    eps_x = mixed_regime_eps(grid, eps0=args.eps0)

    print("=" * 72)
    print("MIXED-REGIME STUDY")
    print(f"method           = {args.method}")
    print(f"reference method = {args.reference_method}")
    print(f"transport        = {args.transport}")
    print(f"Tfinal           = {args.T:.3f}")
    print(f"eps0             = {args.eps0:.3e}")
    print(f"dt               = {args.dt:.3e}")
    print("=" * 72)

    dt_ref = args.dt / args.refine_factor
    t0 = time.perf_counter()
    fref = evolve(f0, args.T, dt_ref, grid, eps_x, args.reference_method, args.transport, args.eta_model)
    tref = time.perf_counter() - t0
    print(f"Reference done in {tref:.3f} s")
    print("-" * 72)

    t1 = time.perf_counter()
    f = evolve(f0, args.T, args.dt, grid, eps_x, args.method, args.transport, args.eta_model)
    elapsed = time.perf_counter() - t1

    errs = compare_solution(f, fref, grid)
    print_run_summary(args.method, args.transport, args.dt, elapsed, errs)

    rho, u, T = moments(f, grid)
    rho_ref, u_ref, T_ref = moments(fref, grid)

    print("Ranges and max-errors vs reference:")
    print(f"rho range        = [{rho.min():.6e}, {rho.max():.6e}]")
    print(f"u range          = [{u.min():.6e}, {u.max():.6e}]")
    print(f"T range          = [{T.min():.6e}, {T.max():.6e}]")
    print(f"max |rho-rho_ref| = {np.max(np.abs(rho - rho_ref)):.6e}")
    print(f"max |u-u_ref|     = {np.max(np.abs(u - u_ref)):.6e}")
    print(f"max |T-T_ref|     = {np.max(np.abs(T - T_ref)):.6e}")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Hu–Shu ExpRK2 BGK experiments reconstructed from the paper.")

    p.add_argument("--case", choices=["timeorder", "table", "mixed"], default="timeorder")
    p.add_argument("--method", choices=["exprk2", "ssp_rk2", "rk4"], default="exprk2")
    p.add_argument("--reference-method", choices=["exprk2", "ssp_rk2", "rk4"], default="exprk2")
    p.add_argument("--transport", choices=["upwind", "muscl2"], default="muscl2")
    p.add_argument("--eta-model", choices=["constant1", "rho", "rho_sqrtT"], default="constant1")

    p.add_argument("--Nx", type=int, default=160, help="Number of x cells.")
    p.add_argument("--Nv", type=int, default=150, help="Number of v cells.")
    p.add_argument("--xmax", type=float, default=2.0, help="Domain length in x.")
    p.add_argument("--vmax", type=float, default=15.0, help="Velocity truncation.")
    p.add_argument("--T", type=float, default=0.1, help="Final time.")
    p.add_argument("--dt", type=float, default=2.5e-4, help="Main dt for table/mixed.")
    p.add_argument("--dt-list", dest="dt_list", type=str, default="4e-4,2e-4,1e-4,5e-5",
                   help="Comma-separated dt list for timeorder.")
    p.add_argument("--refine-factor", type=int, default=4,
                   help="Reference dt = min(dt_list)/refine_factor or dt/refine_factor.")

    p.add_argument("--eps", type=float, default=1e-2, help="Uniform epsilon for timeorder.")
    p.add_argument("--eps-list", type=str, default="1e0,1e-1,1e-2,1e-3,1e-4",
                   help="Comma-separated epsilon list for table.")
    p.add_argument("--eps0", type=float, default=1e-5, help="Minimum epsilon for mixed regime.")
    p.add_argument("--init", choices=["accuracy", "paper", "mixed", "simple"], default="accuracy")

    return p


def main() -> None:
    args = build_parser().parse_args()

    print("Running BGK experiment with:")
    print(f"  case             = {args.case}")
    print(f"  method           = {args.method}")
    print(f"  reference_method = {args.reference_method}")
    print(f"  transport        = {args.transport}")
    print(f"  eta_model        = {args.eta_model}")
    print(f"  Nx, Nv           = {args.Nx}, {args.Nv}")
    print(f"  xmax, vmax       = {args.xmax}, {args.vmax}")
    print(f"  T                = {args.T}")
    print(f"  dt               = {args.dt}")
    print()

    if args.case == "timeorder":
        run_timeorder(args)
    elif args.case == "table":
        run_table(args)
    elif args.case == "mixed":
        run_mixed(args)
    else:
        raise ValueError(f"Unknown case: {args.case}")


if __name__ == "__main__":
    main()