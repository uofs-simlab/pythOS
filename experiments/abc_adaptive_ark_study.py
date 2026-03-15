#!/usr/bin/env python3
"""
experiments/abc_adaptive_ark_study.py

Adaptive time-stepping study for:
  - Strang-3 splitting (2nd order)
  - Yoshida-3 splitting (4th order composition)
  - ARK / embedded explicit RK baseline (Heun–Euler fallback)

Writes per-run CSV logs and (optional) dt-vs-t plots.

Key fix:
  ARK CSV output is now *strictly rectangular* (same number of columns each row),
  so pandas/csv parsers don't choke and ARK plots are no longer blank.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np


# -----------------------------
# Utilities
# -----------------------------

def safe_norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(x, ord=2))


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def adapt_dt(dt: float, err: float, tol: float, p: int, safety: float = 0.9,
             fac_min: float = 0.2, fac_max: float = 5.0) -> float:
    """
    Standard PI-free controller:
        dt_new = dt * safety * (tol/err)^(1/(p+1))
    For embedded RK: p is the order of the *high* method.
    """
    if err <= 0.0 or not np.isfinite(err):
        fac = fac_max
    else:
        fac = safety * (tol / err) ** (1.0 / (p + 1.0))
        fac = clamp(fac, fac_min, fac_max)
    return dt * fac


# -----------------------------
# ABC test problem
# -----------------------------

def make_abc_matrices(Bmult: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Simple 3x3 linear system y' = (A+B+C) y with tunable stiffness via Bmult.
    Designed so that different splittings/orderings show measurable differences.
    """
    # Mild "rotation/shear"
    A = np.array([
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.2],
        [0.0, -0.2, 0.0],
    ], dtype=float)

    # Stiff dissipative part scaled by Bmult
    B = Bmult * np.array([
        [-1.0, 0.0, 0.0],
        [0.0, -5.0, 0.0],
        [0.0, 0.0, -10.0],
    ], dtype=float)

    # Small coupling / drift
    C = np.array([
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0],
        [-0.5, 0.0, 0.0],
    ], dtype=float)

    return A, B, C


def exact_solution(A: np.ndarray, B: np.ndarray, C: np.ndarray, y0: np.ndarray, T: float) -> np.ndarray:
    """
    For the small linear system, use matrix exponential via eigendecomposition.
    """
    M = A + B + C
    w, V = np.linalg.eig(M)
    Vinv = np.linalg.inv(V)
    expw = np.diag(np.exp(w * T))
    yT = V @ expw @ (Vinv @ y0)
    # if small imaginary noise due to eig, drop it
    yT = np.real_if_close(yT, tol=1e8)
    return np.array(yT, dtype=float)


def expm_action(M: np.ndarray, y: np.ndarray, h: float) -> np.ndarray:
    w, V = np.linalg.eig(M)
    Vinv = np.linalg.inv(V)
    expw = np.diag(np.exp(w * h))
    out = V @ expw @ (Vinv @ y)
    out = np.real_if_close(out, tol=1e8)
    return np.array(out, dtype=float)


def flow_A(A: np.ndarray, y: np.ndarray, h: float) -> np.ndarray:
    return expm_action(A, y, h)


def flow_B(B: np.ndarray, y: np.ndarray, h: float) -> np.ndarray:
    return expm_action(B, y, h)


def flow_C(C: np.ndarray, y: np.ndarray, h: float) -> np.ndarray:
    return expm_action(C, y, h)


def apply_flow(letter: str, A: np.ndarray, B: np.ndarray, C: np.ndarray, y: np.ndarray, h: float) -> np.ndarray:
    if letter == "A":
        return flow_A(A, y, h)
    if letter == "B":
        return flow_B(B, y, h)
    if letter == "C":
        return flow_C(C, y, h)
    raise ValueError(f"Unknown operator letter: {letter}")


def strang_step(ordering: str, A: np.ndarray, B: np.ndarray, C: np.ndarray, y: np.ndarray, h: float) -> np.ndarray:
    """
    "Strang-3" for three operators:
      X(h/2) Y(h/2) Z(h) Y(h/2) X(h/2)   with ordering = XYZ
    """
    X, Y, Z = ordering
    y1 = apply_flow(X, A, B, C, y, 0.5 * h)
    y2 = apply_flow(Y, A, B, C, y1, 0.5 * h)
    y3 = apply_flow(Z, A, B, C, y2, 1.0 * h)
    y4 = apply_flow(Y, A, B, C, y3, 0.5 * h)
    y5 = apply_flow(X, A, B, C, y4, 0.5 * h)
    return y5


def yoshida4_step(ordering: str, A: np.ndarray, B: np.ndarray, C: np.ndarray, y: np.ndarray, h: float) -> np.ndarray:
    """
    Yoshida 4th-order composition of Strang:
      S(w1 h) S(w0 h) S(w1 h)
    with:
      w1 = 1/(2-2^(1/3)), w0 = -2^(1/3)/(2-2^(1/3))
    """
    w1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
    w0 = - (2.0 ** (1.0 / 3.0)) / (2.0 - 2.0 ** (1.0 / 3.0))
    y1 = strang_step(ordering, A, B, C, y, w1 * h)
    y2 = strang_step(ordering, A, B, C, y1, w0 * h)
    y3 = strang_step(ordering, A, B, C, y2, w1 * h)
    return y3


# -----------------------------
# Adaptive log row
# -----------------------------

@dataclass
class AdaptiveLogRow:
    t: float
    dt: float
    accepted: int
    err_est: float
    y_norm_ratio: float
    step_ratio: float


def integrate_adaptive_splitting(
    method: str,
    ordering: str,
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    y0: np.ndarray,
    T: float,
    dt0: float,
    tol: float,
    out_path: str,
) -> Tuple[int, int, float, float, float, float]:
    """
    Adaptive splitting using step-doubling error estimate:
      err ~ || y(h) - y(h/2)∘y(h/2) ||
    """
    y_exact_T = exact_solution(A, B, C, y0, T)
    exact_norm = max(safe_norm(y_exact_T), 1e-300)
    y0_norm = max(safe_norm(y0), 1e-300)

    if method == "Strang3":
        stepper = strang_step
        p = 2
    elif method == "Yoshida3":
        stepper = yoshida4_step
        p = 4
    else:
        raise ValueError("Unknown splitting method")

    t = 0.0
    dt = float(dt0)
    y = y0.copy()

    accept = 0
    reject = 0
    max_y_ratio = 1.0
    max_step_ratio = 1.0
    rows: List[AdaptiveLogRow] = []

    while t < T - 1e-15:
        dt = clamp(dt, 1e-12, 0.5)
        if t + dt > T:
            dt = T - t

        # one step of size dt
        y1 = stepper(ordering, A, B, C, y, dt)

        # two half-steps
        yhalf = stepper(ordering, A, B, C, y, 0.5 * dt)
        y2 = stepper(ordering, A, B, C, yhalf, 0.5 * dt)

        err_est = safe_norm(y1 - y2)
        step_ratio = safe_norm(y2) / max(safe_norm(y), 1e-300)

        if err_est <= tol:
            # accept the more accurate y2
            t += dt
            y = y2
            accept += 1

            y_ratio = safe_norm(y) / y0_norm
            max_y_ratio = max(max_y_ratio, y_ratio)
            max_step_ratio = max(max_step_ratio, step_ratio)

            rows.append(AdaptiveLogRow(t=t, dt=dt, accepted=1, err_est=err_est,
                                      y_norm_ratio=y_ratio, step_ratio=step_ratio))
            dt = adapt_dt(dt, err_est, tol, p=p)

        else:
            reject += 1
            rows.append(AdaptiveLogRow(t=t, dt=dt, accepted=0, err_est=err_est,
                                      y_norm_ratio=safe_norm(y) / y0_norm, step_ratio=step_ratio))
            dt = adapt_dt(dt, err_est, tol, p=p)
            if reject > 20000:
                break

    final_err = safe_norm(y - y_exact_T) / exact_norm
    final_y_ratio = safe_norm(y) / y0_norm

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["t", "dt", "accepted", "err_est", "y_norm_ratio", "step_ratio"])
        for r in rows:
            w.writerow([f"{r.t:.16e}", f"{r.dt:.16e}", int(r.accepted),
                        f"{r.err_est:.16e}", f"{r.y_norm_ratio:.16e}", f"{r.step_ratio:.16e}"])

    return accept, reject, final_err, max_y_ratio, final_y_ratio, max_step_ratio


# -----------------------------
# ARK baseline (embedded explicit RK)
# -----------------------------

@dataclass
class EmbeddedPair:
    A: np.ndarray
    b_high: np.ndarray
    b_low: np.ndarray
    c: np.ndarray
    order_high: int
    name: str


def builtin_heun_euler_pair() -> EmbeddedPair:
    """
    Heun (explicit trapezoid) as high method (order 2),
    Euler as low method (order 1).
    """
    A = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
    ], dtype=float)
    c = np.array([0.0, 1.0], dtype=float)
    b_high = np.array([0.5, 0.5], dtype=float)  # Heun
    b_low = np.array([1.0, 0.0], dtype=float)   # Euler
    return EmbeddedPair(A=A, b_high=b_high, b_low=b_low, c=c, order_high=2, name="Heun-Euler (builtin)")


def extract_embedded_pair_from_additive_rk() -> Optional[EmbeddedPair]:
    """
    Try to pull a usable embedded pair from additive_rk.embedded_pairs.
    Falls back to None if not available / incompatible.
    """
    try:
        import additive_rk  # type: ignore
    except Exception:
        return None

    if not hasattr(additive_rk, "embedded_pairs"):
        return None
    pairs = getattr(additive_rk, "embedded_pairs")
    if not isinstance(pairs, dict) or len(pairs) == 0:
        return None

    preferred = ["Heun-Euler", "Bogacki-Shampine", "Dormand-Prince", "RK23", "RK45"]
    key = None
    for k in preferred:
        if k in pairs:
            key = k
            break
    if key is None:
        key = list(pairs.keys())[0]

    pair = pairs[key]

    A = getattr(pair, "A", None)
    c = getattr(pair, "c", None)
    b_hi = getattr(pair, "b_high", None) or getattr(pair, "b", None) or getattr(pair, "bhat", None)
    b_lo = getattr(pair, "b_low", None) or getattr(pair, "btilde", None) or getattr(pair, "b_embedded", None)
    order_hi = getattr(pair, "order_high", None) or getattr(pair, "p_high", None) or getattr(pair, "order", None)

    if A is None or c is None or b_hi is None or b_lo is None:
        return None

    A = np.array(A, dtype=float)
    c = np.array(c, dtype=float)
    b_hi = np.array(b_hi, dtype=float)
    b_lo = np.array(b_lo, dtype=float)

    if order_hi is None:
        # if unknown, default to 2 (safe, conservative controller)
        order_hi = 2

    return EmbeddedPair(A=A, b_high=b_hi, b_low=b_lo, c=c, order_high=int(order_hi), name=str(key))


def run_ark_baseline(
    Bmult: float,
    y0: np.ndarray,
    T: float,
    dt0: float,
    tol: float,
    outdir: str,
) -> Optional[str]:
    """
    Adaptive embedded explicit RK on y' = (A+B+C) y.
    Writes ONLY rectangular numeric CSV rows (no summary row).
    Summary info is printed and additionally written to a separate .txt.
    """
    Aop, Bop, Cop = make_abc_matrices(Bmult)
    M = Aop + Bop + Cop

    pair = extract_embedded_pair_from_additive_rk()
    if pair is None:
        pair = builtin_heun_euler_pair()
        print("[ARK] No embedded_pairs found in additive_rk; using Heun-Euler (builtin).")
    else:
        print(f"[ARK] Using embedded pair key: {pair.name}")

    s = len(pair.c)

    def f(_t: float, y: np.ndarray) -> np.ndarray:
        return M @ y

    def rk_step_embedded(t: float, y: np.ndarray, h: float) -> Tuple[np.ndarray, float, float]:
        K = np.zeros((s, len(y)), dtype=float)
        for i in range(s):
            yi = y.copy()
            for j in range(i):
                aij = pair.A[i, j]
                if aij != 0.0:
                    yi = yi + h * aij * K[j]
            K[i] = f(t + pair.c[i] * h, yi)

        y_hi = y + h * np.sum(pair.b_high[:, None] * K, axis=0)
        y_lo = y + h * np.sum(pair.b_low[:, None] * K, axis=0)

        err = safe_norm(y_hi - y_lo)  # absolute error estimate
        step_ratio = safe_norm(y_hi) / max(safe_norm(y), 1e-300)
        return y_hi, err, float(step_ratio)

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"ark_Bmult{int(Bmult)}.csv")
    summary_path = os.path.join(outdir, f"ark_Bmult{int(Bmult)}_summary.txt")

    y_exact_T = exact_solution(Aop, Bop, Cop, y0, T)
    exact_norm = max(safe_norm(y_exact_T), 1e-300)
    y0_norm = max(safe_norm(y0), 1e-300)

    t = 0.0
    dt = float(dt0)
    y = y0.copy()

    accept = 0
    reject = 0
    max_y_ratio = 1.0
    max_step_ratio = 1.0
    rows: List[AdaptiveLogRow] = []

    while t < T - 1e-15:
        dt = clamp(dt, 1e-12, 0.5)
        if t + dt > T:
            dt = T - t

        y_hi, err_est, step_ratio = rk_step_embedded(t, y, dt)

        if err_est <= tol:
            t += dt
            y = y_hi
            accept += 1

            y_ratio = safe_norm(y) / y0_norm
            max_y_ratio = max(max_y_ratio, y_ratio)
            max_step_ratio = max(max_step_ratio, step_ratio)

            rows.append(AdaptiveLogRow(t=t, dt=dt, accepted=1, err_est=err_est,
                                      y_norm_ratio=y_ratio, step_ratio=step_ratio))
            dt = adapt_dt(dt, err_est, tol, p=pair.order_high)

        else:
            reject += 1
            rows.append(AdaptiveLogRow(t=t, dt=dt, accepted=0, err_est=err_est,
                                      y_norm_ratio=safe_norm(y) / y0_norm, step_ratio=step_ratio))
            dt = adapt_dt(dt, err_est, tol, p=pair.order_high)
            if reject > 20000:
                break

    final_err = safe_norm(y - y_exact_T) / exact_norm
    final_y_ratio = safe_norm(y) / y0_norm

    # Rectangular numeric CSV only
    with open(out_path, "w", newline="") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["t", "dt", "accepted", "err_est", "y_norm_ratio", "step_ratio"])
        for r in rows:
            w.writerow([f"{r.t:.16e}", f"{r.dt:.16e}", int(r.accepted),
                        f"{r.err_est:.16e}", f"{r.y_norm_ratio:.16e}", f"{r.step_ratio:.16e}"])

    # Summary separately (won't break parsing)
    with open(summary_path, "w") as fs:
        fs.write(f"pair={pair.name}\n")
        fs.write(f"final_err={final_err:.16e}\n")
        fs.write(f"accept={accept}\n")
        fs.write(f"reject={reject}\n")
        fs.write(f"max_y_ratio={max_y_ratio:.16e}\n")
        fs.write(f"final_y_ratio={final_y_ratio:.16e}\n")
        fs.write(f"max_step_ratio={max_step_ratio:.16e}\n")

    print(f"[ARK] saved: {out_path}")
    return out_path


# -----------------------------
# Plotting
# -----------------------------

def plot_dt_vs_t(csv_path: str, title: str, out_png: str) -> None:
    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_csv(csv_path)

    # Robust numeric coercion
    df["t"] = pd.to_numeric(df["t"], errors="coerce")
    df["dt"] = pd.to_numeric(df["dt"], errors="coerce")
    df["accepted"] = pd.to_numeric(df["accepted"], errors="coerce")

    df = df.dropna(subset=["t", "dt", "accepted"])
    df_acc = df[df["accepted"] == 1].copy()

    plt.figure()
    if len(df_acc) == 0:
        # fallback: plot all rows if nothing accepted (shouldn't happen, but safe)
        plt.plot(df["t"].values, df["dt"].values, marker="o", linewidth=1)
    else:
        plt.plot(df_acc["t"].values, df_acc["dt"].values, marker="o", linewidth=1)

    plt.yscale("log")
    plt.xlabel("t")
    plt.ylabel("dt")
    plt.title(title)
    plt.tight_layout()

    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=150)
    plt.close()


# -----------------------------
# Main driver
# -----------------------------

def parse_list_floats(s: str) -> List[float]:
    parts = [p.strip() for p in s.split(",") if p.strip()]
    return [float(p) for p in parts]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--Bmults", type=str, default="1,20,50,100")
    ap.add_argument("--T", type=float, default=0.5)
    ap.add_argument("--dt0", type=float, default=1e-2)
    ap.add_argument("--tol", type=float, default=1e-7)
    ap.add_argument("--ordering", type=str, default="CAB")
    ap.add_argument("--outdir", type=str, default="experiments/abc_adaptive_outputs")
    ap.add_argument("--plots", action="store_true")
    ap.add_argument("--do_ark", action="store_true")
    args = ap.parse_args()

    Bmults = parse_list_floats(args.Bmults)
    T = float(args.T)
    dt0 = float(args.dt0)
    tol = float(args.tol)
    ordering = args.ordering.strip().upper()
    outdir = args.outdir

    print("========================================================================")
    print(f"Adaptive study: Bmults={Bmults}  T={T}  dt0={dt0}  tol={tol}  ordering={ordering}")
    print("========================================================================")

    # fixed initial condition for reproducibility
    y0 = np.array([1.0, 0.2, -0.1], dtype=float)

    plot_dir = os.path.join(outdir, "plots")

    for Bmult in Bmults:
        print("\n========================================================================")
        print(f"Adaptive study: Bmult={int(Bmult)}  T={T}  dt0={dt0}  tol={tol}  ordering={ordering}")
        print("========================================================================")

        A, B, C = make_abc_matrices(Bmult)

        # Strang-3 adaptive
        strang_csv = os.path.join(outdir, f"adaptive_Strang3_Bmult{int(Bmult)}_{ordering}.csv")
        acc, rej, ferr, my, fy, ms = integrate_adaptive_splitting(
            method="Strang3",
            ordering=ordering,
            A=A, B=B, C=C,
            y0=y0, T=T, dt0=dt0, tol=tol,
            out_path=strang_csv,
        )
        print(f"[Strang3] accept={acc} reject={rej}  final_err={ferr:.6e}  "
              f"max||y||/||y0||={my:.6e}  final||y||/||y0||={fy:.6e}  max_step={ms:.6e}")
        print(f"  saved: {strang_csv}")

        # Yoshida-3 adaptive (4th order)
        yosh_csv = os.path.join(outdir, f"adaptive_Yoshida3_Bmult{int(Bmult)}_{ordering}.csv")
        acc, rej, ferr, my, fy, ms = integrate_adaptive_splitting(
            method="Yoshida3",
            ordering=ordering,
            A=A, B=B, C=C,
            y0=y0, T=T, dt0=dt0, tol=tol,
            out_path=yosh_csv,
        )
        print(f"[Yoshida3] accept={acc} reject={rej}  final_err={ferr:.6e}  "
              f"max||y||/||y0||={my:.6e}  final||y||/||y0||={fy:.6e}  max_step={ms:.6e}")
        print(f"  saved: {yosh_csv}")

        # ARK baseline
        ark_csv = None
        if args.do_ark:
            ark_csv = run_ark_baseline(Bmult=Bmult, y0=y0, T=T, dt0=dt0, tol=tol, outdir=outdir)

        # Plots
        if args.plots:
            plot_dt_vs_t(
                strang_csv,
                title=f"Adaptive dt vs t | Strang-3 | Bmult={int(Bmult)} | {ordering}",
                out_png=os.path.join(plot_dir, f"dt_vs_t_Strang3_Bmult{int(Bmult)}_{ordering}.png"),
            )
            plot_dt_vs_t(
                yosh_csv,
                title=f"Adaptive dt vs t | Yoshida-3 | Bmult={int(Bmult)} | {ordering}",
                out_png=os.path.join(plot_dir, f"dt_vs_t_Yoshida3_Bmult{int(Bmult)}_{ordering}.png"),
            )
            if ark_csv is not None:
                plot_dt_vs_t(
                    ark_csv,
                    title=f"Adaptive dt vs t | ARK (embedded) | Bmult={int(Bmult)}",
                    out_png=os.path.join(plot_dir, f"dt_vs_t_ARK_Bmult{int(Bmult)}.png"),
                )

    if args.plots:
        print(f"Plots saved to: {plot_dir}")


if __name__ == "__main__":
    main()