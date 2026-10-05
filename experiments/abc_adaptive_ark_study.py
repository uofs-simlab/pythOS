#!/usr/bin/env python3
"""
abc_adaptive_ark_study.py

Adaptive step-size study for:
  - Strang-3 splitting (2nd order)
  - Yoshida-3 splitting (4th order)
  - Embedded RK baseline using pythOS butcher_tableau embedded pairs

Important note:
`additive_rk.py` in pythOS imports and uses tableau objects, but the actual
embedded RK pair definitions live in `butcher_tableau.py` under `embedded_pairs`.
So this script:
  - verifies that `additive_rk.py` is importable from pythOS
  - uses `butcher_tableau.embedded_pairs[...]` to obtain the pair
  - runs the adaptive RK baseline on the full RHS y' = (A+B+C)y
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------
# Repository-local imports
# ---------------------------------------------------------------------
THIS = Path(__file__).resolve()
REPO = THIS.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


# ----------------------------
# Utility helpers
# ----------------------------

def safe_norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(x))


def clamp(x: float, lo: float, hi: float) -> float:
    return float(max(lo, min(hi, x)))


def adapt_dt(dt: float, err_est: float, tol: float, p: int, safety: float = 0.9) -> float:
    """
    Standard adaptive controller:
        dt_new = dt * safety * (tol / err_est)^(1/(p+1))

    where p is the order of the high solution.
    """
    if err_est <= 0.0 or not np.isfinite(err_est):
        return min(2.0 * dt, 10.0 * dt)

    fac = safety * (tol / err_est) ** (1.0 / (p + 1.0))
    fac = clamp(fac, 0.2, 5.0)
    return dt * fac


# ----------------------------
# ABC toy problem definition
# ----------------------------

def make_abc_matrices(Bmult: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Linear test problem:
        y' = (A + B + C) y
    with tunable stiffness in B via Bmult.
    """
    A = np.array([
        [-1.0,  0.2,  0.0],
        [ 0.0, -0.5,  0.1],
        [ 0.0,  0.0, -0.2],
    ], dtype=float)

    B = Bmult * np.array([
        [-2.0,  1.0,  0.0],
        [-1.0, -2.0,  0.0],
        [ 0.0,  0.0, -1.0],
    ], dtype=float)

    C = np.array([
        [-0.1,  0.0,  0.0],
        [ 0.0,  0.0,  0.3],
        [ 0.0, -0.3,  0.0],
    ], dtype=float)

    return A, B, C


def exact_solution(A: np.ndarray, B: np.ndarray, C: np.ndarray, y0: np.ndarray, T: float) -> np.ndarray:
    """
    Exact solution using matrix exponential.
    Uses scipy if available; otherwise falls back to eigendecomposition.
    """
    M = A + B + C
    try:
        from scipy.linalg import expm
        return expm(M * T) @ y0
    except Exception:
        w, V = np.linalg.eig(M)
        Vinv = np.linalg.inv(V)
        expw = np.diag(np.exp(w * T))
        y = (V @ expw @ Vinv) @ y0
        return np.real_if_close(y, tol=1000)


# ----------------------------
# Splitting building blocks
# ----------------------------

def flow_linear(M: np.ndarray, y: np.ndarray, h: float) -> np.ndarray:
    """
    Exact flow for y' = M y over time h.
    """
    try:
        from scipy.linalg import expm
        return expm(M * h) @ y
    except Exception:
        w, V = np.linalg.eig(M)
        Vinv = np.linalg.inv(V)
        expw = np.diag(np.exp(w * h))
        out = (V @ expw @ Vinv) @ y
        return np.real_if_close(out, tol=1000)


def apply_ordering(A: np.ndarray, B: np.ndarray, C: np.ndarray, ordering: str) -> List[np.ndarray]:
    """
    ordering is a string like 'CAB' specifying subflow order.
    """
    mp = {"A": A, "B": B, "C": C}
    return [mp[ch] for ch in ordering]


def strang_step(A: np.ndarray, B: np.ndarray, C: np.ndarray, y: np.ndarray, h: float, ordering: str) -> np.ndarray:
    """
    Strang splitting for 3 operators:
        phi_X(h/2) phi_Y(h/2) phi_Z(h) phi_Y(h/2) phi_X(h/2)
    """
    X, Y, Z = apply_ordering(A, B, C, ordering)
    y1 = flow_linear(X, y, 0.5 * h)
    y2 = flow_linear(Y, y1, 0.5 * h)
    y3 = flow_linear(Z, y2, h)
    y4 = flow_linear(Y, y3, 0.5 * h)
    y5 = flow_linear(X, y4, 0.5 * h)
    return y5


def yoshida_4th_from_strang(A: np.ndarray, B: np.ndarray, C: np.ndarray, y: np.ndarray, h: float, ordering: str) -> np.ndarray:
    """
    4th-order Yoshida composition from Strang:
        S(w1 h) S(w2 h) S(w1 h)
    """
    cbrt2 = 2.0 ** (1.0 / 3.0)
    w1 = 1.0 / (2.0 - cbrt2)
    w2 = -cbrt2 / (2.0 - cbrt2)

    y1 = strang_step(A, B, C, y, w1 * h, ordering)
    y2 = strang_step(A, B, C, y1, w2 * h, ordering)
    y3 = strang_step(A, B, C, y2, w1 * h, ordering)
    return y3


# ----------------------------
# Adaptive CSV logging
# ----------------------------

@dataclass
class AdaptiveLogRow:
    t: float
    dt: float
    accepted: int
    err_est: float
    y_norm_ratio: float
    step_ratio: float


def write_adaptive_csv(path: str, rows: List[AdaptiveLogRow], summary_lines: List[str]) -> None:
    """
    Writes a rectangular CSV with comment-header summary lines.
    """
    with open(path, "w", newline="") as fcsv:
        for line in summary_lines:
            if not line.startswith("#"):
                line = "# " + line
            fcsv.write(line.rstrip() + "\n")

        w = csv.writer(fcsv)
        w.writerow(["t", "dt", "accepted", "err_est", "y_norm_ratio", "step_ratio"])
        for r in rows:
            w.writerow([
                f"{r.t:.16e}",
                f"{r.dt:.16e}",
                int(r.accepted),
                f"{r.err_est:.16e}",
                f"{r.y_norm_ratio:.16e}",
                f"{r.step_ratio:.16e}",
            ])


# ----------------------------
# Adaptive drivers for splitting methods
# ----------------------------

def step_with_two_dts(
    method: str,
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    t: float,
    y: np.ndarray,
    dt: float,
    ordering: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    One dt step and two half steps for local error estimation.
    """
    if method == "Strang3":
        y_big = strang_step(A, B, C, y, dt, ordering)
        y_half = strang_step(A, B, C, y, 0.5 * dt, ordering)
        y_half2 = strang_step(A, B, C, y_half, 0.5 * dt, ordering)
        return y_big, y_half2

    if method == "Yoshida3":
        y_big = yoshida_4th_from_strang(A, B, C, y, dt, ordering)
        y_half = yoshida_4th_from_strang(A, B, C, y, 0.5 * dt, ordering)
        y_half2 = yoshida_4th_from_strang(A, B, C, y_half, 0.5 * dt, ordering)
        return y_big, y_half2

    raise ValueError(f"Unknown method={method}")


def run_adaptive_splitting(
    method: str,
    Bmult: float,
    y0: np.ndarray,
    T: float,
    dt0: float,
    tol: float,
    ordering: str,
    outdir: str,
):
    A, B, C = make_abc_matrices(Bmult)
    y_exact_T = exact_solution(A, B, C, y0, T)
    exact_norm = max(safe_norm(y_exact_T), 1e-300)
    y0_norm = max(safe_norm(y0), 1e-300)

    p = 2 if method == "Strang3" else 4

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

        y_big, y_ref = step_with_two_dts(method, A, B, C, t, y, dt, ordering)
        err_est = safe_norm(y_ref - y_big)
        step_ratio = safe_norm(y_big) / max(safe_norm(y), 1e-300)

        if err_est <= tol:
            t += dt
            y = y_ref
            accept += 1

            y_ratio = safe_norm(y) / y0_norm
            max_y_ratio = max(max_y_ratio, y_ratio)
            max_step_ratio = max(max_step_ratio, step_ratio)

            rows.append(
                AdaptiveLogRow(
                    t=t,
                    dt=dt,
                    accepted=1,
                    err_est=err_est,
                    y_norm_ratio=y_ratio,
                    step_ratio=step_ratio,
                )
            )
            dt = adapt_dt(dt, err_est, tol, p=p)
        else:
            reject += 1
            rows.append(
                AdaptiveLogRow(
                    t=t,
                    dt=dt,
                    accepted=0,
                    err_est=err_est,
                    y_norm_ratio=safe_norm(y) / y0_norm,
                    step_ratio=step_ratio,
                )
            )
            dt = adapt_dt(dt, err_est, tol, p=p)
            if reject > 20000:
                break

    final_err = safe_norm(y - y_exact_T) / exact_norm
    final_y_ratio = safe_norm(y) / y0_norm

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"adaptive_{method}_Bmult{int(Bmult)}_{ordering}.csv")

    summary = [
        f"{method}  Bmult={Bmult}  T={T}  dt0={dt0}  tol={tol}  ordering={ordering}",
        f"accept={accept}  reject={reject}",
        f"final_err={final_err:.16e}",
        f"max_y_norm_ratio={max_y_ratio:.16e}  final_y_norm_ratio={final_y_ratio:.16e}",
        f"max_step_ratio={max_step_ratio:.16e}",
    ]
    write_adaptive_csv(out_path, rows, summary)

    return out_path, accept, reject, final_err, max_y_ratio, final_y_ratio, max_step_ratio


# ----------------------------
# Embedded RK baseline
# ----------------------------

def get_builtin_heun_euler_pair():
    """
    Fallback explicit embedded pair.
    """
    A = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
    ], dtype=float)
    b_high = np.array([0.5, 0.5], dtype=float)
    b_low = np.array([1.0, 0.0], dtype=float)
    c = np.array([0.0, 1.0], dtype=float)
    order_high = 2
    return A, b_high, b_low, c, order_high


def get_embedded_pair_from_pythos(pair_name: str):
    """
    Get an embedded pair from butcher_tableau.py.

    We also import additive_rk.py first to verify the repo wiring is correct.
    """
    import additive_rk
    from butcher_tableau import embedded_pairs

    print("[ARK] Imported additive_rk from:", getattr(additive_rk, "__file__", "unknown"))

    if pair_name not in embedded_pairs:
        available = sorted(list(embedded_pairs.keys()))
        raise ValueError(f"Embedded pair '{pair_name}' not found. Available: {available}")

    pair = embedded_pairs[pair_name]

    A_rk = np.array(pair._a, dtype=float)
    b_hi = np.array(pair._b, dtype=float)
    b_lo = np.array(pair.b_aux, dtype=float)
    c_rk = np.array(pair._c, dtype=float)
    order_high = int(pair.order)

    print(f"[ARK] Using embedded pair from butcher_tableau: {pair_name}")
    return A_rk, b_hi, b_lo, c_rk, order_high


def run_ark_baseline(
    Bmult: float,
    y0: np.ndarray,
    T: float,
    dt0: float,
    tol: float,
    outdir: str,
    pair_name: str = "Dormand-Prince",
) -> Optional[str]:
    """
    Embedded RK baseline on the full RHS:
        y' = (A + B + C) y

    Uses pythOS butcher_tableau embedded pairs if available.
    Falls back to builtin Heun-Euler only if something fails.
    """
    A, B, C = make_abc_matrices(Bmult)

    try:
        A_rk, b_hi, b_lo, c_rk, order_high = get_embedded_pair_from_pythos(pair_name)
    except Exception as e:
        print(f"[ARK] Could not use pythOS embedded pair: {e}")
        print("[ARK] Falling back to Heun-Euler (builtin).")
        A_rk, b_hi, b_lo, c_rk, order_high = get_builtin_heun_euler_pair()

    s = len(c_rk)

    def f(t: float, y: np.ndarray) -> np.ndarray:
        return (A + B + C) @ y

    def rk_step_embedded(t: float, y: np.ndarray, h: float) -> Tuple[np.ndarray, float, float]:
        K = np.zeros((s, len(y)), dtype=float)

        for i in range(s):
            yi = y.copy()
            for j in range(i):
                aij = A_rk[i, j]
                if aij != 0.0:
                    yi = yi + h * aij * K[j]
            K[i] = f(t + c_rk[i] * h, yi)

        y_hi = y + h * np.sum(b_hi[:, None] * K, axis=0)
        y_lo = y + h * np.sum(b_lo[:, None] * K, axis=0)
        err = safe_norm(y_hi - y_lo)
        step_ratio = safe_norm(y_hi) / max(safe_norm(y), 1e-300)
        return y_hi, err, float(step_ratio)

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"ark_Bmult{int(Bmult)}.csv")

    y_exact_T = exact_solution(A, B, C, y0, T)
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

            rows.append(
                AdaptiveLogRow(
                    t=t,
                    dt=dt,
                    accepted=1,
                    err_est=err_est,
                    y_norm_ratio=y_ratio,
                    step_ratio=step_ratio,
                )
            )
            dt = adapt_dt(dt, err_est, tol, p=order_high)
        else:
            reject += 1
            rows.append(
                AdaptiveLogRow(
                    t=t,
                    dt=dt,
                    accepted=0,
                    err_est=err_est,
                    y_norm_ratio=safe_norm(y) / y0_norm,
                    step_ratio=step_ratio,
                )
            )
            dt = adapt_dt(dt, err_est, tol, p=order_high)
            if reject > 20000:
                break

    final_err = safe_norm(y - y_exact_T) / exact_norm
    final_y_ratio = safe_norm(y) / y0_norm

    summary = [
        f"ARK embedded baseline  Bmult={Bmult}  T={T}  dt0={dt0}  tol={tol}",
        f"pair_name={pair_name}",
        f"accept={accept}  reject={reject}",
        f"final_err={final_err:.16e}",
        f"max_y_norm_ratio={max_y_ratio:.16e}  final_y_norm_ratio={final_y_ratio:.16e}",
        f"max_step_ratio={max_step_ratio:.16e}",
    ]
    write_adaptive_csv(out_path, rows, summary)

    print(f"[ARK] saved: {out_path}")
    return out_path


# ----------------------------
# Plotting
# ----------------------------

def load_dt_series(csv_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ts: List[float] = []
    dts: List[float] = []
    acc: List[int] = []

    with open(csv_path, "r") as f:
        reader = csv.reader(line for line in f if not line.lstrip().startswith("#"))
        header = next(reader, None)
        if header is None:
            return np.array([]), np.array([]), np.array([])

        for row in reader:
            if len(row) < 3:
                continue
            try:
                t = float(row[0])
                dt = float(row[1])
                accepted = int(row[2])
            except Exception:
                continue
            ts.append(t)
            dts.append(dt)
            acc.append(accepted)

    return np.array(ts), np.array(dts), np.array(acc)


def make_dt_plot(csv_path: str, title: str, out_png: str) -> None:
    import matplotlib.pyplot as plt

    t, dt, accepted = load_dt_series(csv_path)
    if len(t) == 0:
        print(f"[plots] no data in {csv_path}")
        return

    mask = accepted == 1
    tA = t[mask]
    dtA = dt[mask]

    plt.figure()
    plt.plot(tA, dtA, marker="o", linewidth=1.5)
    plt.yscale("log")
    plt.xlabel("t")
    plt.ylabel("dt")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


# ----------------------------
# CLI
# ----------------------------

def parse_list_floats(s: str) -> List[float]:
    out: List[float] = []
    for part in s.split(","):
        part = part.strip()
        if part:
            out.append(float(part))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--Bmults", type=str, default="1,20,50,100", help="comma-separated list")
    ap.add_argument("--T", type=float, default=0.5)
    ap.add_argument("--dt0", type=float, default=1e-2)
    ap.add_argument("--tol", type=float, default=1e-7)
    ap.add_argument("--ordering", type=str, default="CAB")
    ap.add_argument("--outdir", type=str, default="experiments/abc_adaptive_outputs")
    ap.add_argument("--do_ark", action="store_true")
    ap.add_argument("--ark_pair", type=str, default="Dormand-Prince")
    ap.add_argument("--plots", action="store_true")
    args = ap.parse_args()

    Bmults = parse_list_floats(args.Bmults)

    print("=" * 72)
    print(
        f"Adaptive study: Bmults={Bmults}  T={args.T}  dt0={args.dt0}  "
        f"tol={args.tol}  ordering={args.ordering}"
    )
    print("=" * 72)

    os.makedirs(args.outdir, exist_ok=True)
    plot_dir = os.path.join(args.outdir, "plots")
    if args.plots:
        os.makedirs(plot_dir, exist_ok=True)

    y0 = np.array([1.0, 0.2, -0.1], dtype=float)

    for Bmult in Bmults:
        print("\n" + "=" * 72)
        print(
            f"Adaptive study: Bmult={int(Bmult)}  T={args.T}  dt0={args.dt0}  "
            f"tol={args.tol}  ordering={args.ordering}"
        )
        print("=" * 72)

        path_strang, aS, rS, eS, maxyS, finalyS, maxstepS = run_adaptive_splitting(
            method="Strang3",
            Bmult=Bmult,
            y0=y0,
            T=args.T,
            dt0=args.dt0,
            tol=args.tol,
            ordering=args.ordering,
            outdir=args.outdir,
        )
        print(
            f"[Strang3] accept={aS} reject={rS}  final_err={eS:.6e}  "
            f"max||y||/||y0||={maxyS:.6e}  final||y||/||y0||={finalyS:.6e}  "
            f"max_step={maxstepS:.6e}"
        )
        print(f"  saved: {path_strang}")

        path_yosh, aY, rY, eY, maxyY, finalyY, maxstepY = run_adaptive_splitting(
            method="Yoshida3",
            Bmult=Bmult,
            y0=y0,
            T=args.T,
            dt0=args.dt0,
            tol=args.tol,
            ordering=args.ordering,
            outdir=args.outdir,
        )
        print(
            f"[Yoshida3] accept={aY} reject={rY}  final_err={eY:.6e}  "
            f"max||y||/||y0||={maxyY:.6e}  final||y||/||y0||={finalyY:.6e}  "
            f"max_step={maxstepY:.6e}"
        )
        print(f"  saved: {path_yosh}")

        ark_path = None
        if args.do_ark:
            ark_path = run_ark_baseline(
                Bmult=Bmult,
                y0=y0,
                T=args.T,
                dt0=args.dt0,
                tol=args.tol,
                outdir=args.outdir,
                pair_name=args.ark_pair,
            )

        if args.plots:
            make_dt_plot(
                path_strang,
                f"Adaptive dt vs t | Strang-3 | Bmult={int(Bmult)} | {args.ordering}",
                os.path.join(plot_dir, f"dt_vs_t_Strang3_Bmult{int(Bmult)}_{args.ordering}.png"),
            )
            make_dt_plot(
                path_yosh,
                f"Adaptive dt vs t | Yoshida-3 | Bmult={int(Bmult)} | {args.ordering}",
                os.path.join(plot_dir, f"dt_vs_t_Yoshida3_Bmult{int(Bmult)}_{args.ordering}.png"),
            )
            if ark_path is not None:
                make_dt_plot(
                    ark_path,
                    f"Adaptive dt vs t | ARK ({args.ark_pair}) | Bmult={int(Bmult)}",
                    os.path.join(plot_dir, f"dt_vs_t_ARK_Bmult{int(Bmult)}.png"),
                )

    if args.plots:
        print(f"Plots saved to: {plot_dir}")


if __name__ == "__main__":
    main()