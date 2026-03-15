#!/usr/bin/env python3
"""
ABC splitting ordering study.

What it does
------------
- Defines a linear 2x2 system y' = (A + B + C) y with tunable stiffness via Bmult.
- Compares orderings of 3-operator splittings for:
    * Strang-3 (2nd order; symmetric)
    * Yoshida-3 (4th order composition built from Strang; has negative substeps)
- For each (method, ordering, dt) computes error vs exact expm solution at final time.
- Reports observed order p from multiple dt values.
- Tracks stability/step-growth metrics:
    * max||y||/||y0||
    * final||y||/||y0||
    * max_step_ratio = max_n ||y_{n+1}||/||y_n||
- Saves per-regime CSV and generates plots:
    * error vs ordering
    * max_step_ratio vs ordering

Run examples
------------
python3 experiments/abc_splitting_study.py
python3 experiments/abc_splitting_study.py --Bmults 1,20,50,100
python3 experiments/abc_splitting_study.py --dts 1e-2,5e-3,2.5e-3 --T 0.5
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# --- Make repo-root imports robust when running from repo root ---
THIS = Path(__file__).resolve()
REPO = THIS.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

# Optional SciPy expm (preferred). Fallback provided if SciPy isn't installed.
try:
    from scipy.linalg import expm  # type: ignore
except Exception:
    expm = None


def expm_fallback(M: np.ndarray) -> np.ndarray:
    """
    Small-matrix expm fallback using eigendecomposition.
    Works well for 2x2 (our use case).
    """
    w, V = np.linalg.eig(M)
    Vinv = np.linalg.inv(V)
    return (V @ np.diag(np.exp(w)) @ Vinv).astype(complex)


def mat_expm(M: np.ndarray) -> np.ndarray:
    if expm is not None:
        return expm(M)
    E = expm_fallback(M)
    # If imaginary noise is tiny, cast back to real.
    if np.max(np.abs(E.imag)) < 1e-12:
        return E.real
    return E


def norm2(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def safe_ratio(a: float, b: float) -> float:
    if b == 0.0:
        return float("inf") if a != 0.0 else 1.0
    return a / b


@dataclass
class RunMetrics:
    err: float
    max_norm_ratio: float
    final_norm_ratio: float
    max_step_ratio: float
    stable: bool


# -------------------------------
# Problem definition (A, B, C)
# -------------------------------
def make_abc_matrices(Bmult: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    A, B, C are chosen so ordering matters (non-commuting operators),
    and stiffness is introduced by scaling B.

    You can tweak these matrices if you want a stronger/weaker ordering effect.
    """
    A = np.array([[0.0, 1.0],
                  [-1.0, 0.0]], dtype=float)            # rotation-like
    B = np.array([[-1.0, 0.0],
                  [0.0, -0.2]], dtype=float) * Bmult    # damping/stiffness
    C = np.array([[0.0, 0.6],
                  [0.0, 0.0]], dtype=float)             # shear-like
    return A, B, C


def exact_solution(y0: np.ndarray, T: float, A: np.ndarray, B: np.ndarray, C: np.ndarray) -> np.ndarray:
    L = A + B + C
    return mat_expm(T * L) @ y0


# -------------------------------
# Splitting building blocks
# -------------------------------
def flow_linear(L: np.ndarray, h: float, y: np.ndarray) -> np.ndarray:
    return mat_expm(h * L) @ y


def strang_step(ordering: str, h: float, y: np.ndarray, ops: Dict[str, np.ndarray]) -> np.ndarray:
    """
    Strang splitting for 3 operators using an explicit ordering like 'ABC'.

    For 'ABC', we use:
        exp(h/2 A) exp(h/2 B) exp(h C) exp(h/2 B) exp(h/2 A)
    """
    a, b, c = ordering[0], ordering[1], ordering[2]
    y1 = flow_linear(ops[a], 0.5 * h, y)
    y2 = flow_linear(ops[b], 0.5 * h, y1)
    y3 = flow_linear(ops[c], 1.0 * h, y2)
    y4 = flow_linear(ops[b], 0.5 * h, y3)
    y5 = flow_linear(ops[a], 0.5 * h, y4)
    return y5


def yoshida4_step(ordering: str, h: float, y: np.ndarray, ops: Dict[str, np.ndarray]) -> np.ndarray:
    """
    4th-order Yoshida composition built from Strang:
        S4(h) = S2(w1 h) S2(w0 h) S2(w1 h)
    where
        w1 = 1/(2 - 2^(1/3)),  w0 = -2^(1/3)/(2 - 2^(1/3))
    Note w0 is negative => potential stability issues in stiff regimes.
    """
    cbrt2 = 2.0 ** (1.0 / 3.0)
    w1 = 1.0 / (2.0 - cbrt2)
    w0 = -cbrt2 / (2.0 - cbrt2)

    y1 = strang_step(ordering, w1 * h, y, ops)
    y2 = strang_step(ordering, w0 * h, y1, ops)
    y3 = strang_step(ordering, w1 * h, y2, ops)
    return y3


def simulate(method_name: str, ordering: str, dt: float, T: float, y0: np.ndarray, ops: Dict[str, np.ndarray]) -> RunMetrics:
    steps = int(round(T / dt))
    if steps <= 0:
        raise ValueError("dt too large for given T")
    # Make dt so N*dt exactly hits T (avoid drift)
    dt = T / steps

    y = y0.copy()
    n0 = norm2(y0)
    max_norm_ratio = 1.0
    max_step_ratio = 1.0
    stable = True

    stepper = strang_step if method_name == "Strang-3" else yoshida4_step

    for _ in range(steps):
        y_new = stepper(ordering, dt, y, ops)

        ny = norm2(y)
        ny_new = norm2(y_new)
        if not np.all(np.isfinite(y_new)):
            stable = False
            break

        # step ratio ||y_{n+1}|| / ||y_n||
        if ny > 0:
            max_step_ratio = max(max_step_ratio, ny_new / ny)
        else:
            max_step_ratio = max(max_step_ratio, float("inf"))

        # global ratio ||y||/||y0||
        if n0 > 0:
            max_norm_ratio = max(max_norm_ratio, ny_new / n0)

        # crude stability guard: if it explodes massively, declare unstable
        if ny_new > 1e12 * max(1.0, n0):
            stable = False
            break

        y = y_new

    yT = y
    # If unstable, set err=inf
    return RunMetrics(
        err=float("inf"),
        max_norm_ratio=max_norm_ratio,
        final_norm_ratio=safe_ratio(norm2(yT), n0),
        max_step_ratio=max_step_ratio,
        stable=stable,
    )


def compute_error(y_approx: np.ndarray, y_exact: np.ndarray) -> float:
    return norm2(y_approx - y_exact)


def run_one(method_name: str, ordering: str, dt: float, T: float, y0: np.ndarray, ops: Dict[str, np.ndarray]) -> Tuple[np.ndarray, RunMetrics]:
    steps = int(round(T / dt))
    steps = max(1, steps)
    dt = T / steps

    y = y0.copy()
    n0 = norm2(y0)
    max_norm_ratio = 1.0
    max_step_ratio = 1.0
    stable = True

    stepper = strang_step if method_name == "Strang-3" else yoshida4_step

    for _ in range(steps):
        y_new = stepper(ordering, dt, y, ops)

        if not np.all(np.isfinite(y_new)):
            stable = False
            break

        ny = norm2(y)
        ny_new = norm2(y_new)
        if ny > 0:
            max_step_ratio = max(max_step_ratio, ny_new / ny)
        else:
            max_step_ratio = max(max_step_ratio, float("inf"))

        if n0 > 0:
            max_norm_ratio = max(max_norm_ratio, ny_new / n0)

        if ny_new > 1e12 * max(1.0, n0):
            stable = False
            break

        y = y_new

    metrics = RunMetrics(
        err=float("inf"),
        max_norm_ratio=max_norm_ratio,
        final_norm_ratio=safe_ratio(norm2(y), n0),
        max_step_ratio=max_step_ratio,
        stable=stable,
    )
    return y, metrics


def estimate_order(errs: List[float], dts: List[float]) -> float:
    """
    Estimate observed order from (dt1,dt2,dt3...) and errs at final time.
    Uses least-squares fit on log(err)=p log(dt)+const, ignoring non-finite.
    """
    xs, ys = [], []
    for e, dt in zip(errs, dts):
        if np.isfinite(e) and e > 0 and dt > 0:
            xs.append(math.log(dt))
            ys.append(math.log(e))
    if len(xs) < 2:
        return float("nan")
    x = np.array(xs)
    y = np.array(ys)
    A = np.vstack([x, np.ones_like(x)]).T
    p, _ = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(p)


def all_orderings() -> List[str]:
    return ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]


# -------------------------------
# Plotting
# -------------------------------
def make_plots_for_csv(csv_path: Path, outdir: Path) -> None:
    import matplotlib.pyplot as plt

    rows = []
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)

    if not rows:
        print(f"[plot] No rows in {csv_path}")
        return

    # infer Bmult and method from filename
    stem = csv_path.stem  # abc_study_Bmult20
    Bmult_label = stem.replace("abc_study_", "")

    # group by method
    by_method: Dict[str, List[dict]] = {}
    for r in rows:
        by_method.setdefault(r["method"], []).append(r)

    for method, mr in by_method.items():
        # sort by fixed ordering list
        ords = all_orderings()
        mr_map = {r["ordering"]: r for r in mr}
        errs = [float(mr_map[o]["err_dtmin"]) if o in mr_map else float("nan") for o in ords]
        msteps = [float(mr_map[o]["max_step_ratio"]) if o in mr_map else float("nan") for o in ords]
        stable = [mr_map[o]["status"] if o in mr_map else "?" for o in ords]

        # Error bar plot
        plt.figure()
        x = np.arange(len(ords))
        plt.bar(x, errs)
        plt.xticks(x, ords)
        plt.yscale("log")
        plt.xlabel("Ordering")
        plt.ylabel("Error at dt_min (log scale)")
        plt.title(f"Error vs ordering — {Bmult_label} — {method}")
        for i, st in enumerate(stable):
            if st != "ok":
                plt.text(i, errs[i], st, rotation=90, va="bottom", ha="center")
        pth = outdir / f"err_vs_order_{Bmult_label}_{method}.png"
        plt.tight_layout()
        plt.savefig(pth, dpi=200)
        plt.close()

        # max_step_ratio plot
        plt.figure()
        plt.bar(x, msteps)
        plt.xticks(x, ords)
        plt.yscale("log")
        plt.xlabel("Ordering")
        plt.ylabel("max_step_ratio = max ||y_{n+1}||/||y_n|| (log scale)")
        plt.title(f"max_step_ratio vs ordering — {Bmult_label} — {method}")
        pth = outdir / f"maxstep_vs_order_{Bmult_label}_{method}.png"
        plt.tight_layout()
        plt.savefig(pth, dpi=200)
        plt.close()


# -------------------------------
# Main experiment
# -------------------------------
def run_regime(Bmult: float, T: float, dts: List[float], outdir: Path) -> Path:
    A, B, C = make_abc_matrices(Bmult)
    ops = {"A": A, "B": B, "C": C}

    y0 = np.array([1.0, 0.25], dtype=float)
    y_exact = exact_solution(y0, T, A, B, C)
    exact_norm = norm2(y_exact)

    print("\n" + "=" * 72)
    print(f"Regime: Bmult={Bmult:g}   T={T:g}   dts={dts}")
    print(f"Exact ||y(T)|| = {exact_norm:.12g}   ||y0|| = {norm2(y0):.12g}")
    print("=" * 72 + "\n")

    methods = ["Strang-3", "Yoshida-3"]
    ords = all_orderings()

    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / f"abc_study_Bmult{int(Bmult) if float(Bmult).is_integer() else Bmult}.csv"

    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Bmult", "T", "method", "ordering",
            "err_dtmin", "p_est",
            "max_norm_ratio", "final_norm_ratio",
            "max_step_ratio",
            "status"
        ])

        for method in methods:
            print(f"Method: {method}")
            print("-" * 72)

            # compute per ordering errors at each dt
            results: Dict[str, Dict[str, float]] = {}

            for ordering in ords:
                errs = []
                max_norms = []
                final_norms = []
                max_steps = []
                statuses = []

                for dt in dts:
                    y_num, m = run_one(method, ordering, dt, T, y0, ops)
                    if not m.stable:
                        errs.append(float("inf"))
                        statuses.append("UNSTABLE")
                    else:
                        errs.append(compute_error(y_num, y_exact))
                        statuses.append("ok")
                    max_norms.append(m.max_norm_ratio)
                    final_norms.append(m.final_norm_ratio)
                    max_steps.append(m.max_step_ratio)

                # dt_min is smallest dt
                dt_min_idx = int(np.argmin(np.array(dts)))
                err_dtmin = float(errs[dt_min_idx])

                # order estimate from errs vs dt
                p_est = estimate_order(errs, dts)

                # use dt_min metrics (consistent with err_dtmin)
                max_norm_ratio = float(max_norms[dt_min_idx])
                final_norm_ratio = float(final_norms[dt_min_idx])
                max_step_ratio = float(max_steps[dt_min_idx])

                status = "ok"
                if not np.isfinite(err_dtmin):
                    status = "UNSTABLE"

                results[ordering] = {
                    "err_dtmin": err_dtmin,
                    "p_est": p_est,
                    "max_norm_ratio": max_norm_ratio,
                    "final_norm_ratio": final_norm_ratio,
                    "max_step_ratio": max_step_ratio,
                    "status": status
                }

                print(
                    f"{ordering:3s}  "
                    f"err(dt_min)={err_dtmin:.6e}  "
                    f"p≈{p_est:.2f}  "
                    f"max||y||/||y0||={max_norm_ratio:.6e}  "
                    f"final||y||/||y0||={final_norm_ratio:.6e}  "
                    f"max_step={max_step_ratio:.6e}  "
                    f"{status}"
                )

                writer.writerow([
                    f"{Bmult:g}", f"{T:g}", method, ordering,
                    f"{err_dtmin:.16e}", f"{p_est:.8g}",
                    f"{max_norm_ratio:.16e}", f"{final_norm_ratio:.16e}",
                    f"{max_step_ratio:.16e}",
                    status
                ])

            # best/worst based on err_dtmin among stable
            stable_errs = [(o, results[o]["err_dtmin"]) for o in ords if np.isfinite(results[o]["err_dtmin"])]
            if stable_errs:
                best_o, best_e = min(stable_errs, key=lambda x: x[1])
                worst_o, worst_e = max(stable_errs, key=lambda x: x[1])
                spread = worst_e / best_e if best_e > 0 else float("inf")
                print("\nBest ordering:", best_o, "err=", f"{best_e:.3e}")
                print("Worst ordering:", worst_o, "err=", f"{worst_e:.3e}")
                print(f"Spread (worst/best): {spread:.2f}x\n")
            else:
                print("\nAll orderings unstable for this method/regime.\n")

    print(f"CSV saved: {csv_path}")
    return csv_path


def parse_csv_floats(s: str) -> List[float]:
    out = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        out.append(float(part))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=0.5)
    ap.add_argument("--dts", type=str, default="1e-2,5e-3,2.5e-3")
    ap.add_argument("--Bmults", type=str, default="1,20")
    ap.add_argument("--outdir", type=str, default="experiments/abc_outputs")
    ap.add_argument("--plots", action="store_true", help="Generate plots after running.")
    ap.add_argument("--plot_only", action="store_true", help="Only generate plots from existing CSVs in outdir.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    dts = parse_csv_floats(args.dts)
    Bmults = parse_csv_floats(args.Bmults)

    if args.plot_only:
        if not outdir.exists():
            print(f"[plot_only] outdir not found: {outdir}")
            return
        csvs = sorted(outdir.glob("abc_study_Bmult*.csv"))
        if not csvs:
            print(f"[plot_only] no CSVs found in {outdir}")
            return
        for c in csvs:
            make_plots_for_csv(c, outdir)
        print(f"[plot_only] Plots saved to: {outdir}")
        return

    csv_paths = []
    for Bmult in Bmults:
        csv_paths.append(run_regime(Bmult, args.T, dts, outdir))

    if args.plots:
        for c in csv_paths:
            make_plots_for_csv(c, outdir)
        print(f"Plots saved to: {outdir}")


if __name__ == "__main__":
    main()