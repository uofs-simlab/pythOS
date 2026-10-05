#!/usr/bin/env python3
"""
Make plots from experiments/abc_outputs/abc_study_Bmult*.csv

Outputs (per Bmult):
  - err_vs_order_<Bmult>_<method>.png        (log-scale error by ordering)
  - maxstep_vs_order_<Bmult>_<method>.png    (max_step by ordering)

No pandas required (uses csv + matplotlib).
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Helpers to be robust to header names
# -----------------------------
def _first_present(row: Dict[str, str], keys: List[str]) -> Optional[str]:
    for k in keys:
        if k in row and row[k] != "":
            return row[k]
    return None


def _as_float(x: Optional[str]) -> float:
    if x is None:
        return float("nan")
    x = x.strip()
    if x.lower() == "inf":
        return float("inf")
    if x.lower() == "nan":
        return float("nan")
    return float(x)


def _as_bool(x: Optional[str]) -> bool:
    if x is None:
        return False
    s = x.strip().lower()
    return s in ("1", "true", "t", "yes", "y", "ok", "stable")


@dataclass
class Row:
    bmult: float
    method: str
    ordering: str
    err_dt_min: float
    p_obs: float
    max_norm_ratio: float
    final_norm_ratio: float
    max_step: float
    stable: bool


def read_csv_rows(path: str) -> List[Row]:
    rows: List[Row] = []
    bmult_guess = _infer_bmult_from_filename(path)

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        # print("CSV headers:", headers)

        for r in reader:
            method = _first_present(r, ["method", "Method"])
            ordering = _first_present(r, ["ordering", "Ordering", "perm", "Permutation"])
            if method is None or ordering is None:
                # Skip malformed lines
                continue

            # Error at smallest dt (dt_min)
            err_dt_min = _as_float(
                _first_present(r, ["err_dt_min", "err(dt_min)", "err_dtmin", "err_min", "err"])
            )
            # Observed order
            p_obs = _as_float(_first_present(r, ["p", "p_obs", "observed_p", "p≈"]))

            # Norm ratios
            max_norm_ratio = _as_float(
                _first_present(r, ["max_norm_ratio", "max||y||/||y0||", "max_ratio"])
            )
            final_norm_ratio = _as_float(
                _first_present(r, ["final_norm_ratio", "final||y||/||y0||", "final_ratio"])
            )

            # New metric
            max_step = _as_float(_first_present(r, ["max_step", "max_step_ratio"]))

            # Stability flag
            stable_str = _first_present(r, ["stable", "ok", "status"])
            if stable_str is None:
                # Sometimes there's an "UNSTABLE" status column
                status = _first_present(r, ["Status", "status", "label"])
                if status is not None and status.strip().upper() == "UNSTABLE":
                    stable = False
                else:
                    # fall back: unstable if err is inf or max_norm_ratio huge
                    stable = np.isfinite(err_dt_min) and (not (np.isfinite(max_norm_ratio) and max_norm_ratio > 1e3))
            else:
                # "ok" is stable; "UNSTABLE" is not
                s = stable_str.strip().lower()
                if "unstable" in s:
                    stable = False
                elif s in ("ok", "stable", "true", "1", "yes", "y"):
                    stable = True
                else:
                    stable = _as_bool(stable_str)

            rows.append(
                Row(
                    bmult=bmult_guess,
                    method=method.strip(),
                    ordering=ordering.strip(),
                    err_dt_min=err_dt_min,
                    p_obs=p_obs,
                    max_norm_ratio=max_norm_ratio,
                    final_norm_ratio=final_norm_ratio,
                    max_step=max_step,
                    stable=stable,
                )
            )

    return rows


def _infer_bmult_from_filename(path: str) -> float:
    base = os.path.basename(path)
    # expected like abc_study_Bmult50.csv
    for token in base.replace(".csv", "").split("_"):
        if token.lower().startswith("bmult"):
            try:
                return float(token[len("Bmult"):])
            except Exception:
                pass
    return float("nan")


def ensure_dir(d: str) -> None:
    os.makedirs(d, exist_ok=True)


def sanitize(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in s)


# -----------------------------
# Plotting
# -----------------------------
def plot_err_vs_order(rows: List[Row], outdir: str, tag: str) -> None:
    """
    One plot per method: x=ordering, y=err_dt_min (log scale)
    Stable: circle marker, Unstable: x marker
    """
    ensure_dir(outdir)
    methods = sorted(set(r.method for r in rows))
    for method in methods:
        sub = [r for r in rows if r.method == method]

        # sort by err (stable first), then ordering name
        sub.sort(key=lambda r: (not r.stable, r.ordering))

        xlabels = [r.ordering for r in sub]
        x = np.arange(len(sub))
        y = np.array([r.err_dt_min for r in sub], dtype=float)

        # Replace non-finite with a large sentinel for plotting
        finite = np.isfinite(y) & (y > 0)
        if np.any(finite):
            y_max = float(np.max(y[finite]))
            y_plot = y.copy()
            y_plot[~finite] = y_max * 10.0
        else:
            y_plot = np.ones_like(y)

        fig = plt.figure()
        ax = plt.gca()
        ax.set_yscale("log")

        stable_idx = [i for i, r in enumerate(sub) if r.stable and np.isfinite(r.err_dt_min) and r.err_dt_min > 0]
        unstable_idx = [i for i, r in enumerate(sub) if not r.stable or not (np.isfinite(r.err_dt_min) and r.err_dt_min > 0)]

        if stable_idx:
            ax.plot(x[stable_idx], y_plot[stable_idx], linestyle="None", marker="o", label="stable")
        if unstable_idx:
            ax.plot(x[unstable_idx], y_plot[unstable_idx], linestyle="None", marker="x", label="unstable / inf")

        # annotate unstable points with max_step (if available)
        for i in unstable_idx:
            ms = sub[i].max_step
            if np.isfinite(ms):
                ax.annotate(f"{ms:.2g}", (x[i], y_plot[i]), textcoords="offset points", xytext=(0, 8), ha="center")

        ax.set_xticks(x)
        ax.set_xticklabels(xlabels, rotation=0)
        ax.set_xlabel("Ordering")
        ax.set_ylabel("abs_err at smallest dt (log scale)")
        ax.set_title(f"{tag}  |  {method}  |  error vs ordering")
        ax.legend()

        fig.tight_layout()
        outpath = os.path.join(outdir, f"err_vs_order_{sanitize(tag)}_{sanitize(method)}.png")
        fig.savefig(outpath, dpi=200)
        plt.close(fig)


def plot_maxstep_vs_order(rows: List[Row], outdir: str, tag: str) -> None:
    """
    One plot per method: x=ordering, y=max_step (linear)
    """
    ensure_dir(outdir)
    methods = sorted(set(r.method for r in rows))
    for method in methods:
        sub = [r for r in rows if r.method == method]
        sub.sort(key=lambda r: (not r.stable, r.ordering))

        xlabels = [r.ordering for r in sub]
        x = np.arange(len(sub))
        y = np.array([r.max_step for r in sub], dtype=float)

        fig = plt.figure()
        ax = plt.gca()

        stable_idx = [i for i, r in enumerate(sub) if r.stable and np.isfinite(r.max_step)]
        unstable_idx = [i for i, r in enumerate(sub) if (not r.stable) and np.isfinite(r.max_step)]

        if stable_idx:
            ax.plot(x[stable_idx], y[stable_idx], linestyle="None", marker="o", label="stable")
        if unstable_idx:
            ax.plot(x[unstable_idx], y[unstable_idx], linestyle="None", marker="x", label="unstable")

        ax.axhline(1.0, linewidth=1)

        ax.set_xticks(x)
        ax.set_xticklabels(xlabels, rotation=0)
        ax.set_xlabel("Ordering")
        ax.set_ylabel("max_step = max ||y_{n+1}||/||y_n||")
        ax.set_title(f"{tag}  |  {method}  |  max_step vs ordering")
        ax.legend()

        fig.tight_layout()
        outpath = os.path.join(outdir, f"maxstep_vs_order_{sanitize(tag)}_{sanitize(method)}.png")
        fig.savefig(outpath, dpi=200)
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--csv_glob",
        default="experiments/abc_outputs/abc_study_Bmult*.csv",
        help="Glob for input CSVs",
    )
    ap.add_argument(
        "--outdir",
        default="experiments/abc_outputs",
        help="Where to save PNGs",
    )
    args = ap.parse_args()

    paths = sorted(glob.glob(args.csv_glob))
    if not paths:
        raise SystemExit(f"No CSVs matched: {args.csv_glob}")

    for p in paths:
        rows = read_csv_rows(p)
        if not rows:
            print(f"Skipping (no rows): {p}")
            continue
        bmult = rows[0].bmult
        tag = f"Bmult={bmult:g}"
        plot_err_vs_order(rows, args.outdir, tag)
        plot_maxstep_vs_order(rows, args.outdir, tag)
        print(f"Saved plots for {tag} -> {args.outdir}")


if __name__ == "__main__":
    main()