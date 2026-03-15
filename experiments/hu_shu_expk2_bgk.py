"""
experiments/hu_shu_expk2_bgk.py

Reproduce key Hu–Shu (arXiv:1807.03728) BGK experiments using pythOS as the driver.

What this script does:
1) Implements Hu–Shu ExpRK2 step (paper eq. (2.6)) for BGK.
2) Hooks it into pythOS fractional_step via methods={(1,):'ANALYTIC'} and alpha=[[1.0]].
3) Provides two runnable experiments:
   A) Mixed-regime BGK (eps = eps(x) from eq. (7.6)) at tf=0.5, Nx=40, Nv=150, vmax=15.
      -> Produces macroscopic (rho,u,T) for Hu–Shu and an SSP-RK2 reference.
   B) Accuracy sweep (Table-1-style): error vs refined run (dx/2, dt/2) at tf=0.1,
      for several Nx and eps values.
      -> Writes CSV.

Run examples (from pythOS root with PYTHONPATH set):
  python3 experiments/hu_shu_expk2_bgk.py --case mixed
  python3 experiments/hu_shu_expk2_bgk.py --case table
  python3 experiments/hu_shu_expk2_bgk.py --case timeorder

Notes:
- Spatial discretization: simple first-order upwind for transport (periodic in x).
  This may not match the paper’s high-order WENO exactly, but it gives a correct,
  reproducible baseline wiring of Hu–Shu into pythOS. You can swap transport_rhs_upwind
  with your preferred reconstruction later.
- Positivity limiter: not implemented here. This is Step-1 reproduction of the *method*
  and baseline results workflow. Add limiter after the baseline matches qualitatively.
"""

import argparse
import csv
import os
import sys
from dataclasses import dataclass

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from fractional_step import fractional_step


# ----------------------------
# Grids
# ----------------------------
@dataclass
class Grid1D:
    x: np.ndarray
    dx: float
    v: np.ndarray
    dv: float
    vmax: float


def make_grid(Nx: int, Nv: int, x0: float = 0.0, x1: float = 2.0, vmax: float = 15.0) -> Grid1D:
    x = np.linspace(x0, x1, Nx, endpoint=False)
    dx = (x1 - x0) / Nx
    v = np.linspace(-vmax, vmax, Nv)
    dv = (2.0 * vmax) / (Nv - 1)
    return Grid1D(x=x, dx=dx, v=v, dv=dv, vmax=vmax)


# ----------------------------
# BGK building blocks
# ----------------------------
def moments(f, v, dv, rho_floor=1e-12):
    """
    Robust moments for Maxwellian construction.
    - clips negative f when forming moments (optional but stabilizes Maxwellian)
    - guards against tiny rho causing u overflow
    """
    fpos = np.maximum(f, 0.0)  # stabilize moments if slight negativity appears
    rho = np.sum(fpos, axis=1) * dv

    u = np.zeros_like(rho)
    T = np.ones_like(rho)

    mask = rho > rho_floor
    if np.any(mask):
        rho_m = rho[mask]
        num_u = np.sum(fpos[mask, :] * v[None, :], axis=1) * dv
        u_m = num_u / rho_m

        # variance form: T = (1/rho) ∫ (v-u)^2 f dv
        diff = v[None, :] - u_m[:, None]
        var = (np.sum(fpos[mask, :] * diff * diff, axis=1) * dv) / rho_m
        T_m = np.maximum(var, 1e-14)

        u[mask] = u_m
        T[mask] = T_m

    rho = np.maximum(rho, rho_floor)
    return rho, u, T


def maxwellian(rho: np.ndarray, u: np.ndarray, T: np.ndarray, v: np.ndarray):
    """
    1D Maxwellian:
      M = rho / sqrt(2*pi*T) * exp(-(v-u)^2/(2T))
    """
    coef = rho[:, None] / np.sqrt(2.0 * np.pi * T[:, None])
    arg = -((v[None, :] - u[:, None]) ** 2) / (2.0 * T[:, None])
    return coef * np.exp(arg)


def exp_bgk_homogeneous(f: np.ndarray, dt: float, eps: float, v: np.ndarray, dv: float):
    """
    Exact homogeneous BGK relaxation with frozen moments at input f:
      f_t = (1/eps) (M[f] - f)
    Solution: f(dt) = e^{-dt/eps} f + (1-e^{-dt/eps}) M
    """
    rho, u, T = moments(f, v, dv)
    M = maxwellian(rho, u, T, v)
    lam = dt / eps
    a = np.exp(-lam)
    return a * f + (1.0 - a) * M


def exp_bgk(f: np.ndarray, dt: float, eps_x, grid: Grid1D):
    """
    BGK relaxation with eps possibly x-dependent.
    eps_x: scalar or array (Nx,)
    """
    if np.isscalar(eps_x):
        return exp_bgk_homogeneous(f, dt, float(eps_x), grid.v, grid.dv)

    out = np.empty_like(f)
    for i in range(f.shape[0]):
        out[i, :] = exp_bgk_homogeneous(f[i:i + 1, :], dt, float(eps_x[i]), grid.v, grid.dv)[0, :]
    return out


# ----------------------------
# Transport: upwind FD, periodic in x
# ----------------------------
def transport_rhs_upwind(f: np.ndarray, grid: Grid1D):
    """
    T(f) = - v * d_x f
    Upwind in x for each velocity v.
    """
    fL = np.roll(f, 1, axis=0)
    fR = np.roll(f, -1, axis=0)
    rhs = np.empty_like(f)

    for j, vj in enumerate(grid.v):
        if vj >= 0.0:
            rhs[:, j] = -vj * (f[:, j] - fL[:, j]) / grid.dx
        else:
            rhs[:, j] = -vj * (fR[:, j] - f[:, j]) / grid.dx
    return rhs


def minmod(a, b):
    """
    minmod limiter for arrays a,b (same shape).
    """
    out = np.zeros_like(a)
    s = np.sign(a) + np.sign(b)
    mask = (s == 2) | (s == -2)  # same sign and nonzero
    out[mask] = np.sign(a[mask]) * np.minimum(np.abs(a[mask]), np.abs(b[mask]))
    return out


def minmod3(a, b, c):
    out = np.zeros_like(a)
    s = np.sign(a) + np.sign(b) + np.sign(c)
    mask = (np.abs(s) == 3)  # all same sign
    out[mask] = np.sign(a[mask]) * np.minimum(np.minimum(np.abs(a[mask]), np.abs(b[mask])), np.abs(c[mask]))
    return out


def slope_mc(df_back, df_forw):
    return minmod3(2.0*df_back, 0.5*(df_back + df_forw), 2.0*df_forw)


def transport_rhs_muscl2(f: np.ndarray, grid: Grid1D):
    """
    2nd-order MUSCL finite-volume transport for:
        f_t + d_x (v f) = 0
    periodic in x.

    For each velocity v_j:
      - reconstruct interface states with limited slopes
      - upwind flux based on sign of v_j
      - compute flux divergence

    Returns RHS = - d_x (v f) in conservative form.
    """
    Nx, Nv = f.shape
    rhs = np.zeros_like(f)

    # periodic shifts
    f_im1 = np.roll(f, 1, axis=0)
    f_ip1 = np.roll(f, -1, axis=0)

    # cell-centered slopes (limited)
    df_back = f - f_im1
    df_forw = f_ip1 - f
    slope = slope_mc(df_back, df_forw)

    # left/right states at interfaces i+1/2
    # f_{i+1/2}^- from cell i, f_{i+1/2}^+ from cell i+1
    fL = f + 0.5 * slope  # from i
    slope_ip1 = np.roll(slope, -1, axis=0)
    fR = f_ip1 - 0.5 * slope_ip1  # from i+1

    # flux at interfaces i+1/2: F = v * upwind_state
    # For v>=0 use left state; for v<0 use right state
    v = grid.v
    for j, vj in enumerate(v):
        if vj >= 0.0:
            F_iphalf = vj * fL[:, j]
        else:
            F_iphalf = vj * fR[:, j]

        # flux divergence: (F_{i+1/2} - F_{i-1/2})/dx
        F_imhalf = np.roll(F_iphalf, 1, axis=0)
        rhs[:, j] = -(F_iphalf - F_imhalf) / grid.dx

    return rhs


# ----------------------------
# Hu–Shu ExpRK2 step (paper eq. (2.6))
# ----------------------------
def hu_shu_expk2_step(f: np.ndarray, t: float, dt: float, grid: Grid1D, eps_x, transport_rhs,
                      a0=1 / 3, a1=1 / 3, a2=1 / 3, w=1 / 2, b1=1.0, b2=1.0):
    """
    Hu–Shu ExpRK2 scheme for BGK.

    Parameters (recommended in paper):
      w=1/2, b1=b2=1, a0=a1=a2=1/3
    """

    def expQ(fin, dtsub):
        return exp_bgk(fin, dtsub, eps_x, grid)

    # f^(0) = exp(a0 dt Q/eps) f^n
    f0 = expQ(f, a0 * dt)

    # f^(1) = exp(a1 dt Q/eps) ( f^(0) + b1 dt T(f^(0)) )
    Tf0 = transport_rhs(f0, grid)
    f1 = expQ(f0 + b1 * dt * Tf0, a1 * dt)

    # f^(2) = f^(1) + b2 dt T(f^(1))
    Tf1 = transport_rhs(f1, grid)
    f2 = f1 + b2 * dt * Tf1

    # f^{n+1} = exp(a2 dt Q/eps) ( w f^(2) + (1-w) exp((1-a2)dt Q/eps) f^n )
    fn_relax = expQ(f, (1.0 - a2) * dt)
    mix = w * f2 + (1.0 - w) * fn_relax
    f_next = expQ(mix, a2 * dt)
    return f_next


def make_hu_shu_operator(grid: Grid1D, eps_x, transport_name: str):
    transport_rhs = transport_rhs_upwind if transport_name == "upwind1" else transport_rhs_muscl2

    def op(t, dt, y):
        return hu_shu_expk2_step(y, t, dt, grid, eps_x, transport_rhs)

    return op


# ----------------------------
# Reference SSP-RK2 (Heun) on full RHS
# ----------------------------
def rhs_full(f: np.ndarray, t: float, grid: Grid1D, eps_x):
    """
    Full kinetic RHS for BGK:
      f_t = -v d_x f + (1/eps)(M[f] - f)
    """
    Tpart = transport_rhs_upwind(f, grid)

    # collision part
    rho, u, Temp = moments(f, grid.v, grid.dv)
    M = maxwellian(rho, u, Temp, grid.v)

    if np.isscalar(eps_x):
        Cpart = (M - f) / float(eps_x)
    else:
        Cpart = (M - f) / eps_x[:, None]

    return Tpart + Cpart


def ssp_rk2_solve(f0: np.ndarray, t0: float, tf: float, dt: float, grid: Grid1D, eps_x):
    """
    SSP-RK2 / Heun:
      f1 = f^n + dt * L(f^n)
      f^{n+1} = 0.5 f^n + 0.5 ( f1 + dt * L(f1) )
    """
    f = f0.copy()
    t = t0
    while t < tf - 1e-12:
        dtn = min(dt, tf - t)
        k1 = rhs_full(f, t, grid, eps_x)
        f1 = f + dtn * k1
        k2 = rhs_full(f1, t + dtn, grid, eps_x)
        f = 0.5 * f + 0.5 * (f1 + dtn * k2)
        t += dtn
    return f


# ----------------------------
# Paper setups
# ----------------------------
def eps_mixed_paper(x: np.ndarray):
    """
    Paper eq. (7.6):
      eps(x) = eps0 + (tanh(1 - 11(x-1)) + tanh(1 + 11(x-1))), eps0=1e-5
    """
    eps0 = 1e-5
    return eps0 + (np.tanh(1.0 - 11.0 * (x - 1.0)) + np.tanh(1.0 + 11.0 * (x - 1.0)))


def initial_inconsistent_maxwellians(grid: Grid1D):
    """
    Paper accuracy-test IC:
      f(0,x,v)=0.5 M(rho,u,T) + 0.3 M(rho,-0.5u,T)
      rho=1+0.2 sin(pi x), u=1, T=1/(1+0.2 sin(pi x))
    """
    x = grid.x
    rho = 1.0 + 0.2 * np.sin(np.pi * x)
    u = np.ones_like(x)
    T = 1.0 / (1.0 + 0.2 * np.sin(np.pi * x))
    M1 = maxwellian(rho, u, T, grid.v)
    M2 = maxwellian(rho, -0.5 * u, T, grid.v)
    return 0.5 * M1 + 0.3 * M2


# ----------------------------
# Experiment A: Mixed regime (Figure-2-style)
# ----------------------------
def run_mixed_regime(outdir: str, transport_name: str):
    os.makedirs(outdir, exist_ok=True)

    # paper-like defaults
    Nx, Nv, vmax = 40, 150, 15.0
    grid = make_grid(Nx=Nx, Nv=Nv, vmax=vmax)
    eps_x = eps_mixed_paper(grid.x)

    # dt guideline used in paper accuracy test; used here as a reasonable baseline too
    dt = 0.5 * grid.dx / grid.vmax

    t0, tf = 0.0, 0.5

    f0 = initial_inconsistent_maxwellians(grid)

    # Hu–Shu via pythOS driver
    hu_op = make_hu_shu_operator(grid, eps_x, transport_name)
    f_hu = fractional_step(
        functions=[hu_op],
        delta_t=dt,
        initial_y=f0,
        initial_t=t0,
        final_t=tf,
        alpha=[[1.0]],
        methods={(1,): 'ANALYTIC'}
    )

    # SSP-RK2 reference: smaller dt
    dt_ref = dt / 4.0
    hu_op_ref = make_hu_shu_operator(grid, eps_x, transport_name)
    f_ref = fractional_step(
        functions=[hu_op_ref],
        delta_t=dt_ref,
        initial_y=f0,
        initial_t=t0,
        final_t=tf,
        alpha=[[1.0]],
        methods={(1,): "ANALYTIC"},
    )

    # Macros
    rho_h, u_h, T_h = moments(f_hu, grid.v, grid.dv)
    rho_r, u_r, T_r = moments(f_ref, grid.v, grid.dv)

    # --- Metrics vs refined reference (sanity checks) ---
    def linf(a, b):
        return float(np.max(np.abs(a - b)))

    def l2(a, b, dx):
        return float(np.sqrt(np.sum((a - b) ** 2) * dx))

    print("Differences vs refined Hu–Shu reference:")
    print("  rho: Linf =", linf(rho_h, rho_r), "  L2 =", l2(rho_h, rho_r, grid.dx))
    print("    u: Linf =", linf(u_h, u_r), "  L2 =", l2(u_h, u_r, grid.dx))
    print("    T: Linf =", linf(T_h, T_r), "  L2 =", l2(T_h, T_r, grid.dx))

    # Save
    np.savez(
        os.path.join(outdir, "mixed_regime_bgk.npz"),
        x=grid.x, v=grid.v,
        rho_hu=rho_h, u_hu=u_h, T_hu=T_h,
        rho_ref=rho_r, u_ref=u_r, T_ref=T_r,
        eps_x=eps_x
    )

    # quick console sanity
    print("Mixed regime done.")
    print("Hu–Shu rho range:", rho_h.min(), rho_h.max())
    print("Ref   rho range:", rho_r.min(), rho_r.max())
    print("Hu–Shu min(f):", float(np.min(f_hu)))
    print("Ref   min(f):", float(np.min(f_ref)))
    print("Saved:", os.path.join(outdir, "mixed_regime_bgk.npz"))

    # --- Plots ---
    try:
        import matplotlib.pyplot as plt

        # eps plot
        plt.figure()
        plt.plot(grid.x, eps_x)
        plt.xlabel("x")
        plt.ylabel("epsilon(x)")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "mixed_eps.png"), dpi=200)
        plt.close()

        for name, a, b in [("rho", rho_h, rho_r), ("u", u_h, u_r), ("T", T_h, T_r)]:
            plt.figure()
            plt.plot(grid.x, a, label="Hu–Shu ExpRK2")
            plt.plot(grid.x, b, "--", label="Ref (Hu–Shu, smaller dt)")
            plt.xlabel("x")
            plt.ylabel(name)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f"mixed_{name}.png"), dpi=200)
            plt.close()

        print("Plots saved to:", outdir)

    except Exception as e:
        print("Plotting skipped.")
        print("Reason:", repr(e))


# ----------------------------
# Experiment B: Table-1-style accuracy sweep
# ----------------------------
def l2_xv_error(f: np.ndarray, g: np.ndarray, dx: float, dv: float):
    """
    L2 over x,v: sqrt( sum_{i,j} (f-g)^2 * dx * dv )
    """
    return float(np.sqrt(np.sum((f - g) ** 2) * dx * dv))


def run_accuracy_table(outdir: str, transport_name: str):
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, "table1_like_errors.csv")

    # Paper parameters
    Nv, vmax = 150, 15.0
    tf = 0.1

    # Typical Nx list (edit to match exactly what you want)
    Nx_list = [10, 20, 40, 80]
    eps_list = [1.0, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10]

    rows = []
    header = ["Nx", "dt"] + [f"eps={e:g}" for e in eps_list]

    for Nx in Nx_list:
        grid = make_grid(Nx=Nx, Nv=Nv, vmax=vmax)
        dt = 0.5 * grid.dx / grid.vmax  # paper choice for Table 1

        # refined grid for reference: Nx_ref = 2Nx, dt_ref = dt/2, same Nv
        grid_ref = make_grid(Nx=2 * Nx, Nv=Nv, vmax=vmax)
        dt_ref = dt / 2.0

        # initial data defined on each grid
        f0 = initial_inconsistent_maxwellians(grid)
        f0_ref = initial_inconsistent_maxwellians(grid_ref)

        errs = []
        for eps in eps_list:
            # Hu–Shu coarse
            hu_op = make_hu_shu_operator(grid, eps, transport_name)
            f_coarse = fractional_step(
                functions=[hu_op],
                delta_t=dt,
                initial_y=f0,
                initial_t=0.0,
                final_t=tf,
                alpha=[[1.0]],
                methods={(1,): 'ANALYTIC'}
            )

            # Hu–Shu refined
            hu_op_ref = make_hu_shu_operator(grid_ref, eps, transport_name)
            f_fine = fractional_step(
                functions=[hu_op_ref],
                delta_t=dt_ref,
                initial_y=f0_ref,
                initial_t=0.0,
                final_t=tf,
                alpha=[[1.0]],
                methods={(1,): 'ANALYTIC'}
            )

            # Compare by restricting fine to coarse grid in x (take every other point)
            f_fine_restricted = f_fine[::2, :]
            err = l2_xv_error(f_coarse, f_fine_restricted, grid.dx, grid.dv)
            errs.append(err)

            print(f"Nx={Nx:4d} eps={eps:>8g}  err={err:.6e}")

        row = [Nx, dt] + errs
        rows.append(row)

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)

    print("Accuracy table CSV saved:", csv_path)


def run_time_order(outdir: str, transport_name: str):
    """
    Time-order test at fixed spatial grid (isolates temporal accuracy).
    We keep Nx, Nv fixed and vary dt. Reference = same method with very small dt.
    """
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, "time_order.csv")

    # Fix spatial resolution (pick something not too small)
    Nx, Nv, vmax = 160, 150, 15.0
    grid = make_grid(Nx=Nx, Nv=Nv, vmax=vmax)

    # Use the same IC as before (paper's inconsistent Maxwellian mix)
    f0 = initial_inconsistent_maxwellians(grid)

    # Choose eps values to test time-order robustness across regimes
    eps_list = [1.0, 1e-2, 1e-6, 1e-10]

    # Final time (shorter makes this fast + keeps solution smooth)
    tf = 0.1

    # Base dt tied to dx, then refine in time only
    dt0 = 0.5 * grid.dx / grid.vmax
    dt_list = [dt0, dt0 / 2.0, dt0 / 4.0, dt0 / 8.0]

    # Reference dt: much smaller than smallest dt in dt_list
    dt_ref = dt0 / 512.0

    def l2_xv(f, g):
        return float(np.sqrt(np.sum((f - g) ** 2) * grid.dx * grid.dv))

    rows = []
    header = ["eps", "dt", "L2_xv_error", "observed_order_vs_prev"]
    rows.append(header)

    for eps in eps_list:
        # Reference solution at dt_ref using same Hu–Shu method
        hu_op_ref = make_hu_shu_operator(grid, eps, transport_name)
        f_ref = fractional_step(
            functions=[hu_op_ref],
            delta_t=dt_ref,
            initial_y=f0,
            initial_t=0.0,
            final_t=tf,
            alpha=[[1.0]],
            methods={(1,): "ANALYTIC"},
        )

        prev_err = None
        prev_dt = None

        print(f"\nTime-order test for eps={eps:g} (Nx={Nx}, Nv={Nv}, tf={tf})")
        for dt in dt_list:
            hu_op = make_hu_shu_operator(grid, eps)
            f_dt = fractional_step(
                functions=[hu_op],
                delta_t=dt,
                initial_y=f0,
                initial_t=0.0,
                final_t=tf,
                alpha=[[1.0]],
                methods={(1,): "ANALYTIC"},
            )

            err = l2_xv(f_dt, f_ref)

            if prev_err is None:
                order = ""
            else:
                # observed order between successive dt refinements:
                # err ~ C * dt^p => p ≈ log(err_prev/err)/log(dt_prev/dt)
                order_val = np.log(prev_err / err) / np.log(prev_dt / dt)
                order = f"{order_val:.3f}"

            rows.append([f"{eps:g}", f"{dt:.16e}", f"{err:.16e}", order])

            if order == "":
                print(f"  dt={dt:.3e}  err={err:.6e}")
            else:
                print(f"  dt={dt:.3e}  err={err:.6e}  observed p≈{order}")

            prev_err = err
            prev_dt = dt

    # write CSV
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        for r in rows:
            w.writerow(r)

    print("\nTime-order CSV saved:", csv_path)


# ----------------------------
# Main CLI
# ----------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=["mixed", "table", "timeorder"], required=True)
    parser.add_argument("--outdir", default="experiments/hu_shu_outputs")
    parser.add_argument("--transport", choices=["upwind1", "muscl2"], default="upwind1")
    args = parser.parse_args()

    if args.case == "mixed":
        run_mixed_regime(args.outdir, args.transport)
    elif args.case == "table":
        run_accuracy_table(args.outdir, args.transport)
    else:
        run_time_order(args.outdir, args.transport)


if __name__ == "__main__":
    main()
