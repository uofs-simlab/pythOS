"""
Reusable adaptive driver for operator splitting methods.

The core idea is step-doubling:

- One full step of size h
- Two half-steps of size h/2
- Use the difference as an error estimate
- Accept/reject the step and adapt h

This is intentionally generic so it can wrap ABC, BGK, or other split steppers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple, Any

import math
import time
import numpy as np


ArrayLike = Any


def _to_numpy(x: ArrayLike) -> np.ndarray:
    """
    Convert state to numpy array for error calculations.

    Supports:
    - numpy arrays
    - scalars
    - lists/tuples
    """
    if isinstance(x, np.ndarray):
        return x
    return np.asarray(x, dtype=float)


def l2_norm(x: ArrayLike) -> float:
    arr = _to_numpy(x)
    return float(np.sqrt(np.sum(arr * arr)))


def linf_norm(x: ArrayLike) -> float:
    arr = _to_numpy(x)
    return float(np.max(np.abs(arr)))


def weighted_rms_error(
    y_high: ArrayLike,
    y_low: ArrayLike,
    atol: float,
    rtol: float,
) -> float:
    """
    Weighted RMS error often used in adaptive time stepping.

    err = sqrt( mean( ((y_high - y_low) / (atol + rtol*max(|y_high|,|y_low|)))^2 ) )

    Returns a dimensionless scalar. Accept if err <= 1.
    """
    yh = _to_numpy(y_high)
    yl = _to_numpy(y_low)

    scale = atol + rtol * np.maximum(np.abs(yh), np.abs(yl))
    diff = (yh - yl) / scale
    return float(np.sqrt(np.mean(diff * diff)))


@dataclass
class AdaptiveOptions:
    """
    Configuration for adaptive step control.
    """
    order: int
    atol: float = 1.0e-8
    rtol: float = 1.0e-6
    dt_min: float = 1.0e-12
    dt_max: float = 1.0
    safety: float = 0.9
    growth_max: float = 2.0
    shrink_min: float = 0.2
    max_reject: int = 20
    use_weighted_rms: bool = True
    accept_on_nan_error: bool = False


@dataclass
class StepRecord:
    """
    One accepted or rejected attempt.
    """
    t: float
    dt_try: float
    dt_suggested: float
    accepted: bool
    err_est: float
    cpu_seconds: float
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class AdaptiveResult:
    """
    Full adaptive integration result.
    """
    t_values: List[float]
    y_values: List[ArrayLike]
    records: List[StepRecord]
    n_accept: int
    n_reject: int
    wall_seconds: float
    success: bool
    message: str

    @property
    def final_time(self) -> float:
        return self.t_values[-1]

    @property
    def final_state(self) -> ArrayLike:
        return self.y_values[-1]

    @property
    def dt_min_used(self) -> float:
        if not self.records:
            return 0.0
        return min(r.dt_try for r in self.records if r.accepted)

    @property
    def dt_max_used(self) -> float:
        if not self.records:
            return 0.0
        return max(r.dt_try for r in self.records if r.accepted)

    @property
    def dt_avg_used(self) -> float:
        accepted = [r.dt_try for r in self.records if r.accepted]
        if not accepted:
            return 0.0
        return float(sum(accepted) / len(accepted))


def _default_error_estimate(
    y_two_half: ArrayLike,
    y_full: ArrayLike,
    atol: float,
    rtol: float,
    use_weighted_rms: bool,
    error_norm_fn: Optional[Callable[[ArrayLike], float]] = None,
) -> float:
    """
    Returns a dimensionless error estimate.

    If weighted RMS is used, accept if err <= 1.
    Otherwise, uses user-supplied norm on the raw difference and scales by atol+rtol.
    """
    if use_weighted_rms:
        return weighted_rms_error(y_two_half, y_full, atol=atol, rtol=rtol)

    diff = _to_numpy(y_two_half) - _to_numpy(y_full)
    if error_norm_fn is None:
        err_raw = l2_norm(diff)
    else:
        err_raw = float(error_norm_fn(diff))

    scale = atol + rtol * max(l2_norm(y_two_half), l2_norm(y_full), 1.0)
    return err_raw / scale


def _propose_new_dt(err: float, dt: float, opts: AdaptiveOptions) -> float:
    """
    Standard adaptive controller:

        dt_new = safety * dt * err^(-1/(p+1))

    Because the local error estimator from step-doubling behaves like O(h^(p+1)).
    """
    if not math.isfinite(err):
        factor = opts.shrink_min
    elif err <= 0.0:
        factor = opts.growth_max
    else:
        exponent = -1.0 / (opts.order + 1.0)
        factor = opts.safety * (err ** exponent)
        factor = min(opts.growth_max, max(opts.shrink_min, factor))

    dt_new = dt * factor
    dt_new = min(opts.dt_max, max(opts.dt_min, dt_new))
    return dt_new


def adaptive_integrate(
    stepper: Callable[[ArrayLike, float, float], ArrayLike],
    y0: ArrayLike,
    t0: float,
    tfinal: float,
    dt0: float,
    opts: AdaptiveOptions,
    error_norm_fn: Optional[Callable[[ArrayLike], float]] = None,
    post_accept_callback: Optional[Callable[[float, ArrayLike], None]] = None,
    diagnostics_fn: Optional[Callable[[ArrayLike], Dict[str, Any]]] = None,
) -> AdaptiveResult:
    """
    Adaptive integrator using step-doubling around a user-supplied one-step method.

    Parameters
    ----------
    stepper
        Function stepper(y, t, dt) -> y_next
        This should perform one step of the splitting method of nominal order `opts.order`.
    y0
        Initial state.
    t0
        Initial time.
    tfinal
        Final time.
    dt0
        Initial time step guess.
    opts
        Adaptive options.
    error_norm_fn
        Optional norm function for raw differences if not using weighted RMS.
    post_accept_callback
        Optional hook called after each accepted step.
    diagnostics_fn
        Optional function diagnostics_fn(y) -> dict to store per-step diagnostics.

    Returns
    -------
    AdaptiveResult
    """
    if tfinal <= t0:
        raise ValueError("tfinal must be greater than t0.")
    if dt0 <= 0.0:
        raise ValueError("dt0 must be positive.")
    if opts.order < 1:
        raise ValueError("opts.order must be at least 1.")

    start_all = time.perf_counter()

    t = float(t0)
    y = y0
    dt = min(max(dt0, opts.dt_min), opts.dt_max)

    t_values: List[float] = [t]
    y_values: List[ArrayLike] = [y]
    records: List[StepRecord] = []

    n_accept = 0
    n_reject = 0

    while t < tfinal:
        dt = min(dt, tfinal - t)
        reject_count_this_step = 0

        while True:
            tic = time.perf_counter()

            # One full step
            y_full = stepper(y, t, dt)

            # Two half-steps
            y_half = stepper(y, t, 0.5 * dt)
            y_two_half = stepper(y_half, t + 0.5 * dt, 0.5 * dt)

            err = _default_error_estimate(
                y_two_half=y_two_half,
                y_full=y_full,
                atol=opts.atol,
                rtol=opts.rtol,
                use_weighted_rms=opts.use_weighted_rms,
                error_norm_fn=error_norm_fn,
            )

            cpu = time.perf_counter() - tic
            dt_new = _propose_new_dt(err, dt, opts)

            if not math.isfinite(err) and not opts.accept_on_nan_error:
                accepted = False
            else:
                accepted = (err <= 1.0)

            extra = {}
            if diagnostics_fn is not None:
                try:
                    extra = diagnostics_fn(y_two_half if accepted else y_full)
                except Exception as exc:
                    extra = {"diagnostics_error": str(exc)}

            records.append(
                StepRecord(
                    t=t,
                    dt_try=dt,
                    dt_suggested=dt_new,
                    accepted=accepted,
                    err_est=err,
                    cpu_seconds=cpu,
                    extra=extra,
                )
            )

            if accepted:
                t = t + dt
                y = y_two_half  # higher quality solution than single full step
                t_values.append(t)
                y_values.append(y)
                n_accept += 1

                if post_accept_callback is not None:
                    post_accept_callback(t, y)

                dt = dt_new
                break

            n_reject += 1
            reject_count_this_step += 1
            dt = dt_new

            if dt <= opts.dt_min + 1.0e-30:
                wall = time.perf_counter() - start_all
                return AdaptiveResult(
                    t_values=t_values,
                    y_values=y_values,
                    records=records,
                    n_accept=n_accept,
                    n_reject=n_reject,
                    wall_seconds=wall,
                    success=False,
                    message="Step size reached dt_min during rejection.",
                )

            if reject_count_this_step >= opts.max_reject:
                wall = time.perf_counter() - start_all
                return AdaptiveResult(
                    t_values=t_values,
                    y_values=y_values,
                    records=records,
                    n_accept=n_accept,
                    n_reject=n_reject,
                    wall_seconds=wall,
                    success=False,
                    message="Exceeded maximum rejects for a single advance.",
                )

    wall = time.perf_counter() - start_all
    return AdaptiveResult(
        t_values=t_values,
        y_values=y_values,
        records=records,
        n_accept=n_accept,
        n_reject=n_reject,
        wall_seconds=wall,
        success=True,
        message="Adaptive integration completed successfully.",
    )