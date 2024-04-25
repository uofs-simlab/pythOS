import math as m

import numpy as np
from scipy.integrate import solve_ivp
from butcher_tableau import EmbeddedTableau, embedded_pairs
try:
    from firedrake import Form, Function, Constant, dx, inner, sin, solve
    from firedrake import *
except:
    Form = type(None)
    Function = type(None)
try:
    from sundials_exode import SundialsExode
except:
    SundialsExode = type(None)
def exode(t_out, A, u, tol=(1e-3, 1e-6), task1 = False, method='RK45', **params):
    def fun(t, x):
        ret = A(x)
        for j in range(u.shape[0]-1):
            ret += t**j / m.factorial(j) * u[j+1]
        return ret
    if isinstance(u, list):
        out = [None for _ in range(len(t_out))]
    else:
        out = np.zeros(( len(t_out), u.shape[1]), u.dtype)
    t0 = 0

    if isinstance(u[0], Function):
        tf = A.arguments()[0]
        t0 = Constant(t0)
        F = A
        for j in range(len(u) - 1):
            if u[j+1] != 0:
                F += inner(t0 ** j / m.factorial(j) * u[j+1], tf) * dx
        fun = F

    if isinstance(tol, tuple):
        rtol = tol[0]
        atol = tol[1]
    else:
        rtol = tol
        atol = tol

    if method in embedded_pairs:
        method = embedded_pairs[method]

    y0 = u[0]
    if isinstance(u, list):
        y_s = Function(y0)
    else:
        y0 = u[0].copy()
    t_steps = sorted(t_out)
    
    for i in range(len(t_steps)):
        if isinstance(method, SundialsExode):
            tf = t_steps[i] - t0
            yi = y0.copy()
            y0 = method.solve(fun, y0, t0, tf)
        elif isinstance(method, EmbeddedTableau):
            if isinstance(y0, Function):
                tf = t_steps[i] - t0.values()[0]
                if isinstance(tf, np.complex128):
                    tf = complex(tf)
                else:
                    tf = float(tf)
            else:
                tf = t_steps[i] - t0
            y0 = method.y_step(fun, y0, t0, tf, rtol, **params)
            
        else:
            result = solve_ivp(fun, [t0, t_steps[i]], y0, method=method, rtol = rtol, atol=atol, t_eval = [t0, t_steps[i]], **params)
            if not result.success:
                print('Failed')
                y0 = np.nan * y0
            else:
                y0 = result.y[:,-1]
        out[i] = (y0 * (1/ t_steps[i])).copy()
        if isinstance(y0, Function):
            out[i] = Function(y0)
            out[i].assign(y0 / t_steps[i])
            y0.assign(y_s)
            t0.assign(t_steps[i])
        else:
            t0 = t_steps[i]
    return out, None

