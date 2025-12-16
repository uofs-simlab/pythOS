# Functions to test various integrators from within pythOS
import fractional_step as fs
from additive_rk import ark_solve
from gark_methods import gark_solve
from multirate import multirate_solve
from multirate_infinitesimal import multirate_infinitesimal_solve
import numpy as np
import timeit
try:
    from firedrake import *
except:
    Function = type(None)
    
def calc_error(result, reference):
    if isinstance(result,Function):
        return errornorm(result,reference)
    return np.linalg.norm((result-reference)/reference)

def test_fractional_step(f_list,t0,tf,y0,splitting_method,rk_methods,dt0,reference_solution,num_time_repeats=0,N=8,err_max=0.02,err_min=2e-12,dt_min=1e-10,**kwargs):
    err = []
    dts = []
    times = []
    dt = dt0 * 2
    while len(dts) < N:
        dt = dt / 2
        print(dt)
        if (dt < dt_min):
            break
        try:
            result = fs.fractional_step(f_list,dt,y0,t0,tf,splitting_method,rk_methods,**kwargs)
        except:
            continue
        error = calc_error(result,reference_solution)
        if error < err_max:
            err.append(error)
            dts.append(dt)
            if num_time_repeats > 0:
                time = min(timeit.repeat(lambda:  fs.fractional_step(f_list,dt,y0,t0,tf,splitting_method,rk_methods,**kwargs), number=1,repeat=num_time_repeats))
                times.append(time)
        print(dt, error)
        if error < err_min:
            break
    return (err, dts, times)

def test_additive_rk(f_list,t0,tf,y0,tableaus,dt0,reference_solution,num_time_repeats=0,N=8,err_max=0.02,err_min=2e-12,dt_min=1e-10,**kwargs):
    err = []
    dts = []
    times = []
    dt = dt0 * 2
    while len(dts) < N:
        dt = dt / 2
        if (dt < dt_min):
            break
        result = ark_solve(f_list,dt,y0,t0,tf,tableaus,**kwargs)
        error = calc_error(result,reference_solution)
        if error < err_max:
            err.append(error)
            dts.append(dt)
            if num_time_repeats > 0:
                time = min(timeit.repeat(lambda: ark_solve(f_list,dt,y0,t0,tf,tableaus,**kwargs), number=1,repeat=num_time_repeats))
                times.append(time)
        if error < err_min:
            break
    return (err, dts, times)

def test_gark(f_list,t0,tf,y0,A,b,dt0,reference_solution,num_time_repeats=0,N=8,err_max=0.02,err_min=2e-12,dt_min=1e-10,**kwargs):
    err = []
    dts = []
    times = []
    dt = dt0 * 2
    while len(dts) < N:
        dt = dt / 2
        if (dt < dt_min):
            break
        result = gark_solve(f_list,dt,y0,t0,tf,A,b,**kwargs)
        error = calc_error(result,reference_solution)
        if error < err_max:
            err.append(error)
            dts.append(dt)
            if num_time_repeats > 0:
                time = min(timeit.repeat(lambda: gark_solve(f_list,dt,y0,t0,tf,A,b,**kwargs), number=1,repeat=num_time_repeats))
                times.append(time)
        if error < err_min:
            break
    return (err, dts, times)


def test_multirate(fs,ff,t0,tf,y0,multirate_method,M,dt0,reference_solution,num_time_repeats=0,N=8,err_max=0.02,err_min=2e-12,dt_min=1e-10,**kwargs):
    err = []
    dts = []
    times = []
    dt = dt0 * 2
    while len(dts) < N:
        dt = dt / 2
        if (dt < dt_min):
            break
        try:
            result = multirate_solve(y0,t0,dt,tf,multirate_method,M,ff,fs,**kwargs)
        except:
            continue
        error = calc_error(result,reference_solution)
        if error < err_max:
            err.append(error)
            dts.append(dt)
            if num_time_repeats > 0:
                time = min(timeit.repeat(lambda: multirate_solve(y0,t0,dt,tf,multirate_method,M,ff,fs,**kwargs), number=1,repeat=num_time_repeats))
                times.append(time)
        if error < err_min:
            break
    return (err, dts, times)

def test_multirate_infinitesimal(fs,ff,t0,tf,y0,multirate_method,dt0,reference_solution,fe=None,num_time_repeats=0,N=8,err_max=0.02,err_min=2e-12,dt_min=1e-10,**kwargs):
    err = []
    dts = []
    times = []
    dt = dt0 * 2
    while len(dts) < N:
        dt = dt / 2
        print(dt)
        if (dt < dt_min):
            break
        try:
            result = multirate_infinitesimal_solve(y0,t0,dt,tf,multirate_method,fs,ff,fe=fe,**kwargs)
        except:
            continue
        error = calc_error(result,reference_solution)
        if error < err_max:
            err.append(error)
            dts.append(dt)
            if num_time_repeats > 0:
                time = min(timeit.repeat(lambda: multirate_infinitesimal_solve(y0,t0,dt,tf,multirate_method,fs,ff,fe=fe,**kwargs), number=1,repeat=num_time_repeats))
                times.append(time)
        print(error)
        if error < err_min:
            break
    return (err, dts, times)


