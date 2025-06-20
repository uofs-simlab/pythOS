# Stiff Brusselator work precision plotting
import numpy as np
from math import pi

import sys
sys.path.insert(1, sys.path[0]+"/../")

a = 0.6
b = 2
ep = 1e-3

a_u = 1e-2
a_v = 1e-2
a_w = 1e-2
p_u = 1e-3
p_v = 1e-3
p_w = 1e-3

N = 201

x0 = 0
x1 = 1

dx = 1/(N-1)

x = np.array([x0 + (i) * dx for i in range(N)])

u0 = a + 0.1 * np.sin(pi*x)
v0 = b/a + 0.1 * np.sin(pi*x)
w0 = b + 0.1*np.sin(pi*x)
t0 = 0
tf = 10

def fi(t, y):
    u = y[:N]
    v = y[N:2*N]
    w = y[2*N:]
    dxx_u = np.zeros(N)
    dxx_u[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / (dx**2)
    dxx_v = np.zeros(N)
    dxx_v[1:-1] = (v[2:] - 2 * v[1:-1] + v[:-2]) / (dx**2)
    dxx_w = np.zeros(N)
    dxx_w[1:-1] = (w[2:] - 2* w[1:-1] + w[:-2]) / (dx**2)

    return np.concatenate((a_u * dxx_u, a_v*dxx_v, a_w*dxx_w))

def fe(t, y):
    u = y[:N]
    v = y[N:2*N]
    w = y[2*N:]
    
    dx_u = np.zeros(N)
    dx_u[1:-1] = (u[2:] - u[:-2]) / (dx*2)
    dx_v = np.zeros(N)
    dx_v[1:-1] = (v[2:] - v[:-2]) / (dx*2)
    dx_w = np.zeros(N)
    dx_w[1:-1] = (w[2:] -  w[:-2]) / (dx*2)
    
    return np.concatenate((p_u * dx_u, p_v*dx_v, p_w*dx_w))

def ff(t, y):
    u = y[:N]
    v = y[N:2*N]
    w = y[2*N:]

    dudt = a - (w+1)*u + u**2*v
    dudt[0] = 0
    dudt[-1] = 0

    dvdt = w*u - u**2*v
    dvdt[0] = 0
    dvdt[-1] = 0

    dwdt = (b-w)/ep - w*u
    dwdt[0] = 0
    dwdt[-1] = 0

    return np.concatenate((dudt, dvdt, dwdt))

y0 = np.concatenate((u0, v0, w0))
t0 = 0
tf = 3
dt = 0.1

if __name__ == "__main__":
    import fractional_step as fs
    import additive_rk as ark
    import matplotlib.pyplot as plt
    import timeit
    f = lambda t, y: fe(t, y) + fi(t, y) + ff(t, y)

    ref_sol = fs.fractional_step([f], 1, y0, t0, tf, "Godunov", {(0,): "ADAPTIVE"})

    methods = [
    ('Godunov-2', lambda dt: fs.fractional_step([lambda t, y: fe(t, y) + fi(t,y), ff], dt, y0, t0, tf, 'Godunov', {(1,): "FE", (2,): "RK4"})),
    ("RK4", lambda dt: fs.fractional_step([lambda t, y: fe(t,y) + ff(t,y) + fi(t,y)], dt, y0, t0, tf, 'Godunov', {(1,): "RK4"})),
    ("BE", lambda dt: fs.fractional_step([lambda t,y: fe(t,y) + ff(t,y) + fi(t,y)], dt, y0, t0, tf, 'Godunov', {(1,): "BE"})),
    ("IMEX", lambda dt: ark.ark_solve([ff, lambda t, y: fe(t,y) + fi(t,y)], dt, y0, t0, tf, [fs.Tableau(np.array([0,1]), np.array([[0,0],[1,0]]), np.array([1,0])), fs.Tableau(np.array([0,1]), np.array([[0,0],[0,1]]),np.array([0,1]))]))
]

    for method_t in methods:
        errors = []
        times = []

        dt = tf
        while len(times) < 4:
            dt = dt / 2
            try:
                result = method_t[1](dt)
                time = min(timeit.repeat(lambda : method_t[1](dt),number=1,repeat=1))

                error = np.linalg.norm((result - ref_sol) / ref_sol)

                if error < 0.15:
                    errors.append(error)
                    times.append(time)
                    print(time, error)
            except Exception as e:
                continue
        plt.loglog(times, errors, label=method_t[0],marker='o')

    plt.xlabel('CPU Time')
    plt.ylabel('Error')
    plt.legend()
    plt.show()
        
    

