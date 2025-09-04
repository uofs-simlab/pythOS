# Stiff Brusselator work precision plotting
import numpy as np
from math import pi

import sys
sys.path.insert(1, sys.path[0]+"/../")

# define constants
a = 0.6
b = 2
ep = 1e-3

a_u = 1e-2
a_v = 1e-2
a_w = 1e-2
p_u = 1e-3
p_v = 1e-3
p_w = 1e-3

# define spatial discretization
N = 201

x0 = 0
x1 = 1

dx = 1/(N-1)

x = np.array([x0 + (i) * dx for i in range(N)])

# define initial conditions
u0 = a + 0.1 * np.sin(pi*x)
v0 = b/a + 0.1 * np.sin(pi*x)
w0 = b + 0.1*np.sin(pi*x)
t0 = 0
tf = 10

# fi: diffusion term
def fi(t, y):
    # divide input vector into component variables
    u = y[:N]
    v = y[N:2*N]
    w = y[2*N:]
    # compute second derivatives with zero boundary conditions
    dxx_u = np.zeros(N)
    dxx_u[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / (dx**2)
    dxx_v = np.zeros(N)
    dxx_v[1:-1] = (v[2:] - 2 * v[1:-1] + v[:-2]) / (dx**2)
    dxx_w = np.zeros(N)
    dxx_w[1:-1] = (w[2:] - 2* w[1:-1] + w[:-2]) / (dx**2)

    # combine resulting derivatives
    return np.concatenate((a_u * dxx_u, a_v*dxx_v, a_w*dxx_w))

# fe: advection term
def fe(t, y):
    # divide input vector into component variables
    u = y[:N]
    v = y[N:2*N]
    w = y[2*N:]

    # compute first derivatives with zero boundary conditions
    dx_u = np.zeros(N)
    dx_u[1:-1] = (u[2:] - u[:-2]) / (dx*2)
    dx_v = np.zeros(N)
    dx_v[1:-1] = (v[2:] - v[:-2]) / (dx*2)
    dx_w = np.zeros(N)
    dx_w[1:-1] = (w[2:] -  w[:-2]) / (dx*2)

    # combine resulting derivatives
    return np.concatenate((p_u * dx_u, p_v*dx_v, p_w*dx_w))

# ff: reaction term
def ff(t, y):
    # divide input vector into component variables
    u = y[:N]
    v = y[N:2*N]
    w = y[2*N:]

    # compute reaction terms
    dudt = a - (w+1)*u + u**2*v
    dudt[0] = 0
    dudt[-1] = 0

    dvdt = w*u - u**2*v
    dvdt[0] = 0
    dvdt[-1] = 0

    dwdt = (b-w)/ep - w*u
    dwdt[0] = 0
    dwdt[-1] = 0

    # combine results
    return np.concatenate((dudt, dvdt, dwdt))


# fs: for multirate solvers without implicit/explicit division at slow scale
# advection + diffusion
def fs(t, y):
    return fi(t, y) + fe(t, y)

y0 = np.concatenate((u0, v0, w0))
t0 = 0
tf = 3
dt = 0.1

# f1 (for GARK solver): advection
def f1(t, y):
    return fe(t, y)

# f2 (for GARK solver): diffusion + reaction
def f2(t, y):
    return fi(t, y) + ff(t, y)

