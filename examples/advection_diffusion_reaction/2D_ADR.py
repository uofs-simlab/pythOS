from firedrake import *

# Create the domain - square mesh
mesh = UnitSquareMesh(40, 40)
mesh.init()
U = FunctionSpace(mesh, 'CG', 2)
u = Function(U)

tu = TestFunction(U)

# Define the initial condition and constants
x, y = SpatialCoordinate(mesh)

u.project(256 * (x * y * (1-x) * (1-y))**2 + 0.3)
epsilon = 1/100
alpha = -10
gamma = 100

# Define the operators
# F1: diffusion
F1 = -epsilon * inner(grad(u), grad(tu)) * dx
# F1a: -epsilon * u_xx (for 4-splitting)
F1a = -epsilon * inner(grad(u)[0], grad(tu)[0]) * dx
# F1b: -epsilon * u_yy (for 4-splitting)
F1b = -epsilon * inner(grad(u)[1], grad(tu)[1]) * dx
# F2: advection
F2 = -alpha * (grad(u)[0] * tu + grad(u)[1] * tu) * dx

# F3: reaction
F3 = gamma * u * (u - 1/2) * (1-u) * tu * dx

u0 = Function(u)

import fractional_step as fs
t = Constant(0)

import time as timeos
import cmath as m
print(u((0,0)), u((0.5,0.5)))
# Create reference solution for order comparison
ref_func = fs.fractional_step([F1+ F2+ F3], 0.01, u, Constant(0), 0.1, "Godunov", methods={(0,): "ADAPTIVE"}, ivp_methods={1:('Dormand-Prince', 1e-14,1e-14)})
ref_func = Function(ref_func)
alphas = ["Godunov-N",
    "Strang-N",
    "PP3_4A-3",
    "Yoshida-3",
]

# Run order convergence study for 3-splitting, outputting a table of results
for alpha in alphas:

    prev_err = 1.
    prev_dt = 2.
    print(alpha)
    print("{:<15} {:<15} {:<10} {:<10}".format("dt", "Error", "Order", "Time"))

    if (alpha == 'Godunov-N'):
        methods={(0,): "FE"}
    elif alpha == 'Strang-N' or alpha == 'OSNNS2P2-N':
        methods={(0,): "Heun"}
    elif alpha == 'PP3_4A-3' or alpha == 'OSNNS4P3-N':
        methods={(0,): "RK3"}
    else:
        methods = {(0,): "RK4"}
    for N_step in [200, 400,
                   800, 1600,
            3200, 6400,
            12800]:
        dt = .1 / N_step
        tic = timeos.perf_counter()
        u.project(256 * (x * y * (1-x) * (1-y))**2 + 0.3)
        result = fs.fractional_step([F1, F2, F3], dt, u, Constant(0), 0.1, alpha, methods, fname=f"3split_{alpha}_{N_step}.h5", save_steps=1)
        toc = timeos.perf_counter()
        error = errornorm(u, ref_func)
        print("{:<15.8f} {:<15.3e} {:<10.3f} {:<10.4f}".format(dt, error, (m.log(prev_err/error, prev_dt/dt)).real, toc-tic))
        prev_err = error
        prev_dt = dt
        u.assign(u0)

alphas = ["Strang-N",
    "OSNNS4P3-N",
    "OSNNS2P2-N"
]


# Run order convergence study for 4-splitting, outputting a table of results
for alpha in alphas:

    prev_err = 1.
    prev_dt = 2.
    print(alpha)
    print("{:<15} {:<15} {:<10} {:<10}".format("dt", "Error", "Order", "Time"))

    if (alpha == 'Godunov-N'):
        methods={(0,): "FE"}
    elif alpha == 'Strang-N' or alpha == 'OSNNS2P2-N':
        methods={(0,): "Heun"}
    elif alpha == 'PP3_4A-3' or alpha == 'OSNNS4P3-N':
        methods={(0,): "RK3"}
    else:
        methods = {(0,): "RK4"}
    for N_step in [200, 400,
                   800, 1600,
            3200, 6400,
            12800]:
        dt = .1 / N_step
        tic = timeos.perf_counter()
        u.assign(u0)
        result = fs.fractional_step([F1a, F1b, F2, F3], dt, u, Constant(0), 0.1, alpha, methods, fname=f"3split_{alpha}_{N_step}.h5", save_steps=1)
        toc = timeos.perf_counter()
        error = errornorm(u, ref_func)
        print("{:<15.8f} {:<15.3e} {:<10.3f} {:<10.4f}".format(dt, error, (m.log(prev_err/error, prev_dt/dt)).real, toc-tic))
        prev_err = error
        prev_dt = dt
        u.assign(u0)

