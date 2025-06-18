# Niederer benchmark
from __future__ import division
from scipy import sparse

import fractional_step as fs
import numpy as np
import scipy.linalg as la
import math
from math import *
from numpy import *
import ttp

# Monodomain parameters
chi = 1400.0  # SurfaceAreaToVolumeRatio
C = 1.0   # Capacitance
sigma_l = 0.17 * 0.62 / (0.17 + 0.62) * 10
sigma_t = 0.019 * 0.24 / (0.019 + 0.24) * 10 

xn = 2
yn = 0.7
zn = 0.3
dx = 0.05
dy = dx
dz = dx
Nx = int(round((xn + dx) / dx))
Ny = int(round((yn + dy) / dy))
Nz = int(round((zn + dz) / dz))
N = Nx * Ny * Nz
xi = np.linspace(0, xn, num=Nx)
yi = np.linspace(0, yn, num=Ny)
zi = np.linspace(0, zn, num=Nz)

stim_x = sum(xi <=.15)
stim_y = sum(yi <=0.15)
stim_z = sum(zi <=0.15)

# Create the Laplacian Matrix, M using central differences
# Central difference in x:
firstcol = np.zeros(Nx)
firstcol[0] = -2
firstcol[1] = 1
Mx = la.toeplitz(firstcol)
Mx[0,1] = 2
Mx[-1,-2] = 2
Mx = 1/(dx**2)*Mx

Mx = sparse.kron(Mx, sparse.identity(Nz*Ny))
firstcol = np.zeros(Ny)
firstcol[0] = -2
firstcol[1] = 1
My = la.toeplitz(firstcol)
My[0,1] = 2
My[-1,-2] = 2
My = 1/(dy**2)*My
My = sparse.kron(sparse.identity(Nx), sparse.kron(My, sparse.identity(Nz)))

firstcol = np.zeros(Nz)
firstcol[0] = -2
firstcol[1] = 1
Mz = la.toeplitz(firstcol)
Mz[0,1] = 2
Mz[-1,-2] = 2
Mz = 1/(dz**2)*Mz

Mz = sparse.kron(sparse.identity(Nx*Ny), Mz)
M = sigma_l * Mx + sigma_t * My + sigma_t * Mz

# operators:
# First operator: dV/dt = 1/(chi*C)*sigma*M*V
# Diffusion
def fD(t, y):
    V = y[0:N]
    dVdt = 1/(chi*C)*(M @ V)
    dWdt = np.zeros(18*N)
    dydt = np.append(dVdt, dWdt)
    return dydt

# Second operator: Reaction
stim_amplitude = -50000
nVars = 19
def f_TTP(t,y):
    #print(t)
    dydt = np.zeros(nVars*N)

    i_Stim = np.zeros((Nx, Ny, Nz))
    i_Stim[:stim_x,:stim_y,:stim_z] = stim_amplitude / chi / C
    i_Stim = i_Stim.flatten()
    states,constants = ttp.initConsts()
    constants[5] = 0
    constants[6] = 40
    constants[7] = 2
    constants[8] = i_Stim
    for var in range(nVars):
        states[var] = y[var*N:(var+1)*N]
    rate = ttp.computeRates(t, states, constants)
    for var in range(nVars):
        dydt[var*N:(var+1)*N] = rate[var]
    return dydt

# initial conditions
y0 = np.zeros(nVars*N)

states, _ = ttp.initConsts()
for var in range(nVars):
    y0[var*N:(var+1)*N] = states[var]

t0=0.
tf=4
