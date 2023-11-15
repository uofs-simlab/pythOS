# file to run 1D FHN monodomain problem using Chaste default parameter values

from __future__ import division
import sys
from scipy import sparse
#set path to the Python_OS_Code directory
sys.path.insert(1, sys.path[0]+'/../Python_OS_Code')
from timing_script import time_test
import timeit

import operator_splitting as os
import numpy as np
import scipy.linalg as la
import math

#define the parameters
# Monodomain parameters
chi = 1400.0  # SurfaceAreaToVolumeRatio
C = 1.0   # Capacitance
sigma = 1.75 # Effective Conductivity = lambda/(1+lambda)*sigma_i where lambda = sigma_e/sigma_i


# Cell model parameters 1D FHN
Vrest = -85.
Vthreshold = -70.0
Vpeak = 40.
K = 4.0e-4
L = 0.63
B = 0.013

# dV/dt = -K*(V-Vrest)*(W+(V-Vthreshold)*(V-Vpeak))-i_stim
# dW/dt = L*(V-Vrest) - B*W
# IC: V = voltage + 100*[cos(pi*x/0.2)]^2

# simulation parameters
x1 = 0.
xn = 1.0
Nx = 131
dx=float(xn-x1)/(Nx-1)

# initial condition, may also use a 1 dimensional numpy array here
x = np.linspace(x1, xn, num=Nx)
V0= Vrest + 100*(np.cos(np.pi*x[0:int(round(0.1/dx))]/0.2))**2
V0 = np.append(V0,Vrest*np.ones(Nx-len(V0)))
W0 = np.zeros(Nx)
y0 = np.append(V0,W0)


# DE:
#

# Create the 1D Laplacian Matrix, M using central differences
# Central difference in x:
firstcol = np.zeros(Nx)
firstcol[0] = -2
firstcol[1] = 1
M = la.toeplitz(firstcol)
M[0,1] = 2
M[-1,-2] = 2
M = 1/(dx**2)*M

# operators:
# First operator: dV/dt = 1/(chi*C)*sigma*M*V
# Diffusion
# y = [V;W]
def f1(t, y):
    V = y[0:Nx]
    dVdt = 1/(chi*C)*sigma*np.matmul(M,V)
    dWdt = np.zeros(Nx)
    dydt = np.append(dVdt, dWdt)
    return dydt


# Second operator:
# dV/dt = -K*(V-Vrest)*(W+(V-Vthreshold)*(V-Vpeak))-i_stim
# dW/dt = L*(V-Vrest) - B*W
# Reaction
def f2(t, y):
    V = y[0:Nx]
    W = y[Nx:]
    dVdt = -K*(V-Vrest)*(W+(V-Vthreshold)*(V-Vpeak))
    dWdt = L*(V-Vrest) - B*W
    dydt = np.append(dVdt, dWdt)
    return dydt



# list the operators in the order you wish to use them
f_list=[f1,f2]




# time period and step
t0=0.
tf=1.5
#delta_t=0.2


# splitting method
#OS_method = raw_input("Which OS method would you like to use? (eg:Godunov,SM2,R3, etc.)")
#Intg_method_1 = raw_input("Which integration method would you like to use for operator 1 (eg: FE, BE, RK3, SD2O2,etc.)? ")
#Intg_method_2 = raw_input("Which integration method would you like to use for operator 2 (eg: FE, BE, RK3, SD2O2,etc.)? ")

beta = 0.5

alphas= ['Godunov'] 
Intg_method_1 = 'BE'
Intg_method_2 = 'FE'
# runge--kutta methods

#ivp_method = 'RK45' # default
#rtol = 1e-3  # default
#atol = 1e-6  # default

methods={(0,):'',   #default,
         (1,):Intg_method_1,   #operator 1 default
         (2,):Intg_method_2,  #operator2 default
         # may also define runge--kutta methods for each individual arc using
         # keys (i, j) where i is the operator number and j is the the arc
         # number
         # they do not need to all be specified either
    }

for n in [10,15,20,30,35,50,70,100,130,150,200,300,400,500]: #np.linspace(4,50,24):#
     delta_t = 1.5/n
     print('delta_t = ', delta_t, 'time inverval [', str(t0), ' , ', str(tf), ']')
     for alpha in alphas:
        result = os.operator_splitting(f_list, delta_t, y0, t0, tf, alpha, methods,
                                        fname="./" + str(
                                            alpha) + "+" + str(methods[(1,)]) + "+" + str(
                                            methods[(2,)]) + " on interval " + str(t0) + " to " + str(
                                            tf) + "-dt" + str(delta_t) + ".csv", save_steps=(tf - t0) / delta_t)

    
        print(alpha, ' with ', Intg_method_1, ' and ', Intg_method_2, ', delta t =', delta_t, result)


