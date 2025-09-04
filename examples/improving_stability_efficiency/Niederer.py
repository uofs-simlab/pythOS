# Niederer benchmark
from __future__ import division
import sys
from scipy import sparse
#set path to the Python_OS_Code directory
sys.path.insert(1, sys.path[0]+'/../../')

import fractional_step as fs
import numpy as np
import scipy.linalg as la
import math

#define the parameters
# Monodomain parameters
chi = 1400.0  # SurfaceAreaToVolumeRatio
C = 1.0   # Capacitance
sigma = 1.75 # Effective Conductivity = lambda/(1+lambda)*sigma_i where lambda = sigma_e/sigma_i
sigma_l = 0.17 * 0.62 / (0.17 + 0.62) * 10
sigma_t = 0.019 * 0.24 / (0.019 + 0.24) * 10 

Vrest = -85.
Vthreshold = -70.0
Vpeak = 40.
K = 4.0e-4
L = 0.63
B = 0.013

from numpy import log, exp, sqrt

def pow(y, n):
    return y ** n

x1 = 0
xn = 2
y1 = 0
yn = 0.7
z1 = 0
zn = 0.3
dx = 0.05
dy = dx
dz = dx
Nx = int(round((xn -x1 + dx) / dx))
Ny = int(round((yn - y1 + dy) / dy))
Nz = int(round((zn - z1 + dz) / dz))
print(Nx, Ny, Nz)
N = Nx * Ny * Nz
xi = np.linspace(x1, xn, num=Nx)
yi = np.linspace(y1, yn, num=Ny)
zi = np.linspace(z1, zn, num=Nz)

stim_x = sum(xi <=.15)
stim_y = sum(yi <=0.15)
stim_z = sum(zi <=0.15)

V0 = -85.23 * np.ones(N)
Xr10 = 0.00621 * np.ones(N)
Xr20 = 0.4712 * np.ones(N)
Xs0 = 0.0095 * np.ones(N)
m0 = 0.00172 * np.ones(N)
h0 = 0.7444 * np.ones(N)
j0 = 0.7045 * np.ones(N)
d0 = 3.373e-5 * np.ones(N)
f0 = 0.7888 * np.ones(N)
f20 = 0.9755 * np.ones(N)
fCass0 = 0.9953 * np.ones(N)
s0 = 0.999998 * np.ones(N)
r0 = 2.42e-8 * np.ones(N)
R_prime0 = 0.9073 * np.ones(N)
Ca_i0 = 0.000126 * np.ones(N)
Ca_Sr0 = 3.64 * np.ones(N)
Ca_ss0 = 0.00036 * np.ones(N)
Na_i0 = 8.604 * np.ones(N)
K_i0 = 136.89 * np.ones(N)

y0 = np.concatenate((V0, K_i0, Na_i0, Ca_i0, Xr10, Xr20, Xs0, m0,h0,j0,Ca_ss0,d0,f0,f20,fCass0,s0,r0,Ca_Sr0,R_prime0))
#y0 = np.concatenate((V0, Xr10, Xr20, Xs0, m0, h0, j0, d0,f0,f20,fCass0,s0,r0,R_prime0,Ca_i0,Ca_Sr0,Ca_ss0,Na_i0,K_i0))
ans = [2.31955735e+01, 3.98881690e-01, 1.20461456e-01, 1.36553037e-02,
                9.99728227e-01, 5.61883611e-08, 6.51708166e-08, 9.87265864e-01,
                7.78533609e-01, 8.96946333e-01, 9.95177252e-01, 4.40475973e-01,
                7.77571418e-01, 9.08708925e-01, 1.26411299e-04, 3.63986901e+00,
                3.15789877e-02, 8.61245039e+00, 1.36893465e+02]
(V0_i, Xr1_i, Xr2_i, Xs_i, m_i, h_i, j_i, d_i,f_i,f2_i,fCass_i,s_i,r_i,R_prime_i,Ca_i_i,Ca_Sr_i,Ca_ss_i,Na_i_i,K_i_i) = ans

#ans = np.array([V0_i,K_i_i,Na_i_i,Ca_i_i, Xr1_i, Xr2_i, Xs_i, m_i,h_i,j_i,Ca_ss_i,d_i,f_i,f2_i,fCass_i,s_i,r_i,Ca_Sr_i,R_prime_i])

# Create the Laplacian Matrix, M using central differences
# Central difference in x:
firstcol = np.zeros(Nx)
firstcol[0] = -2
firstcol[1] = 1
Dxx = la.toeplitz(firstcol)
Dxx[0,1] = 2
Dxx[-1,-2] = 2
Dxx = 1/(dx**2)*Dxx

Dxx = sparse.kron(Dxx, sparse.identity(Nz*Ny))
firstcol = np.zeros(Ny)
firstcol[0] = -2
firstcol[1] = 1
Dyy = la.toeplitz(firstcol)
Dyy[0,1] = 2
Dyy[-1,-2] = 2
Dyy = 1/(dy**2)*Dyy
Dyy = sparse.kron(sparse.identity(Nx), sparse.kron(Dyy, sparse.identity(Nz)))

firstcol = np.zeros(Nz)
firstcol[0] = -2
firstcol[1] = 1
Dzz = la.toeplitz(firstcol)
Dzz[0,1] = 2
Dzz[-1,-2] = 2
Dzz = 1/(dz**2)*Dzz

Dzz = sparse.kron(sparse.identity(Nx*Ny), Dzz)

# operators:
# First operator: dV/dt = 1/(chi*C)*sigma*M*V
# Diffusion
from scipy.sparse.linalg import expm
D = sigma_l * Dxx + sigma_t * Dyy + sigma_t * Dzz
def f1(t, y):
    V = y[0:N]
    dVdt = 1/(chi*Cm)*(D @ V)
    dWdt = np.zeros(18*N)
    dydt = np.append(dVdt, dWdt)
    return dydt

stim_amplitude = -50000
stim_period = 40
stim_start = 0
stim_duration = 2
# Second operator:
# Reaction - TTP model
def f2(t, y):
#    print(t)
    V = y[0:N]
    K_i = y[N:2*N]
    Na_i = y[2*N:3*N]
    Ca_i = y[3*N:4*N]
    Xr1 = y[4*N:5*N]
    Xr2 = y[5*N:6*N]
    Xs = y[6*N:7*N]
    m = y[7*N:8*N]
    h = y[8*N:9*N]
    j = y[9*N:10*N]
    Ca_ss = y[10*N:11*N]
    d = y[11*N:12*N]
    f = y[12*N:13*N]
    f2 = y[13*N:14*N]
    fCass = y[14*N:15*N]
    s = y[15*N:16*N]
    r = y[16*N:17*N]
    Ca_Sr = y[17*N:18*N]
    R_prime = y[18*N:]
    Xr1 = y[N:2*N]
    Xr2 = y[2*N:3*N]
    Xs = y[3*N:4*N]
    m = y[4*N:5*N]
    h = y[5*N:6*N]
    j = y[6*N:7*N]
    d = y[7*N:8*N]
    f = y[8*N:9*N]
    f2 = y[9*N:10*N]
    fCass = y[10*N:11*N]
    s = y[11*N:12*N]
    r = y[12*N:13*N]
    R_prime = y[13*N:14*N]
    Ca_i = y[14*N:15*N]
    Ca_Sr = y[15*N:16*N]
    Ca_ss = y[16*N:17*N]
    Na_i = y[17*N:18*N]
    K_i = y[18*N:]
    print(t)

    i_Stim = np.zeros((Nx, Ny, Nz))
    if t < stim_duration:
        i_Stim[:stim_x,:stim_y,:stim_z] = stim_amplitude / chi / C
    i_Stim = i_Stim.flatten()
        
    Cm = 0.185
    

    dXr1 = ((1 + np.exp((-26-V)/7))**-1 - Xr1) / ((450 / (1+np.exp((-45-V)/10)))*(6 / (1+np.exp((V+30)/11.5))))
    dXr2 = ((1 + np.exp((V+88)/24))**-1 - Xr2) / ((3 / (1+np.exp((-60-V)/20)))*(1.12 / (1+np.exp((V-60)/20))))

    dXs = ((1 + np.exp((-5-V)/14))**-1 - Xs) / ((1400 / (1+np.exp((5-V)/6))**0.5)*(1 / (1+np.exp((V-35)/15))) + 80)
    dm = ((1 + np.exp((-56.86-V)/9.03))**-2 - m) / ((1 / (1+np.exp((-60-V)/5)))*(0.1 / (1+np.exp((V+35)/5)) + 0.1 / (1+np.exp((V-50)/200))))

    tau_h = np.zeros(N)
    V0 = V[V < -40]
    V1 = V[V >=-40]
    tau_h[V < -40] += 0.057*np.exp(-(V0+80)/6.8) # alpha
    tau_h[V < -40] += 2.7*np.exp(0.079*V0) + 310000*np.exp(0.3485*V0) # beta
    tau_h[V >= -40] = 0.77 / 0.13 / (1 + np.exp((V1 + 10.66)/-11.1)) # beta
    dh = ((1 + np.exp((V+71.55)/7.43))**-2 - h) * tau_h

    tau_j = np.zeros(N)
    tau_j[V < -40] += (-25428 * np.exp(0.2444*V0) - 6.948e-6 * np.exp(-0.04391 * V0)) * (V0 + 37.78) / (1 + np.exp(0.311 * (V0 + 79.23))) # alpha
    tau_j[V < -40] += 0.02424 * np.exp(-0.01052 * V0) / (1 + np.exp(-0.1378 * (V0 + 40.14))) # beta
    tau_j[V >= -40] += 0.6 * np.exp(0.057 * V1) / (1 + np.exp(-0.1 * (V1 + 32))) # beta
    dj = ((1 + np.exp((V + 71.55) / 7.43))**-2 - j) * tau_j

    dd = ((1 + np.exp((-8-V) / 7.5))**-1 - d) / ((1.4 / (1 + np.exp((-35-V)/13)) + 0.25) * (1.4 / (1 + np.exp((V + 5)/5))) + (1 + np.exp((50 - V) / 20))**-1)

    df = ((1 + np.exp((V + 20)/7))**-1 - f) / (1102.5 * np.exp(-(V+27)**2/225) + 200 / (1 + np.exp((13 - V)/10)) + 180 / (1 + np.exp((V + 30) / 10)) + 20)

    df2 = (0.67 / (1 + np.exp((V + 35)/7)) + 0.33 - f2) / (562 * np.exp(-(V + 27)**2 / 240) + 31 / (1 + np.exp((25 - V)/10)) + 80 / (1 + np.exp((V + 30) / 10)))
    
    dfCass = (0.6 / (1 + (Ca_ss/0.05)**2) + 0.4 - fCass) / (80 / (1 + (Ca_ss / 0.05) **2) + 2)

    ds = ((1 + np.exp((V+20)/5))**-1 -s) / (85 * np.exp(-(V+45)**2/320) + 5 / (1 + np.exp((V - 20) / 5)) + 3)

    dr = ((1 + np.exp((20-V)/6))**-1 - r) / (9.5 * np.exp(-(V + 40)**2 / 1800) + 0.8)

    k2_prime = 0.045
    k4 = 0.005
    max_sr = 2.5
    min_sr = 1
    EC = 1.5
    kcasr = max_sr - ((max_sr - min_sr) / (1 + (EC / Ca_Sr)**2))
    k2 = k2_prime * kcasr
    dR_prime = -k2 * Ca_ss * R_prime + k4 * (1 - R_prime)

    Buf_c = 0.2
    K_buf_c = 0.001
    Ca_i_bufc = (1 + Buf_c*K_buf_c / (Ca_i + K_buf_c)**2)**-1
    
    V_leak = 0.00036
    i_leak = V_leak * (Ca_Sr - Ca_i)
    
    V_sr = 0.001094
    V_c = 0.016404

    Vmax_up = 0.006375
    K_up = 0.00025
    i_up = Vmax_up / (1 + K_up**2 / Ca_i**2)

    V_xfer = 0.0038
    i_xfer = V_xfer * (Ca_ss - Ca_i)

    R = 8314.472
    T = 310
    F = 96485.3415
    Ca_o = 2
    g_bca = 0.000592
    E_Ca = 0.5 * R * T / F * np.log(Ca_o / Ca_i)
    i_b_Ca = g_bca * (V - E_Ca)

    g_pCa = 0.1238
    K_pCa = 0.0005
    i_p_Ca = g_pCa * Ca_i / (Ca_i + K_pCa)

    K_NaCa = 1000
    K_sat = 0.1
    alpha = 2.5
    gamma = 0.35
    Km_Ca = 1.38
    Km_Nai = 87.5
    Na_o = 140
    i_NaCa = K_NaCa * (np.exp(gamma * V * F / (R*T)) * Na_i **3 * Ca_o - (np.exp((gamma-1)*V*F/(R*T)) * Na_o**3 * Ca_i * alpha)) / (Km_Nai**3 + Na_o**3) / (Km_Ca + Ca_o) / (1 + K_sat * np.exp((gamma - 1) * V * F / (R*T)))

    dCa_i = Ca_i_bufc * ((i_leak - i_up) * V_sr / V_c + i_xfer - (i_b_Ca + i_p_Ca - 2 * i_NaCa) * Cm / (2 * V_c *F))

    Buf_sr = 10 
    K_buf_sr = 0.3
    Ca_sr_bufsr = (1 + Buf_sr * K_buf_sr / (Ca_Sr + K_buf_sr)**2)**-1
    V_rel = 0.102
    k1_prime = 0.15
    k3 = 0.06
    k1 = k1_prime / kcasr
    O = k1 * Ca_ss**2 * R_prime / (k3 + k1 * Ca_ss**2)

    i_rel = V_rel * O * (Ca_Sr - Ca_ss) 

    dCa_Sr = Ca_sr_bufsr * (i_up - (i_rel + i_leak))
    
    Buf_ss =0.4
    K_buf_ss =0.00025
    Ca_ss_bufss = (1 + (Buf_ss * K_buf_ss) / (Ca_ss + K_buf_ss)**2)**-1
    g_CaL = 0.0000398
    
    i_CaL = g_CaL * d * f * f2 * fCass * 4 *(V-15) * F**2 / (R*T) * (0.25 * Ca_ss * np.exp(2*(V-15)*F / (R*T)) - Ca_o) / (np.exp(2*(V-15)*F/(R*T)) - 1)
    V_ss = 0.00005468

    dCa_ss = Ca_ss_bufss * (-i_CaL * Cm / (2 * V_ss*F) + i_rel * V_sr / V_ss - i_xfer *V_c / V_ss)

    g_Na = 14.838

    E_Na = R*T/F * np.log(Na_o / Na_i)
    i_Na = g_Na * m**3 * h * j * (V - E_Na)
    g_bna = 0.00029
    i_b_Na = g_bna * (V - E_Na)

    P_NaK = 2.724
    K_mk = 1
    K_mNa = 40
    i_NaK = (P_NaK * K_o) / (K_o + K_mk) * Na_i / (Na_i + K_mNa) / (1 + 0.1245 * np.exp(-0.1 * V*F/(R*T)) + 0.0353 * np.exp(-V * F / (R*T)))
    dNa_i = -(i_Na + i_b_Na + 3*i_NaK + 3*i_NaCa) * Cm / (V_c * F)
    
    E_K = R*T/F * np.log(K_o / K_i)
    alpha_K1 = 0.1 / (1 + np.exp(0.06 * (V - E_K - 200)))
    beta_K1 = (3 * np.exp(0.0002 * (V - E_K + 100)) + np.exp(0.1*(V-E_K-10))) / (1 + np.exp(-0.5 * (V-E_K)))
    g_K1 = 5.405
    i_K1 = g_K1 * (alpha_K1 / (alpha_K1 + beta_K1)) * (K_o / 5.4)**0.5 * (V - E_K)
    g_to = 0.294


    i_t0 = g_to * r * s * (V - E_K)
    g_Ks = 0.392
    P_kna = 0.03

    E_Ks = R * T / F * np.log((K_o + P_kna * Na_o) / (K_i + P_kna * Na_i))
    i_Ks = g_Ks * Xs**2 * (V-E_Ks)
    g_pK = 0.0146
    i_p_K = g_pK * (V - E_K) / (1 + np.exp((25-V)/5.98))
    i_Kr = g_Kr * (K_o/5.4)**0.5 * Xr1*Xr2 *(V-E_K)
    dK_i = -(i_K1 + i_t0 + i_Kr + i_Ks + i_p_K + i_Stim - 2*i_NaK) * Cm / (V_c * F)

    dVdt = -(i_K1 + i_t0 + i_Kr + i_Ks + i_CaL + i_NaK + i_Na + i_b_Na + i_NaCa + i_b_Ca + i_p_K + i_p_Ca + i_Stim)

    # concatenate result for output
    dydt = np.concatenate((dVdt, dXr1, dXr2, dXs, dm, dh, dj, dd, df, df2, dfCass, ds, dr, dR_prime, dCa_i, dCa_Sr, dCa_ss, dNa_i, dK_i))
    #print(dydt[0])
    return dydt

g_Kr = 0.153
K_o = 5.4
Cm = C


t0=0.
tf=40 


neg_exp = 'FE'
exp = 'SD2O3'
imp = 'RK3'

experiment = 'RD'
alpha = 'OS43(7)'
dt = 2/323
if  'DR' in experiment:
    f_list = [f1, f2]
else:
    f_list = [f2, f1]
methods = {}
if experiment[:2] == 'DR':
    methods = {(1,): imp, (2,): exp}
else:
    methods = {(1,): exp, (2,): imp}

if 'neg' in experiment:
    if alpha == 'R3':
        methods[(1,3)] = neg_exp
        methods[(2,2)] = neg_exp
    elif alpha == 'AKS3':
        methods[(1,2)] = neg_exp
        methods[(2,2)] = neg_exp
    else:
        methods[(1,3)] = neg_exp
        methods[(2,2)] = neg_exp

fname = f"N_{alpha}_{experiment}_{neg_exp}_{dt}.csv"
try:
    fs.fractional_step(f_list, dt, y0, t0, tf, alpha, methods, fname = fname, save_steps = tf//2)
except Exception as e:
    print(e)
