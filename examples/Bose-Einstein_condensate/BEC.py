import numpy as np
from numpy import exp, sqrt, cos, sin
from scipy import sparse
import scipy
import scipy.linalg as la
import os
import sys


xmax = 40.96
xmin = -40.96
ymax = 1.72
ymin = -1.72
dx = 0.04

R0 = 30  # length of condensate
omega = 1.0  # Rabi (linear) coupling
gamma = 0.01 # repulsion parameter

dt = 0.0004 # (larger timesteps will lead to instability)
t_start = 0
t_end = 1

deltat = 0.1
step = 0.1

if len(sys.argv) > 1:  # input from file
    infile = sys.argv[1]

    with open(infile, 'r') as file:
        content = file.read()
        words = content.split()
        R0 = float(words[1])
        xmax = float(words[3])
        xmin = float(words[5])
        ymax = float(words[7])
        ymin = float(words[9])
        dx = float(words[11])
        dt = float(words[13])
        t_start = float(words[15])
        t_end = float(words[17])
        deltat = float(words[19]) # total time interval saved in each file
        step = float(words[21])   # time interval between saves in each file
        omega = float(words[23])
        #k = float(words[25])
        gamma = float(words[27])

save_steps = int(deltat/step + 0.5) # num saves in each file

g1 = R0**2/2  # nonlinear couplings
g2 = R0**2/2
g12 = (1 + gamma)*g1

print(f'R0 = {R0}, gamma = {gamma}, omega = {omega}')
print(f'xmin = {xmin:.2f}, xmax = {xmax:.2f}, ymin = {ymin:.2f}, ymax = {ymax:.2f}')

dy = dx
Nx = int((xmax - xmin) / dx) + 1  # Nx = jxmax + 1
Ny = int((ymax - ymin) / dy) + 1  # Ny = jymax + 1

x = np.array([xmin + (i) * dx for i in range(Nx)])
y = np.array([ymin + (j) * dy for j in range(Ny)])


V = (1/2) * pow(x,2)  # harmonic trap potential (x-dir)
V = np.reshape(V, (Nx,1))
V = np.kron(np.ones((1,Ny)), V)
V = V.ravel()


#-----------------------------------------
def make_psi_0(interface=True, perturbation=True): 
    Vx = np.reshape((1/2) * pow(x,2), (Nx,1))
    
    # Thomas-Fermi approx:
    psi_0 = 1 - 2 * Vx / R0**2  
    for i in range(Nx):
        if psi_0[i] < 0:
            psi_0[i] = 0
    psi_0 = sqrt(psi_0)

    if interface: # interface: 
        psi_01 = psi_0 / sqrt(1 + exp( 2 * sqrt(2*g12/(g1+g2)-1) * R0 * np.reshape(x, (Nx,1)) ))
        psi_01 = np.kron(np.ones((1,Ny)), psi_01)
        psi_02 = psi_0 / sqrt(1 + exp(-2 * sqrt(2*g12/(g1+g2)-1) * R0 * np.reshape(x, (Nx,1)) ))
        psi_02 = np.kron(np.ones((1,Ny)), psi_02)
        
    else: # no interface:
        psi_01 = psi_0/sqrt(2)
        psi_02 = psi_0/sqrt(2)
        psi_01 = np.kron(np.ones((1,Ny)), psi_01)
        psi_02 = np.kron(np.ones((1,Ny)), psi_02)

    
    if perturbation:
        # initial perturbation (use linear interpolation):
        perturb = 0.1 * cos(2*np.pi*y/(ymax-ymin)) / dx;
        p = np.ceil(abs(perturb)).astype(int)
        
        for jy in range(Ny):  
            if perturb[jy] >= 0:            
                psi_01[0:Nx-p[jy], jy] = psi_01[p[jy]-1:Nx-1, jy] - (perturb[jy] - np.floor(perturb[jy])) * (psi_01[p[jy]-1:Nx-1, jy] - psi_01[p[jy]:Nx, jy])
                psi_02[0:Nx-p[jy], jy] = psi_02[p[jy]-1:Nx-1, jy] - (perturb[jy] - np.floor(perturb[jy])) * (psi_02[p[jy]-1:Nx-1, jy] - psi_02[p[jy]:Nx, jy])
            else:            
                psi_01[p[jy]:Nx, jy] = psi_01[0:Nx-p[jy], jy] + (perturb[jy] - np.floor(perturb[jy])) * (psi_01[1:Nx-p[jy]+1, jy] - psi_01[0:Nx-p[jy], jy])
                psi_02[p[jy]:Nx, jy] = psi_02[0:Nx-p[jy], jy] + (perturb[jy] - np.floor(perturb[jy])) * (psi_02[1:Nx-p[jy]+1, jy] - psi_02[0:Nx-p[jy], jy])

                
    # (flatten to 1D)
    psi_01 = psi_01.ravel()
    psi_02 = psi_02.ravel()
    psi_0 = np.concatenate((psi_01, psi_02)) # psi = psi1, psi2
    
    return psi_0
    


def get_psi_0(fname, t):
    # Function to get psi_0 (init condition) from file:
    if t < step/2:

        psi_0 = make_psi_0()

        #fname = f"droplet_dx{dx}.csv"
        #np.savetxt(fname, psi_0, delimiter=',')
        #print(f"Initial data saved to {fname}.")
            
        return psi_0

    else:
        # get psi from csv file (python data):
        base_name, ext = os.path.splitext(fname)
        count = int((t-step/2)/deltat) # t=n*deltat, count = n-1
        fname = base_name + "_" + str(count) + ".csv"
        try:
            data = np.loadtxt(fname, delimiter=',', dtype=np.complex128)
        except FileNotFoundError:
            print(f"Error: file {fname} not found. Try using t_start = 0.")
            sys.exit()

        print(f'Retrieving data from {fname}.')
        
        t_arr = np.real(data[...,0]) # different times saved in the file
        ind = int((t - count*deltat) / step + 0.5) # index of desired timestep
        
        psi = np.zeros(2*Nx*Ny) + 0j
        print(len(t_arr.shape))
        if len(t_arr.shape) == 0:
            if t - 1e-06 < t_arr < t + 1e-06:
                psi = data[1:]
            else:
                print(f"Error. t = {t}, t_arr[ind] = {t_arr}")
        else:
            if t - 1e-06 < t_arr[ind] < t + 1e-06:
                psi = data[ind, 1:]
            else:
                print(f"Error. t = {t}, t_arr[ind] = {t_arr[ind]}")
    
    return psi



#---------------------------------
def make_exp_factors_2D(dt0):
    # Function to make the 2D factors for (spectral) time propagation:
    
    N = int((Nx-1)/2)
    M = int((Ny-1)/2)
    
    exp_factor_x = np.array([exp(-1j * 2 * dt0 * (np.pi * p / (xmax-xmin))**2) for p in range(0,N)])
    exp_factor_x = np.concatenate((exp_factor_x, np.array([exp(-1j * 2 * dt0 * (np.pi * p / (xmax-xmin))**2) for p in range(-N,0)])))
    exp_factor_y = np.array([exp(-1j * 2 * dt0 * (np.pi * q / (ymax-ymin))**2) for q in range(0,M)])
    exp_factor_y = np.concatenate((exp_factor_y, np.array([exp(-1j * 2 * dt0 * (np.pi * q / (ymax-ymin))**2) for q in range(-M,0)])))
    
    exp_factor = np.kron(np.reshape(exp_factor_x, (Nx-1,1)), exp_factor_y)
    
    return exp_factor


# Get time-propagation factors:
exp_dt_2 = make_exp_factors_2D(dt/2)  # Strang uses dt and dt/2, Godunov dt
exp_dt = make_exp_factors_2D(dt)




def get_exp_factors_2D (dt0):
    if dt-1e-08 < dt0 < dt+1e-08:
        return exp_dt
    elif 0.5*dt-1e-08 < dt0 < 0.5*dt+1e-08:
        return exp_dt_2
    else:
        exp_factor = make_exp_factors_2D(dt0)
        return exp_factor



#---------------------------------------------------
def f_diff(dt,psi):  # Laplacian (diffusion) operator 
    psi = np.reshape(psi, (2,Nx,Ny))
    psi_L = 1j * np.zeros((2,Nx,Ny))
    psi_L[0,:Nx-1,:Ny-1] = np.fft.fft2(psi[0,:Nx-1,:Ny-1]) # FFT psi1
    psi_L[1,:Nx-1,:Ny-1] = np.fft.fft2(psi[1,:Nx-1,:Ny-1]) # FFT psi2

    exp_factor = get_exp_factors_2D(dt)
    psi_L[0,:Nx-1,:Ny-1] = exp_factor * psi_L[0,:Nx-1,:Ny-1]  # mult by time-propagation factors
    psi_L[1,:Nx-1,:Ny-1] = exp_factor * psi_L[1,:Nx-1,:Ny-1]
    
    psi_L[0,:Nx-1,:Ny-1] = np.fft.ifft2(psi_L[0,:Nx-1,:Ny-1]) # inverse FFT psi1
    psi_L[1,:Nx-1,:Ny-1] = np.fft.ifft2(psi_L[1,:Nx-1,:Ny-1]) # inverse FFT psi2
    
    psi_L[0,Nx-1,:] = psi_L[0,0,:]  # periodic boundaries
    psi_L[0,:,Ny-1] = psi_L[0,:,0]
    psi_L[1,Nx-1,:] = psi_L[1,0,:]
    psi_L[1,:,Ny-1] = psi_L[1,:,0]
    
    psi_L = psi_L.ravel()
    return psi_L    

def f_nl(dt,psi):  # non-linear coupling
    psi_nl = 1j * np.zeros(2*Nx*Ny)
    psi_nl[:Nx*Ny] = psi[:Nx*Ny] * exp(-1j * (V + (g1*(abs(psi[:Nx*Ny]))**2 + g12*(abs(psi[Nx*Ny:]))**2)) * dt)
    psi_nl[Nx*Ny:] = psi[Nx*Ny:] * exp(-1j * (V + (g12*(abs(psi[:Nx*Ny]))**2 + g2*(abs(psi[Nx*Ny:]))**2)) * dt)
    
    return psi_nl

def f_Rabi(dt,psi):  # Rabi (linear) coupling
    psi_Rabi = 1j * np.zeros(2*Nx*Ny)
    psi_Rabi[:Nx*Ny] = cos(omega*dt) * psi[:Nx*Ny]  +  1j * sin(omega*dt) * psi[Nx*Ny:]
    psi_Rabi[Nx*Ny:] = 1j * sin(omega*dt) * psi[:Nx*Ny]  +  cos(omega*dt) * psi[Nx*Ny:]
    
    return psi_Rabi


#-------------------------------------------------------------------
def get_files(fname, alpha, methods, operators, dt0, ivp_methods={}):
    
    base_name, ext = os.path.splitext(fname)
    psi = get_psi_0(fname, t_start)
    psi_0 = psi
    
    t0 = t_start
    tf = t_start + deltat
    count = int((t_start+step/2)/deltat)
    elapsed_time = 0

    t = t0
    count = int(t0/deltat + 0.5)
    dtcount = int(t0/dt + 0.5)
    
    print(f"\nElapsed time: {elapsed_time}")
    
    while tf <= t_end + 1e-06:
        start_time = time.perf_counter() # record time
        
        result = fs.fractional_step(operators, dt0, psi_0, t0, tf, alpha, methods, fname=fname, save_steps=save_steps, ivp_methods=ivp_methods)
        
        psi_0 = result

        end_time = time.perf_counter() # record time
        elapsed_time += end_time - start_time
        print(f't = {tf}, elapsed time = {elapsed_time}', flush=True)

        fname_mod = base_name + "_" + str(count) + ext  # save to file
        psi_t = np.insert(psi_0, 0, tf) # insert t at start
        np.savetxt(fname_mod, psi_t, delimiter=',')
        
        t0 = tf
        tf += deltat
        count += 1

    
    return elapsed_time, fname_mod


    

#---------------------- Use pythOS: ----------------------------
import fractional_step as fs
import time

if R0 == 30:
    Rstr = ''
else:
    Rstr = f'_R{int(R0)}'

fname = f"BEC" + Rstr + f"_omega{omega}/BEC_omega{omega}.csv"
print(f'Outputting to file {fname}')

alpha =  "Strang"

methods = {(0,): "ANALYTIC", (1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"}

operators = [f_diff, f_nl, f_Rabi]

#result = fs.fractional_step([f_diff, f_nl, f_Rabi], dt, psi_0, t0, tf, alpha, methods, fname=fname, save_steps=save_steps)

elapsed_time, fname_mod = get_files(fname, alpha, methods, operators, dt)
