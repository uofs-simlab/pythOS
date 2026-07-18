import numpy as np
from numpy import exp, sqrt, cos, sin, cosh
from scipy import sparse
import scipy.linalg as la
import os
import sys

taskID = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))

dim = 2   # 2D
R0 = 30   # length of condensate in trap, also controls g
omega = 1  # Rabi (linear) coupling
gamma = 0.01  # repulsion parameter
xmax = 40.96
xmin = -40.96
ymax = 1.72
ymin = -1.72
dx = 0.04
dy = dx
dt = 0.0004
t_start = 0
t_end = 1
deltat = 0.1  # time interval saved in each file
step = 0.1    # time step between saves


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
    dy = dx
    dt = float(words[13])
    t_start = float(words[15])
    t_end = float(words[17])
    deltat = float(words[19])
    step = float(words[21])
    omega = float(words[23])
    gamma = float(words[25])
    dim = int(words[27])


g1 = R0**2/2  # nonlinear couplings
g2 = R0**2/2
g12 = (1+gamma)*g1

save_steps = int(deltat/step + 0.5) # num saves in each file

multi = False  # for finding error vs total time, multiple dt values
if t_end - t_start - 1e-06 < deltat < t_end - t_start + 1e-06:
    multi = True

directory = f"BEC_omega{omega}/"
directory_timing = directory + f"dx_{dx}_dt_{dt}/"
directory_data = directory + f"dx_{dx}/"
os.makedirs(directory_timing, exist_ok=True) # make folders if don't exist
os.makedirs(directory_data, exist_ok=True)
timing_file = directory_timing + f"timing_{taskID}.txt" # detailed timings
data_file = f"CPU_vs_error_omega{omega}_dx{dx}_t{t_end}.txt" # error vs total time

if dim == 1:
    ymax = 0
    ymin = 0
    g1 = -abs(g1) # attractive condensate
    g2 = -abs(g2)
    g12 = -abs(g12)
elif dim != 2:
    print("Warning: system should be 1D or 2D. Setting to 2D.")
    dim = 2

print(f"{dim}D system:")
print(f"R0 = {R0}, omega = {omega}, gamma = {gamma}")
print(f"dx = {dx}")
print(f"dt = {dt}")
print(f"t_start = {t_start}, t_end = {t_end}")
#-----------------------------------------------------

Nx = int((xmax - xmin) / dx) + 1  # Nx = jxmax + 1
Ny = int((ymax - ymin) / dy) + 1  # Ny = jymax + 1

x = np.array([xmin + (i) * dx for i in range(Nx)])
y = np.array([ymin + (j) * dy for j in range(Ny)])

V = (1/2) * pow(x,2)  # harmonic trap potential (x-dir)
if dim == 2:
    V = np.kron(np.ones((1,Ny)), np.reshape(V, (Nx,1)))
    V = V.ravel()



def make_psi_0(interface=True, perturbation=True, bright_soliton=False): 
    Vx = (1/2) * pow(x,2)
    
    # Thomas-Fermi approx:
    psi_0 = 1 - 2 * Vx / R0**2  
    for i in range(Nx):
        if psi_0[i] < 0:
            psi_0[i] = 0
    psi_0 = sqrt(psi_0)
    
    if interface: # interface:
        psi_01 = psi_0 / sqrt(1 + exp( 2 * sqrt(2*g12/(g1+g2)-1) * R0 * x ))
        psi_02 = psi_0 / sqrt(1 + exp(-2 * sqrt(2*g12/(g1+g2)-1) * R0 * x ))
    else: # no interface:
        if bright_soliton:
            psi_01 = psi_0 * exp(5*1j*x) / cosh(x/2) # bright soliton
            psi_02 = psi_0 * exp(-5*1j*x) / cosh(x/2)
        else:
            psi_01 = psi_0/sqrt(2)
            psi_02 = psi_0/sqrt(2)

    if dim == 2:
        psi_01 = np.kron(np.ones((1,Ny)), np.reshape(psi_01, (Nx,1))) 
        psi_02 = np.kron(np.ones((1,Ny)), np.reshape(psi_02, (Nx,1)))
    

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



def get_psi_0(fname, t, dt0):
    
    if t < step/2:

        if dim == 1:
            psi_0 = make_psi_0(False,False,True)
        elif dim == 2:
            psi_0 = make_psi_0()
    
        # fname = f"droplet_dx{dx}.csv"
        # np.savetxt(fname, psi_0, delimiter=',')
        # print(f"Initial data saved to {fname}.")
            
        return psi_0
    
    else:
        # Get psi from csv file (python data):
        base_name, ext = os.path.splitext(fname)
        count = int((t-step/2)/deltat) # t=n*deltat, count = n-1

        if multi:

            fname_mod = base_name + "_t" + str(t_start) + "_dt" + str(dt0) + ext
            try:
                data = np.loadtxt(fname_mod, delimiter=',', dtype=np.complex128)
            except FileNotFoundError:
                print(f"Error: file {fname_mod} not found.")
                sys.exit()
        else:
            
            fname_mod = base_name + "_" + str(count) + ".csv"
            try:
                data = np.loadtxt(fname_mod, delimiter=',', dtype=np.complex128)
            except FileNotFoundError:
                print(f"Error: file {fname_mod} not found.")
                sys.exit()
                
        t_arr = np.real(data[...,0]) # different times saved in file
        ind = int((t - count*deltat) / step + 0.5) # index of desired time
            
        psi = np.zeros(2*Nx*Ny) + 0j
            
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



# -----------------------------------------------------
def make_exp_factors(dt0):
    # Function to make the full-dimensional factors for (spectral) time propagation:
    
    N = int((Nx-1)/2) # N = jxmax/2
    M = int((Ny-1)/2) # M = jymax/2
    
    exp_factor_x = np.array([exp(-1j * 2 * dt0 * (np.pi * p / (xmax-xmin))**2) for p in range(0,N)])
    exp_factor_x = np.concatenate((exp_factor_x, np.array([exp(-1j * 2 * dt0 * (np.pi * p / (xmax-xmin))**2) for p in range(-N,0)])))
    
    if dim == 2:
        exp_factor_y = np.array([exp(-1j * 2 * dt0 * (np.pi * q / (ymax-ymin))**2) for q in range(0,M)])
        exp_factor_y = np.concatenate((exp_factor_y, np.array([exp(-1j * 2 * dt0 * (np.pi * q / (ymax-ymin))**2) for q in range(-M,0)])))
    
        exp_factor = np.kron(np.reshape(exp_factor_x, (Nx-1,1)), exp_factor_y)

    elif dim == 1:
        exp_factor = np.reshape(exp_factor_x, (Nx-1,1))

    return exp_factor


def make_exp_factors_xy(dt0):
    # Function to make the x- and y-factors for (spectral) time propagation:
    
    N = int((Nx-1)/2) # N = jxmax/2
    M = int((Ny-1)/2) # M = jymax/2
    
    exp_factor_x = np.array([exp(-1j * 2 * dt0 * (np.pi * p / (xmax-xmin))**2) for p in range(0,N)])
    exp_factor_x = np.concatenate((exp_factor_x, np.array([exp(-1j * 2 * dt0 * (np.pi * p / (xmax-xmin))**2) for p in range(-N,0)])))
    exp_factor_y = np.array([exp(-1j * 2 * dt0 * (np.pi * q / (ymax-ymin))**2) for q in range(0,M)])
    exp_factor_y = np.concatenate((exp_factor_y, np.array([exp(-1j * 2 * dt0 * (np.pi * q / (ymax-ymin))**2) for q in range(-M,0)])))

    if dim == 2:
        exp_factor_x = np.kron( np.reshape(exp_factor_x, (Nx-1,1)), np.ones(Ny) )
        exp_factor_y = np.kron( np.ones((Nx,1)), exp_factor_y )
    
    return exp_factor_x, exp_factor_y


# Get time-propagation factors:
exp_dt_2 = make_exp_factors(dt/2)  # Strang
exp_x_dt_2, exp_y_dt_2 = make_exp_factors_xy(dt/2)
exp_dt = make_exp_factors(dt)
exp_x_dt, exp_y_dt = make_exp_factors_xy(dt)

exp_c1 = make_exp_factors((1-1j)/2*dt) # complex Lie Trotter
exp_c2 = make_exp_factors((1+1j)/2*dt)

exp_7_24 = make_exp_factors((7/24)*dt)  # Ruth method
exp_n1_24 = make_exp_factors(-(1/24)*dt)
exp_2_3 = make_exp_factors((2/3)*dt)
exp_n2_3 = make_exp_factors(-(2/3)*dt)
exp_3_4 = make_exp_factors((3/4)*dt)

AKS3_1 = 0.268330095673069  # AKS3 method
AKS3_2 = -0.187991620228223
AKS3_3 = 0.919661524555154
exp_AKS3_1 = make_exp_factors(AKS3_1*dt)
exp_AKS3_2 = make_exp_factors(AKS3_2*dt)
exp_AKS3_3 = make_exp_factors(AKS3_3*dt)

OS2_1 = 0.2148701498521859  # OS2(4,3)[7] method
OS2_2 = 0.5114860522253670
OS2_3 = 0.6686906878883930
OS2_4 = -0.5014273889798119
OS2_5 = -0.0419569080414939
OS2_6 = 0.9899413367544449
OS2_7 = 0.1583960703009149
exp_OS2_1 = make_exp_factors(OS2_1)
exp_OS2_2 = make_exp_factors(OS2_2)
exp_OS2_3 = make_exp_factors(OS2_3)
exp_OS2_4 = make_exp_factors(OS2_4)
exp_OS2_5 = make_exp_factors(OS2_5)
exp_OS2_6 = make_exp_factors(OS2_6)
exp_OS2_7 = make_exp_factors(OS2_7)




def get_exp_factors_xy (dt0):
    if dt-1e-08 < dt0 < dt+1e-08: # dt
        return exp_x_dt, exp_y_dt
    elif 0.5*dt-1e-08 < dt0 < 0.5*dt+1e-08: # dt/2
        return exp_x_dt_2, exp_y_dt_2
    else:
        exp_factor_x, exp_factor_y = make_exp_factors_xy(dt0)
        return exp_factor_x, exp_factor_y


def get_exp_factors (dt0):
    if np.abs(dt0.imag) < 1e-08:
        if dt-1e-08 < dt0 < dt+1e-08: # dt
            return exp_dt
        elif 0.5*dt-1e-08 < dt0 < 0.5*dt+1e-08: # dt/2
            return exp_dt_2
        elif (7/24)*dt-1e-08 < dt0 < (7/24)*dt+1e-08: # 7/24*dt (Ruth)
            return exp_7_24
        elif -(1/24)*dt-1e-08 < dt0 < -(1/24)*dt+1e-08: # -1/24*dt
            return exp_n1_24
        elif (2/3)*dt-1e-08 < dt0 < (2/3)*dt+1e-08: # 2/3*dt
            return exp_2_3
        elif -(2/3)*dt-1e-08 < dt0 < -(2/3)*dt+1e-08: # -2/3*dt
            return exp_n2_3
        elif (3/4)*dt-1e-08 < dt0 < (3/4)*dt+1e-08: # 3/4*dt
            return exp_3_4
        elif AKS3_1*dt-1e-08 < dt0 < AKS3_1*dt+1e-08: # (AKS3)
            return exp_AKS3_1
        elif AKS3_2*dt-1e-08 < dt0 < AKS3_2*dt+1e-08:
            return exp_AKS3_2
        elif AKS3_3*dt-1e-08 < dt0 < AKS3_3*dt+1e-08:
            return exp_AKS3_3
        elif OS2_1*dt-1e-08 < dt0 < OS2_1*dt+1e-08: # (OS2(3,4)[7])
            return exp_OS2_1
        elif OS2_2*dt-1e-08 < dt0 < OS2_2*dt+1e-08:
            return exp_OS2_2
        elif OS2_3*dt-1e-08 < dt0 < OS2_3*dt+1e-08:
            return exp_OS2_3
        elif OS2_4*dt-1e-08 < dt0 < OS2_4*dt+1e-08:
            return exp_OS2_4
        elif OS2_5*dt-1e-08 < dt0 < OS2_5*dt+1e-08:
            return exp_OS2_5
        elif OS2_6*dt-1e-08 < dt0 < OS2_6*dt+1e-08:
            return exp_OS2_6
        elif OS2_7*dt-1e-08 < dt0 < OS2_7*dt+1e-08:
            return exp_OS2_7
        else:
            exp_factor = make_exp_factors(dt0)
            return exp_factor
    else:
        if -0.5*dt-1e-08 < dt0.imag < -0.5*dt+1e-08: # (1-i)/2*dt
            return exp_c1
        elif 0.5*dt-1e-08 < dt0.imag < 0.5*dt+1e-08: # (1+i)/2*dt
            return exp_c2
        else:
            exp_factor = make_exp_factors(dt0)
            return exp_factor



# --------------------------------------------------
def f_diff(dt0, psi):  # Laplacian (diffusion) operator
    psi = np.reshape(psi, (2,Nx,Ny))
    psi_L = 1j * np.zeros((2,Nx,Ny))
    ny = max(1,Ny-1)
    psi_L[0,:-1,:ny] = np.fft.fft2(psi[0,:-1,:ny]) # FFT psi1
    psi_L[1,:-1,:ny] = np.fft.fft2(psi[1,:-1,:ny]) # FFT psi2
    
    exp_factor = get_exp_factors(dt0)  # mult by time-propagation factors
    
    psi_L[0,:Nx-1,:ny] = exp_factor * psi_L[0,:Nx-1,:ny]
    psi_L[1,:Nx-1,:ny] = exp_factor * psi_L[1,:Nx-1,:ny]
    
    psi_L[0,:Nx-1,:ny] = np.fft.ifft2(psi_L[0,:Nx-1,:ny]) # inverse FFT psi1
    psi_L[1,:Nx-1,:ny] = np.fft.ifft2(psi_L[1,:Nx-1,:ny]) # inverse FFT psi2
    
    psi_L[0,Nx-1,:] = psi_L[0,0,:]  # periodic boundaries
    psi_L[0,:,Ny-1] = psi_L[0,:,0]
    psi_L[1,Nx-1,:] = psi_L[1,0,:]
    psi_L[1,:,Ny-1] = psi_L[1,:,0]
    
    psi_L = psi_L.ravel()
    return psi_L  


def f_xx(dt0, psi):  # x-part of Laplacian (diffusion)
    psi = np.reshape(psi, (2,Nx,Ny))
    psi_xx = 1j * np.zeros((2,Nx,Ny))
    psi_xx[0,:-1,:] = np.fft.fft(psi[0,:-1,:], axis=0) # FFT psi1
    psi_xx[1,:-1,:] = np.fft.fft(psi[1,:-1,:], axis=0) # FFT psi2
    
    exp_factor_x, _ = get_exp_factors_xy(dt0) # mult by time-propagation factors
    psi_xx[0,:-1,:] = exp_factor_x * psi_xx[0,:-1,:]
    psi_xx[1,:-1,:] = exp_factor_x * psi_xx[1,:-1,:]
    
    psi_xx[0,:-1,:] = np.fft.ifft(psi_xx[0,:-1,:], axis=0) # inverse FFT psi1
    psi_xx[1,:-1,:] = np.fft.ifft(psi_xx[1,:-1,:], axis=0) # inverse FFT psi2
    
    psi_xx[0,-1,:] = psi_xx[0,0,:]  # periodic boundaries
    psi_xx[1,-1,:] = psi_xx[1,0,:]
    
    psi_xx = psi_xx.ravel()
    return psi_xx  


def f_yy(dt0, psi):  # y-part of Laplacian (diffusion)
    psi = np.reshape(psi, (2,Nx,Ny))
    psi_yy = 1j * np.zeros((2,Nx,Ny))
    psi_yy[0,:,:-1] = np.fft.fft(psi[0,:,:-1], axis=1) # FFT psi1
    psi_yy[1,:,:-1] = np.fft.fft(psi[1,:,:-1], axis=1) # FFT psi2
    
    _, exp_factor_y = get_exp_factors_xy(dt0) # mult by time-propagation factors
    psi_yy[0,:,:-1] = exp_factor_y * psi_yy[0,:,:-1]
    psi_yy[1,:,:-1] = exp_factor_y * psi_yy[1,:,:-1]
    
    psi_yy[0,:,:-1] = np.fft.ifft(psi_yy[0,:,:-1], axis=1) # inverse FFT psi1
    psi_yy[1,:,:-1] = np.fft.ifft(psi_yy[1,:,:-1], axis=1) # inverse FFT psi2
    
    psi_yy[0,:,-1] = psi_yy[0,:,0]  # periodic boundaries
    psi_yy[1,:,-1] = psi_yy[1,:,0]
    
    psi_yy = psi_yy.ravel()
    return psi_yy   


def f_diffR(dt0, psi):  # Laplacian (diffusion) AND Rabi coupling
    psi = np.reshape(psi, (2,Nx,Ny))
    psi_diffR = 1j * np.zeros((2,Nx,Ny))
    ny = max(1,Ny-1)
    psi_diffR[0,:Nx-1,:ny] = np.fft.fft2(psi[0,:Nx-1,:ny]) # FFT psi1
    psi_diffR[1,:Nx-1,:ny] = np.fft.fft2(psi[1,:Nx-1,:ny]) # FFT psi2

    exp_factor = get_exp_factors(dt0) # checked Apr-22, 2026
     
    c1 = 0.5 * (psi_diffR[0,:-1,:ny] - psi_diffR[1,:-1,:ny]) # find coeff's
    c2 = 0.5 * (psi_diffR[0,:-1,:ny] + psi_diffR[1,:-1,:ny])
    
    c1 = c1 * exp_factor
    c2 = c2 * exp_factor
    psi_diffR[0,:-1,:ny] = c1 * exp(-1j*omega*dt0)  +  c2 * exp(1j*omega*dt0)
    psi_diffR[1,:-1,:ny] = -c1 * exp(-1j*omega*dt0)  +  c2 * exp(1j*omega*dt0)

    
    psi_diffR[0,:-1,:ny] = np.fft.ifft2(psi_diffR[0,:-1,:ny]) # inverse FFT psi1
    psi_diffR[1,:-1,:ny] = np.fft.ifft2(psi_diffR[1,:-1,:ny]) # inverse FFT psi2
    
    psi_diffR[0,Nx-1,:] = psi_diffR[0,0,:]  # periodic boundaries
    psi_diffR[0,:,Ny-1] = psi_diffR[0,:,0]
    psi_diffR[1,Nx-1,:] = psi_diffR[1,0,:]
    psi_diffR[1,:,Ny-1] = psi_diffR[1,:,0]
    
    psi_diffR = psi_diffR.ravel()
    return psi_diffR  


def f_nl(dt0, psi):  # non-linear operator
    psi_nl = 1j * np.zeros(2*Nx*Ny)
    psi_nl[:Nx*Ny] = psi[:Nx*Ny] * exp(-1j * (V + (g1*(abs(psi[:Nx*Ny]))**2 + g12*(abs(psi[Nx*Ny:]))**2)) * dt0)
    psi_nl[Nx*Ny:] = psi[Nx*Ny:] * exp(-1j * (V + (g12*(abs(psi[:Nx*Ny]))**2 + g2*(abs(psi[Nx*Ny:]))**2)) * dt0)
    
    return psi_nl


def f_Rabi(dt0, psi):  # Rabi (linear) coupling
    psi_Rabi = 1j * np.zeros(2*Nx*Ny)
    
    psi_Rabi[:Nx*Ny] = cos(omega*dt0) * psi[:Nx*Ny]  +  1j * sin(omega*dt0) * psi[Nx*Ny:]
    psi_Rabi[Nx*Ny:] = 1j * sin(omega*dt0) * psi[:Nx*Ny]  +  cos(omega*dt0) * psi[Nx*Ny:]
    
    return psi_Rabi


def laplacian_2D(Nx, Ny, dx=1.0, dy=1.0):
    """
    Create the Laplacian Matrix for Nx x Ny grid, using central differences
    Returns a sparse matrix of size (Nx*Ny x Nx*Ny).
    """
    # Central difference in x:
    firstcol = np.zeros(Nx)
    firstcol[0] = -2
    firstcol[1] = 1
    Dxx = la.toeplitz(firstcol)
    #Dxx[0,1] = 2  # Neumann BC
    #Dxx[-1,-2] = 2
    Dxx[0,-1] = 1  # PBC
    Dxx[-1,0] = 1
    Dxx = 1/(dx**2) * sparse.csr_matrix(Dxx)
    
    # Central difference in y:
    firstcol = np.zeros(Ny)
    firstcol[0] = -2
    firstcol[1] = 1
    Dyy = la.toeplitz(firstcol)
    # Dyy[0,1] = 2  # Neumann BC
    # Dyy[-1,-2] = 2
    Dyy[0,-1] = 1  # PBC
    Dyy[-1,0] = 1
    Dyy = 1/(dy**2) * sparse.csr_matrix(Dyy)
    
    # 2D Laplacian using Kronecker product
    L = sparse.kron(Dxx, sparse.identity(Ny)) + sparse.kron(sparse.identity(Nx), Dyy)
    
    return L
  

# L = laplacian_2D(Nx,Ny,dx,dy)


def f_L(t,psi):
    psi_L = 1j * np.zeros(2*Nx*Ny)
    psi_L[:Nx*Ny] = (L @ psi[:Nx*Ny])
    psi_L[Nx*Ny:] = (L @ psi[Nx*Ny:])
    
    return (1j/2) * psi_L    

def f_nl_dpsi(t,psi):
    dpsi_nl = 1j * np.zeros(2*Nx*Ny)
    dpsi_nl[:Nx*Ny] = psi[:Nx*Ny] * (V + g1*(abs(psi[:Nx*Ny]))**2 + g12*(abs(psi[Nx*Ny:]))**2)
    dpsi_nl[Nx*Ny:] = psi[Nx*Ny:] * (V + g12*(abs(psi[:Nx*Ny]))**2 + g2*(abs(psi[Nx*Ny:]))**2)
    
    return -1j * dpsi_nl

def f_Rabi_dpsi(t,psi):
    dpsi_Rabi = 1j * np.zeros(2*Nx*Ny)
    dpsi_Rabi[:Nx*Ny] = -omega * psi[Nx*Ny:]
    dpsi_Rabi[Nx*Ny:] = -omega * psi[:Nx*Ny]
    
    return -1j * dpsi_Rabi

def f_full(t,psi):
    dpsi = f_L(t,psi) + f_nl_dpsi(t,psi) + f_Rabi_dpsi(t,psi)

    return dpsi


# -------------------------------------------------

def get_pdata(fname, dt0):
    # Function to get all data from csv file (python data)
    base_name, ext = os.path.splitext(fname)
    fname_mod = base_name + "_t" + str(t_end) + "_dt" + str(dt0) + ext
    
    try:
        data = np.loadtxt(fname_mod, delimiter=',', dtype=np.complex128)
    except FileNotFoundError:
        return None
    
    return data


def get_ref(t):
    # Function to get csv file with reference solution (python data)
    if 0.04 - 1e-06 < dx < 0.04 + 1e-06:
        folder = directory + "dx_0.04_dt_2.5e-05/"
        #step0 = 0.01 # check these...
        deltat0 = 0.1
    elif 0.016 - 1e-06 < dx < 0.016 + 1e-06:
        folder = directory + "dx_0.016_dt_5e-06/"
        #step0 = 0.1
        deltat0 = 0.1
    elif 0.2 - 1e-06 < dx < 0.2 + 1e-06:
        folder = directory + "dx_0.2_dt_0.0001/"
        #step0 = 0.1
        deltat0 = 0.1
    elif 0.1 - 1e-06 < dx < 0.1 + 1e-06:
        folder = directory + "dx_0.1_dt_0.0001/"
        #step0 = 1
        deltat0 = 1
    else:
        print(f"Error: reference solution not specified for dx = {dx}")
        
    base_name = folder + "3-split/BEC_3-split_large-diff" + my_str
    if 0.1 - 1e-06 < dx < 0.1 + 1e-06 and 1-1e-06 < omega < 1+1e-06:
        base_name = folder + "2-split/BEC_2-split_R3_diffR_nl"
    count = int((t-1e-08)/deltat0) # if t=n*delta, count = n-1
    fname = base_name + "_" + str(count) + ".csv"
    
    try:
        refdata = np.loadtxt(fname, delimiter=',', dtype=np.complex128)
        if 0.04 - 1e-06 < dx < 0.04 + 1e-06 and 1-1e-06 < omega < 1+1e-06:
            refdata = refdata[10,:]
    except FileNotFoundError:
        return None
    
    return refdata



def get_timing_files(fname, alpha, methods, operators, dt0, ivp_methods={}):
    
    base_name, ext = os.path.splitext(fname)
    psi_0 = get_psi_0(fname, t_start, dt0)
    
    t0 = t_start
    tf = t_start + deltat
    count = int((t_start+step/2)/deltat)
    elapsed_time = 0

    if not multi: # if not running multi
        with open(timing_file, "a") as file:
            file.write("\n" + str(elapsed_time))
    print(f"\nElapsed time: {elapsed_time}")
    
    while tf <= t_end + 1e-06:
        start_time = time.perf_counter() # record time
        
        result = fs.fractional_step(operators, dt0, psi_0, t0, tf, alpha, methods, fname=fname, save_steps=save_steps, ivp_methods=ivp_methods)
        
        end_time = time.perf_counter() # record time
        elapsed_time += end_time - start_time
        if not t_end - 1e-06 < deltat < t_end + 1e-06:
            with open(timing_file, "a") as file:
                file.write("\n" + str(elapsed_time))
        print(elapsed_time)
        
        psi_0 = result

        if multi:
            fname_mod = base_name + "_t" + str(t_end) + "_dt" + str(dt0) + ext
        else:
            fname_mod = base_name + "_" + str(count) + ext  # rename file
        os.replace(fname, fname_mod) # modify file name
        if save_steps == 1:
            data = np.loadtxt(fname_mod, delimiter=',', dtype=np.complex128)
            data = data[1:,...] # remove initial data (half of file)
            np.savetxt(fname_mod, data, delimiter=',')
        
        t0 = tf
        tf += deltat
        count += 1

    return elapsed_time



def get_timing_files_multi(fname, alpha, methods, operators, dt0, ivp_methods={}):
    
    for i in range(len(dt0)): # get timings and files
        print(f"\ndt = {dt0[i]}")
        elapsed_time = get_timing_files(fname, alpha, methods, operators, dt0[i], ivp_methods)
        with open(data_file, "a") as file: # save total elapsed time to file
            file.write(f"{elapsed_time:.2f}" + "  ")

    with open(data_file, "a") as file: # linebreak
        file.write("\n")

    for i in range(len(dt0)): # get total errors
        error = get_error(fname, dt0[i])
        with open(data_file, "a") as file: # save errors to file
            file.write(f"{error:.2e}" + "  ")

    with open(data_file, "a") as file: # linebreak
        file.write("\n")



def get_error(fname, dt0):
    # Function to get total error at final time:
    base_name, ext = os.path.splitext(fname)

    ref = get_ref(t_end) 
    pdata = get_pdata(fname, dt0)

    refnorm = sum(abs(ref[1:])**2)
    error = np.sqrt(sum(abs(ref[1:] - pdata[1:])**2) / refnorm)

    return error





# --------------------- Use pythOS: ----------------------- 
import fractional_step as fs
import time


if 0.1 - 1e-06 < dx < 0.1 + 1e-06: # dx = 0.1
    
    dt_Strang = [0.0002, 0.0006, 0.001, 0.0014, 0.0018, 0.002]
    dt_Ruth = [0.001, 0.0014, 0.0018, 0.002, 0.0022]

elif 0.04 - 1e-06 < dx < 0.04 + 1e-06: # dx = 0.04

    dt_Strang = [0.0001, 0.0002, 0.0004, 0.00044, 0.00046]
    dt_Ruth = [0.0002, 0.0004, 0.00045, 0.00047]
    dt_OS2 = [0.0002, 0.00036, 0.00038, 0.0004]

    if dim == 1:
        dt_Strang = [0.001, 0.005, 0.01, 0.02, 0.04, 0.08]
        dt_Ruth = dt_Strang
        dt_OS2 = dt_Strang

my_str = ''
if dim == 1:
    my_str += "_1D"


    
# # -------- 1-splitting: -------- (SUNDIALS - adaptive)
# if multi:
#     fname = directory_data + f"BEC_1-split_Sundials" + my_str + ".csv"
# else:
#     os.makedirs(directory_data + "1-split", exist_ok=True)
#     fname = directory_data + f"1-split/BEC_Sundials_CV_BDF" + my_str + ".csv"
# base_name, ext = os.path.splitext(fname)
# print("\nTesting 1-splitting: SUNDIALS CV_BDF")
# print(f"Filename: {fname}")

# alpha = [[1]]

# methods = {(1,): "ADAPTIVE"}
# #ivp_methods = {(1,): ("CV_ADAMS", 1e-03, 1e-02)}
# ivp_methods = {(1,): ("CV_BDF", 1e-03, 1e-04)}

# operators = [f_full]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"Sundials" + "  ")
#     ivp_methods = {(1,): ("CV_BDF", 1e-01, 1e-02)}
#     get_timing_files_multi(fname, alpha, methods, operators, [dt], ivp_methods)
#     ivp_methods = {(1,): ("CV_BDF", 1e-02, 1e-03)}
#     get_timing_files_multi(fname, alpha, methods, operators, [dt], ivp_methods)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt, ivp_methods)



# ----------------- 1st-order ---------------------#
    
# -------- 2-splitting: -------- (Lie Trotter)
if multi:
    fname = directory_data + f"BEC_2-split_LT" + my_str + ".csv"
else:
    os.makedirs(directory_timing + "2-split", exist_ok=True)
    fname = directory_timing + f"2-split/BEC_2-split_LT" + my_str + ".csv"
print("\nTesting 2-splitting: Lie Trotter method")
print(f"Filename: {fname}")

alpha = [[1,1]] # 'Godunov'

methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

operators = [f_diffR, f_nl]

if multi:
    with open(data_file, "a") as file: # write method name
        file.write(f"LT" + "  ")
    get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
else:
    get_timing_files(fname, alpha, methods, operators, dt)




#----------------- 2nd-order ---------------------#

# -------- 2-splitting: -------- (f_nl middle step)
if multi:
    fname = directory_data + f"BEC_2-split_Strang-nl" + my_str + ".csv"
else:
    os.makedirs(directory_timing + "2-split", exist_ok=True)
    fname = directory_timing + f"2-split/BEC_2-split_large-nl" + my_str + ".csv"
print("\nTesting 2-splitting: large step for nonlinear part")
print(f"Filename: {fname}")

alpha = "Strang-N" 

methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

operators = [f_diffR, f_nl]

if multi:
    with open(data_file, "a") as file: # write method name
        file.write(f"Str-2-nl" + "  ")
    get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
else:
    get_timing_files(fname, alpha, methods, operators, dt)




# -------- 2-splitting: -------- (f_diffR middle step)
if multi:
    fname = directory_data + f"BEC_2-split_Strang-diffR" + my_str + ".csv"
else:
    os.makedirs(directory_timing + "2-split", exist_ok=True)
    fname = directory_timing + f"2-split/BEC_2-split_large-diffR" + my_str + ".csv"
print("\nTesting 2-splitting: large step for diffusion (Laplacian) with Rabi coupling")
print(f"Filename: {fname}")

alpha = "Strang-N" 

methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

operators = [f_nl, f_diffR]

if multi:
    with open(data_file, "a") as file: # write method name
        file.write(f"Str-2-diffR" + "  ")
    get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
else:
    get_timing_files(fname, alpha, methods, operators, dt)




# -------- 3-splitting: -------- (f_Rabi middle step)
if multi:
    fname = directory_data + f"BEC_3-split_Strang-Rabi" + my_str + ".csv"
else:
    os.makedirs(directory_timing + "3-split", exist_ok=True)
    fname = directory_timing + f"3-split/BEC_3-split_large-Rabi" + my_str + ".csv"
print("\nTesting 3-splitting: large step for Rabi coupling")
print(f"Filename: {fname}")

alpha = "Strang-N" 

methods = {(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"}

operators = [f_diff, f_nl, f_Rabi]

if multi:
    with open(data_file, "a") as file: # write method name
        file.write(f"Str-3-R" + "  ")
    get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
else:
    get_timing_files(fname, alpha, methods, operators, dt)




# -------- 3-splitting: -------- (f_nl middle step)
if multi:
    fname = directory_data + f"BEC_3-split_Strang-nl" + my_str + ".csv"
else:
    os.makedirs(directory_timing + "3-split", exist_ok=True)
    fname = directory_timing + f"3-split/BEC_3-split_large-nl" + my_str + ".csv"
print("\nTesting 3-splitting: large step for nonlinear part")
print(f"Filename: {fname}")

alpha = "Strang-N" 

methods = {(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"}

operators = [f_diff, f_Rabi, f_nl]

if multi:
    with open(data_file, "a") as file: # write method name
        file.write(f"Str-3-nl" + "  ")
    get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
else:
    get_timing_files(fname, alpha, methods, operators, dt)




# -------- 3-splitting: -------- (f_diff middle step)
if multi:
    fname = directory_data + f"BEC_3-split_Strang-diff" + my_str + ".csv"
else:
    os.makedirs(directory_timing + "3-split", exist_ok=True)
    fname = directory_timing + f"3-split/BEC_3-split_large-diff" + my_str + ".csv"
print("\nTesting 3-splitting: large step for diffusion (Laplacian)")
print(f"Filename: {fname}")

alpha = "Strang-N" 

methods = {(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC"}

operators = [f_Rabi, f_nl, f_diff]

if multi:
    with open(data_file, "a") as file: # write method name
        file.write(f"Str-3-diff" + "  ")
    get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
else:
    get_timing_files(fname, alpha, methods, operators, dt)




# # -------- 4-splitting: -------- (f_Rabi middle step)
# if multi:
#     fname = directory_data + f"BEC_4-split_Strang-Rabi" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "4-split", exist_ok=True)
#     fname = directory_timing + f"4-split/BEC_4-split_large-Rabi" + my_str + ".csv"
# print("\nTesting 4-splitting: large step for Rabi coupling")
# print(f"Filename: {fname}")

# alpha = "Strang-N" 

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC", (4,): "ANALYTIC"}

# operators = [f_xx, f_yy, f_nl, f_Rabi]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"Str-4-R" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)

    



# # -------- 4-splitting: -------- (f_nl middle step)
# if multi:
#     fname = directory_data + f"BEC_4-split_Strang-nl" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "4-split", exist_ok=True)
#     fname = directory_timing + f"4-split/BEC_4-split_large-nl" + my_str + ".csv"
# print("\nTesting 4-splitting: large step for nonlinear part")
# print(f"Filename: {fname}")

# alpha = "Strang-N" 

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC", (4,): "ANALYTIC"}

# operators = [f_xx, f_yy, f_Rabi, f_nl]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"Str-4-nl" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)





# # -------- 4-splitting: -------- (f_yy middle step)
# if multi:
#     fname = directory_data + f"BEC_4-split_Strang-yy" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "4-split", exist_ok=True)
#     fname = directory_timing + f"4-split/BEC_4-split_large-yy" + my_str + ".csv"
# print("\nTesting 4-splitting: large step for y-part of diffusion (Laplacian)")
# print(f"Filename: {fname}")

# alpha = "Strang-N" 

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC", (4,): "ANALYTIC"}

# operators = [f_Rabi, f_nl, f_xx, f_yy]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"Str-4-y" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)





# # -------- 4-splitting: -------- (f_xx middle step)
# if multi:
#     fname = directory_data + f"BEC_4-split_Strang-xx" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "4-split", exist_ok=True)
#     fname = directory_timing + f"4-split/BEC_4-split_large-xx" + my_str + ".csv"
# print("\nTesting 4-splitting: large step for x-part of diffusion (Laplacian)")
# print(f"Filename: {fname}")

# alpha = "Strang-N" 

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC", (3,): "ANALYTIC", (4,): "ANALYTIC"}

# operators = [f_Rabi, f_nl, f_yy, f_xx]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"Str-4-x" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_Strang)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)





# #----------------- 3rd-order ---------------------#
    
# # -------- 2-splitting: -------- (Ruth method, A=diffR, B=nl)
# if multi:
#     fname = directory_data + f"BEC_2-split_Ruth-1" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "2-split", exist_ok=True)
#     fname = directory_timing + f"2-split/BEC_2-split_R3_diffR_nl" + my_str + ".csv"
# print("\nTesting 2-splitting: Ruth method (3rd-order), A=diffR, B=nl")
# print(f"Filename: {fname}")

# alpha = "R3"  #'SM2' = ABBA, 'Strang' = ABAB, 'Strang-N' is ABCBA, etc 

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

# operators = [f_diffR, f_nl]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"Ruth1" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_Ruth)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)




# # -------- 2-splitting: -------- (Ruth method, A=nl, B=diffR)
# if multi:
#     fname = directory_data + f"BEC_2-split_Ruth-2" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "2-split", exist_ok=True)
#     fname = directory_timing + f"2-split/BEC_2-split_R3_nl_diffR" + my_str + ".csv"
# print("\nTesting 2-splitting: Ruth method (3rd-order), A=nl, B=diffR")
# print(f"Filename: {fname}")

# alpha = "R3"

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

# operators = [f_nl, f_diffR]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"Ruth2" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_Ruth)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)




# # -------- 2-splitting: -------- (AKS3 method, A=diffR, B=nl)
# if multi:
#     fname = directory_data + f"BEC_2-split_AKS3-1" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "2-split", exist_ok=True)
#     fname = directory_timing + f"2-split/BEC_2-split_AKS3-1" + my_str + ".csv"
# print("\nTesting 2-splitting: AKS3 method (3rd-order), A=diffR, B=nl")
# print(f"Filename: {fname}")

# alpha = "AKS3" # palidromic method, p. 102 Wei's thesis

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

# operators = [f_diffR, f_nl]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"AKS3-1" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_Ruth)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)




# # -------- 2-splitting: -------- (AKS3 method, A=nl, B=diffR)
# if multi:
#     fname = directory_data + f"BEC_2-split_AKS3-2" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "2-split", exist_ok=True)
#     fname = directory_timing + f"2-split/BEC_2-split_AKS3-2" + my_str + ".csv"
# print("\nTesting 2-splitting: AKS3 method (3rd-order), A=nl, B=diffR")
# print(f"Filename: {fname}")

# alpha = "AKS3" # palidromic method, p. 102 Wei's thesis

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

# operators = [f_nl, f_diffR]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"AKS3-2" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_Ruth)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)




# # -------- 2-splitting: -------- (OS2(3,4)[7] method, A=diffR, B=nl)
# if multi:
#     fname = directory_data + f"BEC_2-split_OS2-1" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "2-split", exist_ok=True)
#     fname = directory_timing + f"2-split/BEC_2-split_OS2-1" + my_str + ".csv"
# print("\nTesting 2-splitting: OS2(3,4)[7] method (3rd-order), A=diffR, B=nl")
# print(f"Filename: {fname}")

# alpha = [[0,OS2_1], [OS2_2,OS2_3], [OS2_4,OS2_5], [OS2_6,OS2_7]] # p. ?? Wei's thesis

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

# operators = [f_diffR, f_nl]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"OS2-1" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_OS2)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)




# # -------- 2-splitting: -------- (OS2(3,4)[7] method, A=nl, B=diffR)
# if multi:
#     fname = directory_data + f"BEC_2-split_OS2-2" + my_str + ".csv"
# else:
#     os.makedirs(directory_timing + "2-split", exist_ok=True)
#     fname = directory_timing + f"2-split/BEC_2-split_OS2-2" + my_str + ".csv"
# print("\nTesting 2-splitting: OS2(3,4)[7] method (3rd-order), A=nl, B=diffR")
# print(f"Filename: {fname}")

# alpha = [[0,OS2_1], [OS2_2,OS2_3], [OS2_4,OS2_5], [OS2_6,OS2_7]] # p. ?? Wei's thesis

# methods = {(1,): "ANALYTIC", (2,): "ANALYTIC"}

# operators = [f_nl, f_diffR]

# if multi:
#     with open(data_file, "a") as file: # write method name
#         file.write(f"OS2-2" + "  ")
#     get_timing_files_multi(fname, alpha, methods, operators, dt_OS2)
# else:
#     get_timing_files(fname, alpha, methods, operators, dt)

    


# -------------------------------------
if multi:
    with open(data_file, "a") as file: # linebreak
        file.write("\n")

else:
    with open(timing_file, "a") as file: # linebreak
        file.write("\n")
