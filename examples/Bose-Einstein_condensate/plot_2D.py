import numpy as np
from matplotlib import pyplot as plt
import sys


xmax = 40.96
xmin = -40.96
ymax = 1.72
ymin = -1.72
dx = 0.04
dy = dx

omega = 1.0

if len(sys.argv) > 1:  # system parameters from file
    infile = sys.argv[1]

    with open(infile, 'r') as file:
        content = file.read()
        words = content.split()
        #R0 = float(words[1])
        xmax = float(words[3])
        xmin = float(words[5])
        ymax = float(words[7])
        ymin = float(words[9])
        dx = float(words[11])
        dy = dx
        #dt = float(words[13])
        #t_start = float(words[15])
        #t_end = float(words[17])
        omega = float(words[23])

fname = f"BEC_omega{omega}/BEC_omega{omega}_9.csv"

if len(sys.argv) > 2: # use filename from command line
    fname = sys.argv[2]

# if len(sys.argv) > 3:
#     nt = float(sys.argv[3])  # which plot to view

try:
    pdata = np.loadtxt(fname, delimiter=',', dtype=np.complex128)
except FileNotFoundError:
    print("Error: no file given to plot or file not found.")
    sys.exit()


Nx = int((xmax - xmin) / dx) + 1
Ny = int((ymax - ymin) / dy) + 1

x = np.array([xmin + (i) * dx for i in range(Nx)])
y = np.array([ymin + (j) * dy for j in range(Ny)])

X, Y = np.meshgrid(x, y)


Nt = 0 # number of saves (different times)

if len(pdata.shape) == 1:
    Nt = 1
    t = np.real(pdata[0])
    pdata = pdata[1:]
    pdata = np.reshape(pdata, (2,Nx,Ny))
    Z1 = np.transpose(pdata[0, :, :]) # psi1
    Z2 = np.transpose(pdata[1, :, :]) # psi2
else:
    Nt = pdata.shape[0]
    t = np.real(pdata[:,0])
    pdata = pdata[:,1:]
    pdata = np.reshape(pdata, (Nt,2,Nx,Ny))
    # if len(sys.argv) > 3 and nt < Nt + 1e-06:
    #     Nt = int(nt)
    # elif nt > Nt + 1e-06:
    #     print(f"Number of plot data is {Nt}.")
    t = t[Nt-1]
    Z1 = np.transpose(pdata[Nt-1, 0, :, :]) # psi1
    Z2 = np.transpose(pdata[Nt-1, 1, :, :]) # psi2

phase_rel = (np.angle(Z1) - np.angle(Z2)) % (2*np.pi)
for i in range(phase_rel.shape[0]):
    for j in range(phase_rel.shape[1]):
        if phase_rel[i,j] > np.pi:
            phase_rel[i,j] += -2*np.pi
    
print(f"t = {t}")
print(f"omega = {omega}")


plt.rcParams.update({'font.size':22})

# component 1 (density):
fig1, ax1 = plt.subplots(figsize=(7.4,4.8))
im1 = ax1.pcolormesh(X, Y, abs(Z1)**2, cmap='Spectral_r')
cbar1 = fig1.colorbar(im1, ax=ax1, aspect=8, shrink=0.8, anchor=(0,0.05))

# component 2 (density):
fig2, ax2 = plt.subplots(figsize=(7.4,4.8))
im2 = ax2.pcolormesh(X, Y, abs(Z2)**2, cmap='Spectral_r')
cbar2 = fig2.colorbar(im2, ax=ax2, aspect=8, shrink=0.8, anchor=(0,0.05))

# relative phase:
fig3, ax3 = plt.subplots(figsize=(7.2,4.8))
im3 = ax3.pcolormesh(X, Y, phase_rel, cmap='seismic')
cbar3 = fig3.colorbar(im3, ax=ax3, aspect=8, shrink=0.8, anchor=(0,0.05))

ax1.set_aspect('equal', adjustable='box')
ax1.set_ylabel('x')
ax1.set_xlabel('y')
ax1.set_xlim(-4,4)
ax1.set_ylim(-2,2)
ax1.set_title(f'$t = ${t:.1f}', style='italic', loc='left')

ax2.set_aspect('equal', adjustable='box')
ax2.set_ylabel('x')
ax2.set_xlabel('y')
ax2.set_xlim(-4,4)
ax2.set_ylim(-2,2)
ax2.set_title(f'$t = ${t:.1f}', style='italic', loc='left')

ax3.set_aspect('equal', adjustable='box')
ax3.set_ylabel('x')
ax3.set_xlabel('y')
ax3.set_xlim(-4,4)
ax3.set_ylim(-2,2)
ax3.set_title(f'$t = ${t:.1f}', style='italic', loc='left')

ax1.set_yticks([-2,-1,0,1,2])
ax2.set_yticks([-2,-1,0,1,2])
ax3.set_yticks([-2,-1,0,1,2])

im1.set_clim(0,1)
cbar1.ax.locator_params(nbins=2)

im2.set_clim(0,1)
cbar2.ax.locator_params(nbins=2)

val = np.pi
im3.set_clim(-val, val)
#clabels = [r'$-\frac{\pi}{2}$', r'$0$', r'$\frac{\pi}{2}$']
clabels = [r'$-\pi$', r'$0$', r'$\pi$']
cbar3.set_ticks([-val, 0, val])
cbar3.set_ticklabels(clabels)

# Optional: Adjust label properties (e.g., rotation, padding)
# By default, vertical colorbar labels are rotated. To make it horizontal:
cbar1.ax.xaxis.set_label_position('top')
cbar2.ax.xaxis.set_label_position('top')
cbar3.ax.xaxis.set_label_position('top')
cbar1.ax.set_xlabel('$n_1$', rotation=0, labelpad=10)
cbar2.ax.set_xlabel('$n_2$', rotation=0, labelpad=10)
cbar3.ax.set_xlabel(r'$\theta$', rotation=0, labelpad=10)


ax1.set_xlabel('y', labelpad=4)
for tick in ax1.xaxis.get_major_ticks():
    tick.set_pad(10) # Increase the padding (shift down)

ax1.set_ylabel('x', labelpad=4)
for tick in ax1.yaxis.get_major_ticks():
    tick.set_pad(10) # Increase the padding (shift left)

ax2.set_xlabel('y', labelpad=4)
for tick in ax2.xaxis.get_major_ticks():
    tick.set_pad(10) # Increase the padding (shift down)

ax2.set_ylabel('x', labelpad=4)
for tick in ax2.yaxis.get_major_ticks():
    tick.set_pad(10) # Increase the padding (shift left)

ax3.set_xlabel('y', labelpad=4)
for tick in ax3.xaxis.get_major_ticks():
    tick.set_pad(10) # Increase the padding (shift down)

ax3.set_ylabel('x', labelpad=4)
for tick in ax3.yaxis.get_major_ticks():
    tick.set_pad(10) # Increase the padding (shift left)


plt.show()

