import numpy as np

# second order MrGARK method (Ex2-Im2) from Arash Sarshar, Steven Roberts,
# and Adrian Sandu. 2019. Design of High-Order Decoupled Multirate GARK Schemes.
A_ff = np.array([[0,0],
                 [2/3,0]])
A_ss = np.array([[1-1/np.sqrt(2), 0],
                 [1/np.sqrt(2), 1-1/np.sqrt(2)]])
def A_sf(lam, M):
    if lam == 1:
        return np.array([[M-M/np.sqrt(2), 0],
                         [1/4, 3/4]])
    else:
        return np.array([[0, 0],
                         [1/4, 3/4]])

def A_fs(lam, M):
    return np.array([[(lam-1) / M, 0], [(3*lam-1)/(3*M), 0]])

b_ff = np.array([1/4,3/4])
b_ss = np.array([1/np.sqrt(2), 1-1/np.sqrt(2)])


# import problem definition
from stiff_brusselator import *

dt = tf/256

# solve the problem
from multirate import multirate_solve, Multirate

method = Multirate(A_ff, A_ss, b_ff, b_ss, A_fs, A_sf)
result = multirate_solve(y0, t0, dt, tf, method, 32, fs, ff)

print(result)
