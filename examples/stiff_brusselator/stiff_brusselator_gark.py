import numpy as np

# second order, stability-decoupled IMEX GARK method from Example 6 of Adrian Sandu and Michael Günther. 2015. A Generalized-Structure Approach to Additive Runge–Kutta Methods

beta = 0.5
Aee = np.array([[0,0,0],
                [1/2,0,0],
                [1 - beta, beta, 0]])
Aei = np.array([[0,0],
                [1/2, 0],
                [1/2, 1/2]])
Aie = np.array([[1/4,0,0],
                [1/4,1/2,0]])
Aii = np.array([[1/4,0],
                [1/2,1/4]])

be = np.array([1/4,1/2,1/4])
bi = np.array([1/2,1/2])

A = [[Aee, Aei], [Aie, Aii]]
B = [be, bi]

# import problem definition
from stiff_brusselator import *

dt = tf/32

# solve the problem
from gark_methods import gark_solve
result = gark_solve([f1, f2], dt, y0, t0, tf, A, B)

print(result)
