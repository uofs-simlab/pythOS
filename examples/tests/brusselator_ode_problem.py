import numpy as np
A = 1
B = 3
t0 = 0
tf = 15

# dx/dt = x^2 y
# dydx = -x^2 y
def f1(t, y):
    x = y[0]
    y = y[1]
    return np.array([x*x*y, -x*x*y])

# dx/dt = A - B x^2
# dy/dt = Bx
def f2(t, y):
    x = y[0]
    y = y[1]
    return np.array([A - B*x-x, B*x])

# splitting f2 into two component:
# dx/dt = - B x^2
# dy/dt = Bx
def fi(t, y):
    x = y[0]
    y = y[1]
    return np.array([-B*x-x, B*x])
# dx/dt = A
# dy/dt = 0
def fe(t, y):
    x = y[0]
    y = y[1]
    return np.array([A, 0])

# initial condition
y0 = np.array([1.0, 1.0])
