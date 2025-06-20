import numpy as np
A = 1
B = 3
t0 = 0
tf = 15

def f1(t, y):
    x = y[0]
    y = y[1]
    return np.array([x*x*y, -x*x*y])
def f2(t, y):
    x = y[0]
    y = y[1]
    return np.array([A - B*x-x, B*x])

def fi(t, y):
    x = y[0]
    y = y[1]
    return np.array([-B*x-x, B*x])
def fe(t, y):
    x = y[0]
    y = y[1]
    return np.array([1, 0])
y0 = np.array([1.0, 1.0])
