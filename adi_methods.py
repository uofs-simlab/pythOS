import numpy as np
import scipy.linalg as linalg
from butcher_tableau import constant_matrix
from scipy import sparse

from scipy.sparse.linalg import spsolve

def mcs(functions, initial_y, initial_t, delta_t):
    """
    This function takes a step in time in the differential equation using the
    Modified Craig-Sneyd method
    Inputs
    ------
    functions : list of functions to use
        functions to use to approximate the differential equation in order
        the functions will be numbered 1 to n
        if any element returns np.nan, it will not be considered 
        inputs are (t, y)
    initial_y
        the current value of y to use
    initial_t
        the current value of t to use
    delta_t : float
        the ammount time will increase by
    Return
    ------
        the approximate next value of y after the given delta_t
    """

    k=len(functions)
    alpha= [[1]*k,
            [0]+[1./3]*(k-1),
            [1./3]+[0]*(k-1),
            [1./2-1./3]*k,
            [0]+[1./3]*(k-1)]
    y0=np.array(initial_y, initial_y.dtype)
    for i in range(len(functions)):
        y0+=delta_t*functions[i](initial_t, initial_y)*alpha[0][i]
    yj=y0
    for j in range(1,len(functions)):
        A_matrix=constant_matrix(functions[j], yj, initial_t+delta_t)
        left=sparse.identity(np.size(initial_y))-alpha[1][j]*delta_t*A_matrix
        right=yj+alpha[1][j]*delta_t*(-functions[j](initial_t, initial_y)+functions[j](initial_t+delta_t, yj)-(A_matrix*yj))
        yj=spsolve(left, right)
    for i in range(2,len(alpha)-1):
        for j in range(len(functions)):
            y0+=alpha[i][j]*delta_t*(functions[j](initial_t+delta_t, yj)-functions[j](initial_t, initial_y))
    yj=y0
    for j in range(1,len(functions)):
        A_matrix=constant_matrix(functions[j], yj, initial_t)
        left=sparse.identity(np.size(initial_y))-alpha[-1][j]*delta_t*A_matrix
        right=yj+alpha[-1][j]*delta_t*(-functions[j](initial_t, initial_y)+functions[j](initial_t+delta_t, yj)-A_matrix*yj)
        yj=spsolve(left, right)
    return yj

def hv(functions, initial_y, initial_t, delta_t):
    """
    This function takes a step in time in the differential equation using the
    Nudsdorfer-Verwer method
    Inputs
    ------
    functions : list of functions to use
        functions to use to approximate the differential equation in order
        the functions will be numbered 1 to n
        if any element returns np.nan, it will not be considered 
        inputs are (t, y)
    initial_y
        the current value of y to use
    initial_t
        the current value of t to use
    delta_t : float
        the ammount time will increase by
    Return
    ------
        the approximate next value of y after the given delta_t
    """

    k=len(functions)
    alpha= [[1]*k,
            [0]+[1./3]*(k-1),
            [1./2]*k,
            [0]+[1./3]*(k-1)]
    y0=np.array(initial_y, initial_y.dtype)
    for i in range(len(functions)):
        y0+=delta_t*functions[i](initial_t, initial_y)*alpha[0][i]
    yj=y0
    for j in range(1,len(functions)):
        A_matrix=constant_matrix(functions[j], yj, initial_t+delta_t)
        left=np.identity(np.size(initial_y))-alpha[1][j]*delta_t*A_matrix
        right=yj+alpha[1][j]*delta_t*(-functions[j](initial_t, initial_y)+functions[j](initial_t+delta_t, yj)-np.matmul(A_matrix, yj))
        yj=linalg.solve(left, right)
    yk=np.array(yj)
    for i in range(2,len(alpha)-1):
        for j in range(len(functions)):
            y0+=alpha[i][j]*delta_t*(functions[j](initial_t+delta_t, yj)-functions[j](initial_t, initial_y))
    yj=y0
    for j in range(1,len(functions)):
        A_matrix=constant_matrix(functions[j], yj, initial_t)
        left=np.identity(np.size(initial_y))-alpha[-1][j]*delta_t*A_matrix
        right=yj+alpha[-1][j]*delta_t*(-functions[j](initial_t+delta_t, yk)+functions[j](initial_t+delta_t, yj)-np.matmul(A_matrix, yj))
        yj=linalg.solve(left, right)
    return yj

def dr(functions, initial_y, initial_t, delta_t):
    """
    This function takes a step in time in the differential equation using the
    Douglas-Rachford method
    Inputs
    ------
    functions : list of functions to use
        functions to use to approximate the differential equation in order
        the functions will be numbered 1 to n
        if any element returns np.nan, it will not be considered 
        inputs are (t, y)
    initial_y
        the current value of y to use
    initial_t
        the current value of t to use
    delta_t : float
        the ammount time will increase by
    Return
    ------
        the approximate next value of y after the given delta_t
    """

    A_matrix=constant_matrix(functions[0], initial_y, initial_t+delta_t)
    left=np.identity(np.size(initial_y))-delta_t*A_matrix
    right=initial_y+delta_t*functions[1](initial_t,initial_y)
    try:
        yj=linalg.solve(left,right)
    except:
        return np.array([np.nan])
    A_matrix=constant_matrix(functions[1], yj, initial_t+delta_t)
    left=np.identity(np.size(initial_y))-delta_t*A_matrix
    right=initial_y+delta_t*functions[0](initial_t+delta_t,yj)
    try:
        y=linalg.solve(left,right)
    except:
        return np.array([np.nan])
    return y

methods={'MCS': mcs, 'HV': hv, 'DR': dr}

def adi_step(functions, initial_t, delta_t, initial_y, method):
    """
    This function takes a step in time in the differential equation using an 
    ADI type method
    Inputs
    ------
    functions : list of functions to use
        functions to use to approximate the differential equation in order
        the functions will be numbered 1 to n
        if any element returns np.nan, it will not be considered 
        inputs are (t, y)
    initial_t
        the current value of t to use
    delta_t : float
        the ammount time will increase by
    initial_y
        the current value of y to use
    methods : string
        the method to use
    Return
    ------
    float or array
        type depends on type of initial_y
        the approximate next value of y after the given delta_t
    """

    return methods[method](functions, initial_y, initial_t, delta_t)

