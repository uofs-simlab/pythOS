import numpy as np
from gark_methods import gark_solve

def multirate_solve(y0, t0, dt, tf, method, M, fs, ff, bc=None, solver_parameters={}, fname=None, save_steps = 0):
    """ This function uses a multirate method to solve a differential equation
    -----
    Inputs:
    y0 - the current value of y
         if the functions are Forms, this should be of type Function
    t0 - the current value of t
         if using the finite element version, this should be of type Constant
    dt - the amount time will increase by
    tf - the time to solve until
    method - the method to use (of type Multirate)
    fs - the slow operator of the differential equation, with inputs (t, y)
    ff - the fast operator of the differential equation, with inputs (t, y)
    bcs - optional.  Only used if using the finite element version.  The boundary condition(s) to apply.
    solver_parameters - optional.  Any solver parameters to use (see firedrake documentation or scipy.roots documentation for details)
    fname - optional.  If provided, will save intermediate results to this file.
          - if using the finite element version of the code, this is a HDF5 file
            otherwise it is a csv file.
    save_steps - the number of intermediate steps to save if fname is provided
           if it is not provided, the default is to save every step
           (or after every dt if embedded methods are being used).
    Return
    -----
    the approximate value of y at tf
    """

    
    A_ff = method.A_ff
    A_ss = method.A_ss
    b_f = method.b_f
    b_s = method.b_s
    A_fs = method.A_fs
    A_sf = method.A_sf
    
    ss = A_ss.shape[0]
    sf = A_ff.shape[0]
    Asf = np.zeros((ss, sf*M))
    for lam in range(M):
        Asf[:,lam*sf:(lam+1)*sf] = A_sf(lam+1, M)/M
    Afs = np.zeros((sf*M, ss))
    for lam in range(M):
        Afs[lam*sf:(lam+1)*sf,:] = A_fs(lam+1, M)
    Aff = np.zeros((M*sf, M*sf))

    b_block = (np.ones((sf, 1)) @ b_f[np.newaxis]) / M
    for i in range(M):
    
        for j in range(i):
            Aff[i*sf:(i+1)*sf, j*sf:(j+1)*sf] = b_block
        Aff[i*sf:(i+1)*sf, i*sf:(i+1)*sf] = A_ff/M
    b_new = np.zeros((sf*M))
    for i in range(M):
        b_new[i*sf:(i+1)*sf] = 1/M * b_f

    if not np.all(Asf * Afs.transpose() == 0):
        # method is not decoupled, solve as is.
        return gark_solve([ff, fs], dt, y0, t0, tf, [[Aff, Afs], [Asf, A_ss]], [b_new, b_s], bc=bc, solver_parameters=solver_parameters, fname=fname, save_steps=save_steps)
    order = []
    last_f = 0
    for j in range(ss):
        t_f = 0
        for i in range(M*sf):
            if Asf[j,i] != 0:
                t_f = i+1
        for i in range(last_f, t_f):
            order.append(i)
        last_f = max(last_f, t_f)
        order.append(M*sf + j)
    for i in range(last_f, M*sf):
        order.append(i)
    return gark_solve([ff, fs], dt, y0, t0, tf, [[Aff, Afs], [Asf, A_ss]], [b_new, b_s], bc=bc, solver_parameters=solver_parameters, fname=fname, save_steps=save_steps, order=order)



class Multirate:
    def __init__(self, A_ff, A_ss, b_f, b_s, A_fs, A_sf):
        self.A_ff = A_ff
        self.A_ss = A_ss
        self.b_f = b_f
        self.b_s = b_s
        self.A_fs = A_fs
        self.A_sf = A_sf

mrgark_ex2_im2 = Multirate(np.array([[0, 0], [2/3, 0]]),
                           np.array([[1-1/np.sqrt(2), 0],
                                     [1/np.sqrt(2), 1-1/np.sqrt(2)]]),
                           np.array([1/4, 3/4]),
                           np.array([1/np.sqrt(2), 1-1/np.sqrt(2)]),
                           lambda lam, M: np.array([[(lam-1) / M, 0], [(3*lam-1)/(3*M), 0]]),
                           lambda lam, M: np.array([[M-M/np.sqrt(2) if lam == 1 else 0, 0], [1/4, 3/4]]))
