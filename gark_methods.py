import numpy as np
from additive_rk import ark_solve
from butcher_tableau import Tableau

def gark_convert(A, b, order=None):
    """ Convert the provided GARK structure to an ARK structure """
    # This assumes M = N, and J = I
    tableaus = []
    Ni = sum([x.size for x in b])
    pre = 0
    for i in range(len(b)):
        si = b[i].size
        ai = np.zeros((Ni, Ni))
        ai[:,pre:pre+si] = np.concatenate([A[q][i] for q in range(len(b))])
        bi = np.zeros(Ni) 
        bi[pre:pre+si] = b[i]
        if (order is not None):
            ai = ai[:,order][order,:]
            bi = bi[order]
        ci = ai.sum(axis=1)
        tableaus.append(Tableau(ci, ai, bi))
        pre += si
    return tableaus


def gark_solve(f, dt, y0, t0, tf, A, b, bc=None, solver_parameters={}, fname=None, save_steps = 0, jacobian=None, order=None):
    """ This function uses an GARK method to solve a differential equation
    This is done by converting to an ARK structure
    -----
    Inputs:
    functions - the operators of the differential equation, each with inputs (t, y)
                these may also be finite element Forms as provided by firedrake
    dt - the amount time will increase by
    y0 - the current value of y
         if using the finite element version, this should be of type Function
    t0 - the current value of t
         if using the finite element version, this should be of type Constant
    tf - the time to solve until
    A - the list of lists containing the component arrays a{J(q), i}
    b - the list containing the b{i} arrays
    bc - optional.  Only used if using the finite element version.  The boundary condition(s) to apply.
    solver_parameters - optional.  Only used for the finite element version.  Any solver parameters to use (see firedrake documentation for details)
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

    if len(f) > len(b):
        print("ERROR: not enough tableau provided")
        return

    tableaus = gark_convert(A, b, order)

    return ark_solve(f, dt, y0, t0, tf, tableaus, bc=bc, solver_parameters=solver_parameters, fname=fname, save_steps = save_steps, jacobian=jacobian)


