import numpy as np
from kiops import kiops
from exode import exode
try:
    from firedrake import Function, solve, dx, replace, inner, derivative, Constant
    from firedrake import *
    from kiops_fem import kiops_fem
except:
    Function = type(None)

class EPIMethod:
    def __init__(self, p, g, ab):
        """ Initialize EPI Method object.
        array ab contains the a coefficients, and the b coefficients in the last row
        """
        self.p = p
        self.g = g
        self.a = ab
        self.S = ab.shape[0]

    def step_fem(self, fi, t, y, dt, options=('kiops', 1e-3), params=None, bc=None):
        if dt == 0:
            return y
        if isinstance(fi, tuple):
            bc = fi[1]
            fi = fi[0]
        test_f = fi.arguments()[0]

        def solver(t_out, A, u, task1=False, **kwargs):
            if options[0] == 'kiops':
                return kiops_fem(t_out, A, u, task1 = task1, tol=options[1], **kwargs)
            else:
                return exode(t_out, A, u, task1=task1, tol=options[1], method=options[0], **kwargs)
        Ri = [[Function(y).assign(0) for _ in range(self.a.shape[0])] for _ in range(self.a.shape[0])]

        F = Function(y)
        solve(inner(F, test_f) * dx - fi == 0, F)

        Yi = [Function(y).assign(0) for _ in range(self.S)]
            
        v = Function(y)
        dJ = derivative(fi, y)
        matvec_handle = dt * dJ * v

        for i in range(self.S):
            Yi[i].assign(y)
            vec = [Function(y).assign(0) for _ in range(2)]
            vec[1].assign( dt * float(self.p[0,0] * self.a[i, 0]) * F)
            
            ph, _ = solver([self.g[i, 0]], matvec_handle, vec, task1=True, params=params)
            Yi[i] += ph[0]

            for j in range(1, self.p.shape[0]):
                if self.a[i,j] != 0:
                    vec = [Function(y).assign(0) for _ in range(self.p.shape[1] + 1)]
                    vec[0] = v
                    vec[0].assign(0)
                    for k in range(self.p.shape[1]):
                        vec[k+1].assign ( float(self.p[j,k]) * Ri[j][j-1] * self.a[i, j]*dt / (self.g[i,j]**k))
           
                    ph, _ = solver([self.g[i, j]], matvec_handle, vec, task1=True, params=params)
                    Yi[i] += ph[0]
                t0 = Constant(t)
                t.assign(t0 + dt * i / self.S)
                if bc is not None:
                    bc.apply(Yi[i])
                t.assign(t0)
            solve(inner(Ri[0][i], test_f) * dx - replace(fi, {y: Yi[i]}) == 0, Ri[0][i])
            F2 = Function(y)
            solve(inner(F2, test_f) * dx - 1/dt * replace(matvec_handle, {v: Yi[i] - y}) == 0, F2)
            Ri[0][i].assign(Ri[0][i] - F - F2)
                        
            for j in range(1, self.S):
                Ri[j][i].assign (Ri[j-1][ i] - Ri[j-1][ i-1])
        y.assign(Yi[-1])
        return y

    def step(self, fi, t, y, dt, options=('kiops', 1e-3), **kwargs):
        """ Find the value of y after time dt """
        if isinstance(y, Function):
            return self.step_fem(fi, t, y, dt, options, **kwargs)
        y = np.concatenate((y, [t]))
        f = lambda v: np.concatenate((fi(v[-1], v[:-1]), [1]))
        
        def solver(t_out, A, u, task1=False):
            if options[0] == 'kiops':
                return kiops(t_out, A, u, task1 = task1, tol=options[1])
            else:
                return exode(t_out, A, u, task1=task1, tol=options[1], method=options[0])
        Ri = np.zeros((self.a.shape[0], self.a.shape[0], y.size), y.dtype)

        F = f(y)
        Yi = np.zeros((self.S, y.size), y.dtype)
        An = jac(f, y, t)
        
        matvec_handle = lambda v: (An * dt) @ v

        for i in range(self.S):
            Yi[i] = y
            vec = np.zeros((2, y.size), y.dtype)
            vec[1,:] = F * dt * self.p[0,0] * self.a[i, 0]
            ph, _ = solver([self.g[i, 0]], matvec_handle, vec, task1=True)
            Yi[i] += ph[0,:]
            for j in range(1, self.p.shape[0]):
                if self.a[i,j] != 0:
                    vec = np.zeros((self.p.shape[1]+1, y.size), y.dtype)
                    for k in range(self.p.shape[1]):
                        vec[k+1,:] = self.p[j,k] * Ri[j, j-1] * self.a[i, j]*dt / (self.g[i,j]**k)

                    ph, _ = solver([self.g[i, j]], matvec_handle, vec, task1=True)
                    Yi[i] += ph[0,:]
            
            Ri[0,i] = f(Yi[i]) - F - An @ (Yi[i]-y)
                        
            for j in range(1, self.S):
                Ri[j, i] = (Ri[j-1, i] - Ri[j-1, i-1])
        return Yi[-1][:-1]

    def y_step(self, f, y, t0, dt, mask):
        """Return the adjustment on y for the next step"""
        return (self.step(f, t0, y, dt) - y) / dt

def jac(f, y, t):
    out = np.zeros((y.size, y.size), y.dtype)
    eps = 1e-5
    f0 = f(y)
    for i in range(y.size):
        y[i] += eps
        out[:,i] = (f(y) - f0) / eps
        y[i] -= eps
    return out


epi_methods = {
        'EPI2': EPIMethod(np.array([[1, 0], [0, 1]]), np.array([[0.5, 0], [1, 1]]), np.array([[0, 0], [1, 0/3]])),
        'EPI3': EPIMethod(np.array([[1, 0], [0, 1]]), np.array([[0.5, 0], [1, 1]]), np.array([[1, 0], [1, 2/3]])),
        }
s = 30 ** 0.5
epi_methods['EPI4']= EPIMethod(np.array([[1, 0, 0], [0, 1, 0], [0, -1, 6]]),
        np.array([[1/3, 0, 0], [2/3, 2/3, 0], [1, 1, 1]]),
        np.array([[27 * (s**2 + 18)/(12*(54-3*s**2+2*s**3)), 0,0], [18*s*(s**2+18)/(24*(54 - 3*s**2+2*s**3)), 0, 0], [1, 96*(54-s**2)*(54-3*s**2+2*s**3)**2/ (729 * (s**2 + 18)**3), 384*(54-3*s**2 + 2*s**3)**2 / (162 * (s**2 + 18) **3)]]))


epi_methods['EPI4s3'] = EPIMethod(np.array([[1, 0, 0, 0], [0, 0, 1892, -42336],
    [0, 0, 1458, -34992]]),
    np.array([[1/8, 0, 0], [1/9, 0, 0], [1, 1, 1]]),
    np.array([[1/8, 0, 0],
        [1/9, 0, 0], [1, 1, 1]]))
        
epi_methods['EPI5s3'] = EPIMethod(np.array([[1, 0, 0], [2/3, 2/3, 0], [2/3, 1/2, 1]]), np.array([[0.41657015580651858694, 0, 0], [0.86246743701274574979, 0.5, 0], [1, 0.730416157608327661916, 0.325076967060782773227]]), np.array([[0.41657015580651858694, 0, 0], [0.86246743701274574979, 1.32931146991722972036, 0], [1, 1.15467303405015770322, 0.30931492086655796815]]))

epi_methods['EPI5s4'] = EPIMethod(np.array([[1, 0, 0, 0], [0, 4, 0, 0], [0, -2, 16, 0], [0, 4/3, -16, 64]]), np.array([[1/4,0,0,0],[1/2,1/2,0,0],[3/4,1/2,1/2,0],[1, (195+3315**0.5)/260, (255-3315**0.5)/204, (60+3315**0.5)/104]]), np.array([[1/4,0,0,0],[2/4, -1/108, 0, 0], [3/4, 1/2, 1/4,0], [1, 52/45, 68/75,52/45]]))
