import numpy as np
from kiops import kiops
from exode import exode
try:
    from firedrake import *
    from kiops_fem import kiops_fem
except:
    Function = type(None)

try:
    from sundials_exode import SundialsExode
except:
    SundialsExode = type(None)

class EpiMultistep:
    def __init__(self, A):
        self.A = A
        self.dt = None
        self.yn = []
        self.fn = []
        self.prev_steps = A.shape[1]
    def step_fem(self, fi, t, y, dt, options=('kiops', 1e-10), bc = None, **kwargs):
        if dt == 0:
            return y
        f_in = fi
        if isinstance(fi, tuple):
            bc = fi[1]
            fi = fi[0]
        test_f = fi.arguments()[0]
        def solver(t_out, A, u, task1=False):
            if options[0] == 'kiops':
                return kiops_fem(t_out, A, u, task1 = task1, tol=options[1], **kwargs)
            else:
                return exode(t_out, A, u, task1=task1, tol=options[1], method=options[0], **kwargs)
        F = Function(y)
        solve(inner(F, test_f) * dx - fi == 0, F, bcs=bc)
        if len(self.fn) < self.prev_steps:
            self.yn.insert(0, Function(y))
            self.fn.insert(0, F)
            result = operator_splitting([f_in], dt, y, t, t.values()[0]+dt, 'Godunov', methods={(0,): "EXACT"}, ivp_methods={1: ('Dormand-Prince', 1e-12, 1e-14)})
            return result

        v = Function(y)
        dJ = derivative(fi, y)
        
        if len(fi.coefficients()) > 0:
            matvec_handle = dt * dJ * v
        else:
            Z = Constant(0.0)
            matvec_handle = Z *  test_f * dx
            dJ = 0

        vm = [Function(y).assign(0) for _ in range(self.A.shape[0] + 2)]
        vm[0] = v
        vm[0].assign(0)
        vm[1].assign(F)
        for i in range(self.prev_steps):
            r = Function(y)
            if dJ != 0:
                solve(inner(r, test_f) * dx - (dJ * self.yn[i] - dJ *y) == 0, r)
            else:
                r.assign(0)
            r.assign(self.fn[i] - F - r)

            for k in range(self.A.shape[0]):

                vm[k+2].assign(float(self.A[k, i]) * r + vm[k+2])

        ph, _ = solver([1], matvec_handle, vm)
        if bc is not None:
            if isinstance(bc, list):
                for bc_i in bc:
                    bc_i.apply(ph[0])
            else:
                bc.apply(ph[0])
        self.fn.insert(0, F)
        self.yn.insert(0, Function(y))
        y.assign(y + ph[0] * dt)
        return y
        
                       
    def step(self, fi, t, y, dt, options=('kiops', 1e-10), **kwargs):
        if isinstance(y, Function):
            return self.step_fem(fi, t, y, dt, options, **kwargs)
        y = np.concatenate((y, [t]))
        f = lambda v: np.concatenate((fi(v[-1], v[:-1]), [1]))

        if isinstance(options[0], SundialsExode):
            y = y.copy()
        def solver(t_out, A, u, task1=False):
            if options[0] == 'kiops':
                return kiops(t_out, A, u, task1 = task1, tol=options[1], **kwargs)
            else:
                return exode(t_out, A, u, task1=task1, tol=options[1], method=options[0], **kwargs)
        if len(self.fn) < self.prev_steps:
            self.yn.insert(0, y)
            self.fn.insert(0, f(y))
            result = operator_splitting([lambda t, y: f(y)], dt, y, t, t+dt, 'Godunov', methods={(0,): "EXACT"})
            
            return result
            
        F = f(y)
        An = jac(f, y, t)
        matvec_handle = lambda v: (An * dt) @ v

        vm = np.zeros((self.A.shape[0]+2, y.size), dtype=y.dtype)
        vm[1,:] = F
        for i in range(self.prev_steps):
            r = self.fn[i] - F - (An @ (self.yn[i] - y))
            for k in range(self.A.shape[0]):
                vm[k+2,:] += self.A[k, i] * r
        ph, _ = solver([1], matvec_handle, vm)
        self.fn.insert(0, F)
        self.yn.insert(0, y)
        return y + ph[0] * dt
    def solve(self, f, t, y, dt, n_steps=20, options=('Dormand-Prince', (1e-10, 1e-12)), monolithic=False, **kwargs):
        if not monolithic:
            self.fn = []
            self.yn = []
        if self.dt is None or abs(dt/n_steps - self.dt) > 1e-6:
            self.fn = []
            self.yn = []
        self.dt = dt/n_steps
        for i in range(n_steps):
            y = self.step(f, t, y, self.dt, options, **kwargs)
            if not isinstance(y, Function):
                y = y[:-1]
            if not isinstance(t, Constant):
                t += self.dt
            else:
                t.assign(t + self.dt)
            
        if len(self.fn) > self.prev_steps:
            self.fn = self.fn[:self.prev_steps]
            self.yn = self.yn[:self.prev_steps]
        return y
        
def jac(f, y, t):
    out = np.zeros((y.size, y.size), y.dtype)
    eps = 1e-5
    f0 = f(y)
    for i in range(y.size):
        y[i] += eps
        out[:,i] = (f(y) - f0) / eps
        y[i] -= eps
    return out


epi_methods = {'EPI2': np.array([[]]),
               'EPI3': np.array([[2/3]]),
               'EPI4': np.array([[-3/10, 3/40], [32/5, -11/10]]),
               'EPI5': np.array([[-4/5, 2/5, -4/45],
                                 [12, -9/2, 8/9],
                                 [3, 0, -1/3]
                                 ]),
               'EPI6': np.array([[ -49/60, 351/560, -359/1260, 367/6720],
                                 [92/7, -99/14, 176/63, -1/2],
                                 [485/21, -151/14, 23/9, -31/168]])
}
