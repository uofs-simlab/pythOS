import numpy as np
from scipy import linalg, sparse
from math import sqrt, cos, pi
#from scipy.sparse.linalg.dsolve import linsolve
from scipy.optimize import root
import timeit
try:
    from firedrake import Constant, Function, solve, replace, inner, TestFunction, split, dx, errornorm, DirichletBC, norm, Form
except:
    fem = False
    Function = type(None)

class Tableau:
    def __init__(self, c, a, b,explicit=False, sdirk=False):
        """ Creates a butcher tableau from the given arrays"""
        self._c=c
        self._a=a
        self._b=b
        self._explicit=np.all(np.equal(np.triu(a), 0))
        self._sizeb=b.size
        self._sdirk = not self._explicit and np.all(np.equal(np.triu(a,1),0)) and np.all(np.isclose(np.diag(a), a[0,0]))
        self._implicit=not self._sdirk and not self._explicit
        if self._implicit:
            try:
                self._Ainv = np.linalg.inv(self._a)
            except:
                self._Ainv = None

    def get_b_values(self):
        """Return the b array from the butcher tableau"""
        return self._b
    def get_x_values(self, current_x, step):
        """Return the x values for the next step of the equation"""
        return step*self._c+current_x
    def get_k_values_explicit(self, function, current_y, current_x, step):
        """Return the k values (slopes) for the next step of the solution
           for an explicit Runge--Kutta method"""
        x_values=step*self._c+current_x
        next=np.zeros((np.size(current_y),self._sizeb), current_y.dtype)
        for i in range(np.size(self._b)):
            k=function(x_values[i], current_y+step*np.dot(next, self._a[i]))
            next[:,i]=k
        return next
    
    def get_k_implicit(self, function, current_y, current_x,
                              step, **kwargs):
        """Return the k values (slopes) for the next step of the solution
           for an implicit Runge--Kutta method"""
        Y = np.zeros((current_y.size, self._sizeb), current_y.dtype)
        for i in range(self._sizeb):
            Y[:,i] = current_y
        
        y0 = Y.flatten(order='F')
        
        # scipy.optimize.roots cannot handle complex numbers, so if the input is complex, interpret each complex number as two floating points
        if current_y.dtype == np.complex128:
            y0 = y0.view(np.float64)

        if 'method' not in kwargs:
            kwargs['method'] = 'krylov'
        sol = root(implicit_root_func, y0, args=(current_x, current_y, function, step, self), **kwargs)
        if not sol.success:
            print(sol.message)
            print(np.max(abs(sol.fun)))
        
        sol_x = sol.x
        # if necessary, reinterpret the result as an array of complex numbers
        if current_y.dtype == np.complex128:
            sol_x = sol_x.view(np.complex128)
        
        Y = sol_x.reshape(Y.shape, order='F')
        k_i = np.zeros(Y.shape, current_y.dtype)
        
        for i in range(self._sizeb):
            k_i[:,i] = function(current_x + step * self._c[i], Y[:,i])

        return k_i

    def get_k_values_sdirk(self, function, current_y, current_x, step, **kwargs):

        k=np.zeros((np.size(current_y), self._sizeb), current_y.dtype)
        x=step*self._c+current_x
        if 'method' not in kwargs:
            kwargs['method'] = 'krylov'

        for i in range(self._sizeb):

            y0 = current_y
            
            # scipy.optimize.roots cannot handle complex numbers, so if the input is complex, interpret each complex number as two floating points
            if current_y.dtype == np.complex128:
                y0 = y0.view(np.float64)

            sol = root(dirk_root_function, y0, args = (current_x, current_y, function, step, k, i, self), **kwargs)
            sol_x = sol.x
            
            if not sol.success:
                print(sol.message)
                print(np.max(abs(sol.fun)))
                 
            # if necessary, reinterpret the result as an array of complex numbers
            if k.dtype == np.complex128:
                sol_x = sol_x.view(np.complex128)

            k[:,i]=function(x[i], sol_x)
        return k

    def get_k_values(self, function, current_y, current_x, step, bc=None, params=None):
        """Return the k values (slopes) for the next step of the solution"""
        if isinstance(function, tuple) or isinstance(current_y, Function):
            if isinstance(function, tuple):
                f = function[0]
                test_f = f.arguments()[0]
                if len(function) == 2:
                    if isinstance(function[1], DirichletBC) or isinstance(function[1], list):
                        bcs = function[1]
                    else:
                        f = function
                else:
                    bcs = function[1]
                    f = (function[0], function[2])
            elif isinstance(function, Form):
                f = function
                test_f = function.arguments()[0]
                bcs = bc
            else:
                f = function
                test_f = None
                bcs = bc
            y_s = Function(current_y)
            M = 1
            g = None
            if isinstance(f, tuple):
                M = f[1][0]
                g = f[1][1]
                f = f[0]
            if self._explicit:
                k = self.k_erk_fem(f, current_y, current_x, step, test_f, bcs=bcs, params = params)
            elif self._sdirk:
                k = self.k_dirk_fem(f, current_y, current_x, step, test_f, bcs=bcs, params = params, M=M)
            else:
                k = self.k_irk_fem(f, current_y, current_x, step, test_f, bcs=bcs, params=params, M=M)
            current_y.assign(y_s)
            return k
                
        if self._explicit:
            return self.get_k_values_explicit(function, current_y, current_x,step)
        elif self._sdirk:
            return self.get_k_values_sdirk(function, current_y, current_x, step)
        else:
            return self.get_k_implicit(function, current_y, current_x, step)


    def k_erk_fem(self, F, y0, t, dt, test_f, bcs=None, params={}):
        ki = [None for _ in range(self._sizeb)]
        t_0 = Constant(t)
        for i in range(self._sizeb):
            ki[i] = Function(y0)
            yi = Function(y0)
            t.assign(t_0 + dt * self._c[i])
            for j in range(i):
                if self._a[i,j] != 0:
                    yi += float(self._a[i,j]) * ki[j] * dt
            
            if isinstance(F, Form):
                F2 = inner(ki[i], test_f) * dx - F

                #if bcs is not None:
                #    bcs.apply(yi)
                solve(replace(F2, {y0: yi}) == 0,  ki[i],
                      bcs=bcs,
                      solver_parameters=params)
            else:
                ki[i] = F(t, yi)
        t.assign(t_0)
        return ki
    
    def k_dirk_fem(self, F, y0, t, dt, test_f, bcs=None, M = 1, params={}):
        ki = [None for _ in range(self._sizeb)]
        t0 = Constant(t)
        y_s = Function(y0)
        for i in range(self._sizeb):
            yi = Function(y_s)
            for j in range(i):
                if self._a[i,j] != 0:
                    yi += dt * float(self._a[i,j]) * ki[j]
            F2 = inner((y0 - yi) / (dt * float(self._a[i,i])), M * test_f) * dx - F
            t.assign(dt * self._c[i] + t0)
            solve(F2 == 0, y0, bcs=bcs, solver_parameters=params)
            ki[i] = Function(y_s)
            if self._a[i,i] != 0:
                ki[i].assign((y0 - yi) / dt / float(self._a[i,i]))
            else:
                F2=inner(ki[i], test_f) * dx - F
                solve(F2 == 0, ki[i])
        t.assign(t0)
        return ki

    def k_irk_fem(self, F, y0, t, dt, test_f, bcs=None, M = 1, params={}):
        Vbig = y0.function_space()
        N = len(y0.function_space())
        for i in range(1, self._sizeb):
            Vbig = y0.function_space() * Vbig
        test_b = TestFunction(Vbig)
        yis = Function(Vbig)

        Fnew = 0
        ys = [y0]
        ts = [Constant(t + dt * self._c[j]) for j in range(self._sizeb)]
        tfs = [test_f]
        if N > 1:
            ys = split(y0)
            tfs = split(test_f)
        for i in range(self._sizeb):
            for j in range(self._sizeb):
                rd = {t: ts[j]}
                for kk in range(N):
                    rd[ys[kk]] = split(yis)[j * N + kk]
                    rd[tfs[kk]] = split(test_b)[i * N + kk]
                Fnew += self._a[i, j] * replace(F, rd)
            for kk in range(N):
                if M == 1:
                    Fnew -= inner((split(yis)[i * N + kk] - ys[kk])/dt, M * split(test_b)[i * N + kk]) * dx
                else:
                    tf_i = 0
                    for kj in range(N):
                        if M[kk][kj] != 0:
                            tf_i += M[kk][kj] * split(test_b)[i * N + kj]
                    if tf_i != 0:
                        Fnew -= inner((split(yis)[i * N + kk] - ys[kk])/dt, tf_i)*dx
        new_bcs = []
        t0 = Constant(t)
        if bcs is not None:
            if isinstance(bcs, DirichletBC):
                bcs = [bcs]
            for bc in bcs:
                if N == 1:
                    c = bc.function_space().component
                    if c is not None:
                        Vbi = lambda i: Vbig[i].sub(c)
                    else:
                        Vbi = lambda i: Vbig[i]
                else:
                    s = bc.function_space_index()
                    c = bc.function_space().component
                    if c is not None:
                        Vbi = lambda i: Vbig[s + N * i].sub(c)
                    else:
                        Vbi = lambda i: Vbig[s + N * i]
                for j in range(self._sizeb):
                    t.assign(t0 + dt * self._c[j])
                    if bc.function_arg != 0:
                        new_bcs.append(DirichletBC(Vbi(j), bc.function_arg.copy(deepcopy=True), bc.sub_domain))
                    else:
                        new_bcs.append(DirichletBC(Vbi(j), 0, bc.sub_domain))
        solve(Fnew == 0, yis, bcs=new_bcs, solver_parameters=params)
        y_out = []

        for j in range(self._sizeb):
            y_i = Function(y0)
            if N == 1:
                y_i.assign(yis.sub(j*N + kk))
            else:
                for kk in range(N):
                    y_i.sub(kk).assign(yis.sub(j * N + kk))
            y_out.append(y_i)

        kis = []
        for j in range(self._sizeb):
            if self._Ainv is not None:
                kis.append(Function(y0))
                kis[-1].assign(0)
                for i in range(self._sizeb):
                    kis[-1].assign(kis[-1] + (y_out[i] - y0)/dt * self._Ainv[j, i])
            else:
                t.assign(t0 + dt * self._c[j])
                kis.append(Function(y0))
                F2 = inner(kis[j], test_f) * dx - replace(F, {y0: y_out[j]})
                solve(F2 == 0, kis[j])

        t.assign(t0)
        return kis

    def step_fem(self, F, y0, t, dt, bcs=None, **params):
        y_s = Function(y0)
        M = 1
        g = None
        if dt == 0:
            return y0
        if isinstance(F, tuple):
            M = F[1][0]
            g = F[1][1]
            F = F[0]
        if isinstance(F, Form):
            test_f = F.arguments()[0]
        else:
            test_f = None
        if self._explicit:
            ki = self.k_erk_fem(F, y0, t, dt, test_f, bcs, params)
        elif self._sdirk:
            ki = self.k_dirk_fem(F, y0, t, dt, test_f, bcs, M, params)
        else:
            ki = self.k_irk_fem(F, y0, t, dt, test_f, bcs, M, params)
        y0.assign(y_s)
        for i in range(self._sizeb):
            y_s += ki[i] * float(self._b[i])
        #if self._explicit:
            #if bcs is not None:
            #    bcs.apply(y_s)
        y_s.assign(y_s - y0)
        
        return y_s

    def y_step(self, function, current_y, current_x, step, **kwargs):
        """Return the adjustment on y for the next step"""
        if isinstance(function, tuple) or isinstance(current_y, Function):
            if isinstance(function, tuple):
                if len(function) == 3:
                    return self.step_fem((function[0], function[2]), current_y, current_x, step, bcs=function[1], **kwargs)
                if isinstance(function[1], DirichletBC):
                    return self.step_fem(function[0], current_y, current_x, step, bcs=function[1], **kwargs)
                else:
                    return self.step_fem(function, current_y, current_x, step, **kwargs)
            else:
                return self.step_fem(function, current_y, current_x, step, **kwargs)
        if self._explicit:
            return np.dot(self.get_k_values_explicit(function, current_y, current_x, step), self._b)
        elif self._sdirk:
            return np.dot(self.get_k_values_sdirk(function, current_y, current_x, step,  **kwargs), self._b)
        else:
            return np.dot(self.get_k_implicit(function, current_y, current_x, step, **kwargs), self._b)


class EmbeddedTableau(Tableau):
    """ Embedded RK pair to allow variable step sizes.
        The order is the lesser of the orders of the two methods, and is 
        used to calculate the next step size."""
    def __init__(self, c, a, b, b_aux, order, explicit=False, sdirk = False):
        super().__init__(c, a, b, explicit, sdirk)
        self.b_aux = b_aux
        self.order = order
    
    def y_step(self, function, current_y, current_x, step, rtol=1e-3, atol=1e-6, **kwargs):
        """ Returns the value of y after time step step.
        Step size is controlled to get a solution that meets the tolerance. """
        dt = 0
        tn = step
        while abs(dt - step) > 1e-8:
            if abs(tn) < 1e-10:
                if isinstance(current_y, Function):
                    return current_y.assign(np.nan)
                return current_y * np.nan
            if abs(dt + tn) > abs(step):
                tn = step - dt
            fem = False
            k = self.get_k_values(function, current_y, current_x, tn, **kwargs)
            if isinstance(k, list):
                fem = True
                y1 = Function(current_y)
                y2 = Function(current_y)
                for i in range(self._sizeb):
                    if self._b[i] != 0:
                        y1 += k[i] * tn * float(self._b[i])
                    if self.b_aux[i] != 0:
                        y2 += k[i] * tn * float(self.b_aux[i])
            else:
                y1 = current_y + tn * np.dot(k, self._b)
                y2 = current_y + tn * np.dot(k, self.b_aux)
            accept, err = measure_error(current_y, y1, y2, rtol, atol)
            
            if accept:
                dt += tn
                if fem:
                    current_y.assign(y2)
                    current_x.assign(current_x + tn)
                else:
                    current_y = y2
                    current_x += tn
            tn = compute_time(err, self.order, tn)
        return current_y




def implicit_root_func(Y, current_x, current_y, function, step, tableau):
    """ Function used for determining k_i in the fully implicit methods.  
        This wraps the equation to solve with the logic to interpret complex numbers as two floating point numbers if necessary and to ignore elements where the function returns nan"""
    if current_y.dtype == np.complex128:
        Y = Y.view(np.complex128)
    Y = Y.reshape((current_y.size, tableau._b.size), order='F')
    R = np.array(Y)
    
    for i in range(tableau._b.size):
        for j in range(tableau._b.size):
            R[:,i] -= step * tableau._a[i,j] * function(current_x + step * tableau._c[j], Y[:,j])
        R[:,i] -= current_y
    R = R.flatten(order='F')
    if current_y.dtype == np.complex128:
        R = R.view(np.float64)
    return R
    
    
def dirk_root_function(Y, current_x, current_y, function, step, k_i, i, tableau):
    """ Function used for determining Y_i in the SDIRK methods.  
        This wraps the equation to solve with the logic to interpret complex numbers as two floating point numbers if necessary and to ignore elements where the function returns nan"""
    if current_y.dtype == np.complex128:
        Y = Y.view(np.complex128)
    R = np.array(Y)
    R -= step * np.dot(k_i, tableau._a[i]) + step * tableau._a[i,i] * function(current_x + step * tableau._c[i], Y)
    R -= current_y

    if current_y.dtype == np.complex128:
        R = R.view(np.float64)

    return R
    
def constant_matrix(function, current_y, time):
    """ Assuming the system is linear y_dot=A*y+b, determines the matrix A"""
    # This is the FD approximation of the Jacobian for the LINEAR case only.
    size_y=np.size(current_y)
    a = sparse.lil_matrix((size_y,size_y),dtype=current_y.dtype) 
    y=current_y.copy()
    for i in range(size_y):
        y[i]+=1
        a[:,i]=  np.reshape(function(time, y)- function(time, current_y), (np.size(current_y),1))

        y[i]-=1
    return a.tocsr()

def create_error_function(sol_1, sol_2, rel_tol, abs_tol, k_factor):
    if sol_1.ufl_shape == ():
        return Function(sol_1).interpolate(k_factor * (sol_1 - sol_2) / (abs(sol_1) * rel_tol+ abs_tol))
    if len(split(sol_1)) == 1:
        try:
            return Function(sol_1).interpolate(k_factor * (sol_1 - sol_2) * (abs(sol_1) * rel_tol + abs_tol) ** -1)
        except:
            return Function(sol_1).interpolate(k_factor * (sol_1 - sol_2) / norm((abs(sol_1) * rel_tol + abs_tol)))
    out = Function(sol_1)


    N = len(split(sol_1))
    for i in range(N):
        out.sub(i).assign(create_error_function(sol_1.sub(i), sol_2.sub(i), rel_tol, abs_tol.sub(i), k_factor))
    return out

# utilities for embedded methods
def measure_error(y, sol_1, sol_2, rel_tol, abs_tol, k_factor=1):
    # tolerance is defined as rtol * sol_1 + atol
    # k_factor is used for palindromic or Milne methods
    if k_factor is None:
        k_factor = 1
    
    if isinstance(y, Function):
        if not isinstance(abs_tol, Function):
            abs_tol = Function(y).assign(rel_tol)
        E = create_error_function(sol_1, sol_2, rel_tol, abs_tol, k_factor)
        err = norm(E)
        return err < 1, err
    if isinstance(abs_tol, (int, float)):
        abs_tol = np.full_like(sol_1, abs_tol)
    assert len(abs_tol) == len(sol_1), "Error should have same length as solution"
    tol = rel_tol * abs(sol_1) + abs_tol
    err = (sol_1-sol_2)
    if k_factor is not None:
        err = k_factor*err
    return np.all(abs(err) < abs(tol)), np.linalg.norm(err/tol)


def compute_time(err, order, tn):
    a_min = 0.2
    a_max = 5
    a = 0.9
    if err == 0:
        tn = tn * a_max
    else:
        tn = tn * min(a_max, max(a_min, a * (1 / err)**(1/(order+1))))
    return tn

tableaus={'FE': Tableau(np.array([0]), np.array([[0]]), np.array([1]),explicit=True),
          'RK4': Tableau(np.array([0, 0.5, 0.5, 1]),
                         np.array([[0,0,0,0],
                                   [0.5,0,0,0],
                                   [0,0.5,0,0],
                                   [0,0,1,0]]),
                         np.array([1.0/6,1.0/3,1.0/3,1.0/6]), explicit=True),
          'Heun': Tableau(np.array([0,1]),
                          np.array([[0,0],
                                    [1,0]]),
                          np.array([0.5,0.5]),
                          explicit=True),
          'RK3': Tableau(np.array([0,0.5,1]),
                         np.array([[0,0,0],[0.5,0,0],[-1,2,0]]),
                         np.array([1./6,2./3,1./6]),explicit=True),
          'BE': Tableau(np.array([1]), np.array([[1]]), np.array([1]),sdirk=True),
          'GL4': Tableau(np.array([0.5-sqrt(3)/6, 0.5+sqrt(3)/6]),
                         np.array([[0.25, 0.25-sqrt(3)/6],
                                   [0.25+sqrt(3)/6, 0.25]]),
                         np.array([0.5, 0.5])),
          'GL2': Tableau(np.array([0.5]), np.array([[0.5]]), np.array([1])),
          'CN': Tableau(np.array([0, 1]), np.array([[0,0],[0.5,0.5]]),
                        np.array([0.5, 0.5])),
          'GL6': Tableau(np.array([0.5-sqrt(15)/10, 0.5, 0.5+sqrt(15)/10]),
                         np.array([[5./36, 2./9-sqrt(15)/15,
                                    5./36-sqrt(15)/30],
                                   [5./36+sqrt(15)/24, 2./9,
                                    5./36-sqrt(15)/24],
                                   [5./36+sqrt(15)/30, 2./9+sqrt(15)/15,
                                    5./36]]),
                         np.array([5./18, 4./9, 5./18])),
          'SSP(5,4)': Tableau(np.array([0, 0.39175222700392, 0.58607968896779,
                                        0.47454236302687, 0.93501063100924]),
                              np.array([[0, 0, 0, 0, 0],
                                        [0.39175222700392, 0, 0, 0, 0],
                                        [0.21766909633821, 0.36841059262959,
                                         0, 0, 0],
                                        [0.08269208670950, 0.13995850206999,
                                         0.25189177424738, 0, 0],
                                        [0.06796628370320, 0.11503469844438,
                                         0.20703489864929, 0.54497475021237,
                                         0]]),
                              np.array([0.14681187618661, 0.24848290924556,
                                        0.10425883036650, 0.27443890091960,
                                        0.22600748319395]), explicit=True),
          'SSPRK3': Tableau(np.array([0,1.,0.5]),
                         np.array([[0,0,0],[1.,0,0],[1./4,1./4,0]]),
                         np.array([1./6,1./6,2./3]),explicit=True),
}

tableaus['IM']=tableaus['GL2']
tableaus['IT']=tableaus['CN']
gamma=(2-sqrt(2))/2
tableaus['SD2O2']=Tableau(np.array([gamma, 1]),
                          np.array([[gamma, 0],
                                    [1-gamma, gamma]]),
                          np.array([1-gamma, gamma]), sdirk=True)   # Ascher and Petzold
gamma=2*cos(pi/18)/sqrt(3)
tableaus['SD3O4']=Tableau(np.array([(1+gamma)/2, 1./2, (1-gamma)/2]),
                          np.array([[(1+gamma)/2, 0, 0],
                                    [-gamma/2, (1+gamma)/2, 0],
                                    [(1+gamma), -(1+2*gamma), (1+gamma)/2]]),
                          np.array([1./(6*gamma**2), 1-1./(3*gamma**2), 1./(6*gamma**2)]),sdirk=True)
gamma=(3+sqrt(3))/6
tableaus['SD2O3']=Tableau(np.array([gamma, 1-gamma]),
                          np.array([[gamma, 0],
                                    [1-2*gamma, gamma]]),
                          np.array([0.5, 0.5]),sdirk=True)
gamma= 0.5
tableaus['SDAstable']=Tableau(np.array([gamma, 1-gamma]),
                          np.array([[gamma, 0],
                                    [1-2*gamma, gamma]]),
                          np.array([0.5, 0.5]),sdirk=True)
gamma= (2+sqrt(2))/2
tableaus['SDLstable']=Tableau(np.array([gamma, 1-gamma]),
                          np.array([[gamma, 0],
                                    [1-2*gamma, gamma]]),
                          np.array([0.5, 0.5]),sdirk=True)
gamma= 0.4358665215
b1gamma = -3.0/2.0*(gamma**2) + 4.0*gamma - 1.0/4.0
b2gamma = 3.0/2.0*(gamma**2) - 5.0*gamma + 5.0/4.0
tableaus['SD3O3Lstable']=Tableau(np.array([gamma, (1+gamma)/2, 1]),
                          np.array([[gamma, 0, 0],
                                    [(1-gamma)/2, gamma, 0],
                                    [b1gamma, b2gamma, gamma]]),
                          np.array([b1gamma, b2gamma, gamma]),sdirk=True)

tableaus['SD5O4'] = Tableau(np.array([1/4, 3/4, 11/20, 1/2, 1]),
                            np.array([[1/4, 0, 0, 0, 0],
                                      [1/2, 1/4, 0, 0, 0],
                                      [17/50, -1/25, 1/4, 0, 0],
                                      [371/1360, -137/2720, 15/544, 1/4, 0],
                                      [25/24, -49/48, 125/16, -85/12, 1/4]]),
                            np.array([25/24, -49/48, 125/16, -85/12, 1/4]))
#tableaus: 'FE': forward euler,
#          'RK4': classical Runge--Kutta of order 4,
#          'Heun': Heun's,
#          'RK3': Runge--Kutta order 3,
#          'BE': backward euler,
#          'GL6': Gauss-Legendre of order 6
#          'GL4': Gauss-Legendre of order 4,
#          'GL2', Gauss-Legendre of order 2 aka implicit midpoint ('IM'),
#          'CN': Crank-Nicolson aka implicit trapezoidal ('IT')
#          'SD2O3': 2 stage SDIRK method of order 3
#          'SD3O4': 3 stage SDIRK method of order 4
#          'SSP(5,4)'
#          'SSPRK3': Strong stability preserving RK3

embedded_pairs = {
        'Heun-Euler': EmbeddedTableau(np.array([0, 1]), np.array([[0, 0], [1, 0]]), np.array([1/2, 1/2]), np.array([1, 0]), 1),
        'Bogacki-Shampine': EmbeddedTableau(np.array([0, 1/2, 3/4, 1]), np.array([[0,0,0,0],[1/2,0,0,0], [0, 3/4, 0, 0], [2/9,1/3,4/9,0]]), np.array([2/9,1/3,4/9,0]), np.array([7/24,1/4,1/3,1/8]), 2),
        'Fehlberg': EmbeddedTableau(np.array([0,1/4,3/8,12/13,1,1/2]), np.array([[0,0,0,0,0,0],[1/4,0,0,0,0,0], [3/32, 9/32, 0,0,0,0], [1932/2197, -7200/2197, 7296/2197, 0,0,0], [439/216, -8, 3680/513, -845/4104, 0,0], [-8/27, 2, -3544/2565, 1859/4104, -11/40, 0]]), np.array([16/135, 0, 6656/12825, 28561/56430,-9/50, 2/55]), np.array([25/216, 0, 1408/2565,2197/4104, -1/5, 0]), 4),
        'Cash-Karp': EmbeddedTableau(np.array([0,1/5,3/10,3/5,1,7/8]), np.array([[0,0,0,0,0,0],[1/5,0,0,0,0,0], [3/40,9/40,0,0,0,0],[3/10,-9/10,6/5,0,0,0],[-11/54,5/2,-70/27,35/27,0,0],[1631/55296,175/512,575/13824,44275/110592,253/4096,0]]), np.array([37/378,0,250/621,125/594,0,512/1771]), np.array([2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4]), 4),
        'Dormand-Prince': EmbeddedTableau(np.array([0,1/5,3/10,4/5,8/9,1,1]), np.array([[0,0,0,0,0,0,0],[1/5,0,0,0,0,0,0],[3/40,9/40,0,0,0,0,0],[44/45,-56/15,32/9,0,0,0,0],[19372/6561, -25360/2187, 64448/6561, -212/729, 0,0,0], [9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,0], [35/384,0,500/1113, 125/192, -2187/6784,11/84, 0]]), np.array([35/384,0,500/1113,125/192,-2187/6784,11/84,0]), np.array([5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40]), 4)
        #'MERSON4': EmbeddedTableau()
}


