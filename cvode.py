from ctypes import *
import numpy as np

from sundials import *
if os.name == 'nt':
    sundials = sundialsCDLL("sundials_cvode")
else:
    sundials = sundialsCDLL("libsundials_cvode.so")


sundials.CVodeCreate.argtypes = [c_int, c_void_p]
sundials.CVodeCreate.restype = c_void_p
sundials.CVodeCreate.errcheck = memory_create_check

sundials.CVodeInit.argtypes = [c_void_p, c_void_p, c_double, c_void_p]
sundials.CVodeInit.restype = c_int

sundials.CVodeSStolerances.argtypes = [c_void_p, c_double, c_double]
sundials.CVodeSStolerances.restype = c_int
sundials.CVodeSVtolerances.argtypes = [c_void_p, c_double, c_void_p]

sundials.CVodeSetLinearSolver.argtypes = [c_void_p, c_void_p, c_void_p]
sundials.CVodeSetLinearSolver.restype = c_int

sundials.CVode.argtypes = [c_void_p, c_double, c_void_p, c_void_p, c_int]
sundials.CVode.restype = c_int

sundials.CVodeReInit.argtypes = [c_void_p, c_double, c_void_p]
sundials.CVodeReInit.restype = c_int

sundials.CVodeFree.argtypes = [c_void_p]
sundials.CVodeFree.restype = None

sundials.CVodeSetMaxNumSteps.argtypes = [c_void_p,c_long]
sundials.CVodeSetMaxNumSteps.restype = c_int

sundials.CVodeSetJacTimes.argtypes = [c_void_p, c_void_p, c_void_p]
sundials.CVodeSetJacTimes.restyp = c_int

class CVODE(SundialsSolver):
    def __init__(self, ode_type, y0, f, t0, rtol, atol, jac=None, max_steps=0,**kwargs):
        super().__init__(y0)
        if ode_type == 'CV_ADAMS':
            ode_type = 1
        else:
            ode_type = 2
        u = self.u
        LS = self.LS
        ctx = self.ctx
        cvode_mem = c_void_p(sundials.CVodeCreate(ode_type, ctx))

        f_cv = self.to_sundials_function(f)
        
        CF = CFUNCTYPE(c_void_p, c_double, c_void_p, c_void_p, c_void_p)(f_cv)
        sundials.CVodeInit(cvode_mem, CF, t0, u)
        if isinstance(atol, float):
            sundials.CVodeSStolerances(cvode_mem, rtol, atol)
        else:
            atol_v = nvector.N_VMake_Serial(self.N, atol.ctypes.data_as(POINTER(c_double)), self.ctx)
            sundials.CVodeSVtolerances(cvode_mem, rtol, atol_v)
            self.atol_v = atol_v

        sundials.CVodeSetLinearSolver(cvode_mem, LS, None)

        if jac is not None:
            def jac_c(v, jv, tt, yy, fy, data, tmp):
                y = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(yy), dtype=np.float64, count=self.N)).view(self.dtype)
                vv = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(v), dtype=np.float64, count=self.N)).view(self.dtype)
                result = jac(tt, y, vv).view(np.float64)
                jv = nvector.N_VGetArrayPointer_Serial(Jv)
                for i in range(self.N):
                    jv[i] = result[i]
                return 0
            JF = CFUNCTYPE(c_int, c_void_p, c_void_p, c_double, c_void_p, c_void_p, c_void_p, c_void_p)(jac)
            self.JF = JF
            sundials.CVodeSetJacTimes(cvode_mem, None, JF)
        
        sundials.CVodeSetMaxNumSteps(cvode_mem, max_steps)
        
        self.ctx = ctx
        self.LS = LS
        self.u = u
        self.cvode_mem = cvode_mem
        self.CF = CF

    def solve(self, y0, t0, tf, J = None):
        self.J = J
        yi = nvector.N_VGetArrayPointer_Serial(self.u)
        yn = y0.view(dtype=np.float64)
        for i in range(self.N):
            yi[i] = yn[i]
        sundials.CVodeReInit(self.cvode_mem, t0, self.u)

        t_out = c_double(0.0)
        try:
            sundials.CVode(self.cvode_mem, tf, self.u, byref(t_out), 1)
        except:
            print("CVODE integration failed")
            return y0 * np.nan


        return np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(self.u), dtype=np.float64, count=self.N)).view(dtype=self.dtype)

    def free(self):
        assert (self.cvode_mem is not None)
        super().free(self.cvode_mem, sundials.CVodeFree)
        self.cvode_mem = None
