from ctypes import *
import numpy as np

class sundialsCDLL(CDLL):
    def __getattr__(self, name):
        result = super().__getattr__(name)
        if result.errcheck is None:
            result.errcheck = error_check
        return result
sundials = sundialsCDLL("libsundials_core.so")
nvector = sundialsCDLL("libsundials_nvecserial.so")
linsol = sundialsCDLL("libsundials_sunlinsolspgmr.so")
wrapper = sundialsCDLL("sundials_wrapper.so")

def memory_create_check(result, func, arguments):
    if result is None:
        raise RuntimeError("Failed to create object in {}".format(func.__name__))
    return result

def memory_access_check(result, func, arguments):
    if result is None:
        raise RuntimeError("Failed to access object in {}".format(func.__name__))
    return result

def error_check(result, func, arguments):
    if result is None:
        return result
    if result != 0:
        raise RuntimeError("Error {} in function {}".format(result, func.__name__))
    return result


wrapper.allocate_context.restype = c_void_p
wrapper.allocate_context.errcheck = memory_create_check
wrapper.allocate_inner_stepper.restype = c_void_p
wrapper.allocate_inner_stepper.errcheck = memory_create_check

sundials.SUNContext_Create.argtypes = [c_void_p, c_void_p]
sundials.SUNContext_Create.restype = c_int
sundials.SUNContext_Free.argtypes = [c_void_p]
sundials.SUNContext_Free.restype = None

nvector.N_VMake_Serial.argtypes = [c_int, c_void_p, c_void_p]
nvector.N_VMake_Serial.restype = c_void_p
nvector.N_VMake_Serial.errcheck = memory_create_check
nvector.N_VGetArrayPointer_Serial.argtypes = [c_void_p]
nvector.N_VGetArrayPointer_Serial.restype = POINTER(c_double)
nvector.N_VGetArrayPointer_Serial.errcheck = memory_access_check
nvector.N_VDestroy.argtypes = [c_void_p]
nvector.N_VDestroy.restype = None

linsol.SUNLinSol_SPGMR.argtypes = [c_void_p, c_int, c_int, c_void_p]
linsol.SUNLinSol_SPGMR.restype = c_void_p
linsol.SUNLinSol_SPGMR.errcheck = memory_create_check
linsol.SUNLinSolFree.argtypes = [c_void_p]

class SundialsSolver:
    def __init__(self, y0, linear_solver = True):
        ctx = wrapper.allocate_context()
        sundials.SUNContext_Create(None, ctx)
        self.N = y0.size

        if y0.dtype == complex:
            self.N = self.N * 2

        u = nvector.N_VMake_Serial(self.N, y0.ctypes.data_as(POINTER(c_double)), ctx)
        if linear_solver:
            LS = linsol.SUNLinSol_SPGMR(u, 0, 20, ctx)
        else:
            LS = None
        self.ctx = ctx
        self.u = u
        self.LS = LS
        self.dtype = y0.dtype
        self.J = None
        
    def free(self, solver_mem, solver_free):
        nvector.N_VDestroy(self.u)
        if self.LS is not None:
            linsol.SUNLinSolFree(self.LS)
        solver_free(byref(solver_mem))
        sundials.SUNContext_Free(self.ctx)

    def to_sundials_function(self, f):
        def f_c(t, y, ydot, data):
            yi = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(y), dtype=np.float64, count=self.N))
            yd = nvector.N_VGetArrayPointer_Serial(ydot)
            if self.J is None:
                result = f(t, yi.view(self.dtype)).view(np.float64)
            else:
                result = f(t, yi.view(self.dtype), self.J).view(np.float64)
            for i in range(self.N):
                yd[i] = result[i]
            return 0
        return f_c

            
