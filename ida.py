from sundials import *

ida = sundialsCDLL("libsundials_ida.so.6")

ida.IDACreate.argtypes = [c_void_p]
ida.IDACreate.restype = c_void_p
ida.IDACreate.errcheck = memory_create_check

ida.IDAInit.argtypes = [c_void_p, c_void_p, c_double, c_void_p, c_void_p]
ida.IDASStolerances.argtypes = [c_void_p, c_double, c_double]
ida.IDASVtolerances.argtypes = [c_void_p, c_double, c_void_p]
ida.IDASetLinearSolver.argtypes = [c_void_p, c_void_p, c_void_p]
ida.IDASolve.argtypes = [c_void_p, c_double, POINTER(c_double), c_void_p, c_void_p, c_int]

ida.IDAReInit.argtypes = [c_void_p, c_double, c_void_p, c_void_p]
ida.IDACalcIC.argtypes = [c_void_p, c_int, c_double]

ida.IDAFree.restype = None
class IDA(SundialsSolver):
    def __init__(self, fn, y0, rtol, atol, t0, ydot0 = None, id = None, jac=None, **kwargs):
        super().__init__(y0, linear_solver = True)

        if ydot0 is None:
            ydot0 = np.zeros(y0.size, dtype=y0.dtype)
        self.udot = nvector.N_VMake_Serial(y0.size, ydot0.ctypes.data_as(POINTER(c_double)), self.ctx)

        ida_mem = ida.IDACreate(self.ctx)

        def f_c(t, yy, yp, rr, data):
            yi = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(yy), dtype=np.float64, count=self.N)).view(self.dtype)
            ydi = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(yp), dtype = np.float64, count=self.N)).view(self.dtype)
            res = fn(t, yi, ydi).view(np.float64)
            r = nvector.N_VGetArrayPointer_Serial(rr)
            for i in range(self.N):
                r[i] = res[i]
            return 0
        F = CFUNCTYPE(c_int, c_double, c_void_p, c_void_p, c_void_p, c_void_p)(f_c)
        ida.IDAInit(ida_mem, F, t0, self.u, self.udot)

        if isinstance(atol, float):
            ida.IDASStolerances(ida_mem, rtol, atol)
        else:
            atol_v = nvector.N_VMake_Serial(self.N, atol.ctypes.data_as(POINTER(c_double)), self.ctx)
            ida.IDASVtolerances(ida_mem, rtol, atol_v)
            self.atol_v = atol_v


        ida.IDASetLinearSolver(ida_mem, self.LS, None)

        if id is not None:
            if self.dtype == complex:
                id = id * (1 + 1j)
            ids = nvector.N_VMake_Serial(self.N, id.ctypes.data_as(POINTER(c_double)), self.ctx)
            ida.IDASetId(ida_mem, ids)
            self.id = ids
        else:
            self.id = None

        if jac is not None:
            def jac_c(tt, yy, yp, rr, v, Jv, cj, data, t1, t2):
                y = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(yy), dtype=np.float64, count=self.N)).view(self.dtype)
                yd = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(yp), dtype = np.float64, count=self.N)).view(self.dtype)
                vv = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(v), dtype = np.float64, count=self.N)).view(self.dtype)
                result = jac(tt, y, yd, vv, cj).view(np.float64)
                jv = nvector.N_VGetArrayPointer_Serial(Jv)
                for i in range(self.N):
                    jv[i] = result[i]
                return 0

            JF = CFUNCTYPE(c_int,
                                              c_double, #t
                                              c_void_p, #y
                                              c_void_p, #yp
                                              c_void_p, #r
                                              c_void_p, #v
                                              c_void_p, #Jv
                                              c_double, #cj
                                              c_void_p, #data
                                              c_void_p, #t1
                                              c_void_p)(jac_c)
            self.JF = JF
            ida.IDASetJacTimes(ida_mem, None, JF)

        self.F = F
        self.ida_mem = c_void_p(ida_mem)

    def solve(self, y0, t0, tf, J):
        self.J = J
        yi = nvector.N_VGetArrayPointer_Serial(self.u)
        yn = y0.view(np.float64)
        for i in range(self.N):
            yi[i] = yn[i]


        #ydot_i = nvector.N_VGetArrayPointer_Serial(self.udot)
        #for i in range(self.N):
        #    ydot_i[i] = ydot0[i]


        ida.IDAReInit(self.ida_mem, t0, self.u, self.udot)

        if self.id is not None:
            ida.IDACalcIC(self.ida_mem, 1, tf)
        t_ret = c_double(0.0)

        ida.IDASolve(self.ida_mem, tf, byref(t_ret), self.u, self.udot, 1)
        

        return np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(self.u), dtype=np.float64, count=self.N)).view(self.dtype)

    def free(self):
        assert self.ida_mem is not None
        nvector.N_VDestroy(self.udot)
        super().free(self.ida_mem, ida.IDAFree)
        self.ida_mem = None
