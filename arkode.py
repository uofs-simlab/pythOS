from ctypes import *
import numpy as np
from sundials import *
if os.name == 'nt':
    arkode = sundialsCDLL("sundials_arkode.dll")
else:
    arkode = sundialsCDLL("libsundials_arkode.so")

arkode.ARKStepCreate.argtypes = [c_void_p, c_void_p, c_double, c_void_p, c_void_p]
arkode.ARKStepCreate.restype = c_void_p
arkode.ARKStepCreate.errcheck = memory_create_check
arkode.ARKodeSStolerances.argtypes = [c_void_p, c_double, c_double]
arkode.ARKodeSStolerances.restype = c_int
arkode.ARKodeSVtolerances.argtypes = [c_void_p, c_double, c_void_p]
arkode.ARKodeSVtolerances.restype = c_int
arkode.ARKodeSetLinearSolver.argtypes = [c_void_p, c_void_p, c_void_p]
arkode.ARKodeSetLinearSolver.restype = c_int
arkode.ARKStepSetTableName.argtypes = [c_void_p, c_char_p, c_char_p]
arkode.ARKStepSetTableName.restype = c_int
arkode.ARKodeReset.argtypes = [c_void_p, c_double, c_void_p]
arkode.ARKodeReset.restype = c_int
arkode.ARKodeEvolve.argtypes = [c_void_p, c_double, c_void_p, POINTER(c_double), c_int]
arkode.ARKodeEvolve.restype = c_int
arkode.ARKodeFree.argtypes = [c_void_p]
arkode.ARKodeFree.restype = None
arkode.ARKodeSetMaxNumSteps.argtypes = [c_void_p,c_long]
arkode.ARKodeSetMaxNumSteps.restype = c_int

arkode.ERKStepCreate.argtypes = [c_void_p, c_double, c_void_p, c_void_p]
arkode.ERKStepCreate.restype = c_void_p
arkode.ERKStepCreate.errcheck = memory_create_check
arkode.ERKStepSetTableName.argtypes = [c_void_p, c_char_p]
arkode.ERKStepSetTableName.restype = c_int

arkode.ARKodeCreateMRIStepInnerStepper.argtypes = [c_void_p, c_void_p]
arkode.ARKodeCreateMRIStepInnerStepper.restype = c_int
arkode.MRIStepCreate.argtypes = [c_void_p, c_void_p, c_double, c_void_p, c_void_p, c_void_p]
arkode.MRIStepCreate.restype = c_void_p
arkode.MRIStepCreate.errcheck = memory_create_check
arkode.ARKodeSetFixedStep.argtypes = [c_void_p, c_double]
arkode.ARKodeSetFixedStep.restype = c_int
arkode.MRIStepSStolerances.argtypes = [c_void_p, c_double, c_double]
arkode.MRIStepSetLinearSolver.argtypes = [c_void_p, c_void_p, c_void_p]
arkode.MRIStepReset.argtypes = [c_void_p, c_double, c_void_p]
arkode.MRIStepEvolve.argtypes = [c_void_p, c_double, c_void_p, POINTER(c_double), c_int]
arkode.MRIStepInnerStepper_Free.argtypes = [c_void_p]
arkode.MRIStepInnerStepper_Free.restype = None
arkode.MRIStepFree.argtypes = [c_void_p]
arkode.MRIStepFree.restype = None

class ARKStep(SundialsSolver):
    def __init__(self, ode_type, y0, fe, fi, t0, rtol, atol, jac=None, max_steps=0,**kwargs):
        super().__init__(y0, linear_solver = fi)
        fe_c = self.to_sundials_function(fe)
        fi_c = self.to_sundials_function(fi)

        if fe is not None:
            FE = CFUNCTYPE(c_void_p, c_double, c_void_p, c_void_p, c_void_p)(fe_c)
        else:
            FE = None
        if fi is not None:
            FI = CFUNCTYPE(c_void_p, c_double, c_void_p, c_void_p, c_void_p)(fi_c)
        else:
            FI = None

        arkode_mem = arkode.ARKStepCreate(FE, FI, t0, self.u, self.ctx)

        if isinstance(atol, float):
            arkode.ARKodeSStolerances(arkode_mem, rtol, atol)
        else:
            atol_v = nvector.N_VMake_Serial(self.N, atol.ctypes.data_as(POINTER(c_double)), self.ctx)
            arkode.ARKodeSVtolerances(arkode_mem, rtol, atol_v)
            self.atol_v = atol_v


        if fi is not None:
            arkode.ARKodeSetLinearSolver(arkode_mem, self.LS, None)

        if fi is None:
            itable = "ARKODE_DIRK_NONE"
        else:
            itable = ode_type[0]
        if fe is None:
            etable = "ARKODE_ERK_NONE"
        else:
            etable = ode_type[1]

        itable = itable.encode()
        etable = etable.encode()
        arkode.ARKStepSetTableName(arkode_mem, itable, etable)
        arkode.ARKodeSetMaxNumSteps(arkode_mem, max_steps)

        if jac is not None:
            def jac_c(v, Jv, t, yy, fy, data, tmp):
                y = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(yy), dtype = np.float64, count=self.N)).view(self.dtype)
                vv = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(v), dtype = np.float64, count=self.N)).view(self.dtype)
                result = jac(t, y, vv).view(np.float64)
                jv = nvector.N_VGetArrayPointer_Serial(Jv)
                for i in range(self.N):
                    jv[i] = result[i]
                return 0
            JF = CFUNCTYPE(c_int,
                           c_void_p,
                           c_void_p,
                           c_double,
                           c_void_p,
                           c_void_p,
                           c_void_p,
                           c_void_p)(jac_c)
            self.JF = JF
            arkode.ARKStepSetJacTimes(arkode_mem, None, JF)

        self.arkode_mem = c_void_p(arkode_mem)
        self.FI = FI
        self.FE = FE

    def solve(self, y0, t0, tf, J):
        self.J = J
        yi = nvector.N_VGetArrayPointer_Serial(self.u)
        yn = y0.view(dtype=np.float64)
        for i in range(self.N):
            yi[i] = yn[i]
        arkode.ARKodeReset(self.arkode_mem, t0, self.u)

        t_out = c_double(0.0)
        try:
            arkode.ARKodeEvolve(self.arkode_mem, tf, self.u, byref(t_out), 1)
        except:
            print("ARKODE integration failed")
            return y0 * np.nan
        return np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(self.u), dtype = np.float64, count=self.N)).view(self.dtype)

    def free(self):
        assert self.arkode_mem is not None
        super().free(self.arkode_mem, arkode.ARKodeFree)
        self.arkode_mem = None

class ERKStep(SundialsSolver):
    def __init__(self, ode_type, y0, f, t0, rtol, atol, max_steps=0,**kwargs):
        super().__init__(y0, linear_solver=None)

        f_c = self.to_sundials_function(f)

        F = CFUNCTYPE(c_void_p, c_double, c_void_p, c_void_p, c_void_p)(f_c)
        arkode_mem = arkode.ERKStepCreate(F, t0, self.u, self.ctx)

        if isinstance(atol, float):
            arkode.ARKodeSStolerances(arkode_mem, rtol, atol)
        else:
            atol_v = nvector.N_VMake_Serial(self.N, atol.ctypes.data_as(POINTER(c_double)), self.ctx)
            arkode.ARKodeSVtolerances(arkode_mem, rtol, atol_v)
            self.atol_v = atol_v

        arkode.ERKStepSetTableName(arkode_mem, ode_type.encode())

        arkode.ARKodeSetMaxNumSteps(arkode_mem, max_steps)

        self.arkode_mem = c_void_p(arkode_mem)
        self.F = F


    def solve(self, y0, t0, tf, J):
        self.J = J
        yi = nvector.N_VGetArrayPointer_Serial(self.u)
        yn = y0.view(np.float64)
        for i in range(self.N):
            yi[i] = yn[i]
        arkode.ARKodeReset(self.arkode_mem, t0, self.u)

        t_out = c_double(0.0)
        try:
            arkode.ARKodeEvolve(self.arkode_mem, tf, self.u, byref(t_out), 1)
        except:
            print("ERKStep integration failed")
            return y0 * np.nan
        return np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(self.u), dtype=np.float64, count=self.N)).view(self.dtype)

    def free(self):
        assert self.arkode_mem is not None
        super().free(self.arkode_mem, arkode.ARKodeFree)
        self.arkode_mem = None
        
class MRIStep(SundialsSolver):
    def __init__(self, y0, ff, t0, rtol, atol, fs, dt, rtol_slow=1e-4, atol_slow=1e-9, slow_jac = None, fast_jac=None, 
                 max_steps_slow=0, max_steps_fast=0,**kwargs):
        super().__init__(y0)
        ff_i, ff_e = ff
        fs_i, fs_e = fs

        ffe_c = self.to_sundials_function(ff_e)
        
        # create inner stepper with ARKStepCreate
        ffi_c = self.to_sundials_function(ff_i)

        if ff_e is not None:
            FFE = CFUNCTYPE(c_void_p, c_double, c_void_p, c_void_p, c_void_p)(ffe_c)
        else:
            FFE = None
        if ff_i is not None:
            FFI = CFUNCTYPE(c_void_p, c_double, c_void_p, c_void_p, c_void_p)(ffi_c)
        else:
            FFI = None
        arkode_mem = arkode.ARKStepCreate(FFE, FFI, t0, self.u, self.ctx)
        if isinstance(atol, float):
            arkode.ARKodeSStolerances(arkode_mem, rtol, atol)
        else:
            atol_v = nvector.N_VMake_Serial(self.N, atol.ctypes.data_as(POINTER(c_double)), self.ctx)
            arkode.ARKodeSVtolerances(arkode_mem, rtol, atol_v)
            self.atol_v = atol_v

        LS_F = linsol.SUNLinSol_SPGMR(self.u, 0, 20, self.ctx)
        if ff_i is not None:
            arkode.ARKodeSetLinearSolver(arkode_mem, LS_F, None)
        arkode.ARKodeSetMaxNumSteps(arkode_mem, max_steps_fast)
        
        self.FFI = FFI
        self.FFE = FFE
        self.arkode_inner = c_void_p(arkode_mem)
        self.LS_F = LS_F
        if fast_jac is not None:
            def fast_jac_c(v, Jv, t, yy, fy, data, tmp):
                y = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(yy), dtype = np.float64, count=self.N)).view(self.dtype)
                vv = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(v), dtype = np.float64, count=self.N)).view(self.dtype)
                result = fast_jac(t, y, vv).view(np.float64)
                jv = nvector.N_VGetArrayPointer_Serial(Jv)
                for i in range(self.N):
                    jv[i] = result[i]
                return 0
            JF = CFUNCTYPE(c_int,
                           c_void_p,
                           c_void_p,
                           c_double,
                           c_void_p,
                           c_void_p,
                           c_void_p,
                           c_void_p)(fast_jac_c)
            self.fast_JF = JF
            arkode.ARKStepSetJacTimes(self.arkode_inner, None, self.fast_JF)


        inner_stepper = c_void_p(wrapper.allocate_inner_stepper())

        arkode.ARKodeCreateMRIStepInnerStepper(self.arkode_inner, byref((inner_stepper)))

        self.inner_stepper = inner_stepper
        fse_c = self.to_sundials_function(fs_e)
        fsi_c = self.to_sundials_function(fs_i)

        if fs_e is not None:
            FSE = CFUNCTYPE(c_void_p, c_double, c_void_p, c_void_p, c_void_p)(fse_c)
        else:
            FSE = None
        if fs_i is not None:
            FSI = CFUNCTYPE(c_void_p, c_double, c_void_p, c_void_p, c_void_p)(fsi_c)
        else:
            FSI = None

        mri_mem = arkode.MRIStepCreate(FSE, FSI, t0, self.u, (inner_stepper), self.ctx)


        arkode.ARKodeSetFixedStep(mri_mem, dt)

        if fs_i is not None:
            if isinstance(atol_slow, float):
                arkode.ARKodeSStolerances(mri_mem, rtol_slow, atol_slow)
            else:
                atol_sv = nvector.N_VMake_Serial(self.N, atol_slow.ctypes.data_as(POINTER(c_double)), self.ctx)
                arkode.ARKodeSVtolerances(mri_mem, rtol, atol_sv)
                self.atol_sv = atol_sv

            arkode.ARKodeSetLinearSolver(mri_mem, self.LS, None)
            if slow_jac is not None:
                def slow_jac_c(v, Jv, t, yy, fy, data, tmp):
                    y = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(yy), dtype = np.float64, count=self.N)).view(self.dtype)
                    vv = np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(v), dtype = np.float64, count=self.N)).view(self.dtype)
                    result = slow_jac(t, y, vv).view(np.float64)
                    jv = nvector.N_VGetArrayPointer_Serial(Jv)
                    for i in range(self.N):
                        jv[i] = result[i]
                    return 0
                JF = CFUNCTYPE(c_int,
                               c_void_p,
                               c_void_p,
                               c_double,
                               c_void_p,
                               c_void_p,
                               c_void_p,
                               c_void_p)(slow_jac_c)
                self.slow_JF = JF
                arkode.MRIStepSetJacTimes(mri_mem, None, self.slow_JF)
            
        self.FSI = FSI
        self.FSE = FSE
        self.mri_mem = c_void_p(mri_mem)

        ## optional inputs - controlling method

    def solve(self, y0, t0, tf, J):
        self.J = J
        yi = nvector.N_VGetArrayPointer_Serial(self.u)
        yn = y0.view(np.float64)
        for i in range(self.N):
            yi[i] = y0[i]

        arkode.ARKodeReset(self.mri_mem, t0, self.u)

        t_out = c_double(0.0)
        try:
            arkode.ARKodeEvolve(self.mri_mem, tf, self.u, t_out, 1)
        except:
            print("MRIStep failed")
            return y0 * np.nan

        return np.array(np.fromiter(nvector.N_VGetArrayPointer_Serial(self.u), dtype=np.float64, count=self.N)).view(self.dtype)

    def free(self):
        assert self.ctx is not None
        nvector.N_VDestroy_Serial(self.u)
        arkode.ARKodeFree(byref(self.arkode_inner))
        
        arkode.MRIStepInnerStepper_Free(byref(self.inner_stepper))
        arkode.ARKodeFree(byref(self.mri_mem))
        if self.LS is not None:
            linsol.SUNLinSolFree(self.LS)
        if self.LS_F is not None:
            linsol.SUNLinSolFree(self.LS_F)
        sundials.SUNContext_Free(self.ctx)
        self.ctx = None
