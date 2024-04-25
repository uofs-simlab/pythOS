from cvode import CVODE
from arkode import ARKStep, ERKStep
import numpy as np
class SundialsExode:
    def __init__(self, y0, method, rtol, atol, concatenate=True, **options):
        self.f = lambda t, y: y
        if concatenate:
            self.initial_y = np.concatenate((y0, [0])).copy()
        else:
            self.initial_y = y0.copy()
        def fun(t, y):
            return self.f(t, y)
        if method in ['CV_ADAMS', 'CV_BDF']:
            self.solver = CVODE(method, self.initial_y, fun, 0, rtol, atol, **options)
        elif 'ARKODE' in method:
            self.solver = ERKStep(method, self.initial_y, fun, 0, rtol, atol, **options)
        else:
            self.solver = ARKStep(method, self.initial_y, None, fun, 0, rtol, atol, **options)
    def solve(self, f, y0, t0, tf):
        self.f = f
        return self.solver.solve(y0, t0, tf, J = None)
    def free(self):
        self.solver.free()
        
            
        
