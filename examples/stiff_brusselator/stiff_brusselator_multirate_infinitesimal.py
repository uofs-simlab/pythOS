#mri_esdirk3a
#mri_imex3

# There are a collectioon of multirate infinitesimal methods already defined
# which we simply import and use

# import problem definition
from stiff_brusselator import *
dt = tf/4

# solve the problem
from multirate_infinitesimal import multirate_infinitesimal_solve, mri_imex3, mri_esdirk3a, mri_irk2

# solve using MRI-IRK2
result = multirate_infinitesimal_solve(y0, t0, dt, tf, mri_irk2, fs, ff)
print(result)

# solve using MRI-ESDIRK3A
result = multirate_infinitesimal_solve(y0, t0, dt, tf, mri_esdirk3a, fs, ff)
print(result)

# solve using MRI-IMEX3
result = multirate_infinitesimal_solve(y0, t0, dt, tf, mri_imex3, fi, ff, fe)
print(result)
