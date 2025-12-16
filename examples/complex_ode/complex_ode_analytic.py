# Complexe nonlinear ODE
# du/dt = i*u + 0.1*(1-u^2)u = i*u + 0.1*u -0.1*u^3
# u(0) = 0.1, t \in [0,100]
# u = x+i*y

import fractional_step as fs
import numpy as np

#parameter
a = 1j
ar = abs(a)
b = 0.1
c = -0.1

# time period
t0=0
tf= 100.0

# initial condition
u0 = 0.1+0j
u0_real = np.array([0.1, 0])


# operators:
# f = splitting of complex ODE
def f1(t, u):
    dudt = a * u
    return dudt
# Operator 2
def f2(t,dt,u):
    dudt = u * np.exp(b * dt)
    return dudt
# operator 3
def f3(t,u):
    dudt = c*(u**3)
    #print('dudt = ', dudt)
    return dudt

# g = splitting of as system of real ODE
def g1(t, u):
    x = u[0]
    y = u[1]
    dxdt = -abs(a)*y
    dydt = abs(a)*x
    dudt = np.append(dxdt, dydt)
    return dudt
# Operator 2
def g2(t,dt,u):
    x = u[0]
    y = u[1]
    dxdt = x * np.exp(b * dt)
    dydt = y * np.exp(b * dt)
    dudt = np.append(dxdt, dydt)
    return dudt
# operator 3
def g3(t,u):
    x = u[0]
    y = u[1]
    dxdt = c * (x**3-3*x*(y**2))
    dydt = c * (-y**3+3*(x**2)*y)
    dudt = np.append(dxdt, dydt)
    return dudt


# list the operators in the order you wish to use them
f_list=[f1, f2, f3]
g_list=[g1, g2, g3]


# splitting method
alpha='Strang-N'
# runge--kutta methods
Intg_method_1 = 'RK3'
Intg_method_2 = "ANALYTIC"
Intg_method_3 = Intg_method_1
methods={(0,):'',   #default for all operators,
         (1,):Intg_method_1,
         (2,):Intg_method_2,
         (3,):Intg_method_3,
}

alpha = 'Strang-N'
# solve the complex system
print('Solve complex ODE ', alpha, ' with ', methods) #, ' ivp_method ', ivp_method_f1, ' rtol = ', rtol_f1, ' atol = ', atol_f1)

dt = tf/100000
fname = "complex_solve_"+str(alpha)+"_" + str(t0) + "_to_" + str(tf) + "+" + str(methods[(1,)]) + "+" + str(methods[(2,)]) + "+" + str(methods[(3,)]) + "_dt_" + str(
    dt) + ".csv"
num_of_saved_steps = 100
result = fs.fractional_step(f_list, dt, u0, t0, tf, alpha, methods, fname=fname,save_steps=num_of_saved_steps)
print(result)

# solve the real system
print('Solve system of real ODE ', alpha, ' with ', methods)

dt = tf/100000
gname = "real_solve_"+str(alpha)+"_" + str(t0) + "_to_" + str(tf) + "+" + str(methods[(1,)]) + "+" + str(methods[(2,)]) + "+" + str(methods[(3,)]) + "_dt_" + str(
    dt) + ".csv"
num_of_saved_steps = 100
result = fs.fractional_step(g_list, dt, u0_real, t0, tf, alpha, methods, fname=gname, save_steps=num_of_saved_steps)

print(result)