from firedrake import *
import math

norm_phi = []
time_step = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125]

for dt in time_step:

    N =16

    mesh = RectangleMesh(N, N, 5.0, 5.0) 

    V = FunctionSpace(mesh, "CG", 2)  

    dof = V.dim()
    print("DOF = ", dof)

    phi = Function(V)
    v = TestFunction(V)

    D = 1.0/30
    Sigma_a = 1.0
    A = 10
    omega = pi/10
    b = 5.0

    t = Constant(0.0)

    x, y = SpatialCoordinate(mesh)
    S = A * (omega * (b*x - x**2) * cos(omega* t) + sin(omega* t) * (Sigma_a * (b* x - x**2)+ 2 * D))

    n = FacetNormal(mesh)
    F1 = D * inner(inner(grad(phi), n), v) * (ds(1) + ds(2)) - D * inner(grad(phi), grad(v)) * dx
    F2 = (- inner(phi, v)* Sigma_a + inner(S, v)) * dx

    bc_D = DirichletBC(V, 0.0, [1, 2])  

    phi.sub(0).assign(0)

    import sys
    sys.path.insert(1, sys.path[0]+"/../")
    import operator_splitting as os

    T = 10.0
    #dt = 0.1


    os.operator_splitting([F1,F2], dt, phi, t, T, "Strang", methods={(0,): '',(1,): "SD2O2", (2,): "Heun"} , bc = bc_D)

    norm_phi.append(norm(phi, "L2"))
    print("L2 norm of phi:", norm_phi, ", dt = ", dt)
    

orders = []
orders.append(None)
for i in range(len(norm_phi)-1):
    orders.append(math.log(norm_phi[i]/norm_phi[i+1],10)/math.log(time_step[i]/time_step[i+1],10))

print("{:<15} {:<15} {:<10}".format("dt", "Error", "Order"))
for dt, e, o in zip(time_step, norm_phi, orders):
    if o is None:
        print("{:<15.4f} {:<15.3e}".format(dt, e))
    else:
        print("{:<15.4f} {:<15.3e} {:<10.3f}".format(dt, e, o))