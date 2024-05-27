from scipy.integrate import solve_ivp
from scipy.optimize import root
import numpy as np
try:
    from firedrake import Function, dx, solve, Constant, replace, inner, CheckpointFile
except:
    Function = type(None)
    Constant = type(None)
try:
    from sundials_exode import SundialsExode
except:
    SundialsExode = type(None)
from butcher_tableau import embedded_pairs, EmbeddedTableau
class Multirate:
    def __init__(self, c, gamma, omega=None):
        self.c = c
        self.gamm = gamma
        self.omeg = omega
        self.s = len(c)

    def step(self, yn, dt, fi, ff, tn, fe=None, bcs=None, ivp_options={}, ivp_method = 'RK45', implicit_solve_options={}):
        if isinstance(yn, Function) and isinstance(ivp_method, str) and ivp_method not in embedded_pairs:
            ivp_method = 'Dormand-Prince'
        if ivp_method in embedded_pairs:
            ivp_method = embedded_pairs[ivp_method]
        def zeros(s):
            if isinstance(yn, np.ndarray):
                return np.zeros((s, yn.size))
            else:
                l = [Function(yn).assign(0) for _ in range(s)]
                return l
        def dot(a, b):
            if isinstance(b, np.ndarray):
                return np.dot(a, b)
            else:
                out = Function(b[0])
                out.assign(0)
                for i in range(len(a)):
                    if a[i] != 0:
                        out += b[i] * float(a[i])
                return out
        def eval(f, t, y):
            if isinstance(y, np.ndarray):
                return f(t, y)
            else:
                out = Function(y)
                solve(replace(f, {tn: t, yn: y}) - inner(out, f.arguments()[0]) * dx == 0, out, bcs=bcs)
                return out
        def ivp_step(i, y0):
            if isinstance(y0, np.ndarray):
                func = lambda t, y: (dc) * ff(tn1 + dc*t, y) + g_int(t/dt, i)
                t = 0
                yi=y0
            else:
                t = Constant(0)
                yi = Function(y0)
                test_f = ff.arguments()[0]
                func = dc * replace(ff, {tn: tn1 + dc * t, yn: yi})
                for j in range(len(self.gamm)):
                    func += inner((dot(self.gamm[j,i], FI) * ((t/dt) ** j)), test_f) * dx
                if self.omeg is not None and fe is not None:
                    for j in range(len(self.omeg)):
                        func += inner((dot(self.omeg[j, i], FE) * ((t / dt) ** j)), test_f) * dx
            if isinstance(ivp_method, SundialsExode):
                return ivp_method.solve(func, yi, t, t + dt)
            elif isinstance(ivp_method, EmbeddedTableau):
                return ivp_method.y_step(func, yi, t, dt, **ivp_options)
            else:
                sol = solve_ivp(func, [0, dt], yi, method=ivp_method, t_eval = [0, dt], **ivp_options)
                if sol.success:
                    return sol.y[:, -1]
                print('fast solve failed')
                return np.nan * yi
        def implicit_solve(i):
            if isinstance(y0, np.ndarray):
                def root_fun(Y):
                    result = y0 - Y
                    FI[i] = eval(fi, tn + dt * self.c[i], Y)
                    if self.omeg is not None and fe is not None:
                        FE[i] = eval(fe, tn + dt * self.c[i], Y)
                    for j in range(len(self.gamm)):
                        result += dot(self.gamm[j, i], FI) * (dt ** (j+1))/(j+1)
                    if self.omeg is not None and fe is not None:
                        for j in range(len(self.omeg)):
                            result += dot(self.omeg[j, i], FE) * (dt ** (j+1))/(j+1)
                    return result
                if 'method' not in implicit_solve_options:
                    implicit_solve_options['method'] = 'krylov'
                sol = root(root_fun, y0, **implicit_solve_options)
                if sol.success:
                    return sol.x
                print('implicit solve failed')
                return np.nan * sol.x
            else:
                result = Function(y0)
                result.assign(0)
                for j in range(len(self.gamm)):
                    result += dot(self.gamm[j, i,:i], FI[:i]) * (dt ** (j+1))/(j+1)
                if self.omeg is not None and fe is not None:
                    for j in range(len(self.omeg)):
                        result += dot(self.omeg[j, i,:i], FE[:i]) * (dt ** (j+1))/(j+1)
                Y = Function(y0)
                test_f = fi.arguments()[0]
                F2 = inner((result + Y - y0), test_f) * dx
                for j in range(len(self.gamm)):
                    F2 +=  (dt ** (j+1))/(j+1) * self.gamm[j,i,i]*replace(fi, {tn: tn + dt * self.c[i], yn: Y})
                if self.omeg is not None and fe is not None:
                    for j in range(len(self.omeg)):
                        F2 += self.omeg[j, i, i] * (dt ** (j+1))/(j+1) * replace(fe, {tn: tn + dt * self.c[i], yn: Y})
                solve(F2 == 0, Y, **implicit_solve_options)
                return Y
                
        if self.omeg is not None and fe is not None:
            Y = zeros(self.s)
            FI = zeros(self.s)
            FE = zeros(self.s)
    
            def g_int(t, j):
                result = np.zeros(yn.size)
                for i in range(len(self.gamm)):
                    result += dot(self.gamm[i,j], FI) * (t ** i)
                for i in range(len(self.omeg)):
                    result += dot(self.omeg[i, j], FE) * (t ** i)
                return result
            Y[0] = yn
            FI[0] = eval(fi, tn+dt*self.c[0], Y[0])
            FE[0] = eval(fe, tn + dt * self.c[0], Y[0])
            for i in range(1,self.s):
                y0 = Y[i-1]
                tn1 = tn + self.c[i-1] * dt
                dc = (self.c[i] - self.c[i-1])
                if abs(dc) <= 1e-10 and np.any(self.gamm[:,i,i]):
                    Y[i] = implicit_solve(i)
                else:
                    result = ivp_step(i, y0)
                    Y[i] = result
                FI[i] = eval(fi,tn + dt * self.c[i], Y[i])
                FE[i] = eval(fe,tn + dt * self.c[i], Y[i])
            if isinstance(y0, Function):
                y0.assign(Y[-1])
            return Y[-1]
        else:
            Y = zeros(self.s)
            FI = zeros(self.s)
            def g_int(t, j):
                result = np.zeros(yn.size)
                for i in range(len(self.gamm)):
                    result += dot(self.gamm[i,j], FI) * (t ** i)
                return result
            Y[0] = yn
            FI[0] = eval(fi,tn + dt * self.c[0], Y[0])
            for i in range(1,self.s):
                y0 = Y[i-1]
                tn1 = tn + self.c[i-1] * dt
                dc = (self.c[i] - self.c[i-1])
                if abs(dc) <= 1e-10 and np.any(self.gamm[:,i,i]):
                    Y[i] = implicit_solve(i)
                else:
                    result = ivp_step(i, y0)
                    Y[i] = result
                FI[i] = eval(fi, tn + dt * self.c[i], Y[i])
            if isinstance(yn, Function):
                yn.assign(Y[-1])
                return yn
            return Y[-1]



def multirate_solve(y0, t0, dt, tf, method, fi, ff, fe=None, fname=None, save_steps=None, bcs=None, ivp_options={}, ivp_method = 'RK45', implicit_solve_options={}):
    """ This function uses a multirate infinitesimal method to solve a differential equation
    -----
    Inputs:
    y0 - the current value of y
         if the functions are Forms, this should be of type Function
    t0 - the current value of t
         if using the finite element version, this should be of type Constant
    dt - the amount time will increase by
    tf - the time to solve until
    method - the method to use (of type Multirate)
    fi - the slow operator of the differential equation, with inputs (t, y)
           if there are more than one (for imex-gark), this is the one that will be solved implicitly (with the gamma coefficients)
    ff - the fast operator of the differential equation, with inputs (t, y)
    fe (optional) - the second slow operator, with inputs (t, y)
                these operators may also be finite element Forms from firedrake.
    bcs - optional.  Only used if using the finite element version.  The boundary condition(s) to apply.
    implicit_solve_parameters - optional.  Any solver parameters to use (see firedrake documentation or scipy.roots documentation for details)
    fname - optional.  If provided, will save intermediate results to this file.
          - if using the finite element version of the code, this is a HDF5 file
            otherwise it is a csv file.
    save_steps - the number of intermediate steps to save if fname is provided
           if it is not provided, the default is to save every step
           (or after every dt if embedded methods are being used).
    ivp_options (optional) - a dictionary of keyword arguments to pass to the adaptive solver
    ivp_method (optional) - the adaptive solver to use. Default is RK45, or Dormand-Prince for the finite element version.
    Return
    -----
    the approximate value of y at tf
    """
    
    t = t0
    if isinstance(t, Constant):
        t = float(t.values()[0])
    if fname is not None:
        if isinstance(y0, Function):
            f = CheckpointFile(fname, 'w')
            f.save_mesh(y0.function_space().mesh())
            f.save_function(initial_y, idx=0)
            f.create_group('times')
            f.set_attr('/times', '0', t)
            count_save = 1
        else:
            f = open(fname, 'wb')
            np.savetxt(f, [[t0] + [x for x in y0]], delimiter=',')
        saved = t

    if save_steps is not None:
        save_interval = (tf - t) / save_steps
    else:
        save_interval = dt
    if 'rtol' not in ivp_options:
        ivp_options['rtol'] = 1e-10
    if 'atol' not in ivp_options:
        ivp_options['atol'] = 1e-12
    if ivp_method in ['CV_ADAMS', 'CV_BDF'] or 'ARKODE' in ivp_method or 'ARKODE' in ivp_method[1]:
        ivp_method = SundialsExode(y0, ivp_method, concatenate=False, **ivp_options)
    while t < tf:
        y0 = method.step(y0, dt, fi, ff, t0, fe=fe, bcs=bcs, ivp_options=ivp_options, ivp_method=ivp_method, implicit_solve_options=implicit_solve_options)
        t += dt
        if fname is not None and t - saved - save_interval > -1e-8:
            if isinstance(y0, Function):
                f.save_function(y0, idx=count_save)
                f.set_attr('/times', str(count_save), t)
                f.set_attr('/times', 'last_idx', count_save)
                count_save += 1
            else:
                np.savetxt(f, [[t] + [x for x in y0]], delimiter=',')
            saved += ((t - saved + 1e-8) // save_interval) * save_interval
        if dt > tf - t:
            dt = tf - t
        if isinstance(t0, Constant):
            t0.assign(t)
        else:
            t0 = t
    if fname is not None and abs(t - saved) > 1e-8:
        if isinstance(y0, Function):
            f.save_function(y0, idx=count_save)
            f.set_attr('/times', str(count_save), t)
            f.set_attr('/times', 'last_idx', count_save)
            count_save += 1
        else:
            np.savetxt(f, [[t] + [x for x in y0]], delimiter=',')
    if fname is not None:
        f.close()
    if isinstance(ivp_method, SundialsExode):
        ivp_method.free()
    return y0

cs = np.array([0, 1/3, 3/4, 1])
gamm = np.array([[[0,0,0,0],
                  [1/3,0,0,0],
                  [-3/16-1/3, 15/16, 0, 0],
                  [1/6+3/16, 3/10-15/16, 8/15, 0]]])

# Knoth and Wolke (1997) order 3
mri_kw3 = Multirate(cs, gamm)

cs = np.array([0, 1/3, 2/3, 1])
d=-1/2
gamm = np.array([[[0,0,0,0],
                  [1/3,0,0,0],
                  [(-6*d-7)/12, (6*d+11)/12,0,0],
                  [0, (6*d-5)/12, (3-2*d)/4, 0]],
                 [[0,0,0,0],
                  [0,0,0,0],
                  [(2*d+1)/2, -(2*d+1)/2, 0,0],
                  [1/2, -(2*d+1)/2, d, 0]]])

c2 = 1/2
cs = np.array([0, c2, 1])
gamm = np.array([[[0,0,0],
                  [c2, 0, 0],
                  [-(2 * c2**2 - 2*c2 + 1)/(2*c2), 1/(2*c2), 0]]])
# Sandu (2019) mri-gark-erk22a
mri_erk2a = Multirate(cs, gamm)

c2 = 1
cs = np.array([0, c2, 1])
gamm = np.array([[[0,0,0],
                  [c2, 0, 0],
                  [-(2 * c2**2 - 2*c2 + 1)/(2*c2), 1/(2*c2), 0]]])
# Sandu (2019) mri-gark-erk22b
mri_erk2b = Multirate(cs, gamm)

# Sandu (2019) mri-gark-erk33 (order 3)
mri_erk3 = Multirate(cs, gamm)


cs = np.array([0, 1/5, 2/5, 3/5, 4/5, 1])
gamm = np.array([[[0,0,0,0,0,0],
                  [1/5,0,0,0,0,0],
                  [-53/16, 281/80,0,0,0,0],
                  [-36562993/71394880, 34903117/17848720, -88770499/71394880,0,0,0],
                  [-7631593/71394880, -166232021/35697440, 6068517/1519040, 8644289/8924360, 0, 0],
                  [277061/303808, -209323/1139280, -1360217/1139280, -148789/56964, 147889/45120, 0]],
                 [[0,0,0,0,0,0],
                  [0,0,0,0,0,0],
                  [503/80, -503/80,0,0,0,0],
                  [-1365537/35697440, 4963773/7139488, -1465833/2231090, 0,0,0],
                  [66974357/35697440, 21445367/7139488, -3, -8388609/4462180, 0,0],
                  [-18227/7520, 2, 1, 5, -41933/7520, 0]]])
# Sandu (2019) mri-gark-erk45a
mri_erk4 = Multirate(cs, gamm)

cs = np.array([0, 1, 1, 1])
gamm = np.array([[[0,0,0,0],
                  [1,0,0,0],
                  [-1/2,0,1/2,0],
                  [0,0,0,0]]])
# Sandu (2019) mri-gark-irk21a
mri_irk2 = Multirate(cs, gamm)

g = 0.435866521508458999416019
cs = np.array([0, 1/3, 1/3, 2/3, 2/3, 1, 1, 1])

gamm = np.array([[[0,0,0,0,0,0,0,0],
                  [1/3,0,0,0,0,0,0,0],
                  [-g, 0, g, 0, 0, 0, 0, 0],
                  [(3-10*g)/(24*g - 6), 0, (5-18*g)/(6-24*g), 0, 0, 0, 0, 0],
                  [(-24*g*g+6*g+1)/(6-24*g), 0, (-48*g*g+12*g+1)/(24*g-6), 0, g, 0, 0, 0],
                  [(3-16*g)/(12-48*g), 0,(48*g*g - 21*g + 2)/(12*g-3),0,(3-16*g)/4, 0, 0, 0],
                  [-g, 0, 0, 0, 0, 0, g, 0],
                  [0,0,0,0,0,0,0,0]],
                 ])
# Sandu (2019) mri-gark-esdirk34a
mri_esdirk3a = Multirate(cs, gamm)

g = 0.435866521508458999416019
cs = np.array([0, g, g, (6 * g*g - 9 *g + 2)/(6*g**2 - 12 * g + 3), (6*g**2-9*g+2)/(6*g**2-12*g+3),1,1,1])
gamm = np.array([[[0,0,0,0,0,0,0,0],
                  [1/3,0,0,0,0,0,0,0],
                  [-g,0,g,0,0,0,0,0],
                  [(3-10*g)/(24*g - 6), 0,(5-18*g)/(6-24*g),0,0,0,0,0],
                  [(-24*g**2+6*g+1)/(6-24*g), 0, (-48*g**2+12*g+1)/(24*g-6), 0, g, 0, 0, 0],
                  [(3-16*g)/(12-48*g), 0, (48*g**2 - 21*g + 2)/(12*g-3), 0, (3-16*g)/4, 0, 0, 0],
                  [-g, 0, 0, 0, 0, 0, g, 0],
                  [0,0,0,0,0,0,0,0]]])
# Sandu (2019) mri-gark-sdirk33a
mri_sdirk3 = Multirate(cs, gamm)

cs = np.zeros((8,))
cs[1] = 0.4358665215084589994160194511935568425
cs[2] = cs[1]
cs[3] = 0.7179332607542294997080097255967784213
cs[4] = cs[3]
cs[5] = 1
cs[6] =  1
cs[7] = 1

gamm = np.array([[[0,0,0,0,0,0,0,0],
                  [.4358665215084589994160194511935568425,0,0,0,0,0,0,0],
                  [-.4358665215084589994160194511935568425,0,.4358665215084589994160194511935568425,0,0,0,0,0],
                  [-0.4103336962288525014599513720161078937,0,.6924004354746230017519416464193294724,0,0,0,0,0],
                  [0.4103336962288525014599513720161078937,0,-.8462002177373115008759708232096647362,0,.4358665215084589994160194511935568425,0,0,0],
                  [.4358665215084589994160194511935568425,0,.9264299099302395700444874096601015328,0,-1.080229692192928069168516586450436797,0,0,0],
                  [-.4358665215084589994160194511935568425,0,0,0,0,0,.4358665215084589994160194511935568425,0],
                  [0,0,0,0,0,0,0,0]]])

omeg = np.array([[[0,0,0,0,0,0,0,0],
                  [.4358665215084589994160194511935568425,0,0,0,0,0,0,0],
                  [0,0,0,0,0,0,0,0],
                  [-.5688715801234400928465032925317932021,0,.8509383193692105931384935669350147809,0,0,0,0,0],
                  [.454283944643608855878770886900124654,0,-.454283944643608855878770886900124654,0,0,0,0,0],
                  [-.4271371821005074011706645050390732474,0,.1562747733103380821014660497037023496,0,.5529291480359398193611887297385924765,0,0,0],
                  [0,0,0,0,0,0,0,0],
                  [.105858296071879638722377459477184953,0,.655567501140070250975288954324730635,0,-1.197292318720408889113685864995472431,0,.4358665215084589994160194511935568425,0]]])
# Chinomona and Reynolds (2021) imex-mri-gark3a
mri_imex3 = Multirate(cs, gamm, omeg)

#cs = np.array([0,1/2,1/2,5/8,5/8,3/4,3/4,7/8,7/8,1,1,1])
cs = np.array([0, .5, .5, .625, .625, .75, .75, .875, .875, 1, 1, 1, 1])
gamm = np.zeros((2,12,12))
gamm[0,1,0] = 1/2
gamm[0,2,0] = -1/4
gamm[0,2,2] = 1/4
gamm[0,4,4] = 1/4
gamm[0,6,6] = 1/4
gamm[0,8,8] = 1/4
gamm[0,10,10] = 1/4
gamm[0,3,0] = -3.97728124810848818306703385146227889
gamm[0,3,2] = 4.10228124810848818306703385146227889
gamm[0,4,0] = -0.0690538874140169123272414708480937406
gamm[0,4,2] = -.180946112585983087672758529151906259
gamm[0,5,0] = -1.76176766375792052886337896482241241
gamm[0,5,2] = 2.69452469837729861015533815079146138
gamm[0,5,4] = -0.807757034619378081291959185969048978
gamm[0,6,0] = 0.555872179155396948730508100958808496
gamm[0,6,2] = -0.679914050157999501395850152788348695
gamm[0,6,4] = -.125958128997397447334657948170459801
gamm[0,7,4] = .125958128997397447334657948170459801
gamm[0,7,0] = -5.84017602872495595444642665754106511
gamm[0,7,2] = 8.17445668429191508919127080571071637
gamm[0,7,6] = -2.33523878456435658207950209634011106
gamm[0,8,0] = -1.9067926451678118080947593050360523
gamm[0,8,2] = -1.54705781138512393363298457924938844
gamm[0,9,2] = 1.54705781138512393363298457924938844
gamm[0,8,4] = 4.12988801314935030595449173802031322
gamm[0,9,4] = -4.12988801314935030595449173802031322
gamm[0,8,6] = -0.926037556596414564226747853734872477
gamm[0,9,6] = 0.926037556596414564226747853734872477
gamm[0,9,0] = 3.33702815168872605455765278252966252
gamm[0,9,8] = -1.55523550652091424646289347749361021
gamm[0,10,0] = -.821293629221007618720524112312446752
gamm[0,10,2] = 0.328610356068599988551677264268969646
gamm[0,10,4] = 0.678001812102026694142641232421139516
gamm[0,10,6] = -0.342779287862800022896645471462060708
gamm[0,10,8] = -.0925392510868190410771489129156017025

gamm[1,3,0] = 8.0925392510868190410771489129156017025
gamm[1,3,2] = -8.0925392510868190410771489129156017025
gamm[1,5,0] = 3.91164310234387488238124087134101229
gamm[1,5,2] = -5.02715717158263104496515924327911025
gamm[1,5,4] = 1.11551406923875616258391837193809796
gamm[1,7,0] = 10.8186076991391180114318371131645132
gamm[1,7,2] = -14.9890852682678311755908413058447354
gamm[1,7,6] = 4.17047756912871316415900419268022213
gamm[1,9,0] = -2.61047101304182849292578695498722043
gamm[1,9,8] = 2.61047101304182849292578695498722043

gamm = np.array([[[0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0.5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-0.25,	0,	0.25,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-3.977281248108488183067033851462278892,	0,	4.102281248108488183067033851462278892,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-0.06905388741401691232724147084809374064,	0,	-0.1809461125859830876727585291519062594,	0,	0.25,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-1.761767663757920528863378964822412405,	0,	2.694524698377298610155338150791461384,	0,	-0.8077570346193780812919591859690489783,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0.5558721791553969487305081009588084962,	0,	-0.6799140501579995013958501527883486949,	0,	-0.1259581289973974473346579481704598013,	0,	0.25,	0,	0,	0,	0,	0,	0],
                  [-5.840176028724955954446426657541065113,	0,	8.174456684291915089191270805710716374,	0,	0.1259581289973974473346579481704598013,	0,	-2.335238784564356582079502096340111063,	0,	0,	0,	0,	0,	0],
                  [-1.906792645167811808094759305036052304,	0,	-1.547057811385123933632984579249388443,	0,	4.129888013149350305954491738020313225,	0,	-0.9260375565964145642267478537348724775,	0,	0.25,	0,	0,	0,	0],
                  [3.337028151688726054557652782529662519,	0,	1.547057811385123933632984579249388443,	0,	-4.129888013149350305954491738020313225,	0,	0.9260375565964145642267478537348724775,	0,	-1.555235506520914246462893477493610215,	0,	0,	0,	0],
                  [-0.8212936292210076187205241123124467518,	0,	0.328610356068599988551677264268969646,	0,	0.6780018121020266941426412324211395162,	0,	-0.3427792878628000228966454714620607079,	0,	-0.0925392510868190410771489129156017025,	0,	0.25,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0]],
                 [[0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [8.704562496216976366134067702924557783,	0,	-8.704562496216976366134067702924557783,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [3.911643102343874882381240871341012292,	0,	-5.027157171582631044965159243279110249,	0,	1.115514069238756162583918371938097957,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [10.81860769913911801143183711316451323,	0,	-14.98908526826783117559084130584473536,	0,	0,	0,	4.170477569128713164159004192680222125,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-2.61047101304182849292578695498722043,	0,	0,	0,	0,	0,	0,	0,	2.61047101304182849292578695498722043,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0]]])

omeg = np.zeros((2,12,12))
omeg[0,1,0] = 1/2
omeg[0,3,0] = -1.91716534363662868878172216064946905
omeg[0,3,2] = 2.04216534363662868878172216064946905
omeg[0,4,0] = -0.404751031801105942697915907046990469
omeg[0,4,2] = 0.404751031801105942697915907046990469
omeg[0,5,0] = 11.4514660224922163666569802860263173
omeg[0,5,2] = -30.2107574752650427144064781557395061
omeg[0,5,4] = 18.8842914527728263477494978697131888
omeg[0,6,0] = -0.709033564760261450684711672946330144
omeg[0,6,2] = 1.03030720858751876652616190884004718
omeg[0,6,4] = -0.321273643827257315841450235893717036
omeg[0,7,4] = 0.321273643827257315841450235893717036
omeg[0,7,0] = -29.9954871645582843984091068494419927
omeg[0,7,2] = 37.605982774991801805364896856243857
omeg[0,7,6] = -7.80676925426077472279724024269558129
omeg[0,8,0] = 3.10466505427296211633876939184912422
omeg[0,8,2] = -2.43032501975716229713206592741556636
omeg[0,9,2] = 2.43032501975716229713206592741556636
omeg[0,8,4] = -1.90547930115152463521920165948384213
omeg[0,9,4] = 1.90547930115152463521920165948384213
omeg[0,8,6] = 1.23113926663572481601249819505028427
omeg[0,9,6] = -1.23113926663572481601249819505028427
omeg[0,9,0] = -2.42442954775204786987587591435551401
omeg[0,9,8] = -.555235506520914246462893477493610215
omeg[0,10,0] = -.010441350444797485902945189451653542
omeg[0,10,2] = .0726030361465507450515210450548814161
omeg[0,10,4] = -.128827595167726095223945409857642431
omeg[0,10,6] = 0.112935535009382356613944010712215408
omeg[0,10,8] = -.0462696255434095205385744564578008512
omeg[0,11,8] = -.0462696255434095205385744564578008512
omeg[0,11,0] = -.81085227877621013281757892286079321
omeg[0,11,2] = .25600731992204924350015621921408823
omeg[0,11,4] = .806829407269752789366586642278781947
omeg[0,11,6] = -.455714822872182379510589482174276116
omeg[0,11,10] = 1/4

omeg[1,3,0] = 4.0843306872732573775634443212989381
omeg[1,3,2] = -4.0843306872732573775634443212989381
omeg[1,5,0] = -21.8434299813822208479181287579586536
omeg[1,5,2] = 59.6120128869278735434171244973850312
omeg[1,5,4] = -37.7685829055456526954989957394263776
omeg[1,7,0] = 61.6590414586370916981876370447766458
omeg[1,7,2] = -77.2725799671586411437821175301678084
omeg[1,7,6] = 15.6135385085215494455944804853911626
omeg[1,9,0] = -1.11047101304182849292578695498722043
omeg[1,9,8] = 1.11047101304182849292578695498722043

omeg = np.array([[[0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0.5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-1.91716534363662868878172216064946905,	0,	2.04216534363662868878172216064946905,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-0.4047510318011059426979159070469904691,	0,	0.4047510318011059426979159070469904691,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [11.45146602249221636665698028602631728,	0,	-30.21075747526504271440647815573950607,	0,	18.88429145277282634774949786971318879,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-0.7090335647602614506847116729463301439,	0,	1.03030720858751876652616190884004718,	0,	-0.3212736438272573158414502358937170357,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-29.99548716455828439840910684944199275,	0,	37.60598277499180180536489685624385701,	0,	0.3212736438272573158414502358937170357,	0,	-7.806769254260774722797240242695581295,	0,	0,	0,	0,	0,	0],
                  [3.104665054272962116338769391849124223,	0,	-2.430325019757162297132065927415566359,	0,	-1.905479301151524635219201659483842131,	0,	1.231139266635724816012498195050284266,	0,	0,	0,	0,	0,	0],
                  [-2.424429547752047869875875914355514008,	0,	2.430325019757162297132065927415566359,	0,	1.905479301151524635219201659483842131,	0,	-1.231139266635724816012498195050284266,	0,	-0.555235506520914246462893477493610215,	0,	0,	0,	0],
                  [-0.01044135044479748590294518945165354204,	0,	0.07260303614655074505152104505488141613,	0,	-0.1288275951677260952239454098576424313,	0,	0.1129355350093823566139440107122154084,	0,	-0.04626962554340952053857445645780085125,	0,	0,	0,	0],
                  [-0.8108522787762101328175789228607932098,	0,	0.2560073199220492435001562192140882299,	0,	0.8068294072697527893665866422787819475,	0,	-0.4557148228721823795105894821742761164,	0,	-0.04626962554340952053857445645780085125,	0,	0.25,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0]],
                 [[0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [4.084330687273257377563444321298938099,	0,	-4.084330687273257377563444321298938099,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-21.84342998138222084791812875795865363,	0,	59.61201288692787354341712449738503121,	0,	-37.76858290554565269549899573942637758,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [61.65904145863709169818763704477664579,	0,	-77.27257996715864114378211753016780838,	0,	0,	0,	15.61353850852154944559448048539116259,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [-1.11047101304182849292578695498722043,	0,	0,	0,	0,	0,	0,	0,	1.11047101304182849292578695498722043,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0],
                  [0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0]]])


# Chinomona and Reynolds imex-mri-gark4
mri_imex4 = Multirate(cs, gamm, omeg)


