import numpy as np
import math as m
import cmath as cm
from butcher_tableau import Tableau,tableaus, embedded_pairs, EmbeddedTableau, measure_error, compute_time
from adi_methods import adi_step
from scipy.integrate import solve_ivp
import sys
import time as timeos
import csv

try:
    from firedrake import Function, Constant, replace, CheckpointFile, Form
    fem = True
except:
    fem = False
    Function = type(None)
    Constant = type(None)

epi_methods = {}
EPIMethod = type(None)
def time_step(function, delta_t, y, initial_t,
                    tableau, index, **kwargs):
    """
    This function takes a step in time in the differential equation.
    Inputs
    ------
    function : callable function(x)
        the function to use as y prime
    delta_t : float
        the amount the time will increase by   (i.e. Delta_t * alpha_j^[i])
    y : array
        the current value of y
    initial_t
        the current value of t
    tableau
        the butcher tableau or EPI method to use, default is forward euler
    Return
    ------
    array or float
        type depends on type of initial_y
        the approximate value of y after the given delta_t
    """
    if isinstance(tableau, Tableau):
        f=tableau.y_step(function, y, initial_t, delta_t, index, **kwargs)
        if isinstance(y, Function):
            y.assign(y + delta_t * f)
        elif np.size(y)!=1:
            y[index]+=delta_t*f[index]  
        else:
            y+=delta_t*f
    else:
        tableau, info = tableau
        y = tableau.step(function, initial_t, y, delta_t, info, **kwargs)
    return y

def time_step_analytic(function, delta_t, y, initial_t, index):
    # this function takes the analytic solution of the DE and returns it.
    # Only used when we want to use analytic solution for each sub-integrators instead of numerically exact values.
    f = function(delta_t,y)   
    if isinstance(y, Function):
        y.assign(f)
        return y
    f = np.array(f)   
    if np.size(y) != 1:
        y[index] = f[index]
    else:
        y = f
    return y

def exact_solution(function, initial_t, delta_t, initial_y, ivp_method, mask,
                   rtol, atol, **kwargs):
    """
    ----------
    This function finds a highly accurate numerical solution for a 
    time step in the ODE 
    ----------
    function : the function to use as y prime
        the calling signature is function(t, y)
    initial_t : float
        the value of t the solver will start with
    delta_t : 
        the value t will increase by
    initial_y : 
        initial state. 
        For problems in the complex domain, pass y0 
        with a complex data type (even if the initial value is purely real).
    ivp_method :
        string or OdeSolver or EmbeddedTableau. Integration method to use.
    rtol, atol : 
        float or array_like. Relative and absolute tolerances.
        If using an EmbeddedTableau as the method, must be float, 
            and only rtol is used
        
    Returns
    -------
    y : Value of the solution at the final t in the subinterval.
    """
    tf = initial_t + delta_t
    t_span = (initial_t, tf)
    y0 = initial_y
    if isinstance(ivp_method, EmbeddedTableau):
        y = ivp_method.y_step(function, initial_y, initial_t, delta_t, mask, rtol, **kwargs)
    else:
        sol = solve_ivp(function, t_span, y0, method = ivp_method, 
                        rtol = rtol, atol = atol, dense_output = True)
        y = sol.y[:,-1]
    return y


def process_os_options(functions, initial_y, initial_t, delta_t, alpha, methods, b=None, ivp_methods={}, epi_options={}, jacobian = None):
    alpha, order = alphas_repo(alpha, b, len(functions))
    
    if order is not None:
        delta_t = 1
    function_list = []
    if not isinstance(initial_y, Function):
        complex_flag = (initial_y.dtype == np.complex128)
    else:
        complex_flag = False
    start_times = [0] * len(functions)
    if jacobian is not None:
        J = jacobian(initial_t, initial_y)
    for j in range(len(alpha)):   # j for stages of the OS
        for i in range(len(functions)):  # i for operators
            if alpha[j][-1] is True:  #index = -1 for last element
                k=len(functions)-i-1
                function=functions[-i-1]
                if k==(len(functions)-1):
                    step=alpha[j][-2]
                else:
                    step=alpha[j][(-i-2)%(len(alpha[j])-2)]
            else:
                k=i
                function=functions[i]
                if i==(len(functions)-1):
                    step=alpha[j][-1]
                else:
                    step=alpha[j][(i)%(len(alpha[j])-1)]   # Get the step = delta_t*alpha_j^i

            # use k,j to determine which RK method to use, if none, use FE.
            if (k+1, j+1) in methods:
                tableau=methods[(k+1, j+1)]
            elif (k+1,) in methods:
                tableau=methods[(k+1,)]
            elif (0,) in methods and methods[(0,)] != 'ANALYTIC':
                tableau = methods[(0,)]
            else:
                tableau='FE'
            

            if tableau in tableaus:
                tableau=tableaus[tableau]   # matching RK method name with RK coefficients
            elif tableau in epi_methods:
                tableau = epi_methods[tableau]
            if isinstance(tableau, EPIMethod):
                if k+1 in epi_options:
                    info = epi_options[k+1]
                elif 0 in epi_options:
                    info = epi_options[0]
                else:
                    info = ('kiops', 1e-10)
                tableau = (tableau, info)
            exact_flag = (tableau == 'EXACT') 
            if tableau == 'EXACT':
                if k+1 in ivp_methods:
                    tableau = ivp_methods[k+1]
                else:
                    tableau = ('RK45', 1e-3, 1e-6)
                if tableau[0] in embedded_pairs:
                    tableau = (embedded_pairs[tableau[0]], tableau[1], tableau[2])

            if not isinstance(initial_y, Function):
                if jacobian is not None:
                    mask=np.logical_not(np.isnan(function(initial_t, initial_y, J))) #check if f(t_0,y_0) is defined
                else:
                    mask = np.logical_not(np.isnan(function(initial_t, initial_y)))
            else:
                mask = None
            if isinstance(step, complex):
                complex_flag=True
            if step!=0:
                function_list.append((k, step * delta_t, function, tableau, mask, exact_flag, start_times[k]))
            if isinstance(step, tuple):
                if isinstance(start_times[k], tuple):
                    start_times[k] = (start_times[k][0] + step[0] * delta_t, start_times[k][1] + step[1] * delta_t)
                else:
                    start_times[k] = (start_times[k] + step[0] * delta_t, start_times[k] + step[1] * delta_t)
            else:
                start_times[k] += step * delta_t
    return function_list, complex_flag, order


def ave_Godunov_step(functions,initial_t,y,delta_t):
    tableau1=tableaus['FE']
    tableau2=tableaus['BE']
    
    mask1=np.logical_not(np.isnan(functions[0](initial_t, y)))
    mask2=np.logical_not(np.isnan(functions[1](initial_t, y)))
    y1=np.array(y)
    y1=time_step(functions[0], delta_t, y1, initial_t, tableau1, mask1)
    y1=time_step(functions[1], delta_t, y1, initial_t, tableau2, mask2)

    y2=np.array(y)
    y2=time_step(functions[1], delta_t, y2, initial_t, tableau1, mask1)
    y2=time_step(functions[0], delta_t, y2, initial_t, tableau2, mask2)

    return (y1+y2)/2

def operator_splitting_inner(functions, delta_t, initial_y, initial_t, final_t, 
                             alpha, methods, b=None, fname=None, save_steps=0, ivp_methods={}, epi_options={}, os_tol = 1e-3, solver_parameters={}, jacobian=None, bc=None):
    y=initial_y
    #process it
    function_list, complex_flag, order = process_os_options(functions, initial_y, initial_t, delta_t, alpha, methods, b, ivp_methods, epi_options, jacobian)
    if complex_flag and not isinstance(initial_y, Function):
        y=np.array(initial_y, np.complex128)

    t = initial_t
    if isinstance(t, Constant):
        t = t.values()[0]
        if isinstance(t, np.complex128):
            t = complex(t)
        else:
            t = float(t)
    if fname is not None:
        if isinstance(initial_y, Function):
            f = CheckpointFile(fname, 'w')
            f.save_mesh(initial_y.function_space().mesh())
            f.save_function(initial_y, idx=0)
            f.create_group('times')
            f.set_attr('/times', '0', t)
            count_save = 1
        else:
            f=open(fname, 'wb')
            np.savetxt(f, [[initial_t] + [x for x in y]], delimiter=',')
        saved = t

    if save_steps != 0:
        save_interval = (final_t - t) / save_steps
    else:
        save_interval = delta_t

    tic = timeos.perf_counter()
    scale = 1
    i = 0
    if order is not None:
        scale = delta_t
    
    while abs(t - final_t) > 1e-8:
        if (t + delta_t - final_t).real > 0:
            delta_tn = final_t - t
            if order is None:
                scale = delta_tn / delta_t
            else:
                scale = delta_tn
            delta_t = delta_tn
        if order is not None:
            if not isinstance(y, Function):
                y1 = np.array(y)
            else:
                ys = Function(y)
                y1 = y
        else:
            y1 = y
        split = False
        if jacobian is not None:
            J = jacobian(t, y)
        for line in function_list:
            (k, step, function_i, tableau, mask, exact_flag, start_time)=line     # Here, if exact solution needed, use y=exact_solution(), if RK method needed use y=time_step() (as indicated by exact_flag)
            if not isinstance(initial_y, Function) and jacobian is not None:
                function = lambda t, y: function_i(t, y, J)
            else:
                function = function_i
            if isinstance(step, tuple):
                if not split:
                    if isinstance(initial_y, Function):
                        y2 = Function(y1)
                        if isinstance(function, Form):
                            function2 = replace(function, {y: y2})
                        elif len(function) == 2:
                            function2 = (replace(function[0], {y: y2}), function[1])
                        else:
                            function2 = (replace(function[0], {y: y2}), function[1], function[2])
                    else:
                        y2 = np.array(y1)
                        function2 = function
                split = True
                step2 = step[1] * scale
                step = step[0]
                if isinstance(start_time, tuple):
                    start_time2 = start_time[1] * scale
                    start_time = start_time[0]
                else:
                    start_time2 = start_time * scale
            step = step * scale
            start_time = start_time * scale

            ti = t + start_time
            if isinstance(initial_t, Constant):
                initial_t.assign(ti)
                ti = initial_t
            #Calculate y for each step, taking y from previous step as initial condition
            if (k+1,) in methods and methods[(k+1,)] == 'ANALYTIC':
                y1 = time_step_analytic(function, step, y1, ti, mask)
                if split:
                    if isinstance(initial_t, Constant):
                        t0 = Constant(initial_t)
                        initial_t.assign(t + start_time2)
                        ti = initial_t
                    else:
                        ti = t + start_time2
                    if (k+1) in solver_parameters:
                        params = solver_parameters[k+1]
                    else:
                        params = None
                    y2 = time_step_analytic(function2, step2, y2, ti, mask, params = params)
                    if isinstance(initial_t, Constant):
                        initial_t.assign(t0)
            elif exact_flag:
                ivp_method, rtol, atol = tableau
                if k+1 in solver_parameters:
                    params = solver_parameters[k+1]
                else:
                    params = None
                y1 = exact_solution(function, ti, step, y1, ivp_method, mask,
                                    rtol, atol, params=params)
                if split:
                    if isinstance(initial_t, Constant):
                        t0 = Constant(initial_t)
                        initial_t.assign(t + start_time2)
                        ti = initial_t
                        
                    else:
                        ti = t + start_time2
                    y2 = exact_solution(function2, ti, step2, y2, ivp_method, mask,
                                        rtol, atol, params=params)
                    if isinstance(initial_t, Constant):
                        initial_t.assign(t0)
            else:
                if k+1 in solver_parameters:
                    params = solver_parameters[k+1]
                else:
                    params = None
                y1=time_step(function, step, y1, ti, tableau, mask, params = params)
                if split:
                    if isinstance(initial_t, Constant):
                        t0 = Constant(initial_t)
                        initial_t.assign(t + start_time2)
                        ti = initial_t
                    else:
                        ti = t + start_time2
                    y2 = time_step(function2, step2, y2, ti, tableau, mask, params=params)
                    if isinstance(initial_t, Constant):
                        initial_t.assign(t0)
        if bc is not None:
            bc.apply(y1)
        if order is not None:
            accept, err = measure_error(y, y1, y2, os_tol)
        else:
            accept = True
        if accept:
            t += delta_t
            y = y1
        elif isinstance(y, Function):
            y.assign(ys)
        if order is not None:
            delta_t = compute_time(err, order, delta_t, os_tol)
            scale = delta_t
        if fname is not None and t - saved - save_interval > -1e-8:
            if isinstance(initial_y, Function):
                f.save_function(y, idx=count_save)
                f.set_attr('/times', str(count_save), t)
                f.set_attr('/times', 'last_idx', count_save)
                count_save += 1
            else:

                np.savetxt(f, [[t]+[x for x in y]], delimiter=',') 
            saved += ((t - saved + 1e-8) // save_interval) * save_interval
        if not isinstance(y, Function) and not np.all(np.isfinite(y)):
            return np.nan*y


    toc = timeos.perf_counter()
    print(f"{toc - tic:0.4f}")
    

    if fname is not None:
        f.close()
    return y

def operator_splitting(
        functions, delta_t, initial_y, initial_t, final_t, alpha,
        methods={},b=None, fname=None, save_steps=0, ivp_methods={},
        epi_options={}, os_tol = 1e-3, solver_parameters={}, jacobian=None,
        bc=None):
    """
    This function uses operator splitting with n functions
    to approximate a differential equation
    Inputs
    ------
    functions : list of functions to use
        functions to use to approximate the differential equation in order
        the functions will be numbered 1 to n
        if any element returns np.nan, that element will not be integrated for 
        that time step
        inputs are (t, y)
        These may also be finite element Forms as provided by firedrake, or
        tuples of (Form, boundary condition)
    delta_t : float
        the amount the time will increase by
    initial_y
        the current value of y to use
        if using the finite element version, this should be of type Function
    initial_t
        the current value of t to use
        if using the finite element version, this should be of type Constant
    final_t
        the time to solve until
    alpha
        a string with the name of the operator splitting method to use,
        or a list of lists representing each substep to take
        at each location
    b value needs to be determined if OS22b method is called to use a particular 2nd-order 2-stage OS method
    methods
        a dictionary providing the methods to use at each step.
        if tuple (i, j) is provided, where i is the number of the 
        operator and j is the stage, starting from 1, it is used.
        else if (i) is provided, where i is the number of the 
        operator, it is used.
        may provide an overall default with the key (0,)
        if none are provided, it defaults to Forward Euler
    fname : string
        the file to save intermediate results to
        If this is not provided, the intermediate results will not be saved
    save_steps : integer
        the number of intermediate steps to save if fname is provided
        if it is not provided, the default is to save all steps to the file
        (or after every delta_t if embedded splitting methods are being used)
    The remaining parameters are 
    Return
    ------
    float or array
        type depends on type of initial_y
        the approximate next value of y after the given delta_t
    """
    if isinstance(initial_y, Function):
        initial_y = initial_y
    elif np.size(initial_y)==1:
        if isinstance(initial_y, complex):
            initial_y=np.array([0], np.complex128)+initial_y
        else:
            initial_y=np.array([0], np.float64)+initial_y
    else:
        
        if initial_y.dtype==np.complex128:
            initial_y=np.array(initial_y, np.complex128)
        else:
            initial_y=np.array(initial_y, np.float64)
            
    if alpha in adi_list:
        time=initial_t
        y=initial_y
        while time<final_t:
            y=adi_step(functions, time, delta_t, y, alpha)
            time+=delta_t
    elif alpha=='average Godunov':
        time=initial_t
        y=initial_y
        while time<final_t:
            y=ave_Godunov_step(functions, time, y, delta_t)
            time+=delta_t
    else:
        y=operator_splitting_inner(functions, delta_t, initial_y, initial_t, 
                                   final_t, alpha, methods,b, 
                                   fname, save_steps, ivp_methods, epi_options,
                                   os_tol, solver_parameters, jacobian, bc)

    if not isinstance(y, Function) and np.all(np.isreal(y)):
        y=y.real
    if not isinstance(y, Function) and np.size(y)==1:
        return y[0]
    else:
        return y

def alphas_repo(alpha,b, N):
    """
    This function contains a dictionary of potential OS methods. It can accept
    a user-defined parameter and apply it to a specified OS method. 

    Parameters
    ----------
    alpha : string
        user-defined splitting method
    b : None, int or float
        Parameter for splitting method

    Returns
    -------
    the coefficients corresponding to the splitting method (an array)

    """
    alphas={}
    embedded_methods = {}
    alphas['Godunov']=[[1, 1]]
    alphas['Godunov-3']=[[1,1,1]]
    alphas['SM2']=[[0.5, 0.5],
                   [0.5, 0.5, True]]
    alphas['Strang'] = [[0.5,1],
                        [0.5,0]]
    alphas['Strang-3'] = [[0.5, 0.5, 1],
                         [0,0.5,0],
                      [0.5,0,0]]
    alphas['StrangBCA-3'] = [[0, 0.5, 0.5],
                          [1.0, 0, 0.5],
                          [0, 0.5, 0]]
    alphas['StrangBAC-3'] = [[0, 0.5, 0],
                          [0.5, 0, 1.0],
                          [0.5, 0.5, 0]]
    alphas['StrangCBA-3'] = [[0, 0, 0.5],
                             [0, 0.5, 0],
                             [1.0, 0.5, 0.5]]
    alphas['StrangCAB-3'] = [[0, 0, 0.5],
                             [0.5, 1.0, 0],
                             [0.5, 0, 0.5]]
    alphas['StrangACB-3'] = [[0.5, 0, 0.5],
                            [0., 1.0, 0.5],
                            [0.5, 0, 0.]]
    alphas['Strang-3split'] = [[0.5, 0.5, 0.5],
                               [0.0,0.0,0.5],
                                [0, 0.5, 0],
                                [0.5, 0, 0]]
    alphas['StrangACB-3split'] = [[0.5, 0, 0.5],
                             [0., 0.5, 0.0],
                             [0., 0.5, 0.0],
                             [0.0,0.0,0.5],
                             [0.5, 0, 0.]]

    alphas['Best22'] = [[1.0-m.sqrt(2.0)/2.0, m.sqrt(2.0)/2.0],
                        [m.sqrt(2.0)/2.0, 1.0-m.sqrt(2.0)/2.0]]
    alphas['R3']=[[7./24, 2./3],
               [3./4, -2./3],
               [-1./24, 1.]]
    theta=1./(2-(2**(1.0/3)))
    
    alphas['Y4']=[[theta/2, theta],
	 	             [(1-theta)/2, 1-2*theta],
	 	             [(1-theta)/2, theta],
	 	              [theta/2, 0]]

    alphas['Y4-3']=[[0,0, theta/2],
                [0, theta/2, 0],
                [theta, theta/2, (1-theta)/2],
                [0, (1-2*theta)/2, 0],
                [1-2*theta, (1-2*theta)/2, (1-theta)/2],
                [0, theta/2, 0],
                [theta, theta/2, theta/2]]

    alphas['C3']=[[(1+1j/m.sqrt(3))/4, (1+1j/m.sqrt(3))/2],
               [1./2, (1-1j/m.sqrt(3))/2],
               [(1-1j/m.sqrt(3))/4, 0]]

    a11=0.0935003487263305760
    a12=-0.0690943698810950380
    a13=0.4755940211547644620
    a21=0.439051727817158558
    a22=-0.136536314071511211
    a23=0.394969172508705306
    alphas['M4']=[[a11, a21],
              [a12, a22],
              [a13, a23],
              [a13, a22],
              [a12, a21],
              [a11, 0]]
    alphas['S3']=[[1./6, 1./6],
              [1./6,1./6],
              [1./6,1./6],
              [-1./3, -1./3, True],
              [1./6, 1./6, True],
              [1./6, 1./6],
              [1./6, 1./6],
              [1./6,1./6],
              [1./6, 1./6, True]]
    a11=0.0792036964311957
    a21=0.3531729060497740
    a31=-0.0420650803577195
    a21=0.209515106613362
    a22=-0.143851773179818
    alphas['B4']=[[a11, a21],
              [a12, a22],
              [a13, 0.5-(a21+a22)],
              [1-2*(a11+a12+a13), 0.5-(a21+a22)],
              [a13, a22],
              [a12, a21],
              [a11, 0]
              ]
    alphas['AKS3']=[[0.268330095781759925, 0.919661523017399857],
                [-0.187991618799159782, -0.187991618799159782],
                [0.919661523017399857, 0.268330095781759925]]

    w1=1./(2-2**(1./3)*cm.exp(2j*m.pi/3))
    w0=1-2*w1
    alphas['C4']=[[w1/2, w1],
              [(w0+w1)/2, w0],
              [(w0+w1)/2, w1],
              [w1/2, 0]]
    alphas['C4-3']=[[0, 0, w1/2],
                [0, w1/2, 0],
                [w1, w1/2, (w1+w0)/2],
                [0, w0/2, 0],
                [w0, w0/2, (w1+w0)/2],
                [0, w1/2, 0],
                [w1, w1/2, w1/2]]
    alphas['AKS3C']=[[0, 1./4+1j*m.sqrt(3)/12],
                 [1./2+1j*m.sqrt(3)/6, 1./2],
                 [1./2-1j*m.sqrt(3)/6, 1./4-1j*m.sqrt(3)/12]]
    alphas['AKS3P']=[[0.201639688260407656+0.105972321241365172j, 0.387747410753696807+0.100071120693574555j],
                 [0.410612900985895537-0.206043441934939727j, 0.410612900985895537-0.206043441934939727j],
                 [0.387747410753696807+0.100071120693574555j, 0.201639688260407656+0.105972321241365172j]]
    alphas['PR']=[[0.5,0.5, True],
              [0.5,0.5]]
    alphas['PR2']=[[0.5,0.5, True],
               [0.5,0.5]]


    alphas['AKOpt22'] = [[0.2929,0.7071],[0.7071,0.2929]]
    alphas['AK2s3i-3'] = [[0.5, 1-m.sqrt(2)/2, m.sqrt(2)/2],
                       [0, m.sqrt(2)/2, 1-m.sqrt(2)/2], 
                       [0.5, 0, 0]]
    alphas['AK2s3ii-3'] = [[0.316620935432115636, 0.273890572734778059, 0.662265355057626845],
                       [-0.030373607778656857, 0.438287559165397521, 0.066439991053339223], 
                       [0.713752672346541221, 0.287821868099824420, 0.271294653889033932]]
    
    alphas['AK2s5-3'] = [[0.161862914279624, 0.242677859055102, 0.5],
                        [0.338137085720376, 0.514644281889796, 0],
                        [0.338137085720376, 0, 0.5], 
                        [0, 0.242677859055102, 0], 
                        [0.161862914279624, 0, 0]]

    alphas['EmbAK4s5'] = [[0.267171359000977615, -0.361837907604416033],
                        [-0.0338279096695056672, 0.861837907604416033],
                        [0.5333131013370561044, 0.861837907604416033], 
                        [-0.0338279096695056672, -0.361837907604416033,], 
                        [0.267171359000977615, 0]]    

    alphas['EmbAK3s5'] = [[0.267171359000977615, -0.361837907604416033],
                        [-0.0338279096695056672, 0.861837907604416033],
                        [0.5333131013370561044, 0.395088376480991403], 
                        [0.267171359000977615, -0.361837907604416033,], 
                        [-0.0338279096695056672, 0.466749531123424630]]
    alphas['Milne22_Complex_i'] = [[12.0/37.0 - 2.0/37.0j,	25.0/34.0 - 1.0/17.0j],
                  [25.0/37.0 + 2.0/37.0j, 9.0/34.0 + 1.0/17.0j]]

    alphas['OS32_Strang_minLEM-3'] = [[0.5, 0.5, 1.0],
                               [0., 0.5, 0],
                               [0.5,0,0]]   # This is a 3-stage 2nd-order 3 splitting method similar to Strang with LEM =1.1895


    alphas['OS32_7op_minLEM-3'] = [[0.306975546320853,    0.306975546320853,    0.721475263023673,],
                                   [0,     0.414499716702820,                     0,],
                                   [0.693024453679147, 0.278524736976327,    0.278524736976327]] # This is a 3-stage
                                  # 2nd-order 3 splitting method with 7 operations and smalled LEM among 7op OS32


    # 3 splitting complex coefficient OS
    alphas['AKT22_C'] = [[0.5 + 0.5j, 0.5 + 0.5j, 0.5 + 0.5j],
                         [0.5 - 0.5j, 0.5 - 0.5j, 0.5 - 0.5j]]
    alphas['PP3_3C'] = [[0.0442100822731214750-0.0713885293035937610j,  0.0973753110633760580-0.112390152630243038j,0.125415464915697242-0.281916718734615225j],
                        [0.157419072651724312-0.1552628290245811054j, 0.179226865237094561-0.0934263750859694960j, 0.353043498499040389+0.0768951336684972038j],
                        [0.260637333463417766+0.07744172526769638060j, 0.223397823699529381+0.205816527716212534j, 0.059274548196998816+0.354231218126596507j],
                        [0.059274548196998816+0.354231218126596507j, 0.223397823699529381+0.205816527716212534j, 0.260637333463417766+0.07744172526769638060j],
                        [0.353043498499040389+0.0768951336684972038j, 0.179226865237094561-0.0934263750859694960j, 0.157419072651724312-0.1552628290245811054j],
                        [0.125415464915697242-0.281916718734615225j, 0.0973753110633760580-0.112390152630243038j, 0.0442100822731214750-0.0713885293035937610j]]

    
    # 5 splitting complex coefficient OS
    alphas['Godunov-5'] = [[1, 1, 1, 1, 1]]
    alphas['Strang-5'] = [[0.5, 0.5, 0.5, 0.5, 1],
                          [0, 0, 0, 0.5, 0],
                          [0, 0, 0.5, 0, 0],
                          [0, 0.5, 0, 0, 0],
                          [0.5, 0, 0, 0, 0]]

    Y4gamma1 = 1 / (2 - 2 ** (1 / 3))
    Y4gamma2 = 1 - 2 * Y4gamma1
    Y4gamma3 = Y4gamma1
    alphas['Yoshida-5'] = [[Y4gamma1 / 2., Y4gamma1 / 2., Y4gamma1 / 2., Y4gamma1 / 2., Y4gamma1],
                           [0, 0, 0, Y4gamma1 / 2., 0],
                           [0, 0, Y4gamma1 / 2., 0, 0],
                           [0, Y4gamma1 / 2., 0, 0, 0],
                           [(Y4gamma1 + Y4gamma2) / 2., Y4gamma2 / 2., Y4gamma2 / 2., Y4gamma2 / 2., Y4gamma2],
                           [0, 0, 0, Y4gamma2 / 2., 0],
                           [0, 0, Y4gamma2 / 2., 0, 0],
                           [0, Y4gamma2 / 2., 0, 0, 0],
                           [(Y4gamma1 + Y4gamma2) / 2., Y4gamma1 / 2., Y4gamma1 / 2., Y4gamma1 / 2., Y4gamma1],
                           [0, 0, 0, Y4gamma1 / 2., 0],
                           [0, 0, Y4gamma1 / 2., 0, 0],
                           [0, Y4gamma1 / 2., 0, 0, 0],
                           [Y4gamma1 / 2., 0, 0, 0, 0]]

    if alpha == 'OS22b':
        if b == None or b == 1.0: # terminate if b == None or b == 1 
            print('Warning: method not defined')
            sys.exit()
        else:
            alphas['OS22b'] = [[(2*b-1)/(2*b-2), 1-b],
                   [ -1/(2*b - 2), b]]
    if alpha == 'OS33bv1':
        if b == None or b == 0.25 or b == 1 or b == 1.0/3.0:   # terminate for unusable b values
            print('Warning: method not defined.')
            sys.exit()
        else:
            sqrtexp = m.sqrt(144.0*(b**4) + 72.0*(b**3) - 99.0*(b**2) + 30.0*b - 3)
            b3 = b
            b2 = (-12.*(b**2)+15.0*b-3+sqrtexp)/(24.0*b-6)
            b1 = 1-b2-b3
            a3 = -(3.0*b2+3.*b3-1)/(6.0*(b3-1)*b2);
            a2 = (2.*a3*b3-2.*a3+1)/(2.*b1)
            a1 = 1-a2-a3
            alphas['OS33bv1'] = [[a1,b1],[a2,b2],[a3,b3]]


    if alpha == 'OS33bv2':
        if b == None or b == 0.25 or b == 1 or b == 1.0/3.0:   # terminate for unusable b values
            print('Warning: method not defined.')
            sys.exit()
        else:
            sqrtexp = m.sqrt(144.0*(b**4) + 72.0*(b**3) - 99.0*(b**2) + 30.0*b - 3)
            b3 = b
            b2 = (-12.*(b**2)+15.0*b-3-sqrtexp)/(24.0*b-6)
            b1 = 1-b2-b3
            a3 = -(3.0*b2+3.*b3-1)/(6.0*(b3-1)*b2);
            a2 = (2.*a3*b3-2.*a3+1)/(2.*b1)
            a1 = 1-a2-a3
            alphas['OS33bv2'] = [[a1,b1],[a2,b2],[a3,b3]]
            #print(alphas['OS33bv2'])

    alphas['Emb3/2Ra'] = [[1, -1/24], [(-2/3, -12/25), (3/4, 25/24)], [(2/3, 12/25), (7/24, 0)]]
    embedded_methods['Emb3/2Ra'] = 2

    alphas['Godunov-N'] = [[1] * N]
    alphas['Strang-N'] = [ [ 0.5] * (N - 1) + [1],
                           [0.5] * (N - 1) + [0, True]]
    
    
    order = None
    if isinstance(alpha, str) and alpha in alphas: 
        os_scheme = alphas[alpha] # matching OS method name with coefficients
    if isinstance(alpha, str) and alpha in embedded_methods:
        order = embedded_methods[alpha]
    return os_scheme, order



adi_list=['MCS', 'HV', 'DR']
  





