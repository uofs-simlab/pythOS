# Multirate Infinitesimal integrators

The multirate infinitesimal integrators may be used with two or three operators.

The relevant integrator function is `multirate_infinitesimal_solve`. It takes some basic arguments to define the problem and method, and optional arguments may be used to control output and integrator behaviour

## Basic Usage

The `multirate_infinitesimal_solve` takes as inputs:

- y0

    the initial value(s) to use for y (see [Initial Condition](../README.md#initial-condition))

- t0

    the time value to start at
    
    If using the finite element version, this must be of type Constant

- dt

    the time step to advance by, which may be an integer or a float

- tf

    the time value to advance to

- method
    The method to use, which is of type `Multirate_Infinitesimal`
- fi
    the slow operator (see [Operators](../README.md#operators))
    If two slow operators are in use, this is the implicit operator

- ff
    the fast operator (see [Operators](../README.md#operators))

- fe (optional)
    the slow explicit operator if two slow operators are in use (see [Operators](../README.md#operators))

## Method Definition

The `Multirate_Infinitesimal` class takes as inputs:

- c - A 1D numpy array
- gamma - A 3D numpy array defining the coupling
- omega (optional) - A 3D numpy array defining the coupling with the optional second slow operator

For convienience, a number of methods are already defined:
- mri\_kw3 - Knoth and Wolke (1997) order 3
- mri\_erk2a - Sandu (2019) mri-gark-erk22a
- mri\_erk2b - Sandu (2019) mri-gark-erk22b
- mri\_erk3 - Sandu (2019) mri-gark-erk33
- mri\_erk4 - Sandu (2019) mri-gark-erk45a
- mri\_irk2 - Sandu (2019) mri-gark-irk21a
- mri\_esdirk3a - Sandu (2019) mri-gark-esdirk34a
- mri\_sdirk3 - Sandu (2019) mri-gark-sdirk33a
- mri\_imex3 - Chinomona and Reynolds (2021) imex-mri-gark3a
- mri\_imex4 - Chinomona and Reynolds (2021) imex-mri-gark4

## Optional Arguments

### File Output

- fname: a filename to save intermediate results to.  
  
     When using the finite element capabilities of firedrake this is a 
     .h5 file accessible through the `CheckpointFile` interface from 
     `firedrake`. The file has attributes `/times/idx/` that store the time 
     of each entry and the attribute `/time/last_idx` indicates the last 
     valid index. The `Function`s are stored using corresponding indices. 
     Otherwise this is a .csv file containing time in the first entry of each 
     line, and the solution vector in the remaining entries

- save\_steps: the number of intermediate steps to save.  
  
     The default is to save steps at delta\_t time interval if a filename
     is provided.

### Boundary Condition

- bcs: Any boundary conditions to apply when using firedrake. These are used when evaluating the operators.

### Controlling Solver Behaviour

The adaptive solver may be controlled through the `ivp_options` dictionary and the `ivp_method` argument.  The default is RK45 from scipy for the non-finite element problems, and Dormand-Prince from the pythOS collection for problems using Firedrake.  The options available are any of the integrators included as [Adaptive Methods](fractional_step.md#adaptive-methods), excluding the MRI methods.  Only the implicit solver of the ARKStep methods are used.
The `ivp_options` dictionary supplies optional arguments to control the adaptive solver, including the tolerances. By default, the tolerances are set to `rtol=1e-10, atol=1e-12`.  The other options depend on the choice of adaptive method, and are detailed in the [Fractional Step documentation](fractional_step.md#options-to-control-sub-integrator-behaviour)

The implicit solver may be controlled through the `implicit_solve_options` argument.  If the problem is defined through Firedrake, these optional arguments are passed to the `solve` function from Firedrake for any implicit solve stages. This does include any boundary conditions to apply on those stages.  For problems not defined through Firedrake, these optional arguments are passed to `scipy.optimize.root`

For the firedrake finite element integrator, the `solve` function is used both to evaluate the operators, and to do any implicit solves required.  The boundary condition supplied to `multirate_infinitesimal_solve` is used in the evaluation, and any optional arguments are only used for the implicit solve.
