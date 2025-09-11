# Multirate integrators

The multirate integrator may only be used with two operators. 

The relevant integrator function is `multirate_solve`. It takes some basic arguments to define the problem and method, and optional arguments may be used to control output and integrator behaviour

## Basic Usage

The `multirate_solve` function takes as inputs:

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
    The method to use, which is of type `Multirate`
- M
    the ratio of fast steps to slow steps, which must be an integer
- fs
    the slow operator (see [Operators](../README.md#operators))
- ff
    the fast operator (see [Operators](../README.md#operators))

## Method Definition

The method is defined by an object of type `Multirate`

The `Multirate` class takes as inputs:

- A_ff - a 2D numpy array
- A_ss - a 2D numpy array
- b_f - a 1D numpy array
- b_s - a 1D numpy array
- A_fs - a function with arguments (lambda, M) that returns a 2D numpy array (where lambda is the step index and M is the number of fast steps per slow step)
- A_sf - a function with arguments (lambda, M) that returns a 2D numpy array (where lambda is the step index and M is the number of fast steps per slow step)

One pre-defined method is included in the library, which is `mrgark_ex2_im2` - Sarshar, Roberts and Sandu (2019) MrGARK EX2-IM2 2(1)[A]

## Optional Arguments

### Method reordering

- order: a list containing the stages (0 to N-1) in the order that calculation is desired.

If supplied, this will be used. Otherwise, the integrator will attempt to reorder the tableaus into the most explicit (or diagonally implicit) structure possible

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

### Dynamic Linearization

- jacobian: A function to compute a jacobian when using dynamic linearization

    If this option is used, the provided operator functions must take as arguments (t, y, J)

When this is used, the jacobian function is called at the start of each step, and held constant for each substep within.  This jacobian is passed to each operator function, which may be uses to compute the values of each operator

### Boundary Condition

- bc: Any boundary conditions to apply when using firedrake

### Controlling Solver Behaviour

An optional dictionary of additional arguments may be supplied, through the `solver_parameters` argument.

If it is supplied, it is used as optional arguments for `scipy.optimize.root` if the operators are defined as functions, or as optional arguments for firedrake's `solve` function if operators are defined as `Form`s. 

For the firedrake finite element integrator, the `solve` function is applied to determine the new stage value (referred to as y), and then the stage derivatives (k=f(y)). The optional arguments are only applied to determining the stage values y.
