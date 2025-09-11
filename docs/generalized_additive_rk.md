# Generalized Additive Runge Kutta integrators

The generalized additive Runge Kutta integrator may be used with any number of operators. 

The relevant integrator function is `gark_solve`. It takes some basic arguments to define the problem and method, and optional arguments may be used to control output and integrator behaviour

## Basic Usage

The `gark_solve` function takes as inputs:

- the operators to use
         
    A list of operators (see [Operators](../README.md#operators)) in the order you wish to use them.

- dt

    the time step to advance by, which may be an integer or a float

- y0

    the initial value(s) to use for y (see [Initial Condition](../README.md#initial-condition))

- t0

    the time value to start at
    
    If using the finite element version, this must be of type Constant
- tf

    the time value to advance to

- A
	
	a list of lists containing the arrays a{i, j} defining the method, each as a numpy array

- b
	
	a list containing the vectors b{i} defining the method, each as a numpy array

## Method Definition

The method is defined by the list of lists of arrays a{i, j} (labelled A), and the list of vectors b{i} (labelled b)

The integrator takes the information from these arguments to create the additive Runge Kutta tableaus, and solve using the additive Runge Kutta solver.  The width of each array a{i, j} must be the same as the length of the vector b{i}, and the heights of arrays a{i, j} must match for constant j.

Optionally, a user-defined reordering may be supplied to specify the order to calculate the stages through the `order` argument.  If used, this is a list containing the stages in the order calculation is desired.  Otherwise, the integrator will attempt to reorder into the most explicit (or diagonally implicit) structure possible.

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

An optional dictionary of additional arguments may be supplied through the `solver_parameters` argument.

If it is supplied, it is used as optional arguments for `scipy.optimize.root` if the operators are defined as functions, or as optional arguments for firedrake's `solve` function if operators are defined as `Form`s. 

For the firedrake finite element integrator, the `solve` function is applied to determine the new stage value (referred to as y), and then the stage derivatives (k=f(y)). The optional arguments are only applied to determining the stage values y.
